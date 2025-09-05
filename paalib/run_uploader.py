"""
Sample archiving and upload functionality for AmpliconSuite-pipeline
"""
import os
import tarfile
import logging
import subprocess


def archive_and_upload_sample(AA_outdir, 
                             AC_outdir, 
                             cnv_bed,
                             cnvkit_output_directory,
                             run_metadata_filename,
                             sample_metadata_filename,
                             amplified_interval_bed,
                             finish_flag_filename,
                             project_uuid,
                             project_key,
                             username,
                             outpre,
                             server='prod'):
    """
    Archive sample files and upload to AmpliconRepository
    
    Args:
        AA_outdir: AmpliconArchitect output directory (can be None)
        AC_outdir: AmpliconClassifier output directory (can be None) 
        cnv_bed: CNV bed file path (can be None)
        cnvkit_output_directory: CNVkit output directory (can be None)
        run_metadata_filename: Run metadata JSON file
        sample_metadata_filename: Sample metadata JSON file
        amplified_interval_bed: Amplified intervals bed file (can be None)
        finish_flag_filename: Finish flag file
        project_uuid: Project UUID for upload (ignored if server='None')
        project_key: Project key for upload (ignored if server='None')
        username: Username for upload (ignored if server='None')
        outpre: Sample name prefix for archive filename
        server: Server endpoint ('local', 'dev', 'prod')
        
    Returns:
        bool: True if successful, False otherwise
    """
    
    # Determine server URL
    server_urls = {
        'local': 'http://localhost:8000',
        'dev': 'https://dev.ampliconrepository.org',
        'prod': 'https://ampliconrepository.org',
    }
    
    if server not in server_urls:
        logging.error("Invalid server option: {}. Must be one of: {}".format(server, list(server_urls.keys())))
        return False
    
    server_url = server_urls[server]
    
    # Create temporary archive
    archive_path = None
    try:
        archive_path = _create_sample_archive(
            AA_outdir=AA_outdir,
            AC_outdir=AC_outdir,
            cnv_bed=cnv_bed,
            cnvkit_output_directory=cnvkit_output_directory,
            run_metadata_filename=run_metadata_filename,
            sample_metadata_filename=sample_metadata_filename,
            amplified_interval_bed=amplified_interval_bed,
            finish_flag_filename=finish_flag_filename,
            output_prefix=outpre
        )
        
        if not archive_path:
            logging.error("Failed to create sample archive")
            return False
        
        # Upload the archive
        success = _upload_archive(
            archive_path=archive_path,
            server_url=server_url,
            project_uuid=project_uuid,
            project_key=project_key,
            username=username
        )
        
        return success
        
    except Exception as e:
        logging.error("Error during archive and upload: {}".format(str(e)))
        return False
        
    finally:
        # Clean up temporary archive only if we uploaded it
        if archive_path and os.path.exists(archive_path):
            try:
                os.remove(archive_path)
                logging.info("Cleaned up temporary archive: {}".format(archive_path))
            except OSError as e:
                logging.warning("Could not remove temporary archive {}: {}".format(archive_path, str(e)))


def _create_sample_archive(AA_outdir,
                          AC_outdir, 
                          cnv_bed,
                          cnvkit_output_directory,
                          run_metadata_filename,
                          sample_metadata_filename,
                          amplified_interval_bed,
                          finish_flag_filename,
                          output_prefix):
    """
    Create a tar.gz archive of sample files, excluding large sequencing files
    
    Args:
        output_prefix: Full path prefix for archive (e.g., "/path/to/sample_name")
    
    Returns:
        str: Path to created archive, or None if failed
    """
    
    # Extensions to exclude
    excluded_extensions = {'.bam', '.cram', '.bai', '.crai', '.fq', '.fastq', 
                          '.fasta', '.fai', '.fa', '.cnr.gz'}
    
    # Create archive file path
    archive_path = "{}.tar.gz".format(output_prefix)
    
    # Make sure output directory exists
    output_dir = os.path.dirname(archive_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        with tarfile.open(archive_path, 'w:gz') as tar:
            
            # Always include metadata files
            for filepath in [run_metadata_filename, sample_metadata_filename, finish_flag_filename]:
                if filepath and os.path.exists(filepath):
                    arcname = os.path.basename(filepath)
                    tar.add(filepath, arcname=arcname)
                    logging.debug("Added to archive: {}".format(arcname))
            
            # Add amplified interval bed if it exists
            if amplified_interval_bed and os.path.exists(amplified_interval_bed):
                arcname = os.path.basename(amplified_interval_bed)
                tar.add(amplified_interval_bed, arcname=arcname)
                logging.debug("Added to archive: {}".format(arcname))
            
            # Handle CNV files - prefer cnvkit_output_directory over cnv_bed
            if cnvkit_output_directory and os.path.exists(cnvkit_output_directory):
                _add_directory_to_archive(tar, cnvkit_output_directory, "cnvkit_output", excluded_extensions)
            elif cnv_bed and os.path.exists(cnv_bed):
                arcname = os.path.basename(cnv_bed)
                tar.add(cnv_bed, arcname=arcname)
                logging.debug("Added to archive: {}".format(arcname))
            
            # Add AA output directory if it exists
            if AA_outdir and os.path.exists(AA_outdir):
                _add_directory_to_archive(tar, AA_outdir, "AA_results", excluded_extensions)
            
            # Add AC output directory if it exists  
            if AC_outdir and os.path.exists(AC_outdir):
                _add_directory_to_archive(tar, AC_outdir, "AC_results", excluded_extensions)
        
        logging.info("Created sample archive: {}".format(archive_path))
        return archive_path
        
    except Exception as e:
        logging.error("Failed to create archive: {}".format(str(e)))
        if os.path.exists(archive_path):
            try:
                os.remove(archive_path)
            except OSError:
                pass
        return None


def _add_directory_to_archive(tar, 
                             directory, 
                             archive_dirname,
                             excluded_extensions):
    """
    Add a directory to the archive, excluding files with certain extensions
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            filepath = os.path.join(root, file)
            
            # Check if file should be excluded
            if _should_exclude_file(file, excluded_extensions):
                logging.debug("Excluding file: {}".format(file))
                continue
            
            # Create archive path preserving directory structure
            rel_path = os.path.relpath(filepath, directory)
            arcname = os.path.join(archive_dirname, rel_path)
            
            tar.add(filepath, arcname=arcname)
            logging.debug("Added to archive: {}".format(arcname))


def _should_exclude_file(filename, excluded_extensions):
    """
    Check if a file should be excluded based on its extension
    """
    filename_lower = filename.lower()
    
    for ext in excluded_extensions:
        if filename_lower.endswith(ext):
            return True
        # Handle special cases like .fq.gz, .fastq.gz
        if ext in ['.fq', '.fastq'] and (filename_lower.endswith('{}.gz'.format(ext)) or filename_lower.endswith('{}.bz2'.format(ext))):
            return True
    
    return False


def _upload_archive(archive_path,
                   server_url, 
                   project_uuid,
                   project_key,
                   username):
    """
    Upload archive to AmpliconRepository using curl
    
    Returns:
        bool: True if upload successful, False otherwise
    """
    
    upload_url = "{}/add_samples_to_project_api/".format(server_url)
    
    # Construct curl command
    curl_cmd = [
        'curl', '-X', 'POST',
        upload_url,
        '-F', 'project_uuid={}'.format(project_uuid),
        '-F', 'project_key={}'.format(project_key), 
        '-F', 'username={}'.format(username),
        '-F', 'file=@{}'.format(archive_path)
    ]
    
    try:
        logging.info("Uploading archive to {}...".format(server_url))
        logging.debug("Upload command: {}".format(' '.join(curl_cmd)))
        
        result = subprocess.run(curl_cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            logging.info("Archive uploaded successfully")
            logging.debug("Upload response: {}".format(result.stdout))
            return True
        else:
            logging.error("Upload failed with return code {}".format(result.returncode))
            logging.error("Upload error output: {}".format(result.stderr))
            return False
            
    except subprocess.TimeoutExpired:
        logging.error("Upload timed out after 5 minutes")
        return False
    except Exception as e:
        logging.error("Upload failed with exception: {}".format(str(e)))
        return False