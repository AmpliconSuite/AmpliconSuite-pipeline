"""
Data repository download functionality for AmpliconSuite-pipeline
"""
import os
import sys
import tarfile


def download_file(url, destination_folder):
    import urllib.request  # here because python2 not work with it
    filename = os.path.join(destination_folder, url.split("/")[-1])
    try:
        response = urllib.request.urlopen(url)
        file_size = int(response.headers.get('Content-Length', 0))
        response.close()
        file_size = round(file_size / (1024**3), 2)
        if file_size > 0.1:
            print("\nDownloading " + url + " ... (" + str(file_size) + "GB)")
        else:
            print("\nDownloading " + url + " ...")

        urllib.request.urlretrieve(url, filename)
        print("File downloaded and saved to: " + str(filename))
    except Exception as e:
        sys.stderr.write("Failed to download file. Error: " + str(e) + "\n")


def extract_tar_gz(file_path, destination_folder):
    if not file_path.endswith('.tar.gz'):
        sys.stderr.write("Cannot extract file " + file_path)
        sys.exit(1)

    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(destination_folder)

    os.remove(file_path)


def handle_repo_download(args, AA_REPO):
    """Handle data repository download and exit"""
    # Import download functions from your existing module

    data_repo_base_url = "https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/"
    
    for ref in args.download_repo:
        print("Downloading " + ref)
        ref_base_url = data_repo_base_url + ref
        md5file = ref_base_url + "_md5sum.txt"
        ref_file = ref_base_url + ".tar.gz"
        
        if os.path.exists(AA_REPO + ref):
            print("An AA data repo directory already exists for " + ref + " and it will be replaced!")
        
        download_file(md5file, AA_REPO)
        download_file(ref_file, AA_REPO)
        print("Extracting...\n")
        extract_tar_gz(AA_REPO + ref + ".tar.gz", AA_REPO)
    
    print("Finished")
    sys.exit(0)