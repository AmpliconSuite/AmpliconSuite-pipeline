import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import argparse
import math

def load_cnv_bed_file(file_path):
    """Load CNV bed file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', header=None, 
                       names=['chrom', 'start', 'end', 'package', 'CN']) 


def load_centromere_file(aa_data_repo, ref_genome='GRCh38'):
    """Load centromere coordinates from AA_DATA_REPO."""
    file_path = os.path.join(aa_data_repo, ref_genome, f'{ref_genome}_centromere.bed')
    
    if os.path.exists(file_path):
        try:
            centromeres_raw = pd.read_csv(file_path, sep='\t', header=None)
            centromeres = centromeres_raw.iloc[:, :3].copy()
            centromeres.columns = ['chrom', 'start', 'end']
            centromeres['start'] = pd.to_numeric(centromeres['start'], errors='coerce')
            centromeres['end'] = pd.to_numeric(centromeres['end'], errors='coerce')
            centromeres = centromeres.dropna()
            return centromeres
        except Exception as e:
            print(f"Error loading centromere file {file_path}: {e}")
            return None
    else:
        print(f"Warning: No centromere file found in {aa_data_repo}/{ref_genome}/")
        return None

def highlight_centromere_regions(ax, chrom, centromeres, debug=False):
    """Add vertical light grey bars to highlight centromere regions."""
    if centromeres is None:
        return

    chrom_variations = [chrom]
    chrom_variations.append(chrom[3:])  #Remove 'chr' prefix
    chrom_centromeres = pd.DataFrame()
    for chrom_var in chrom_variations:
        chrom_centromeres = centromeres[centromeres['chrom'] == chrom_var]
        if not chrom_centromeres.empty:
            break

    for i, centro in chrom_centromeres.iterrows():
        #Add a vertical span covering the centromere region
        ax.axvspan(centro['start'], centro['end'], alpha=0.3, color='lightgrey', zorder=0)

def chrom_sort_key(chrom):
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == 'X':
        return (0, 23)
    elif chrom == 'Y':
        return (0, 24)
    elif chrom.isdigit():
        return (0, int(chrom))
    else:
        return (1, chrom)
    
def closest_divisors(n):
    a = round(math.sqrt(n))
    while n % a > 0:
        a -= 1
    return a, n // a

def plot_segments_as_lines(ax, chrom_data):
    """
    Plot each CNV segment as an individual horizontal line segment.
    """
    chrom_data = chrom_data.sort_values(by='start')
    for i, row in chrom_data.iterrows():
        start = row['start']
        end = row['end']
        cn = row['CN']
        ax.plot([start, end], [cn, cn], color='black', linewidth=1.5, solid_capstyle='butt', zorder=1)
    return ax
    
def plot_cnv_distribution_chromosomes(df, sample_name, output_file, centromeres=None, log_base=2):
    """
    Plot CNV profiles in a grid of subplots — one per chromosome, log-transformed CN values.
    Highlights centromere regions with light grey bars.
    """
    #Sort chromosome labels
    df['chrom'] = df['chrom'].astype(str)
    chromosomes = sorted(df['chrom'].unique(), key=chrom_sort_key)

    #Force CN numeric and log-transform
    df['CN'] = pd.to_numeric(df['CN'], errors='coerce')
    df['log_CN'] = np.log(df['CN']+1) / np.log(log_base)

    #automate setting nrows and ncoles based on number of chromosomes
    n1, n2 = closest_divisors(len(chromosomes))
    fig, axes = plt.subplots(nrows=min(n1, n2), ncols=max(n1, n2), figsize=(18, 12), sharey=False)
    axes = axes.flatten()
    for i, chrom in enumerate(chromosomes):   
        ax = axes[i]
        chrom_data = df[df['chrom'] == chrom].sort_values(by='start')
        
        #highlight centromere regions
        highlight_centromere_regions(ax, chrom, centromeres, debug=(i==0))
        #plot CNV segments
        plot_segments_as_lines(ax, chrom_data)
        
        ax.set_title(f'{chrom}', fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_xlim(left=0)
        if not chrom_data.empty:
            if (max(8, chrom_data['CN'].max() * 1.1) < 50): #if max CN is less than 50
                ax.set_ylim(0, max(8, chrom_data['CN'].max() * 1.1)) #take max between 8 and CN
            else:
                ax.set_ylim(0, 50) #if CN exceeds 50, cut off at 50
        else:
            ax.set_ylim(0, 8)

    #Hide unused subplots
    for j in range(len(chromosomes), len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"{sample_name} Genome-wide CNV Profiles (Centromeres Highlighted)", fontsize=16)
    fig.text(0.5, 0.04, "Genomic Position (bp)", ha="center", fontsize=12)
    fig.text(0.06, 0.5, f"Copy Number", va="center", rotation="vertical", fontsize=12)

    plt.tight_layout(rect=[0.07, 0.07, 1, 0.93])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot CNV distribution from bed files.')
    parser.add_argument('--samples', nargs='+', help='List of sample CNV bed files')
    parser.add_argument('--chromosome', default='chr1', help='Chromosome to plot (default: chr1)')
    parser.add_argument('--output_directory', default='.', help='Output directory for plots')
    parser.add_argument('--output_samples', default='cnv_distribution_samples.png', help='Output file for sample CNV distribution')
    parser.add_argument('--output_chromosomes', default='cnv_distribution_chromosomes.png', help='Output file for chromosome CNV distribution')
    parser.add_argument('--aa_data_repo', default=os.environ.get('AA_DATA_REPO', ''), help='Path to AA_DATA_REPO (default: from environment variable)')
    parser.add_argument('--ref_genome', default='GRCh38', help='Reference genome (default: GRCh38)')
    parser.add_argument('--highlight_centromeres', action='store_true', default=True, help='Highlight centromere regions with grey bars (default: True)')
    parser.add_argument('--no_highlight_centromeres', action='store_true', help='Do not highlight centromere regions (overrides --highlight_centromeres)')
    
    args = parser.parse_args() 
    samples = args.samples
    chromosome = args.chromosome
    output_dir = args.output_directory
    output_samples = args.output_samples
    output_chromosomes = args.output_chromosomes
    aa_data_repo = args.aa_data_repo
    ref_genome = args.ref_genome
    
    # Handle centromere highlighting logic
    highlight_centromeres = args.highlight_centromeres and not args.no_highlight_centromeres
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load centromere coordinates if highlighting centromeres
    centromeres = None
    if highlight_centromeres:
        if not aa_data_repo:
            print("Warning: AA_DATA_REPO not found. Cannot highlight centromeres.")
        else:
            centromeres = load_centromere_file(aa_data_repo, ref_genome)
            if centromeres is not None:
                print(f"Loaded {len(centromeres)} centromere regions for highlighting")
    
    # Plot individual sample chromosome profiles
    for sample in samples:
        cnv_data = load_cnv_bed_file(sample)
        print("CNV data sample:")
        print(cnv_data.head())
        print("Unique chromosomes in CNV data:", sorted(cnv_data['chrom'].unique()))
        match = re.search(r'/([^/]+)_alignment[^/]*\.bed$', sample)
        if match:
            sample_id = match.group(1)
            print(sample_id)
        else:
            print("No match found")
            sample_id = os.path.basename(sample).replace('.bed', '')
        plot_cnv_distribution_chromosomes(cnv_data, sample_id, f'{output_dir}/{sample_id}_{output_chromosomes}', centromeres)

if __name__ == "__main__":
    main()
