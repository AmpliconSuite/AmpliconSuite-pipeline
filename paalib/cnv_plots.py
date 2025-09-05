#packages for plotting
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import os


# Plot CNV distribution from a CNV bed file
def load_cnv_bed_file(file_path):
    """Load CNV bed file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', header=None, 
                       names=['chrom', 'start', 'end', 'package', 'CN']) 


def load_centromere_file(aa_data_repo, ref_genome='GRCh38'):
    """Load centromere coordinates from AA_DATA_REPO."""
    file_path = os.path.join(aa_data_repo, ref_genome, f'{ref_genome}_centromere.bed')
    
    if os.path.exists(file_path):
        centromeres_raw = pd.read_csv(file_path, sep='\t', header=None)
        centromeres = centromeres_raw.iloc[:, :3].copy()
        centromeres.columns = ['chrom', 'start', 'end']
        centromeres['start'] = pd.to_numeric(centromeres['start'], errors='coerce')
        centromeres['end'] = pd.to_numeric(centromeres['end'], errors='coerce')
        centromeres = centromeres.dropna()
        return centromeres
    else:
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
            if max(8, chrom_data['CN'].max() * 1.1) == 8: #if max is 8
                ax.set_ylim(0, 8)
            elif (8 < max(8, chrom_data['CN'].max() * 1.1) < 50): #if max CN is less than 50
                ax.set_ylim(0, chrom_data['CN'].max()+5) #take max CN + 5   
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
    plt.savefig(output_file + ".png", dpi=300, bbox_inches='tight')
    plt.savefig(output_file + ".pdf", bbox_inches='tight')
    plt.close()
