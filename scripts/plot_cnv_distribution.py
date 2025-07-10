import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import argparse

def load_cnv_bed_file(file_path):
    """Load CNV bed file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', header=None, 
                       names=['chrom', 'start', 'end', 'package', 'CN'])    

def plot_cnv_distribution_samples(samples, chromosome, output_file):
    """Plot genome-wide CNV distribution of specfific chromosome between samples."""
    plt.figure(figsize=(10, 6))
    
    for sample in samples:
        cnv_data = load_cnv_bed_file(sample)
        cnv_data = cnv_data[cnv_data['chrom'] == chromosome]
        
        if not cnv_data.empty:
            plt.plot(cnv_data['start'], cnv_data['CN'], label=sample)
    
    plt.title(f'CNV Distribution on {chromosome} between Samples')
    plt.xlabel('Genomic Position')
    plt.ylabel('Copy Number')
    plt.legend()
    plt.grid()
    plt.savefig(output_file)
    plt.close()
    
def plot_cnv_distribution_chromosomes(df, output_file):
    """Plot genome-wide CNV distribution across all chromosomes."""
    plt.figure(figsize=(12, 8))
    
    for chrom in df['chrom'].unique():
        chrom_data = df[df['chrom'] == chrom]
        plt.step(chrom_data['start'], chrom_data['CN'], label=chrom)
    
    plt.title('Genome-wide CNV Distribution across Chromosomes')
    plt.xlabel('Genomic Position')
    plt.ylabel('Copy Number')
    plt.legend()
    plt.grid()
    plt.savefig(output_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot CNV distribution from bed files.')
    parser.add_argument('--samples', nargs='+', help='List of sample CNV bed files')
    parser.add_argument('--chromosome', default='chr1', help='Chromosome to plot (default: chr1)')
    parser.add_argument('--output_directory', default='.', help='Output directory for plots')
    parser.add_argument('--output_samples', default='cnv_distribution_samples.png', help='Output file for sample CNV distribution')
    parser.add_argument('--output_chromosomes', default='cnv_distribution_chromosomes.png', help='Output file for chromosome CNV distribution')
    args = parser.parse_args() 
    samples = args.samples
    chromosome = args.chromosome
    output_dir = args.output_directory
    output_samples = args.output_samples
    output_chromosomes = args.output_chromosomes
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plot_cnv_distribution_samples(samples, chromosome, f'{output_dir}/{output_samples}')
    for sample in samples:
        cnv_data = load_cnv_bed_file(sample)
        print(cnv_data.head())
        match = re.search(r'/([^/]+)_alignment[^/]*\.bed$', sample)
        if match:
            sample_id = match.group(1)
            print(sample_id)
        else:
            print("No match found")
        plot_cnv_distribution_chromosomes(cnv_data, f'{output_dir}/{sample}_{output_chromosomes}')

if __name__ == "__main__":
    main()
    