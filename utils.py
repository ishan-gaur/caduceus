import re
import csv
import numpy as np
import pickle as pkl
from matplotlib import pyplot as plt
from constants import TFBS, BC_10X, GFP_TAG, BC_TRUSEQ1, BC_TRUSEQ2, OUTPUT_DIR

def find_barcodes(s):
    barcode_pattern = r'(GAGT(?:[ATGC]{4}GAGT)+)'
    try:
        st, end = [(match.start(), match.end()) for match in re.finditer(barcode_pattern, s)][0]
    except IndexError:
        return None, -1, -1
    barcodes = s[st:end].split('GAGT')
    assert barcodes[0] == ''
    assert barcodes[-1] == ''
    if len(barcodes) <= 2:
        return None, -1, -1
    return barcodes[1:-1], st, end

def parse_paired_read(start_read, end_read):
    start_read_bcs, start_read_bc_st, start_read_bc_end = find_barcodes(str(start_read))
    end_read_bcs, end_read_bc_st, end_read_bc_end = find_barcodes(str(end_read))

    if start_read_bcs is None or end_read_bcs is None:
        return None, None
    if len(start_read_bcs) != len(end_read_bcs):
        return None, None
    if not all([s == e for s, e in zip(start_read_bcs, end_read_bcs)]):
        return None, None

    merged = start_read[:start_read_bc_st] + end_read[end_read_bc_st:]

    if not BC_10X in merged:
        return None, None
    if not GFP_TAG in merged:
        return None, None

    end_10x = merged.find(BC_10X) + len(BC_10X)
    st_gfp = merged.find(GFP_TAG)
    if merged[end_10x:st_gfp] != start_read[start_read_bc_st:start_read_bc_end]:
        return None, None

    if not all([s in TFBS for s in start_read_bcs]):
        return None, None

    if BC_TRUSEQ1 not in merged or BC_TRUSEQ2 not in merged:
        return None, None

    end_truseq_1 = merged.find(BC_TRUSEQ1) + len(BC_TRUSEQ1)
    st_10x = merged.find(BC_10X)
    umi1 = merged[end_truseq_1:st_10x]

    end_gfp = st_gfp + len(GFP_TAG)
    st_truseq2 = merged.find(BC_TRUSEQ2)
    umi2 = merged[end_gfp:st_truseq2]

    barcode = "_".join(start_read_bcs)
    umi = umi1 + umi2

    return barcode, umi

shard_file = lambda sample: OUTPUT_DIR / f"{sample}_shards.pkl"

def get_shard_starts(sample):
    with open(shard_file(sample), 'rb') as f:
        r1, r2 = pkl.load(f)
    return r1, r2

def save_shard_starts(sample, r1, r2):
    with open(shard_file(sample), 'wb') as f:
        pkl.dump((r1, r2), f)

counts_csv = lambda sample: OUTPUT_DIR / f"{sample}_counts.csv"
counts_hist = lambda sample: OUTPUT_DIR / f"{sample}_hist.png"
def enrichment_csv(sample, pseudocount): 
    return OUTPUT_DIR / f"{sample}_enrichment{'_pc' if pseudocount else ''}.csv"
def enrichment_hist(sample, pseudocount): 
    return OUTPUT_DIR / f"{sample}_enrichment_hist{'_pc' if pseudocount else ''}.png"

def read_counts_csv_to_dict(file_path):
    with open(file_path, mode='r') as infile:
        reader = csv.reader(infile)
        return {rows[0]: int(rows[1]) for rows in reader}

file_format_check_summary = lambda sample: OUTPUT_DIR / f"{sample}_format_check_summary.txt"
enrichment_summary = lambda sample: OUTPUT_DIR / f"{sample}_enrichment.txt"
reads_check_summary = OUTPUT_DIR / f"reads_summary.txt"

def write_variants_to_csv(variants, output_file, log_fold_enrichments, reference_reads):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        header = ['Variant'] + [f'Pos{i}' for i in range(1, 12)] + [f'LFE_{sample}' for sample in log_fold_enrichments] + ['Counts in Library']
        writer.writerow(header)
        
        for variant in variants:
            position_split = variant.split("_")
            position_split += [' '] * (11 - len(position_split))
            row = [variant] + position_split
            row += [log_fold_enrichments[sample].get(variant, 'NA') for sample in log_fold_enrichments]
            row.append(reference_reads.get(variant, 'NA'))
            writer.writerow(row)

def get_sample_lfe(sample_names, lib_name):
    selection_samples = {sample: read_counts_csv_to_dict(counts_csv(sample)) for sample in sample_names}
    reference_reads = read_counts_csv_to_dict(counts_csv(lib_name))

    for sample in selection_samples:
        sample_total = sum(selection_samples[sample].values())
        selection_samples[sample] = {variant: count / sample_total for variant, count in selection_samples[sample].items()}

    reference_total = sum(reference_reads.values())
    reference_reads = {variant: count / reference_total for variant, count in reference_reads.items()}

    sample_lfe = {sample: {variant: np.log2(selection_samples[sample][variant] / reference_reads[variant] if variant in reference_reads else 1) for variant in selection_samples[sample]} for sample in selection_samples}
    return selection_samples, reference_reads, sample_lfe

def variant_len(variant):
    return len(variant.split("_"))

def plot_variant_pairwise_correlation(filtered_variants, sample_lfe, samples):
    fig, axes = plt.subplots(3, 3, figsize=(10, 10))
    for i, sample_a in enumerate(samples):
        for j, sample_b in enumerate(samples):
            axes[i, j].scatter(
                [sample_lfe[sample_a][variant] for variant in filtered_variants],
                [sample_lfe[sample_b][variant] for variant in filtered_variants]
            )
            axes[i, j].set_title(f"{sample_a} vs {sample_b}")
            axes[i, j].set_xlabel(f"{sample_a} log fold enrichment")
            axes[i, j].set_ylabel(f"{sample_b} log fold enrichment")
            # Fit a regression line
            x = np.array([sample_lfe[sample_a][variant] for variant in filtered_variants])
            y = np.array([sample_lfe[sample_b][variant] for variant in filtered_variants])
            coeffs = np.polyfit(x, y, 1)
            poly_eq = np.poly1d(coeffs)
            y_hat = poly_eq(x)
            
            # Plot the regression line
            axes[i, j].plot(x, y_hat, color='red')
            
            # Calculate and display the R^2 value
            ss_res = np.sum((y - y_hat) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            axes[i, j].text(0.05, 0.95, f"$R^2$ = {r2:.2f}", transform=axes[i, j].transAxes, 
                    verticalalignment='top', color='red')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "pairwise_correlation_filtered.png")