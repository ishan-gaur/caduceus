import csv
from tqdm import tqdm
from Bio import SeqIO
import matplotlib.pyplot as plt
from constants import *
from utils import *

SAMPLE = "TEST"
OUTPUT_CSV = OUTPUT_DIR / f"{SAMPLE}_counts.csv"
OUTPUT_HIST = OUTPUT_DIR / f"{SAMPLE}_hist.png"
read_counts = {}

# trackers for monitoring errors and sanity checking the outputs
lens = []
bcs = set()
umis = set()
bcs_not_found_count = 0
bc_mismatch_count = 0
missing_10x_count = 0
missing_gfp_count = 0
umi2_cutoff_count = 0
missing_truseq_count = 0
count_bc_template_mismatch = 0
wrong_library_count = 0
count_success = 0
# bcs_truncated_count = 0 # will never be seen if we check for the bcs being the same first--might need to relax this assumption if not enough hits

n_records = len(SeqIO.index(get_sample_files(SAMPLE)[0], 'fastq'))
r1_reads = SeqIO.parse(get_sample_files(SAMPLE)[0], 'fastq')
r2_reads = SeqIO.parse(get_sample_files(SAMPLE)[1], 'fastq')
for r1_read, r2_read in tqdm(zip(r1_reads, r2_reads), total=n_records):
    # aligning/merging the reads
    start_read = r2_read.seq.reverse_complement()
    end_read = r1_read.seq
    start_read_bcs, start_read_bc_st, start_read_bc_end = find_barcodes(str(start_read))
    end_read_bcs, end_read_bc_st, end_read_bc_end = find_barcodes(str(end_read))

    if start_read_bcs is None or end_read_bcs is None:
        bcs_not_found_count += 1
        continue
    if len(start_read_bcs) != len(end_read_bcs) or not all([s == e for s, e in zip(start_read_bcs, end_read_bcs)]):
        bc_mismatch_count += 1
        continue

    merged = start_read[:start_read_bc_st] + end_read[end_read_bc_st:]

    if not BC_10X in merged:
        missing_10x_count += 1
        continue
    if not GFP_TAG in merged:
        missing_gfp_count += 1
        continue

    end_10x = merged.find(BC_10X) + len(BC_10X)
    st_gfp = merged.find(GFP_TAG)
    if merged[end_10x:st_gfp] != start_read[start_read_bc_st:start_read_bc_end]:
        count_bc_template_mismatch += 1
        continue

    if not all([s in TFBS for s in start_read_bcs]):
        wrong_library_count += 1
        continue

    bcs.update(start_read_bcs)
    lens.append(len(merged))

    if BC_TRUSEQ1 not in merged or BC_TRUSEQ2 not in merged:
        missing_truseq_count += 1
        continue

    end_truseq_1 = merged.find(BC_TRUSEQ1) + len(BC_TRUSEQ1)
    st_10x = merged.find(BC_10X)
    umi1 = merged[end_truseq_1:st_10x]

    end_gfp = st_gfp + len(GFP_TAG)
    st_truseq2 = merged.find(BC_TRUSEQ2)
    umi2 = merged[end_gfp:st_truseq2]

    umi = umi1 + umi2
    umis.add(umi)
    barcode = "_".join(start_read_bcs)
    if barcode not in read_counts:
        read_counts[barcode] = set()
    read_counts[barcode].add(umi)

    count_success += 1

read_counts = {k:len(v) for k, v in read_counts.items()}
csv.writer(OUTPUT_CSV.open("w")).writerows(read_counts.items())

print("reads missing barcodes", bcs_not_found_count)
print("reads with mismatched barcodes", bc_mismatch_count)
print("reads missing the 10x barcode", missing_10x_count)
print("reads missing gfp", missing_gfp_count)
print("barcodes not found between 10x and gfp", count_bc_template_mismatch)
print("reads with TFBS from the wrong library", wrong_library_count)
print("reads missing truseq bcs", missing_truseq_count)
print(count_success)
print(len(umis))
print(min(lens), sum(lens) / len(lens), max(lens))

plt.hist(read_counts.values())
plt.yscale("log")
plt.title(f"Histogram of UMIs per Barcode for {SAMPLE}")
plt.savefig(OUTPUT_HIST)