import csv
from tqdm import tqdm
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from constants import OUTPUT_DIR, get_sample_files
from utils import parse_paired_read, get_shard_starts, counts_csv, counts_hist
import multiprocessing as mp
import os
from fire import Fire
import threading

CHUNK_SIZE = int(1e8)

def process_chunk(chunk_start_r1, chunk_start_r2, chunk_end_r1, chunk_end_r2, r1_file, r2_file):
    read_counts = {}
    with open(r1_file, 'r') as r1, open(r2_file, 'r') as r2:
        r1.seek(chunk_start_r1)
        r2.seek(chunk_start_r2)
        while r1.tell() < chunk_end_r1 and r2.tell() < chunk_end_r2:
            r1_read = Seq([r1.readline() for _ in range(4)][1])
            r2_read = Seq([r2.readline() for _ in range(4)][1])
            start_read = r2_read.reverse_complement()
            end_read = r1_read
            barcode, umi = parse_paired_read(start_read, end_read)
            if barcode is None:
                if None not in read_counts:
                    read_counts[None] = 0
                read_counts[None] += 1
                continue
            if barcode not in read_counts:
                read_counts[barcode] = set()
            read_counts[barcode].add(umi)
    return read_counts

# here output_tsv is actually a .csv file but in write_results_to_csv it is written as a tsv to accomodate the list
def write_results_to_csv(results, output_tsv):
    with output_tsv.open("a") as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        for res in results:
            read_counts = res.get()
            if None in read_counts:
                # n_bad_reads = read_counts[None]
                del read_counts[None]
                # print(f"Bad reads were {n_bad_reads / sum([len(v) for v in read_counts.values()]) * 100:.2f}% of total reads.")
            for barcode, umis in read_counts.items():
                writer.writerow([barcode, list(umis)])
                assert len(barcode) > 0

def merge_duplicates_in_csv(input_tsv):
    output_csv = input_tsv # the input tsv is actually a .csv file but in write_results_to_csv it is written as a tsv to accomodate the list
    merged_counts = {}
    with input_tsv.open("r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            barcode, umis = row
            umis_set = set(eval(umis))  # Convert string representation of list back to set
            if barcode in merged_counts:
                merged_counts[barcode].update(umis_set)
            else:
                merged_counts[barcode] = umis_set

    with output_csv.open("w") as csvfile:
        writer = csv.writer(csvfile)
        for barcode, umis_set in merged_counts.items():
            writer.writerow([barcode, len(umis_set)])

def main(sample, processes=None):
    print("Starting the analysis...")
    r1_file, r2_file = get_sample_files(sample)

    output_csv = counts_csv(sample)
    print(f"Writing results to {output_csv}...")
    if output_csv.exists():
        output_csv.unlink()

    pool = mp.Pool(mp.cpu_count() if processes is None else int(processes))
    starting_points = []

    print("Finding starting points for each chunk...")
    sample_file_size = os.path.getsize(r1_file)
    r1_starts, r2_starts = get_shard_starts(sample)
    starting_points = list(zip(r1_starts, r2_starts))
    starting_points.append((sample_file_size, sample_file_size))  # Add end of file as the last chunk end
    print("Starting points found.")

    results = []
    print("Processing chunks...")
    for i in range(len(starting_points) - 1):
        chunk_start_r1, chunk_start_r2 = starting_points[i]
        chunk_end_r1, chunk_end_r2 = starting_points[i + 1]
        results.append(pool.apply_async(process_chunk, (chunk_start_r1, chunk_start_r2, chunk_end_r1, chunk_end_r2, r1_file, r2_file)))


    def monitor_results():
        with tqdm(total=len(starting_points) - 1) as pbar:
            # while pbar.n < len(starting_points) - 1:
            while len(results) > 0:
                ready_results = [r for r in results if r.ready()]
                if ready_results:
                    write_results_to_csv(ready_results, output_csv)
                for r in ready_results:
                    results.remove(r)
                pbar.n += len(ready_results)
                pbar.refresh()
                # print(pbar.n)
                # print(len(results))

    monitor_thread = threading.Thread(target=monitor_results)
    monitor_thread.start()
    monitor_thread.join()

    pool.close()
    pool.join()
    print("Chunk processing complete.")

    print("Merging duplicates in the CSV...")
    merge_duplicates_in_csv(output_csv)
    print("Duplicates merged.")

    output_hist = counts_hist(sample)
    print(f"Generating histogram and saving to {output_hist}...")
    read_counts = {}
    with output_csv.open("r") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            barcode, umi_count = row
            read_counts[barcode] = int(umi_count)

    plt.hist(read_counts.values())
    plt.yscale("log")
    plt.title(f"Histogram of UMIs per Barcode for {sample}")
    plt.savefig(output_hist)
    print("Histogram saved.")

if __name__ == "__main__":
    Fire(main)
