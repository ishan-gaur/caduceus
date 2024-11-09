import os
from tqdm import tqdm
from warnings import warn
import multiprocessing as mp
from utils import file_format_check_summary
from constants import get_sample_files, SAMPLES
import fire

def process_sample(sample, position, lock):
    r1_file, r2_file = get_sample_files(sample)
    r1_size, r2_size = os.path.getsize(r1_file), os.path.getsize(r2_file)

    unmatched_r1_lines = set()
    unmatched_r2_lines = set()
    total_r1_at_lines = 0
    total_r2_at_lines = 0

    with open(r1_file, 'r') as r1, open(r2_file, 'r') as r2:
        total_size = r1_size + r2_size
        progress_bar = tqdm(total=total_size, unit='B', unit_scale=True, desc=f"Checking {sample}", position=position)

        r1_line, r2_line = r1.readline(), r2.readline()
        while r1_line or r2_line:
            if r1_line and r1_line[0] == '@':
                r1_line = r1_line.split(" ")[0]
                total_r1_at_lines += 1
                if r1_line in unmatched_r2_lines:
                    unmatched_r2_lines.remove(r1_line)
                else:
                    unmatched_r1_lines.add(r1_line)
            if r2_line and r2_line[0] == '@':
                r2_line = r2_line.split(" ")[0]
                total_r2_at_lines += 1
                if r2_line in unmatched_r1_lines:
                    unmatched_r1_lines.remove(r2_line)
                else:
                    unmatched_r2_lines.add(r2_line)

            bytes_read = (len(r1_line) if r1_line else 0) + (len(r2_line) if r2_line else 0)
            progress_bar.update(bytes_read)
            r1_line, r2_line = r1.readline(), r2.readline()

        progress_bar.close()

        summary = []
        with lock:
            if unmatched_r1_lines or unmatched_r2_lines:
                summary.append(f"Files {r1_file} and {r2_file} have unmatched lines")
                summary.append(f"Unmatched lines in r1: {len(unmatched_r1_lines)}")
                summary.append(f"Unmatched lines in r2: {len(unmatched_r2_lines)}")
                summary.append(f"Total '@' lines in r1: {total_r1_at_lines}")
                summary.append(f"Total '@' lines in r2: {total_r2_at_lines}")
            else:
                summary.append(f"Files {r1_file} and {r2_file} are in sync with {total_r1_at_lines} '@' lines in r1 and {total_r2_at_lines} '@' lines in r2")

            for line in summary:
                print(line)

            with open(file_format_check_summary(sample), 'w') as summary_file:
                for line in summary:
                    summary_file.write(line + '\n')

def main(sample=None):
    if sample:
        if sample not in SAMPLES:
            raise ValueError(f"Sample {sample} not found in samples list. Must be one of {SAMPLES}")
            return
        samples_to_process = [sample]
    else:
        samples_to_process = SAMPLES

    with mp.Manager() as manager:
        lock = manager.Lock()
        with mp.Pool(processes=min(len(samples_to_process), mp.cpu_count())) as pool:
            pool.starmap(process_sample, [(sample, i, lock) for i, sample in enumerate(samples_to_process)])

if __name__ == "__main__":
    fire.Fire(main)
