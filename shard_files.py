import os
from tqdm import tqdm
from warnings import warn
import multiprocessing as mp
from utils import save_shard_starts, shard_file
from constants import get_sample_files, SAMPLES

SHARD_RECORD_CT = 100000

def process_sample(sample, position):
    if shard_file(sample).exists():
        print(f"Shards for {sample} already exist. Skipping...")
        return

    r1_file, r2_file = get_sample_files(sample)
    r1_size, r2_size = os.path.getsize(r1_file), os.path.getsize(r2_file)

    with open(r1_file, 'r') as r1, open(r2_file, 'r') as r2:
        r1_shard_start = []
        r2_shard_start = []
        curr_shard_record_ct = SHARD_RECORD_CT

        total_size = r1_size + r2_size
        progress_bar = tqdm(total=total_size, unit='B', unit_scale=True, desc=f"Processing Shards for {sample}", position=position)

        r1_line, r2_line = r1.readline(), r2.readline()
        while r1_line and r2_line:
            if r1_line[0] == '@':
                if r2_line[0] != '@':
                    warn(f"Files {r1_file} and {r2_file} are not in sync")
                    print("R1 line:", r1_line)
                    print("R2 line:", r2_line)
                    print("R1 byte:", r1.tell() - len(r1_line))
                    print("R2 byte:", r2.tell() - len(r2_line))
                    break
                if curr_shard_record_ct == SHARD_RECORD_CT:
                    r1_shard_start.append(r1.tell() - len(r1_line))
                    r2_shard_start.append(r2.tell() - len(r2_line))
                    curr_shard_record_ct = 0
                curr_shard_record_ct += 1
            bytes_read = len(r1_line) + len(r2_line)
            progress_bar.update(bytes_read)
            r1_line, r2_line = r1.readline(), r2.readline()

        progress_bar.close()

        if r1_line or r2_line:
            print(r1_line)
            print(r2_line)
            warn(f"Files {r1_file} and {r2_file} have different number of lines")
        
        save_shard_starts(sample, r1_shard_start, r2_shard_start)

if __name__ == "__main__":
    with mp.Pool(processes=min(len(SAMPLES), mp.cpu_count())) as pool:
        pool.starmap(process_sample, [(sample, i) for i, sample in enumerate(SAMPLES)])