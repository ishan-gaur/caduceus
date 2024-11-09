# Caduceus
Get counts from paired-end reads fast in pure python
1. Run check_paired_files_match.py
2. sanity_check_samples.py this file is specific to my project, but might be useful to model some of your own debugging checks after
3. Once you're confident that your data looks as expected, run shard_files.py
4. This generates the chunk start/end positions for each file
5. It automatically runs on all your SAMPLES, a constant configured in constants.py
6. You'll probably want to go through and change the values there
7. Once the chunks are defined, read_counts.py actually extracts the deduplicated counts
8. You'll probably need to hunt down the places to add the parsing/filtering logic you need in utils.py
