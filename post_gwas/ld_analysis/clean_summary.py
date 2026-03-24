# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import sys

input_file = "results/wild/ld_analysis/ld_decay_summary.tsv"
output_file = "results/wild/ld_analysis/ld_decay_summary_clean.tsv"

with open(input_file, 'r') as f:
    lines = f.readlines()

header = "distance_bin_start\tdistance_bin_end\tmean_r2\tpair_count\tsampling_factor\n"
# Find first header line
found_header = False
clean_lines = []
for line in lines:
    if line.startswith("distance_bin_start"):
        if not found_header:
            # Use this header
            clean_lines.append(header)
            found_header = True
        else:
            # skip duplicate header
            continue
    else:
        clean_lines.append(line)

# If no header found, add it
if not found_header:
    clean_lines.insert(0, header)

# Write cleaned file
with open(output_file, 'w') as fout:
    fout.writelines(clean_lines)

print(f"Cleaned summary written to {output_file}", file=sys.stderr)

# Replace original with cleaned
import shutil
shutil.move(output_file, input_file)