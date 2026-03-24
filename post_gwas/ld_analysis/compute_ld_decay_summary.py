# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import sys
import math

def main():
    ld_file = "results/wild/ld_analysis/wild_ld_decay_full.ld"
    out_file = "results/wild/ld_analysis/ld_decay_summary.tsv"
    
    # Define bin size (1 kb = 1000 bp)
    bin_size = 1000
    # Max distance to consider (maybe 1 Mb)
    max_dist = 1000000
    num_bins = max_dist // bin_size + 1
    
    # Initialize accumulators
    sum_r2 = [0.0] * num_bins
    count = [0] * num_bins
    
    # Process file
    with open(ld_file, 'r') as f:
        line_count = 0
        for line in f:
            # Skip header lines that contain column names
            if line.startswith("CHR_A") or line.startswith(" CHR_A"):
                continue
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            # Try to parse positions and R2
            try:
                bp_a = int(parts[1])
                bp_b = int(parts[4])
                r2 = float(parts[6])
            except ValueError:
                # skip malformed lines
                continue
            
            dist = abs(bp_a - bp_b)
            if dist > max_dist:
                continue
            bin_idx = dist // bin_size
            sum_r2[bin_idx] += r2
            count[bin_idx] += 1
            
            line_count += 1
            if line_count % 10000000 == 0:
                print(f"Processed {line_count} lines", file=sys.stderr)
    
    # Write output
    with open(out_file, 'w') as fout:
        fout.write("distance_bin_start\tdistance_bin_end\tmean_r2\tpair_count\n")
        for i in range(num_bins):
            if count[i] > 0:
                mean_r2 = sum_r2[i] / count[i]
                start = i * bin_size
                end = (i + 1) * bin_size - 1
                fout.write(f"{start}\t{end}\t{mean_r2:.6f}\t{count[i]}\n")
            else:
                start = i * bin_size
                end = (i + 1) * bin_size - 1
                fout.write(f"{start}\t{end}\tNA\t0\n")
    
    print(f"Summary written to {out_file}", file=sys.stderr)

if __name__ == "__main__":
    main()