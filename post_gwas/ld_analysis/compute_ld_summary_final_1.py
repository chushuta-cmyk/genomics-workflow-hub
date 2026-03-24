# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
import sys
import gzip

def main():
    ld_file = "results/cultivated/ld_analysis/cultivated_ld_decay_full.ld"
    out_file = "results/cultivated/ld_analysis/ld_decay_summary.tsv"
    
    bin_size = 1000
    max_dist = 2000000
    num_bins = max_dist // bin_size + 1
    sample_mod = 1000
    
    sum_r2 = [0.0] * num_bins
    count = [0] * num_bins
    
    total_lines = 0
    sampled_lines = 0
    
    with open(ld_file, 'r') as f:
        for line in f:
            total_lines += 1
            if total_lines % sample_mod != 0:
                continue
            if line.startswith("CHR_A") or line.startswith(" CHR_A"):
                continue
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                bp_a = int(parts[1])
                bp_b = int(parts[4])
                r2 = float(parts[6])
            except ValueError:
                continue
            
            dist = abs(bp_a - bp_b)
            if dist > max_dist:
                continue
            bin_idx = dist // bin_size
            sum_r2[bin_idx] += r2
            count[bin_idx] += 1
            sampled_lines += 1
            
            if sampled_lines % 100000 == 0:
                print(f"Processed {sampled_lines} sampled lines (total {total_lines})", file=sys.stderr)
    
    print(f"Finished processing. Total lines read: {total_lines}, sampled: {sampled_lines}", file=sys.stderr)
    
    with open(out_file, 'w') as fout:
        fout.write("distance_bin_start\tdistance_bin_end\tmean_r2\tpair_count\tsampling_factor\n")
        for i in range(num_bins):
            start = i * bin_size
            end = (i + 1) * bin_size - 1
            if count[i] > 0:
                mean_r2 = sum_r2[i] / count[i]
                fout.write(f"{start}\t{end}\t{mean_r2:.6f}\t{count[i]}\t{sample_mod}\n")
            else:
                fout.write(f"{start}\t{end}\tNA\t0\t{sample_mod}\n")
    
    print(f"Summary written to {out_file}", file=sys.stderr)

if __name__ == "__main__":
    main()