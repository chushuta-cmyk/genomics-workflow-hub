#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path

# Test loading one GWAS file
WILD_GWAS_DIR = Path("/data03/karama/soybean_gwas_filtered/workflow/04_gemma_wild")
file_path = WILD_GWAS_DIR / "wild_trait_1.assoc.txt"
print(f"Testing loading: {file_path}")

# Load just first 1000 rows to test
df = pd.read_csv(file_path, sep='\t', nrows=1000)
print(f"Loaded {len(df)} rows")
print(f"Columns: {df.columns.tolist()}")
print(f"First few p-values: {df['p_wald'].head().tolist()}")

# Test plotting
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(df['ps'], -np.log10(df['p_wald']), s=5)
ax.set_xlabel('Position')
ax.set_ylabel('-log10(p)')
ax.set_title('Test Manhattan')
plt.savefig('/tmp/test_manhattan.png')
print("Plot saved to /tmp/test_manhattan.png")