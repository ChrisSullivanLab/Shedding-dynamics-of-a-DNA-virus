#!/usr/bin/env python
# coding: utf-8

# In[5]:


from pathlib import Path
import re

# =========================
# Paths
# =========================
in_fa  = Path("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26/mmu.mirna.cleaned.fa")
out_fa = in_fa.with_name("mmu.mirna.cleaned_1to800.fa")



LOW, HIGH = 1, 800  # numeric range for mmu-miR-<number>

# =========================
# FASTA reader
# =========================
def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

def write_fasta(seqs, path):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")

def extract_mirna_number(name):
    m = re.search(r"^mmu-miR-(\d+)", name)
    return int(m.group(1)) if m else None

# =========================
# Filter
# =========================
mirnas = read_fasta(in_fa)

kept = {}
kept_mir = 0
kept_let = 0

for name, seq in mirnas.items():
    if name.startswith("mmu-let-"):
        kept[name] = seq
        kept_let += 1
        continue

    n = extract_mirna_number(name)
    if n is not None and LOW <= n <= HIGH:
        kept[name] = seq
        kept_mir += 1

# =========================
# Save + show
# =========================
write_fasta(kept, out_fa)

print(f"Input total      : {len(mirnas)}")
print(f"Kept mmu-miR {LOW}-{HIGH}: {kept_mir}")
print(f"Kept mmu-let-*    : {kept_let}")
print(f"Kept total       : {len(kept)}")
print(f"Output file      : {out_fa}")

print("\nFirst 10 kept headers:")
for i, name in enumerate(kept.keys()):
    if i == 10:
        break
    print(" ", name)


# In[6]:


from pathlib import Path
import pandas as pd

# =========================
# User settings (edit these)
# =========================
raw_dir = Path("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26")
barcode_fa = raw_dir / "cutoff_99pct_stock_barcodes_with_rev.fa"

mupyv_mirna_fa = raw_dir / "muPyV-mirna.burke.2018.fa"
mmu_mirna_fa   = raw_dir / "mmu.mirna.cleaned.fa"
#mmu_mirna_fa   = raw_dir / "mmu.mirna.cleaned_1to800.fa" 

seed_len = 6                 # customizable: classic is 6 (positions 2–7)
require_full_seed = True     # keep True for "matches completely"

# Output files (saved in raw_dir = where your FASTA files live)
out_mupyv = raw_dir / f"seedmatch_muPyV_seed{seed_len}_cutoffFile.tsv"
out_mmu   = raw_dir / f"seedmatch_mmu_seed{seed_len}_cutoffFile.tsv"


# =========================
# Helpers
# =========================
def read_fasta(path: Path) -> dict:
    """Return dict: {name: sequence} (sequence is concatenated, no spaces)."""
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).replace(" ", "")
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.replace(" ", ""))
        if name is not None:
            seqs[name] = "".join(buf).replace(" ", "")
    return seqs

def rna_to_dna(seq: str) -> str:
    return seq.replace("U", "T").replace("u", "t")

def rc_dna(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def count_overlapping(haystack: str, needle: str) -> int:
    """Count overlapping occurrences of needle in haystack."""
    if not needle:
        return 0
    count = 0
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1  # +1 allows overlaps
    return count

def miRNA_seed_rc_in_dna(miRNA_rna: str, seed_len: int) -> str:
    """
    Take miRNA positions 2..(1+seed_len) (1-based) => miRNA_rna[1:1+seed_len],
    convert U->T, then return reverse-complement in DNA alphabet.
    """
    if len(miRNA_rna) < 1 + seed_len:
        return ""
    seed_rna = miRNA_rna[1:1+seed_len]
    seed_dna = rna_to_dna(seed_rna)
    return rc_dna(seed_dna)

def make_seedmatch_table(barcodes: dict, mirnas: dict, seed_len: int) -> pd.DataFrame:
    rows = []
    for bc_name, bc_seq in barcodes.items():
        for mir_name, mir_seq_rna in mirnas.items():
            seed_rc = miRNA_seed_rc_in_dna(mir_seq_rna, seed_len)
            if not seed_rc:
                continue
            hits = count_overlapping(bc_seq, seed_rc)
            if require_full_seed and hits == 0:
                continue
            rows.append({
                "barcode_name": bc_name,
                "barcode_seq": bc_seq,
                "matched_miRNA": mir_name,
                "hit_count": hits,
                "seed_len": seed_len,
                "seed_rc_seq_dna": seed_rc,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["barcode_name", "hit_count", "matched_miRNA"], ascending=[True, False, True]).reset_index(drop=True)
    return df


# =========================
# Run
# =========================
barcodes = read_fasta(barcode_fa)

mupyv_mirnas = read_fasta(mupyv_mirna_fa)
mmu_mirnas   = read_fasta(mmu_mirna_fa)

df_mupyv = make_seedmatch_table(barcodes, mupyv_mirnas, seed_len=seed_len)
df_mmu   = make_seedmatch_table(barcodes, mmu_mirnas,   seed_len=seed_len)

# Save
df_mupyv.to_csv(out_mupyv, sep="\t", index=False)
df_mmu.to_csv(out_mmu, sep="\t", index=False)

print("Saved:")
print(" ", out_mupyv)
print(" ", out_mmu)

print("\nFirst 10 lines (muPyV):")
display(df_mupyv.head(10))

print("\nFirst 10 lines (mmu):")
display(df_mmu.head(10))


# In[15]:


from pathlib import Path
import pandas as pd

# =========================
# User settings (edit these)
# =========================
raw_dir = Path("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26")

# with_rev file (targets you're scanning)
barcode_fa_with_rev = raw_dir / "cutoff_99pct_stock_barcodes_with_rev.fa"

# ORIGINAL cutoff file (no -rev)  ✅ you need this
barcode_fa_original = raw_dir / "cutoff_99pct_stock_barcodes.fa"

mupyv_mirna_fa = raw_dir / "muPyV-mirna.burke.2018.fa"
mmu_mirna_fa   = raw_dir / "mmu.mirna.cleaned.fa"
#mmu_mirna_fa   = raw_dir / "mmu.mirna.cleaned_1to800.fa"

seed_len = 7
require_full_seed = True

out_mupyv = raw_dir / f"seedmatch_muPyV_seed{seed_len}_cutoffFile.tsv"
out_mmu   = raw_dir / f"seedmatch_mmu_seed{seed_len}_cutoffFile.tsv"


# =========================
# Helpers
# =========================
def read_fasta(path: Path) -> dict:
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).replace(" ", "")
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.replace(" ", ""))
        if name is not None:
            seqs[name] = "".join(buf).replace(" ", "")
    return seqs

def rna_to_dna(seq: str) -> str:
    return seq.replace("U", "T").replace("u", "t")

def rc_dna(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def count_overlapping(haystack: str, needle: str) -> int:
    if not needle:
        return 0
    count = 0
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count

def miRNA_seed_rc_in_dna(miRNA_rna: str, seed_len: int) -> str:
    if len(miRNA_rna) < 1 + seed_len:
        return ""
    seed_rna = miRNA_rna[1:1+seed_len]  # positions 2..(1+seed_len)
    seed_dna = rna_to_dna(seed_rna)
    return rc_dna(seed_dna)

def base_barcode_name(bc_name: str) -> str:
    """barcode123-rev -> barcode123 ; barcode123 -> barcode123"""
    return bc_name[:-4] if bc_name.endswith("-rev") else bc_name

def make_seedmatch_table(barcodes_with_rev: dict, barcodes_original: dict, mirnas: dict, seed_len: int) -> pd.DataFrame:
    rows = []
    for bc_name, bc_seq in barcodes_with_rev.items():
        orig_name = base_barcode_name(bc_name)
        orig_seq  = barcodes_original.get(orig_name, "")  # empty if not found

        for mir_name, mir_seq_rna in mirnas.items():
            seed_rc = miRNA_seed_rc_in_dna(mir_seq_rna, seed_len)
            if not seed_rc:
                continue
            hits = count_overlapping(bc_seq, seed_rc)
            if require_full_seed and hits == 0:
                continue
            rows.append({
                "barcode_name": bc_name,
                "barcode_seq": bc_seq,

                # ✅ new columns
                "original_barcode_name": orig_name,
                "original_barcode_seq": orig_seq,

                "matched_miRNA": mir_name,
                "hit_count": hits,
                "seed_len": seed_len,
                "seed_rc_seq_dna": seed_rc,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["barcode_name", "hit_count", "matched_miRNA"],
                            ascending=[True, False, True]).reset_index(drop=True)
    return df


# =========================
# Run
# =========================
barcodes_with_rev = read_fasta(barcode_fa_with_rev)
barcodes_original = read_fasta(barcode_fa_original)

mupyv_mirnas = read_fasta(mupyv_mirna_fa)
mmu_mirnas   = read_fasta(mmu_mirna_fa)

df_mupyv = make_seedmatch_table(barcodes_with_rev, barcodes_original, mupyv_mirnas, seed_len=seed_len)
df_mmu   = make_seedmatch_table(barcodes_with_rev, barcodes_original, mmu_mirnas,   seed_len=seed_len)

df_mupyv.to_csv(out_mupyv, sep="\t", index=False)
df_mmu.to_csv(out_mmu, sep="\t", index=False)

print("Saved:")
print(" ", out_mupyv)
print(" ", out_mmu)

print("\nFirst 10 lines (muPyV):")
display(df_mupyv.head(10))

print("\nFirst 10 lines (mmu):")
display(df_mmu.head(10))


# In[14]:


# Number of unique original barcodes (muPyV)
n_unique_mupyv = df_mupyv["original_barcode_name"].nunique()

# Number of unique original barcodes (mmu)
n_unique_mmu = df_mmu["original_barcode_name"].nunique()

print(f"Unique original barcodes (muPyV): {n_unique_mupyv}")
print(f"Unique original barcodes (mmu):   {n_unique_mmu}")

