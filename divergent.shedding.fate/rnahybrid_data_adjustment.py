
import re
import sys
from pathlib import Path
import pandas as pd


RNAHYBRID_TXT = "muPyV.rnahybrid.kidney_mirna.txt"
BARCODES_FASTA = "kidney_barcodes_unique_with_rev.fa"
MIRNA_FASTA    = "muPyV-mirna.burke.2018.fa"
OUT_TSV        = "muPyV_rev_rnahybrid.kidney_mirna.parsed.tsv"
MFE_MAX        = 0.0   # keep hits with mfe <= this (more negative = stronger)
PVAL_MAX       = 1    # keep hits with p-value <= this
SEED_LEN       = 6       # classic seed 2–7


def read_fasta(path):
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
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.replace(" ", ""))
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

def rc_dna(s):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]

def rna_to_dna(s):
    return s.replace("U","T").replace("u","t")

def seed_rc_in_barcode(miRNA_seq, barcode_seq, seed_len=6):
    """
    Check if the reverse complement of miRNA positions 2..(1+seed_len) (1-based)
    exists in the barcode (DNA). Returns (bool, seed_rc_string).
    """
    if len(miRNA_seq) < 1 + seed_len:
        return (False, "")
    seed_rna = miRNA_seq[1:1+seed_len]          # positions 2..7 for seed_len=6
    seed_dna = rna_to_dna(seed_rna)             # U->T
    seed_rc  = rc_dna(seed_dna)                 # reverse complement in DNA alphabet
    return (seed_rc in barcode_seq, seed_rc)

def parse_rnahybrid_blocks(path):
    lines = Path(path).read_text().splitlines()
    starts = [i for i,l in enumerate(lines) if l.startswith("target:")]
    records = []
    for i, si in enumerate(starts):
        ei = starts[i+1] if i+1 < len(starts) else len(lines)
        block = lines[si:ei]

        
        def grab(pattern):
            for ln in block:
                m = re.search(pattern, ln)
                if m: return m.group(1)
            return None

        target   = grab(r"^target:\s*(\S+)")
        mirna    = grab(r"^miRNA\s*:\s*(\S+)")
        mfe_str  = grab(r"^mfe:\s*([-\d\.]+)")
        pval_str = grab(r"^p-value:\s*([\d\.]+)")
        pos_str  = grab(r"^position\s+(\d+)")

        try:
            mfe  = float(mfe_str) if mfe_str is not None else None
            pval = float(pval_str) if pval_str is not None else None
            pos  = int(pos_str) if pos_str is not None else None
        except ValueError:
            mfe = pval = pos = None

        records.append({
            "target": target,
            "mirna": mirna,
            "mfe": mfe,
            "p_value": pval,
            "position_on_target_5prime_1based": pos
        })
    return pd.DataFrame.from_records(records)


if __name__ == "__main__":
    
    df = parse_rnahybrid_blocks(RNAHYBRID_TXT)
    barcodes = read_fasta(BARCODES_FASTA)    # dict: bcID -> DNA seq
    mirnas   = read_fasta(MIRNA_FASTA)       # dict: miRID -> RNA seq

    
    df["barcode_seq"] = df["target"].map(barcodes).fillna("")
    df["mirna_seq_rna"] = df["mirna"].map(mirnas).fillna("")

    
    has_seed_list, seed_rc_list = [], []
    for _, row in df.iterrows():
        has_seed, seed_rc = seed_rc_in_barcode(row["mirna_seq_rna"], row["barcode_seq"], seed_len=SEED_LEN)
        has_seed_list.append(has_seed)
        seed_rc_list.append(seed_rc)
    df["seed2_7_rc_in_barcode"] = has_seed_list
    df["seed_rc_seq"] = seed_rc_list

    # Filters
    df_filt = df.copy()
    if MFE_MAX is not None:
        df_filt = df_filt[df_filt["mfe"].notna() & (df_filt["mfe"] <= MFE_MAX)]
    if PVAL_MAX is not None:
        df_filt = df_filt[df_filt["p_value"].notna() & (df_filt["p_value"] <= PVAL_MAX)]

    
    df_filt.to_csv(OUT_TSV, sep="\t", index=False)

    #quick summary 
    print(f"Total pairs: {len(df)}")
    print(f"Kept after filters (mfe<={MFE_MAX}, p<={PVAL_MAX}): {len(df_filt)}")
    print("Example rows:")
    print(df_filt.head(10).to_string(index=False))
