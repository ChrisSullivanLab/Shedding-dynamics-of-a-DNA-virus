#!/usr/bin/env python3
import sys
from pathlib import Path

def reverse_complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def read_fasta(file_path):
    seqs = {}
    name = None
    buf = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs

def write_fasta(seqs, out_path):
    with open(out_path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 add_reverse_complements.py input.fa")
        sys.exit(1)

    in_path = Path(sys.argv[1])
    out_path = in_path.with_name(in_path.stem + "_with_rev.fa")

    seqs = read_fasta(in_path)
    out_seqs = {}

    for name, seq in seqs.items():
        out_seqs[name] = seq
        out_seqs[f"{name}-rev"] = reverse_complement(seq)

    write_fasta(out_seqs, out_path)
    print(f"✅ Created {out_path} with {len(out_seqs)} total entries")

if __name__ == "__main__":
    main()
