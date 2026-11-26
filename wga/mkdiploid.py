#!/usr/bin/env python3
"""从单倍体 ref 造两个 haplotype（含指定杂合率），用于验证 MSA 的杂合 → N 行为。

用法: mkdiploid.py ref.fa hap1.fa hap2.fa het.bed hetRate [seed]
het.bed 每行: 记录名<TAB>位置(0-based)<TAB>hap1碱基<TAB>hap2碱基
"""
import random
import sys


def read_fa(p):
    recs = []
    name, buf = None, []
    for line in open(p):
        line = line.rstrip()
        if not line:
            continue
        if line[0] == '>':
            if name is not None:
                recs.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def main():
    ref_p, h1_p, h2_p, bed_p = sys.argv[1:5]
    het = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    rng = random.Random(seed)

    recs = read_fa(ref_p)
    h1_out, h2_out, bed = [], [], []
    for name, seq in recs:
        a, b = list(seq), list(seq)
        n = 0
        for i, c in enumerate(seq):
            if c not in 'ACGT':
                continue
            if rng.random() >= het:
                continue
            others = [x for x in 'ACGT' if x != c]
            x, y = rng.sample(others, 2)
            a[i], b[i] = x, y
            bed.append(f"{name}\t{i}\t{x}\t{y}")
            n += 1
        h1_out.append((name, ''.join(a)))
        h2_out.append((name, ''.join(b)))
        print(f"[het] {name}: {n} 位点 ({100.0*n/max(1,len(seq)):.3f}%)", file=sys.stderr)

    for path, rs in ((h1_p, h1_out), (h2_p, h2_out)):
        with open(path, 'w') as f:
            for name, s in rs:
                f.write(f">{name}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + "\n")
    with open(bed_p, 'w') as f:
        f.write("\n".join(bed) + "\n")
    print(f"[het] 写出 {len(bed)} 个杂合位点 → {bed_p}", file=sys.stderr)


if __name__ == '__main__':
    main()
