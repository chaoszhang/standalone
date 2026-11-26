#!/usr/bin/env python3
"""核对 MSA 的杂合判定：het 位点应尽量为 N，非 het 位点应尽量等于 ref。

用法: check_het.py ref.fa msa.fa het.bed
"""
import sys
from collections import Counter


def read_fa(p):
    d = {}
    name, buf = None, []
    for line in open(p):
        line = line.rstrip()
        if not line:
            continue
        if line[0] == '>':
            if name is not None:
                d[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        d[name] = ''.join(buf)
    return d


def main():
    ref = read_fa(sys.argv[1])
    msa = read_fa(sys.argv[2])
    het = {}
    for line in open(sys.argv[3]):
        f = line.split()
        if len(f) >= 4:
            het.setdefault(f[0], set()).add(int(f[1]))

    tot_het = Counter()
    tot_non = Counter()
    sample_non = []
    for name, rseq in ref.items():
        m = msa.get(name)
        if m is None or len(m) != len(rseq):
            print(f"[ERR] {name}: missing or length mismatch", file=sys.stderr)
            continue
        hs = het.get(name, set())
        for i, rc in enumerate(rseq):
            c = m[i]
            if i in hs:
                tot_het[c] += 1
            else:
                tot_non[c] += 1
                if c == 'N' and len(sample_non) < 5:
                    sample_non.append((name, i, rc, c))

    nh = sum(tot_het.values())
    nn = sum(tot_non.values())
    nN_h = tot_het.get('N', 0)
    nN_n = tot_non.get('N', 0)
    print(f"杂合位点 {nh}: N={nN_h} ({100.0*nN_h/max(1,nh):.2f}%)  其他={dict((k,v) for k,v in tot_het.items() if k!='N')}")
    print(f"非杂合位点 {nn}: N={nN_n} ({100.0*nN_n/max(1,nn):.3f}% 假杂合)  '-'={(tot_non.get('-',0))}")
    # 非 het 位点里与 ref 一致的比例（排除 N 和 '-'）
    for name, rseq in ref.items():
        m = msa[name]
        hs = het.get(name, set())
        same = tot = 0
        for i, rc in enumerate(rseq):
            c = m[i]
            if i in hs or c in 'N-':
                continue
            tot += 1
            if c == rc:
                same += 1
        if tot:
            print(f"  {name}: 非杂合已填 {tot} 列，与 ref 一致率 {same/tot:.4f}")
    if sample_non:
        print("  非杂合却判成 N 的样例:", sample_non)


if __name__ == '__main__':
    main()
