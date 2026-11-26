#!/usr/bin/env python3
"""check_reads.py — 用真值读名（t<start>_<f|r>）核对 wga3 读模式 PAF

用法: check_reads.py reads.fa out.paf [tol=20]
报告: 映射数 / 位置正确率(±tol) / 链一致率 / 平均 identity
"""
import sys
from collections import defaultdict

reads, paf = sys.argv[1], sys.argv[2]
tol = int(sys.argv[3]) if len(sys.argv) > 3 else 20

truth = []   # 读序号 -> (起点, 链)
with open(reads) as f:
    for line in f:
        if line[0] == '>':
            name = line[1:].strip()
            body, strand = name.rsplit('_', 1)
            truth.append((int(body[1:]), strand))

hit = defaultdict(int)
tot = ok = okStrand = 0
idSum = 0.0
for line in open(paf):
    c = line.rstrip('\n').split('\t')
    if len(c) < 12:
        continue
    tot += 1
    try:
        ri = int(c[0][1:])
    except ValueError:
        continue
    if ri >= len(truth):
        continue
    tpos, tstrand = truth[ri]
    rstart = int(c[7])
    idSum += float([x for x in c[12:] if x.startswith('id:f=')][0][5:])
    if abs(rstart - tpos) <= tol or abs(rstart + 150 - tpos) <= tol:
        ok += 1
        if (c[4] == '-') == (tstrand == 'r'):
            okStrand += 1
    hit[ri] += 1

print(f"reads_total={len(truth)}  mapped={tot}  ({100.0*tot/max(1,len(truth)):.2f}%)")
if tot:
    print(f"pos_correct(±{tol})={ok}  ({100.0*ok/tot:.2f}%)  strand_ok={okStrand} ({100.0*okStrand/max(1,ok):.2f}%)")
    print(f"mean_identity={idSum/tot:.4f}  multi_mapped_reads={sum(1 for v in hit.values() if v>1)}")
