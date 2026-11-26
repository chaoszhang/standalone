#!/usr/bin/env python3
# uni.py — 计算 PAF 的 ref 侧去重叠覆盖（按记录合并区间）
import sys
from collections import defaultdict

iv = defaultdict(list)
n = 0
for line in open(sys.argv[1]):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12:
        continue
    n += 1
    iv[f[5]].append((int(f[7]), int(f[8])))

tot = 0
for k, v in iv.items():
    v.sort()
    a = b = 0
    for s, e in v:
        if s > b:
            if b > a:
                tot += b - a
            a, b = s, e
        else:
            b = max(b, e)
    if b > a:
        tot += b - a

denom = int(sys.argv[2]) if len(sys.argv) > 2 else 1046932099
print('blocks=%d  ref_union_coverage=%d bp (%.2f%%)' % (n, tot, 100.0 * tot / denom))
