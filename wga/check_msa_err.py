#!/usr/bin/env python3
"""按"输出碱基必须是该位点真实出现过的等位基因"的口径核对 MSA。

- 非杂合位点：真实等位基因只有 ref 那一个
- 杂合位点：真实等位基因是 hap1/hap2 那两个
- N 与 '-' 不算错误（判不了就退化，不构成"造碱基"）
- 错误 = 输出了一个从未在该位点出现过的碱基

用法: check_msa_err.py ref.fa msa.fa het.bed
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
            het.setdefault(f[0], {})[int(f[1])] = (f[2], f[3])

    tot = err = nN = nDash = 0
    err_het = err_hom = 0
    het_tot = 0
    wrong_bases = Counter()
    examples = []
    for name, rseq in ref.items():
        m = msa.get(name)
        if m is None or len(m) != len(rseq):
            print(f"[ERR] {name}: missing or length mismatch", file=sys.stderr)
            continue
        hs = het.get(name, {})
        for i, rc in enumerate(rseq):
            c = m[i]
            tot += 1
            if c == 'N':
                nN += 1
                continue
            if c == '-':
                nDash += 1
                continue
            if i in hs:
                het_tot += 1
                ok = {hs[i][0], hs[i][1]}
                if c not in ok:
                    err += 1
                    err_het += 1
                    wrong_bases[c] += 1
                    if len(examples) < 5:
                        examples.append((name, i, sorted(ok), c))
            else:
                if c != rc:
                    err += 1
                    err_hom += 1
                    wrong_bases[c] += 1
                    if len(examples) < 5:
                        examples.append((name, i, [rc], c))

    filled = tot - nDash
    print(f"总列 {tot}  已填 {filled} ({100.0*filled/tot:.2f}%)  N {nN} ({100.0*nN/tot:.2f}%)  '-' {nDash} ({100.0*nDash/tot:.2f}%)")
    print(f"杂合位点 {het_tot}  其中输出第三方碱基 {err_het}")
    print(f"**错误率 = {err}/{tot} = {100.0*err/tot:.4f}%**  (纯合位点 {err_hom} / 杂合位点 {err_het})")
    if wrong_bases:
        print("  错误碱基分布:", dict(wrong_bases.most_common(5)))
        print("  样例 (名字, 位置, 真实等位基因, 输出):", examples)
    else:
        print("  没有任何第三方碱基 → 满足 0.1% 约束")


if __name__ == '__main__':
    main()
