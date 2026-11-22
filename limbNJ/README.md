# limbNJ — standalone distance-based tree builders

Three small, dependency-free C++20 programs that reconstruct an unrooted
phylogenetic tree from a distance matrix. They share one implementation
(`limbNJ.hpp`) and differ only in how two quantities are computed.

| program | `Lf` estimate | pair selection | main loop |
|---|---|---|---|
| `NJ` | — | brute-force scan of `Q` | O(n³) |
| `limbNJ1.3` | pq-side average | binary heap + lazy deletion | O(n² log n) |
| `limbNJ1.3.1` | pq-side average | per-row min / 2nd-min | O(n²) amortized |

`limbNJ1.3` and `limbNJ1.3.1` are variants of the same additive-style
agglomerative scheme. They maintain, for every active cluster `i`, a limb
length `L_i` and a weight `c_i`; each round merges the pair minimising

```
S(i)   = Σ_k c_k·D(i,k) − D(i,i)
Q(i,j) = (n−2)·D(i,j) − S(i) − S(j)
```

with `D(i,i) = 2·L_i`. After merging `(i,j)` into `i`:

```
ā       = (c_i·L_i + c_j·L_j) / (c_i + c_j)
D(i,p) ← (c_i·D(i,p) + c_j·D(j,p)) / (c_i + c_j)
L_i    ← d(i,w) − L_i + ā          (w = the parent of the cherry)
```

`d(i,w)` is obtained with a three-point formula, `c_i += c_j`, and the cluster
is flattened into its parent (no leaf is ever deleted).

### The `Lf` estimate

`L_i` is initialised from the input matrix. The textbook LimbLength theorem
takes `min` over all pairs `(i,k)`; that is a downward-biased estimator under
noise. `limbNJ1.3` / `limbNJ1.3.1` instead average over the pairs that lie on
**opposite sides** of leaf `j`'s attachment point, which are identified without
knowing the tree:

1. pick `p`; let `q = argmin_k limb(j,p,k)`
2. let `p' = argmin_k limb(j,p',q)`
3. for every leaf `y`, compare `limb(j,y,p')` against `limb(j,y,q)` — the
   smaller one tells which side of the attachment point `y` lies on
4. average the three-point formula over all `P × Q` pairs

All terms in that average equal `Lf(j)` when the matrix is additive, so the
estimator is exact there and less noisy elsewhere.

## Build

```bash
make            # produces NJ, limbNJ1.3, limbNJ1.3.1
```

Only a C++20 compiler is needed — no external libraries.

## Input / output

Input is a plain text matrix:

```
n
d00 d01 ... d0(n-1)
...
d(n-1)0 ... d(n-1)(n-1)
```

Usage:

```bash
./limbNJ1.3.1 --matrix matrix.txt            # splits to stdout
./limbNJ1.3.1 --matrix matrix.txt --out s.txt
```

Output is `n−3` lines, one split per line, as `k v1 v2 ... vk` — the members
of one side of the split, taking the lexicographically smaller side as
canonical. Timings and the split count go to stderr.

`--lf min` switches the `Lf` estimator of `limbNJ1.3` / `limbNJ1.3.1` back to
the textbook `min` (useful for comparison).

## Performance

Exact additive matrices, GCC 13.3 `-O2`, single thread, no I/O included in the
main loop column:

| n | `limbNJ1.3` | `limbNJ1.3.1` | `NJ` |
|---|---|---|---|
| 400 | 0.05 s | 0.03 s | 0.08 s |
| 800 | 0.25 | 0.14 | 0.53 |
| 1600 | 1.37 | **0.73** | 4.17 |

At n = 1600, `limbNJ1.3.1` is about 5.7× faster than `NJ` and `limbNJ1.3`
about 3.0×.

Notes:

- The `Lf` initialisation is O(n³) (O(n²) per leaf). At n = 1600 it costs about
  0.25 s, so it is not yet the bottleneck, but it will overtake the main loop
  at larger sizes.
- The heap variant is slower than the per-row variant in practice, and uses
  about 3× more memory (stale heap entries are not reclaimed).
  **Prefer `limbNJ1.3.1`.**

## Reference

LimbLength theorem and the additive phylogeny construction:
Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*,
Chapter 7.
