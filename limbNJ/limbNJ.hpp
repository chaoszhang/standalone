// limbNJ.hpp —— 三种距离法建树的共用实现
//
//   NJ        经典邻居连接（朴素扫描）                      O(n^3)
//   limbNJ1.3  Lf 用 pq 阵营平均，选对用堆                  O(n^2 log n) + O(n^3) Lf
//   limbNJ1.3.1 Lf 用 pq 阵营平均，选对用每行 min/次小       O(n^2)       + O(n^3) Lf
//
// 输入：第一行 n，随后 n 行每行 n 个数（距离矩阵）
// 输出：n-3 行 split，每行 "k v1 v2 ... vk"（canonical：取字典序小的一侧）
// 计时与规模信息到 stderr。

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

namespace limbNJ {

enum class Mode { A13, A131, NJ };

static const double kInf = 1e300;

struct Result {
    std::vector<std::vector<int>> clades;
    std::vector<std::vector<int>> tree;   // 内部节点 -> 孩子节点（叶子 0..n-1）
    int root = -1;
    double t_lf = 0, t_main = 0;
    bool ok = true;
    std::string why;
};

static inline double at(const std::vector<double>& D, int n, int i, int j) {
    return D[(size_t)i * n + j];
}
static inline double& ref(std::vector<double>& D, int n, int i, int j) {
    return D[(size_t)i * n + j];
}

// 三点式 limb(j; i,k) = (d(j,i)+d(j,k)-d(i,k))/2
static inline double tri(const std::vector<double>& D, int n, int j, int i, int k) {
    return (at(D, n, j, i) + at(D, n, j, k) - at(D, n, i, k)) * 0.5;
}

// ── Lf 估计：pq 阵营平均 ──────────────────────────────────────
// 1. 固定 p，选 q 使 limb(j,p,q) 最小
// 2. 固定 q，选 p' 使 limb(j,p',q) 最小
// 3. 对每片叶 y，比较 limb(j,y,p') 与 limb(j,y,q)，谁小跟谁同阵营
// 4. P x Q 配对，三点式取平均
static double lf_pq(const std::vector<double>& D, int n, int j, int p) {
    int q = -1, p2 = -1;
    double bq = kInf, bp = kInf;
    for (int k = 0; k < n; ++k) {
        if (k == j || k == p) continue;
        double v = tri(D, n, j, p, k);
        if (v < bq) { bq = v; q = k; }
    }
    if (q < 0) return 0.0;
    for (int k = 0; k < n; ++k) {
        if (k == j || k == q) continue;
        double v = tri(D, n, j, k, q);
        if (v < bp) { bp = v; p2 = k; }
    }
    if (p2 < 0) return bq;

    std::vector<int> P, Q;
    P.reserve(n);
    Q.reserve(n);
    for (int y = 0; y < n; ++y) {
        if (y == j) continue;
        if (y == p2) { P.push_back(y); continue; }
        if (y == q) { Q.push_back(y); continue; }
        double a = tri(D, n, j, y, q);
        double b = tri(D, n, j, y, p2);
        if (a < b) P.push_back(y); else Q.push_back(y);
    }
    if (P.empty() || Q.empty()) return bq;

    double s = 0.0;
    long long m = 0;
    for (int a : P) {
        for (int b : Q) { s += tri(D, n, j, a, b); ++m; }
    }
    return s / (double)m;
}

// 朴素 min 版（--lf min，用作对照）
static double lf_min(const std::vector<double>& D, int n, int j) {
    double best = kInf;
    for (int i = 0; i < n; ++i) {
        if (i == j) continue;
        for (int k = i + 1; k < n; ++k) {
            if (k == j) continue;
            best = std::min(best, tri(D, n, j, i, k));
        }
    }
    return best;
}

static std::vector<double> compute_lf(const std::vector<double>& D, int n, bool use_pq) {
    std::vector<double> Lf(n, 0.0);
    for (int j = 0; j < n; ++j) {
        double v;
        if (use_pq) {
            int p = (j == 0 ? 1 : 0);
            v = lf_pq(D, n, j, p);
        } else {
            v = lf_min(D, n, j);
        }
        Lf[j] = std::max(0.0, v);
    }
    return Lf;
}

// ── limbNJ1.3 / limbNJ1.3.1 主循环 ────────────────────────────
// use_heap = true  -> 堆 + stamp 惰性删除        O(n^2 log n)
// use_heap = false -> 每行 min/argmin/2nd        O(n^2) 摊还
static Result solve_a1x(int n, const std::vector<double>& D0,
                        const std::vector<double>& Lf, bool use_heap) {
    Result R;
    std::vector<double> D(D0);
    for (int i = 0; i < n; ++i) ref(D, n, i, i) = 2.0 * Lf[i];

    std::vector<char> alive(n, 1);
    std::vector<double> c(n, 1.0), S(n, 0.0);
    std::vector<std::vector<int>> grp(n);
    std::vector<int> cur(n);                 // clade -> 当前树节点
    for (int i = 0; i < n; ++i) { grp[i].assign(1, i); cur[i] = i; }

    int m = n;
    const int K = n - 2;

    auto recompute_S = [&](int a) {
        double s = 0.0;
        for (int k = 0; k < n; ++k)
            if (alive[k]) s += c[k] * at(D, n, a, k);
        return s - at(D, n, a, a);
    };
    for (int i = 0; i < n; ++i) S[i] = recompute_S(i);

    auto Q = [&](int a, int b) {
        return (double)K * at(D, n, a, b) - S[a] - S[b];
    };

    auto do_merge = [&](int i, int j) {
        double Li = at(D, n, i, i) * 0.5, Lj = at(D, n, j, j) * 0.5;
        double ci = c[i], cj = c[j], cn = ci + cj;
        double abar = (ci * Li + cj * Lj) / cn;

        // d(i,w)：固定 p，只扫 q
        double diw = kInf;
        int p0 = -1;
        for (int k = 0; k < n; ++k)
            if (alive[k] && k != i && k != j) { p0 = k; break; }
        if (p0 >= 0) {
            for (int k = 0; k < n; ++k) {
                if (!alive[k] || k == i || k == j || k == p0) continue;
                diw = std::min(diw,
                    (at(D, n, i, p0) + at(D, n, i, k) - at(D, n, p0, k)) * 0.5);
            }
        }
        double Lnew = (diw > kInf * 0.5) ? abar : diw - Li + abar;

        for (int k = 0; k < n; ++k) {
            if (!alive[k] || k == i || k == j) continue;
            double v = (ci * at(D, n, i, k) + cj * at(D, n, j, k)) / cn;
            ref(D, n, i, k) = ref(D, n, k, i) = v;
        }
        ref(D, n, i, i) = 2.0 * Lnew;
        c[i] = cn;
        for (int x : grp[j]) grp[i].push_back(x);
        // 拓扑：i、j 是 cherry，新父节点 u
        {
            int u = n + (int)R.tree.size();   // 内部节点用 n+idx，避免与叶 id 冲突
            R.tree.push_back({cur[i], cur[j]});
            cur[i] = u;
        }
        alive[j] = 0;
        --m;
        if (m > 1) S[i] = recompute_S(i);   // 其余 S 在 claim4 下严格不变
        R.clades.push_back(grp[i]);
    };

    if (use_heap) {
        using Item = std::tuple<double, int, int, int, int>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        std::vector<int> stamp(n, 0);
        auto push_row = [&](int a) {
            for (int k = 0; k < n; ++k)
                if (alive[k] && k != a)
                    pq.emplace(Q(a, k), a, k, stamp[a], stamp[k]);
        };
        for (int i = 0; i < n; ++i) push_row(i);

        while (m > 2) {
            int i = -1, j = -1;
            while (!pq.empty()) {
                auto [qv, a, b, sa, sb] = pq.top();
                pq.pop();
                if (!alive[a] || !alive[b]) continue;
                if (sa != stamp[a] || sb != stamp[b]) continue;
                i = a; j = b;
                break;
            }
            if (i < 0) { R.ok = false; R.why = "heap-empty"; break; }
            if (i > j) std::swap(i, j);
            do_merge(i, j);
            if (m <= 1) break;
            stamp[i]++;
            push_row(i);
        }
    } else {
        struct Row { double v1 = kInf, v2 = kInf; int a1 = -1, a2 = -1; };
        std::vector<Row> row(n);
        auto rescan = [&](int a, Row& r) {
            r = Row{};
            for (int k = 0; k < n; ++k) {
                if (!alive[k] || k == a) continue;
                double v = Q(a, k);
                if (v < r.v1) { r.v2 = r.v1; r.a2 = r.a1; r.v1 = v; r.a1 = k; }
                else if (v < r.v2) { r.v2 = v; r.a2 = k; }
            }
        };
        for (int a = 0; a < n; ++a) rescan(a, row[a]);

        while (m > 2) {
            int i = -1, j = -1;
            double best = kInf;
            for (int a = 0; a < n; ++a) {
                if (!alive[a] || row[a].a1 < 0) continue;
                int b = row[a].a1;
                double v = row[a].v1;
                int x = std::min(a, b), y = std::max(a, b);
                if (v < best - 1e-12 ||
                    (std::fabs(v - best) <= 1e-12 && (x < i || (x == i && y < j)))) {
                    best = v; i = x; j = y;
                }
            }
            if (i < 0) { R.ok = false; R.why = "row-empty"; break; }
            do_merge(i, j);
            if (m <= 1) break;

            rescan(i, row[i]);
            for (int p = 0; p < n; ++p) {
                if (!alive[p] || p == i) continue;
                bool bad = (row[p].a1 == i || row[p].a1 == j ||
                            row[p].a2 == i || row[p].a2 == j);
                if (bad || row[p].a1 < 0) { rescan(p, row[p]); continue; }
                double v = Q(p, i);
                if (v < row[p].v1 - 1e-12) {
                    row[p].v2 = row[p].v1; row[p].a2 = row[p].a1;
                    row[p].v1 = v; row[p].a1 = i;
                } else if (v < row[p].v2 - 1e-12) {
                    row[p].v2 = v; row[p].a2 = i;
                }
            }
        }
    }

    // 收尾：剩下两个 clade 连成一个根
    {
        int a = -1, b = -1;
        for (int k = 0; k < n; ++k)
            if (alive[k]) { if (a < 0) a = k; else if (b < 0) b = k; }
        if (a >= 0 && b >= 0) {
            R.root = n + (int)R.tree.size();
            R.tree.push_back({cur[a], cur[b]});
        } else if (a >= 0) {
            R.root = cur[a];
        }
    }
    return R;
}

// ── 经典 NJ（朴素扫描）──────────────────────────────────────
static Result solve_nj(int n, const std::vector<double>& D0) {
    Result R;
    std::vector<double> d(D0);
    for (int i = 0; i < n; ++i) ref(d, n, i, i) = 0.0;

    std::vector<char> alive(n, 1);
    std::vector<std::vector<int>> grp(n);
    std::vector<int> cur(n);
    for (int i = 0; i < n; ++i) { grp[i].assign(1, i); cur[i] = i; }

    std::vector<double> r(n, 0.0);
    auto row_sum = [&](int a) {
        double s = 0.0;
        for (int k = 0; k < n; ++k)
            if (alive[k] && k != a) s += at(d, n, a, k);
        return s;
    };
    for (int i = 0; i < n; ++i) r[i] = row_sum(i);

    int m = n;
    while (m > 2) {
        double best = kInf;
        int bi = -1, bj = -1;
        for (int a = 0; a < n; ++a) {
            if (!alive[a]) continue;
            for (int b = a + 1; b < n; ++b) {
                if (!alive[b]) continue;
                double q = (double)(m - 2) * at(d, n, a, b) - r[a] - r[b];
                if (q < best) { best = q; bi = a; bj = b; }
            }
        }
        if (bi < 0) { R.ok = false; R.why = "nj-stuck"; break; }

        int i = bi, j = bj;
        double dij = at(d, n, i, j);
        for (int k = 0; k < n; ++k) {
            if (!alive[k] || k == i || k == j) continue;
            double v = 0.5 * (at(d, n, i, k) + at(d, n, j, k) - dij);
            ref(d, n, i, k) = ref(d, n, k, i) = v;
        }
        ref(d, n, i, i) = 0.0;
        for (int x : grp[j]) grp[i].push_back(x);
        {
            int u = n + (int)R.tree.size();
            R.tree.push_back({cur[i], cur[j]});
            cur[i] = u;
        }
        alive[j] = 0;
        --m;
        R.clades.push_back(grp[i]);
        if (m <= 2) break;
        for (int k = 0; k < n; ++k)
            if (alive[k]) r[k] = row_sum(k);
    }
    {
        int a = -1, b = -1;
        for (int k = 0; k < n; ++k)
            if (alive[k]) { if (a < 0) a = k; else if (b < 0) b = k; }
        if (a >= 0 && b >= 0) {
            R.root = n + (int)R.tree.size();
            R.tree.push_back({cur[a], cur[b]});
        } else if (a >= 0) R.root = cur[a];
    }
    return R;
}

// ── 输出 ───────────────────────────────────────────────────
static std::vector<std::string> to_splits(const Result& R, int n) {
    std::vector<std::string> out;
    for (const auto& g : R.clades) {
        if ((int)g.size() < 2 || (int)g.size() > n - 2) continue;
        std::vector<char> in(n, 0);
        for (int x : g) in[x] = 1;
        std::vector<int> a, b;
        for (int x = 0; x < n; ++x) (in[x] ? a : b).push_back(x);
        const std::vector<int>& s = (a < b) ? a : b;   // 字典序小的一侧
        std::string line = std::to_string(s.size());
        for (int x : s) { line += ' '; line += std::to_string(x); }
        out.push_back(line);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// PHYLIP 风格：第一行 n，随后每行 "name<TAB>d0<TAB>d1..."
static bool read_matrix(const std::string& path, int& n, std::vector<double>& D,
                        std::vector<std::string>& names) {
    std::ifstream in(path);
    if (!in) return false;
    if (!(in >> n) || n < 3) return false;
    std::string line;
    std::getline(in, line);              // 吃掉行尾
    names.assign(n, std::string());
    D.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (!(in >> names[i])) return false;
        for (int j = 0; j < n; ++j)
            if (!(in >> D[(size_t)i * n + j])) return false;
    }
    return true;
}

// 无根 Newick：(A,B);
static std::string to_newick(const Result& R, const std::vector<std::string>& names) {
    int n = (int)names.size();
    std::function<std::string(int)> rec = [&](int u) -> std::string {
        if (u < n) return names[u];
        const std::vector<int>& kids = R.tree[u - n];
        std::string s = "(";
        for (size_t k = 0; k < kids.size(); ++k) {
            if (k) s += ',';
            s += rec(kids[k]);
        }
        s += ")";
        return s;
    };
    if (R.root < 0 || R.tree.empty()) return std::string();
    return rec(R.root) + ";";
}

// ── 统一入口 ───────────────────────────────────────────────
static int run(int argc, char** argv, Mode mode, const char* prog) {
    std::string mat, out_path, lfmode = "pq", fmt = "newick";
    for (int a = 1; a < argc; ++a) {
        std::string s = argv[a];
        if ((s == "--matrix" || s == "-i") && a + 1 < argc) mat = argv[++a];
        else if ((s == "--out" || s == "-o") && a + 1 < argc) out_path = argv[++a];
        else if (s == "--lf" && a + 1 < argc) lfmode = argv[++a];
        else if (s == "--format" && a + 1 < argc) fmt = argv[++a];
        else if (s == "-h" || s == "--help") {
            std::printf(
                "用法: %s --matrix <文件> [--out <文件>] [--lf pq|min]\n"
                "  输入：第一行 n，随后 n 行每行 n 个数（距离矩阵）\n"
                "  输出：n-3 行 split，每行 \"k v1 v2 ... vk\"\n", prog);
            return 0;
        }
    }
    if (mat.empty()) { std::fprintf(stderr, "缺少 --matrix\n"); return 2; }

    int n = 0;
    std::vector<double> D;
    std::vector<std::string> names;
    if (!read_matrix(mat, n, D, names)) {
        std::fprintf(stderr, "读取矩阵失败: %s\n", mat.c_str());
        return 2;
    }

    using clk = std::chrono::steady_clock;
    Result R;
    if (mode == Mode::NJ) {
        auto t0 = clk::now();
        R = solve_nj(n, D);
        auto t1 = clk::now();
        R.t_main = std::chrono::duration<double>(t1 - t0).count();
    } else {
        auto t0 = clk::now();
        std::vector<double> Lf = compute_lf(D, n, lfmode == "pq");
        auto t1 = clk::now();
        R = solve_a1x(n, D, Lf, mode == Mode::A13);
        auto t2 = clk::now();
        R.t_lf = std::chrono::duration<double>(t1 - t0).count();
        R.t_main = std::chrono::duration<double>(t2 - t1).count();
    }

    std::vector<std::string> sp = to_splits(R, n);
    std::string outtext;
    if (fmt == "splits") {
        for (auto& t : sp) { outtext += t; outtext += '\n'; }
    } else {
        outtext = to_newick(R, names);
        if (!outtext.empty()) outtext += '\n';
    }
    std::fprintf(stderr, "%s n=%d  t_lf=%.3fs  t_main=%.3fs  splits=%zu%s%s\n",
                 prog, n, R.t_lf, R.t_main, sp.size(),
                 R.ok ? "" : "  [失败] ", R.why.c_str());

    if (out_path.empty()) {
        std::fputs(outtext.c_str(), stdout);
    } else {
        std::ofstream f(out_path);
        f << outtext;
    }
    return R.ok ? 0 : 1;
}

}  // namespace limbNJ
