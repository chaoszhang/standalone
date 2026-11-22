// gtree_dist.cpp —— 从基因树文件算「平均 patristic 距离」矩阵
//
// 输入：每行一棵 Newick 的 .gtrees 文件
// 输出：PHYLIP 风格矩阵
//         n
//         T0<TAB>d00<TAB>d01...
//
// 编译：g++ -O2 -std=c++20 -o gtree_dist gtree_dist.cpp
// 用法：./gtree_dist -i in.gtrees -o out.mat [-n 1000]

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct Tree {
    vector<vector<pair<int, double>>> adj;
    vector<int> leafNode;   // 叶 id -> 节点 id
};

static unordered_map<string, int> g_name2id;
static vector<string> g_id2name;
static bool g_strip = false;      // true：把 T20_0 / T20_1 归并成 T20

static int getId(const string& nm) {
    string key = nm;
    if (g_strip) {
        size_t p = nm.find('_');
        if (p != string::npos) key = nm.substr(0, p);
    }
    auto it = g_name2id.find(key);
    if (it != g_name2id.end()) return it->second;
    int id = (int)g_id2name.size();
    g_name2id.emplace(key, id);
    g_id2name.push_back(key);
    return id;
}

static double readLen(const string& s, size_t& pos) {
    ++pos;   // 跳过 ':'
    size_t st = pos;
    while (pos < s.size()) {
        char c = s[pos];
        if (isdigit((unsigned char)c) || c == '.' || c == '-' || c == '+' ||
            c == 'e' || c == 'E')
            ++pos;
        else
            break;
    }
    if (pos == st) return 0.0;
    try {
        return stod(s.substr(st, pos - st));
    } catch (...) {
        return 0.0;
    }
}

// 返回子树根节点 id；枝长挂在父边
static int parseSubtree(const string& s, size_t& pos, Tree& T) {
    int u = (int)T.adj.size();
    T.adj.push_back({});

    while (pos < s.size() && isspace((unsigned char)s[pos])) ++pos;

    if (pos < s.size() && s[pos] == '(') {
        ++pos;
        while (true) {
            int v = parseSubtree(s, pos, T);
            double len = 0.0;
            if (pos < s.size() && s[pos] == ':') len = readLen(s, pos);
            T.adj[u].push_back({v, len});
            T.adj[v].push_back({u, len});
            while (pos < s.size() && isspace((unsigned char)s[pos])) ++pos;
            if (pos < s.size() && s[pos] == ',') { ++pos; continue; }
            break;
        }
        if (pos < s.size() && s[pos] == ')') ++pos;
        // 内部节点标签（bootstrap 等）跳过
        while (pos < s.size() && s[pos] != ':' && s[pos] != ',' &&
               s[pos] != ')' && s[pos] != ';')
            ++pos;
    } else {
        size_t st = pos;
        while (pos < s.size() && s[pos] != ':' && s[pos] != ',' &&
               s[pos] != ')' && s[pos] != ';')
            ++pos;
        string nm = s.substr(st, pos - st);
        if (nm.size() >= 2 && nm.front() == '\'' && nm.back() == '\'')
            nm = nm.substr(1, nm.size() - 2);
        int id = getId(nm);
        if ((int)T.leafNode.size() <= id) T.leafNode.resize(id + 1, -1);
        T.leafNode[id] = u;
    }
    return u;
}

int main(int argc, char** argv) {
    string in, out;
    int maxTrees = 0;
    for (int a = 1; a < argc; ++a) {
        string s = argv[a];
        if (s == "-i" && a + 1 < argc) in = argv[++a];
        else if (s == "-o" && a + 1 < argc) out = argv[++a];
        else if (s == "-n" && a + 1 < argc) maxTrees = atoi(argv[++a]);
        else if (s == "-h" || s == "--help") {
            printf("用法: gtree_dist -i in.gtrees -o out.mat [-n maxTrees]\n");
            return 0;
        }
    }
    if (in.empty()) { fprintf(stderr, "缺 -i\n"); return 2; }

    ifstream fin(in);
    if (!fin) { fprintf(stderr, "打不开 %s\n", in.c_str()); return 2; }

    vector<vector<double>> acc;
    long long ntree = 0;
    string line;
    vector<int> seen, stack_;
    vector<double> dist;

    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (maxTrees > 0 && ntree >= maxTrees) break;
        if (line.find('(') == string::npos) continue;

        Tree T;
        T.adj.reserve(1024);
        size_t pos = 0;
        parseSubtree(line, pos, T);

        int n = (int)g_id2name.size();
        if ((int)acc.size() < n) acc.assign(n, vector<double>(n, 0.0));

        int V = (int)T.adj.size();
        if ((int)seen.size() < V) { seen.assign(V, 0); dist.assign(V, 0.0); }
        static int stamp = 0;

        for (int lf = 0; lf < n; ++lf) {
            if (lf >= (int)T.leafNode.size() || T.leafNode[lf] < 0) continue;
            int src = T.leafNode[lf];
            ++stamp;
            stack_.clear();
            stack_.push_back(src);
            seen[src] = stamp;
            dist[src] = 0.0;
            while (!stack_.empty()) {
                int u = stack_.back();
                stack_.pop_back();
                for (auto& e : T.adj[u]) {
                    int v = e.first;
                    if (seen[v] == stamp) continue;
                    seen[v] = stamp;
                    dist[v] = dist[u] + e.second;
                    stack_.push_back(v);
                }
            }
            for (int t = lf + 1; t < n; ++t) {
                if (t >= (int)T.leafNode.size() || T.leafNode[t] < 0) continue;
                acc[lf][t] += dist[T.leafNode[t]];
            }
        }
        ++ntree;
    }

    if (ntree == 0) { fprintf(stderr, "没读到树\n"); return 1; }

    int n = (int)g_id2name.size();
    vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    sort(order.begin(), order.end(),
         [](int a, int b) { return g_id2name[a] < g_id2name[b]; });

    FILE* fo = out.empty() ? stdout : fopen(out.c_str(), "w");
    if (!fo) { fprintf(stderr, "写不了 %s\n", out.c_str()); return 2; }
    fprintf(fo, "%d\n", n);
    for (int i : order) {
        fprintf(fo, "%s", g_id2name[i].c_str());
        for (int j : order) {
            double v;
            if (i == j) v = 0.0;
            else {
                int a = min(i, j), b = max(i, j);
                v = acc[a][b] / (double)ntree;
            }
            fprintf(fo, "\t%.10g", v);
        }
        fprintf(fo, "\n");
    }
    if (fo != stdout) fclose(fo);
    fprintf(stderr, "读了 %lld 棵树，%d 个 taxa\n", ntree, n);
    return 0;
}
