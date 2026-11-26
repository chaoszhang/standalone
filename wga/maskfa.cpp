// maskfa.cpp — 用 PAF 里的比对区间把 FASTA 相应位置屏蔽成 N（自比对去重复用）
// 用法: maskfa <in.fa[.gz]> <out.fa> <blocks.paf> [--min-span 300] [--min-id 0.0] [--soft]
// 说明: ref 侧和 query 侧区间都会屏蔽（自比对时两侧是同一基因组，合并成并集）；
//       --soft 输出小写而不是 N（坐标不变，但 wga2 会把小写当正常碱基，默认硬屏蔽更安全）
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

struct Rec { std::string name; std::string seq; };

static std::vector<Rec> readFasta(const char* path) {
    std::vector<Rec> out;
    const size_t plen = strlen(path);
    const bool gz = plen > 3 && strcmp(path + plen - 3, ".gz") == 0;
    FILE* f = gz ? popen((std::string("gzip -dc '") + path + "'").c_str(), "r") : fopen(path, "rb");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); exit(1); }
    std::vector<char> buf(1 << 20);
    std::string line;
    while (fgets(buf.data(), (int)buf.size(), f)) {
        line.assign(buf.data());
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            out.emplace_back();
            size_t p = 1;
            while (p < line.size() && line[p] != ' ' && line[p] != '\t') ++p;
            out.back().name = line.substr(1, p - 1);
        } else if (!out.empty()) {
            out.back().seq += line;
        }
    }
    if (gz) pclose(f); else fclose(f);
    return out;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <in.fa[.gz]> <out.fa> <blocks.paf> [--min-span N] [--min-id F] [--soft]\n", argv[0]);
        return 1;
    }
    int64_t minSpan = 300;
    double minId = 0.0;
    bool soft = false;
    int64_t dropDiag = -1;     // >=0 时：跳过同一记录上 |rs-qs|<=N 的块（自比对主对角线）
    for (int i = 4; i < argc; ++i) {
        if (!strcmp(argv[i], "--min-span") && i + 1 < argc) minSpan = atoll(argv[++i]);
        else if (!strcmp(argv[i], "--min-id") && i + 1 < argc) minId = atof(argv[++i]);
        else if (!strcmp(argv[i], "--soft")) soft = true;
        else if (!strcmp(argv[i], "--drop-diag") && i + 1 < argc) dropDiag = atoll(argv[++i]);
    }

    std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>> iv;
    FILE* pf = fopen(argv[3], "r");
    if (!pf) { fprintf(stderr, "ERROR: cannot open %s\n", argv[3]); return 1; }
    char line[8192];
    int64_t nBlocks = 0, nKept = 0;
    while (fgets(line, sizeof line, pf)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char qn[512], rn[512], strand[4];
        long long qlen, qs, qe, rlen, rs, re;
        int tok = sscanf(line, "%511s %lld %lld %lld %3s %511s %lld %lld %lld",
                         qn, &qlen, &qs, &qe, strand, rn, &rlen, &rs, &re);
        if (tok < 9) continue;
        ++nBlocks;
        double id = 1.0;
        const char* p = strstr(line, "id:f=");
        if (p) id = atof(p + 5);
        if (id < minId) continue;
        if (qe - qs < minSpan || re - rs < minSpan) continue;
        if (dropDiag >= 0 && !strcmp(qn, rn) && std::llabs((long long)rs - (long long)qs) <= dropDiag) continue;
        ++nKept;
        iv[qn].emplace_back(qs, qe);
        iv[rn].emplace_back(rs, re);
    }
    fclose(pf);

    for (auto& kv : iv) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end());
        std::vector<std::pair<int64_t, int64_t>> m;
        for (auto& p : v) {
            if (!m.empty() && p.first <= m.back().second) m.back().second = std::max(m.back().second, p.second);
            else m.push_back(p);
        }
        v.swap(m);
    }

    auto recs = readFasta(argv[1]);
    FILE* of = fopen(argv[2], "w");
    if (!of) { fprintf(stderr, "ERROR: cannot write %s\n", argv[2]); return 1; }
    int64_t totalBp = 0, maskedBp = 0;
    for (auto& r : recs) {
        auto it = iv.find(r.name);
        int64_t len = (int64_t)r.seq.size();
        totalBp += len;
        std::vector<char> mask((size_t)len, 0);
        if (it != iv.end()) {
            for (auto& p : it->second) {
                int64_t a = std::max<int64_t>(0, p.first), b = std::min<int64_t>(len, p.second);
                for (int64_t i = a; i < b; ++i) mask[(size_t)i] = 1;
            }
        }
        for (int64_t i = 0; i < len; ++i) {
            if (mask[(size_t)i]) {
                ++maskedBp;
                r.seq[(size_t)i] = soft ? (char)tolower((unsigned char)r.seq[(size_t)i]) : 'N';
            }
        }
        fprintf(of, ">%s\n", r.name.c_str());
        for (int64_t i = 0; i < len; i += 60)
            fprintf(of, "%.*s\n", (int)std::min<int64_t>(60, len - i), r.seq.c_str() + i);
        std::string().swap(r.seq);
    }
    fclose(of);
    fprintf(stderr, "[mask] blocks=%lld kept=%lld | 基因组 %lld bp | 屏蔽 %lld bp (%.2f%%)\n",
            (long long)nBlocks, (long long)nKept, (long long)totalBp, (long long)maskedBp,
            totalBp ? 100.0 * (double)maskedBp / (double)totalBp : 0.0);
    return 0;
}
