// mkreads.cpp — 定长读生成器（读名带真值：t<start>_<f|r>），用于校验读模式坐标
// 用法: mkreads ref.fa out.fa readLen sub indel seed [strands=0]
//   strands=1 时奇偶交替输出正/反链（反链读名后缀 r，真值起点仍是正链坐标）
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

static const char* B = "ACGT";
static inline char comp(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
}

int main(int argc, char** argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: %s ref.fa out.fa readLen sub indel seed [strands=0]\n", argv[0]);
        return 1;
    }
    const int L = atoi(argv[3]);
    const double sub = atof(argv[4]), indel = atof(argv[5]);
    std::mt19937_64 rng((unsigned long long)atoll(argv[6]));
    const int strands = argc > 7 ? atoi(argv[7]) : 0;
    std::uniform_real_distribution<double> U(0, 1);

    FILE* fi = fopen(argv[1], "rb");
    if (!fi) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    FILE* fo = fopen(argv[2], "w");
    if (!fo) { fprintf(stderr, "cannot write %s\n", argv[2]); return 1; }

    std::string ref;
    {
        std::vector<char> buf(1 << 20);
        bool first = true;
        while (fgets(buf.data(), (int)buf.size(), fi)) {
            if (buf[0] == '>') { if (!first) break; first = false; continue; }   // 只取第一条记录
            for (char* p = buf.data(); *p; ++p)
                if (*p != '\n' && *p != '\r') ref.push_back(*p);
        }
    }
    fclose(fi);
    const long long n = (long long)ref.size();
    long long nRead = 0;
    std::string s;
    s.reserve((size_t)L * 2);
    for (long long st = 0; st + L <= n; st += L) {
        s.clear();
        for (int i = 0; i < L; ++i) {
            char c = ref[(size_t)(st + i)];
            if (U(rng) < indel) continue;                       // 删除
            char o = c;
            if (U(rng) < sub) { do { o = B[rng() & 3]; } while (o == c); }
            s.push_back(o);
            if (U(rng) < indel) s.push_back(B[rng() & 3]);      // 插入
        }
        const bool rc = strands && ((nRead & 1) == 1);
        if (s.size() < (size_t)L) s.append((size_t)L - s.size(), 'N');
        if (s.size() > (size_t)L) s.resize((size_t)L);
        fprintf(fo, ">t%lld_%c\n", st, rc ? 'r' : 'f');
        if (rc) {
            std::string t((size_t)L, 'N');
            for (int i = 0; i < L; ++i) t[(size_t)i] = comp(s[(size_t)(L - 1 - i)]);
            fprintf(fo, "%s\n", t.c_str());
        } else {
            fprintf(fo, "%s\n", s.c_str());
        }
        ++nRead;
    }
    fclose(fo);
    fprintf(stderr, "mkreads: reads=%lld len=%d sub=%.3f indel=%.3f strands=%d\n", nRead, L, sub, indel, strands);
    return 0;
}
