#define main wga3_main
#include "wga3.cpp"
#undef main

int main() {
    auto recs = readFasta("t1.ref.fa");
    if (recs.empty()) { printf("no ref\n"); return 1; }
    std::string padL, padR;
    { std::mt19937 rng(1); for (int i = 0; i < PAD; ++i) padL.push_back("ACGT"[rng() & 3]); }
    { std::mt19937 rng(2); for (int i = 0; i < PAD; ++i) padR.push_back("ACGT"[rng() & 3]); }
    Genome g;
    buildGenome(recs, padL, padR, 1000, g);
    // 故意插入若干 N，验证 N 的往返
    std::vector<int8_t> orig = g.enc;
    for (size_t i = 1000; i < orig.size(); i += 9973) orig[i] = -1;
    for (size_t i = 0; i < orig.size(); ++i) g.enc[i] = orig[i];
    const size_t n = orig.size();
    packRef(g);
    printf("enc=%zu B -> pk=%zu B (压缩比 %.2fx)\n", n, g.pk.size(), (double)n / (double)g.pk.size());
    std::vector<int8_t> back(n);
    decodeRef(g, 0, (int)n, back.data());
    size_t bad = 0, first = SIZE_MAX;
    for (size_t i = 0; i < n; ++i)
        if (orig[i] != back[i]) { if (first == SIZE_MAX) first = i; ++bad; }
    printf("往返差异 = %zu 字节", bad);
    if (first != SIZE_MAX) printf("（首个 @%zu: 原 %d -> 解 %d）", first, (int)orig[first], (int)back[first]);
    printf("\n");
    // 抽查中间任意起点（奇偶偏移都要覆盖）
    size_t bad2 = 0;
    for (int64_t from = 0; from + 150 <= (int64_t)n; from += 4001) {
        std::vector<int8_t> w(150);
        decodeRef(g, from, 150, w.data());
        for (int k = 0; k < 150; ++k) if (w[(size_t)k] != orig[(size_t)(from + k)]) ++bad2;
    }
    printf("窗口抽查差异 = %zu\n -> %s\n", bad2, (bad == 0 && bad2 == 0) ? "PASS" : "FAIL");
    return (bad == 0 && bad2 == 0) ? 0 : 1;
}
