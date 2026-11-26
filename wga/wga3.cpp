// wga3.cpp — pairwise whole-genome mapping（14-mer seed + 192-bit 远近混合指纹）
//
// 规格：
//   建库：扫 Ref 正链所有 14-mer；code % 7 == 0 者入 multi-map M，value = (位置, 指纹)
//   指纹 192 bit：
//     近指纹 64 bit  = 左侧 16-mer 的 2-bit 编码(32b) | 右侧 16-mer 的 2-bit 编码(32b)
//     远指纹 128 bit = 左 64 bit（左近区再往左 32 bp 内的 3-mer 存在位图）
//                    | 右 64 bit（右近区再往右 32 bp 内的 3-mer 存在位图）
//   打分：近指纹按位点精确匹配数 + 远指纹 popcount(AND)，各自按零假设标准化后加权相加
//         score = 100·(wN·z_near + wF·(z_farL + z_farR))
//
//   比对：
//     0. 字典 D：key = Ref 位置，value = 该位置目前最好的指纹拟合度（只存分数）
//        随机生成 48bp 的 L / R，用于给 Ref/Query 两端补齐上下文（两端指纹不再退化）
//     1. 第一遍：扫 Query 正反链，对完全落在 Query 内部的合格 14-mer 查 M，
//        拟合度高于 D[refPos] 则更新（只记分数，不记 Query 位置 → 省内存）
//     2. 第二遍：重扫，拟合度 == D[refPos] 的配对做锚定延伸（左/右 local alignment）
//     3. 整合：按对角线聚类成块 → 块两端 SW 延伸 → 输出 PAF 风格的 mapping
//
// 读模式（--reads L，短读专用，实测配置）：
//   · 查询按 padL+read+padR 紧凑等距排布（无 N 间隔），记录号 O(1) 反算
//   · 跳过 D（"每 ref 位置只留最佳"在多覆盖下会把 5 条读压成 1 条），
//     改为保留所有过阈候选 → 每条读取 fit 最高的 topK 个候选逐个做全读比对，
//     取 identity 最高的一个。指纹只负责缩小候选，判真靠序列比对。
//   · 阈值可用 --fdr 按实测零假设分位数自动定；[calib] 会打印桶均/每读候选数
//
// 实测（sim 生成，20% 替换 + 2% indel，150 bp 读）：
//   100 Mbp × 1X : mapped 27.56%，位置正确率 100.00%（±20bp），identity 0.802
//   4 Gbp  × 5X  : 133.3M 读，覆盖 80.82%，identity 0.805，8m19s / 130 GB
//                   （-f 450 --topk 8；-f 500 时覆盖 79.85%）
//   长序列模式（50 kb contig，30% 分化）：鸡 1.05 Gb × 鸵鸟 1.23 Gb，ref 覆盖 96.54%
//
// 编译: g++ -O2 -std=c++20 -pthread -o wga3 wga3.cpp
// 自检: ./wga3 --selftest

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <sys/mman.h>
#include <thread>
#include <type_traits>
#include <vector>

// ---------------------------------------------------------------- constants

static constexpr int      SEED_K = 14;
static constexpr uint32_t SEED_MASK = (1u << (2 * SEED_K)) - 1;      // 28 bit
static constexpr uint32_t SEED_MOD = 7;
static constexpr int      NEAR_LEN = 16;                             // 近指纹：左右各 16 bp
static constexpr int      FAR_LEN = 32;                              // 远指纹：再往外 32 bp
static constexpr int      K3 = 3;
static constexpr int      K3N = 1 << (2 * K3);                       // 64 种 3-mer
static constexpr int      PAD = 48;                                  // L/R 补齐长度
static constexpr uint32_t RANGE_N = ((1u << 28) / SEED_MOD) + 1;     // key/7 上界 ≈ 3.84e7

// 192-bit 指纹
struct Fp {
    uint64_t near64;   // 低 32 bit = 左 16-mer 编码，高 32 bit = 右 16-mer 编码
    uint64_t farL;     // 左远：32bp 内 3-mer 存在位图
    uint64_t farR;     // 右远
};

// 运行时标定的零假设参数
static double g_muN = 8.0, g_sdN = 2.6;
static double g_muFL = 15.0, g_sdFL = 3.0;
static double g_muFR = 15.0, g_sdFR = 3.0;
static double g_wN = 1.0, g_wF = 1.0;
// 读模式：短读紧凑打包（padL + read + padR 连续存放，无 N 间隔）
static bool g_readMode = false;
static int64_t g_readStride = 0;   // len + 2*PAD
static int64_t g_readLen = 0;

static inline int rawBaseCode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}
static int8_t g_codeTab[256];
struct CodeTabInit { CodeTabInit() { for (int i = 0; i < 256; ++i) g_codeTab[i] = (int8_t)rawBaseCode((char)i); } };
static CodeTabInit g_codeTabInit;
static inline int baseCode(char c) { return g_codeTab[(unsigned char)c]; }

// 64 bit 里 32 个 2-bit 位点的一致个数
static inline int nearMatch(uint64_t a, uint64_t b) {
    uint64_t x = a ^ b;
    uint64_t z = (x | (x >> 1)) & 0x5555555555555555ull;
    return 32 - __builtin_popcountll(z);
}

// 加权拟合度（内部用 100·z，便于整数比较与阈值表达）
static inline int fpFit(const Fp& a, const Fp& b) {
    const int nr = nearMatch(a.near64, b.near64);
    const int fl = __builtin_popcountll(a.farL & b.farL);
    const int fr = __builtin_popcountll(a.farR & b.farR);
    const double z = g_wN * (nr - g_muN) / g_sdN
                   + g_wF * ((fl - g_muFL) / g_sdFL + (fr - g_muFR) / g_sdFR);
    return (int)std::lround(100.0 * z);
}

// ---------------------------------------------------------------- 3-mer 存在窗口

struct Presence3 {
    std::array<uint16_t, K3N> cnt{};
    uint64_t bits = 0;
    inline void add(int c) { if (cnt[c]++ == 0) bits |= 1ull << c; }
    inline void remove(int c) { if (--cnt[c] == 0) bits &= ~(1ull << c); }
};

// ---------------------------------------------------------------- seed 扫描
// enc 已带 PAD 长度的 L/R 补齐；只对完全落在真实区间 [realLo, realHi) 内的 seed 回调。
// 近指纹用两个滚动 32-bit 寄存器（左 16-mer / 右 16-mer），远指纹用两个 64-bit 存在位图。

// 读模式下判断 seed 起点 s 是否完全落在某个读内。
// 块布局 = padL + read + padR（数组坐标）；减掉 PAD 后读正好落在 rel∈[0,readLen)，
// padR 落 [readLen, readLen+PAD)、padL 落 [readLen+PAD, stride)。
static inline bool readSeedOk(int64_t s) {
    const int64_t rel = (s - PAD) % g_readStride;
    return rel + SEED_K <= g_readLen;
}

// 只判断"某位置是不是合格 seed"（不建指纹）：用于并行建库前的计数遍。
// 与 scanSeeds 的 emit 条件逐位一致（同样的 14-mer 滚动 + bad 计数），但省掉窗口/指纹开销。
template <class Fn>
inline void forEachSeedPos(const std::vector<int8_t>& enc, int64_t lo, int64_t hi, Fn&& fn) {
    const int64_t n = (int64_t)enc.size();
    lo = std::max<int64_t>(0, lo);
    hi = std::min<int64_t>(hi, n - SEED_K + 1);
    if (lo >= hi) return;
    uint32_t mer = 0;
    int bad = 0;
    for (int t = 0; t < SEED_K; ++t) {
        int b = enc[(size_t)(lo + t)];
        if (b < 0) ++bad; else mer = ((mer << 2) | (uint32_t)b) & SEED_MASK;
    }
    for (int64_t s = lo; s < hi; ++s) {
        if (bad == 0 && mer % SEED_MOD == 0) fn(s);
        if (s + 1 >= hi) break;
        int out = enc[(size_t)s], in = enc[(size_t)(s + SEED_K)];
        if (out < 0) --bad;
        if (in < 0) ++bad; else mer = ((mer << 2) | (uint32_t)in) & SEED_MASK;
    }
}

template <class Fn>
inline void scanSeeds(const std::vector<int8_t>& enc, int64_t realLo, int64_t realHi, Fn&& fn) {
    const int64_t n = (int64_t)enc.size();
    if (n < SEED_K) return;
    int64_t lo = std::max<int64_t>(realLo, PAD);
    int64_t hi = std::min<int64_t>(realHi, n - PAD) - SEED_K + 1;   // seed 起点上界
    if (lo >= hi) return;

    auto code3 = [&](int64_t j) -> int {
        int c = 0;
        for (int t = 0; t < K3; ++t) { int b = enc[(size_t)(j + t)]; if (b < 0) return -1; c = (c << 2) | b; }
        return c;
    };
    // 左远窗口：3-mer 起点 [s-48, s-19]（30 个）；右远窗口：[s+30, s+59]
    auto fLlo = [](int64_t s) { return s - NEAR_LEN - FAR_LEN; };
    auto fLhi = [](int64_t s) { return s - NEAR_LEN - K3; };
    auto fRlo = [](int64_t s) { return s + SEED_K + NEAR_LEN; };
    auto fRhi = [](int64_t s) { return s + SEED_K + NEAR_LEN + FAR_LEN - K3; };

    Presence3 pl, pr;
    for (int64_t j = fLlo(lo); j <= fLhi(lo); ++j) { int c = code3(j); if (c >= 0) pl.add(c); }
    for (int64_t j = fRlo(lo); j <= fRhi(lo); ++j) { int c = code3(j); if (c >= 0) pr.add(c); }

    // 滚动寄存器：左 16-mer（[s-16, s)）与右 16-mer（[s+14, s+30)）
    uint32_t nL = 0, nR = 0;
    for (int t = 0; t < NEAR_LEN; ++t) {
        int bl = enc[(size_t)(lo - NEAR_LEN + t)];
        int br = enc[(size_t)(lo + SEED_K + t)];
        if (bl >= 0) nL = (nL << 2) | (uint32_t)bl; else nL <<= 2;
        if (br >= 0) nR = (nR << 2) | (uint32_t)br; else nR <<= 2;
    }
    uint32_t mer = 0; int bad = 0;
    for (int t = 0; t < SEED_K; ++t) {
        int b = enc[(size_t)(lo + t)];
        if (b < 0) ++bad; else mer = ((mer << 2) | (uint32_t)b) & SEED_MASK;
    }

    for (int64_t s = lo; s < hi; ++s) {
        if (bad == 0 && mer % SEED_MOD == 0 && (!g_readMode || readSeedOk(s))) {
            Fp fp;
            fp.near64 = (uint64_t)nL | ((uint64_t)nR << 32);
            fp.farL = pl.bits;
            fp.farR = pr.bits;
            fn((uint64_t)s, mer, fp);
        }
        if (s + 1 >= hi) break;
        // 远窗口滑动
        if (fLlo(s) >= 0) { int c = code3(fLlo(s)); if (c >= 0) pl.remove(c); }
        { int c = code3(fLhi(s + 1)); if (c >= 0) pl.add(c); }
        { int c = code3(fRlo(s)); if (c >= 0) pr.remove(c); }
        { int c = code3(fRhi(s + 1)); if (c >= 0) pr.add(c); }
        // 近指纹滚动：各自左移 2 bit 并补入新碱基
        int lOut = enc[(size_t)(s - NEAR_LEN)], lIn = enc[(size_t)s];
        nL <<= 2; if (lIn >= 0) nL |= (uint32_t)lIn; else nL |= 0;
        (void)lOut;
        int rOut = enc[(size_t)(s + SEED_K)], rIn = enc[(size_t)(s + SEED_K + NEAR_LEN)];
        (void)rOut;
        nR <<= 2; if (rIn >= 0) nR |= (uint32_t)rIn;
        // 14-mer 滚动
        int out = enc[(size_t)s], in = enc[(size_t)(s + SEED_K)];
        if (out < 0) --bad;
        if (in < 0) ++bad; else mer = ((mer << 2) | (uint32_t)in) & SEED_MASK;
    }
}

// ---------------------------------------------------------------- FASTA

// 把大块内存标成 THP（alice 的 transparent_hugepage 是 madvise，
// 不显式请求就是 4 KB 页 → 17 GB 索引要 420 万个页表项）
static void useHugePages(void* p, size_t bytes) {
    if (!p || bytes < (size_t)(4 << 20)) return;
    const uintptr_t pg = 4096;
    const uintptr_t a = (uintptr_t)p;
    const uintptr_t start = (a + pg - 1) & ~(uintptr_t)(pg - 1);
    const size_t len = ((a + bytes) & ~(uintptr_t)(pg - 1)) - start;
    if (len) madvise((void*)start, len, MADV_HUGEPAGE);
}

struct FastaRec { std::string name; std::string seq; };

static std::vector<FastaRec> readFasta(const char* path) {
    std::vector<FastaRec> out;
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
    for (auto& r : out) for (auto& c : r.seq) c = (char)toupper((unsigned char)c);
    return out;
}

// 多记录拼接（记录间 N 间隔），并保留记录边界表
struct RecInfo { std::string name; int64_t off, len; };

struct Genome {
    std::vector<RecInfo> recs;
    bool readMode = false;
    int64_t stride = 0, readLen = 0;
    std::vector<int8_t> enc;      // 已含 L/R 补齐（建库扫描/长序列模式用）
    std::vector<uint8_t> pk;      // ref 的 4 bit/碱基打包副本（one-hot：N=0, A=1, C=2, G=4, T=8），2 碱基/字节
    int64_t total = 0;            // 拼接总长（含间隔）
    int64_t realLen = 0;          // 实际碱基+间隔
    int64_t padLeft = PAD;
};

// ref 建库完成后转成 4 bit/碱基：4 Gbp 从 3.73 GiB 降到 2.0 GiB。
// 只有读模式的比对阶段会用到 packed 副本（每次取 150 列解包），
// 建库扫描仍然走 enc（不压缩），所以这条改动不影响热路径。
static void packRef(Genome& g) {
    const size_t n = g.enc.size();
    const size_t bytes = (n + 1) / 2;
    g.pk.reserve(bytes);
    useHugePages(g.pk.data(), g.pk.capacity());     // 先 madvise 再 resize（THP 顺序坑）
    g.pk.assign(bytes, 0);
    for (size_t i = 0; i < n; ++i) {
        const int8_t b = g.enc[i];
        const uint8_t c = (b < 0) ? (uint8_t)0 : (uint8_t)(1u << b);   // N=0000, A/C/G/T=1/2/4/8
        if (i & 1) g.pk[i >> 1] |= (uint8_t)(c << 4);
        else g.pk[i >> 1] = (uint8_t)(g.pk[i >> 1] | c);
    }
    std::vector<int8_t>().swap(g.enc);
}

// one-hot 4 bit -> 内部 0..3（N/非法码 -> -1）
static const int8_t kOh2Base[16] = { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1 };

static inline void decodeRef(const Genome& g, int64_t from, int n, int8_t* out) {
    int64_t i = from;
    for (int k = 0; k < n; ++k, ++i) {
        const uint8_t byte = g.pk[(size_t)(i >> 1)];
        out[k] = kOh2Base[(i & 1) ? (byte >> 4) : (byte & 0x0F)];
    }
}

// L / R 补齐序列（确定性生成，Ref 与 Query 共用，保证两端指纹可比）
static std::string makePad(unsigned seed, int len) {
    std::mt19937 rng(seed);
    const char* B = "ACGT";
    std::string s((size_t)len, 'A');
    for (auto& c : s) c = B[rng() & 3];
    return s;
}

static void buildGenome(std::vector<FastaRec>& recs, const std::string& padL, const std::string& padR,
                        int64_t gap, Genome& g) {
    // 先算总长
    size_t seqLen = 0;
    for (auto& r : recs) seqLen += r.seq.size();
    g.enc.reserve(seqLen + recs.size() * (size_t)gap + 2 * (size_t)PAD);
    useHugePages(g.enc.data(), g.enc.capacity() * sizeof(int8_t));
    // 整块写（原先逐字节 push_back，4 Gb 参考要 4e9 次）
    const size_t total = seqLen + recs.size() * (size_t)gap + 2 * (size_t)PAD;
    g.enc.resize(total);
    int8_t* dst = g.enc.data();
    for (char c : padL) *dst++ = (int8_t)baseCode(c);
    int64_t off = 0;
    for (auto& r : recs) {
        g.recs.push_back({r.name, off, (int64_t)r.seq.size()});
        const unsigned char* src = (const unsigned char*)r.seq.data();
        for (size_t k = 0; k < r.seq.size(); ++k) *dst++ = (int8_t)g_codeTab[src[k]];
        for (int64_t k = 0; k < gap; ++k) *dst++ = -1;
        off += (int64_t)r.seq.size() + gap;
        std::string().swap(r.seq);
    }
    for (char c : padR) *dst++ = (int8_t)baseCode(c);
    g.realLen = off;
    g.total = (int64_t)g.enc.size();
}

// 读模式：把 FASTA/gz 直接流式打包成 padL+read+padR（不建 133M 条 FastaRec 中间体，
// 省 ~27 GB 峰值内存和 1.3 亿次小对象分配）
static void loadReadsPacked(const char* path, const std::string& padL, const std::string& padR,
                            int64_t readLen, Genome& g) {
    const size_t plen = strlen(path);
    const bool gz = plen > 3 && strcmp(path + plen - 3, ".gz") == 0;
    {
        long long sz = 0;
        if (!gz) {
            FILE* ft = fopen(path, "rb");
            if (ft) { fseek(ft, 0, SEEK_END); sz = ftell(ft); fclose(ft); }
        }
        if (sz > 0) {
            const size_t est = (size_t)(sz / (readLen + 2) + 4096);
            g.enc.reserve(est * (size_t)(readLen + 2 * PAD));
            g.recs.reserve(est);
        }
    }
    if (g.enc.capacity()) useHugePages(g.enc.data(), g.enc.capacity() * sizeof(int8_t));

    FILE* f = gz ? popen((std::string("gzip -dc '") + path + "'").c_str(), "r") : fopen(path, "rb");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); exit(1); }

    g.readMode = true;
    g.readLen = readLen;
    g.stride = readLen + 2 * PAD;
    std::vector<char> buf(1 << 20);
    std::string cur;
    cur.reserve((size_t)readLen + 8);
    bool started = false;
    // 每个读先在本块里拼好再整块 memcpy 进 g.enc：原来逐字节 push_back
    // （2e10 次 push_back + 分支判断）是 load 阶段 85 s 的元凶
    std::vector<int8_t> blk((size_t)g.stride);
    for (int64_t i = 0; i < PAD; ++i) {
        blk[(size_t)i] = (int8_t)baseCode(padL[(size_t)i]);
        blk[(size_t)(g.stride - PAD + i)] = (int8_t)baseCode(padR[(size_t)i]);
    }
    auto flush = [&]() {
        const size_t n = std::min(cur.size(), (size_t)readLen);
        int8_t* dst = blk.data() + PAD;
        const unsigned char* src = (const unsigned char*)cur.data();
        for (size_t k = 0; k < n; ++k) dst[k] = (int8_t)g_codeTab[src[k]];
        for (size_t k = n; k < (size_t)readLen; ++k) dst[k] = -1;
        g.enc.insert(g.enc.end(), blk.begin(), blk.end());
        cur.clear();
    };
    while (fgets(buf.data(), (int)buf.size(), f)) {
        char* p = buf.data();
        if (*p == '>') { if (started) flush(); started = true; continue; }
        if (!started) continue;
        for (; *p; ++p)
            if (*p != '\n' && *p != '\r') cur.push_back((char)toupper((unsigned char)*p));
    }
    if (started) flush();
    if (gz) pclose(f); else fclose(f);

    const int64_t nRec = (int64_t)g.enc.size() / g.stride;
    g.recs.resize((size_t)nRec);
    for (int64_t i = 0; i < nRec; ++i) {
        g.recs[(size_t)i].off = i * g.stride + PAD;
        g.recs[(size_t)i].len = readLen;
    }
    g.realLen = nRec * g.stride;
    g.total = (int64_t)g.enc.size();
    g.padLeft = PAD;
}


// 读模式流式读入：每次打包 maxReads 条读进 Genome（块内本地序号 0..n-1）
struct ReadStream {
    FILE* f = nullptr;
    bool gz = false;
    std::vector<char> buf;
    std::string cur;
    std::string padL, padR;
    int64_t readLen = 0, stride = 0;
    std::vector<int8_t> blk;
    bool started = false;

    bool open(const char* path, const std::string& pl, const std::string& pr, int64_t rl) {
        const size_t plen = strlen(path);
        gz = plen > 3 && strcmp(path + plen - 3, ".gz") == 0;
        f = gz ? popen((std::string("gzip -dc '") + path + "'").c_str(), "r") : fopen(path, "rb");
        if (!f) return false;
        padL = pl; padR = pr; readLen = rl; stride = rl + 2 * PAD;
        buf.resize(1 << 20);
        cur.reserve((size_t)rl + 8);
        blk.resize((size_t)stride);
        for (int64_t i = 0; i < PAD; ++i) {
            blk[(size_t)i] = (int8_t)baseCode(padL[(size_t)i]);
            blk[(size_t)(stride - PAD + i)] = (int8_t)baseCode(padR[(size_t)i]);
        }
        return true;
    }
    ~ReadStream() { if (f) { if (gz) pclose(f); else fclose(f); } }

    void flushOne(Genome& g) {
        const size_t n = std::min(cur.size(), (size_t)readLen);
        int8_t* dst = blk.data() + PAD;
        const unsigned char* sb = (const unsigned char*)cur.data();
        for (size_t k = 0; k < n; ++k) dst[k] = (int8_t)g_codeTab[sb[k]];
        for (size_t k = n; k < (size_t)readLen; ++k) dst[k] = -1;
        g.enc.insert(g.enc.end(), blk.begin(), blk.end());
        cur.clear();
    }

    // 读入至多 maxReads 条；返回条数（0 = 结束）
    int64_t nextChunk(Genome& g, size_t maxReads) {
        std::vector<int8_t>().swap(g.enc);
        std::vector<RecInfo>().swap(g.recs);
        g.readMode = true; g.readLen = readLen; g.stride = stride;
        g.enc.reserve(maxReads * (size_t)stride);
        useHugePages(g.enc.data(), g.enc.capacity());
        g.recs.reserve(maxReads);
        size_t got = 0;
        bool boundary = false;
        while (!boundary && f && fgets(buf.data(), (int)buf.size(), f)) {
            char* p = buf.data();
            if (*p == '>') {
                if (started && !cur.empty()) {
                    flushOne(g);
                    ++got;
                    if (got >= maxReads) boundary = true;   // 本块已满，这条 '>' 属下一块
                }
                started = true;
                continue;
            }
            if (!started) continue;
            for (; *p; ++p)
                if (*p != '\n' && *p != '\r') cur.push_back((char)toupper((unsigned char)*p));
        }
        if (!boundary) {
            if (started && !cur.empty()) { flushOne(g); ++got; }
            started = false;
        }
        const int64_t n = (int64_t)(g.enc.size() / (size_t)stride);
        g.recs.resize((size_t)n);
        for (int64_t i = 0; i < n; ++i) {
            g.recs[(size_t)i].off = i * stride + PAD;
            g.recs[(size_t)i].len = readLen;
        }
        g.realLen = n * stride;
        g.total = (int64_t)g.enc.size();
        g.padLeft = PAD;
        return n;
    }
};

// 读模式：每个读 = padL + read + padR 紧凑排布（无 N 间隔），记录号 O(1) 反算
static void buildReads(std::vector<FastaRec>& recs, const std::string& padL, const std::string& padR,
                       int64_t readLen, Genome& g) {
    g.readMode = true;
    g.readLen = readLen;
    g.stride = readLen + 2 * PAD;
    g.recs.reserve(recs.size());
    g.enc.reserve(recs.size() * (size_t)g.stride);
    for (size_t i = 0; i < recs.size(); ++i) {
        auto& r = recs[i];
        g.recs.push_back({"", (int64_t)i * g.stride + PAD, readLen});
        for (char c : padL) g.enc.push_back((int8_t)baseCode(c));
        const size_t L = std::min(r.seq.size(), (size_t)readLen);
        for (size_t k = 0; k < L; ++k) g.enc.push_back((int8_t)baseCode(r.seq[k]));
        for (size_t k = L; k < (size_t)readLen; ++k) g.enc.push_back(-1);   // 短读用 N 补
        for (char c : padR) g.enc.push_back((int8_t)baseCode(c));
        std::string().swap(r.seq);
    }
    g.realLen = (int64_t)recs.size() * g.stride;   // 打包总长（用于反链坐标翻转）
    g.total = (int64_t)g.enc.size();
    g.padLeft = PAD;
}

static int locateIdx(const Genome& g, int64_t pos) {
    if (g.readMode) {                     // 读模式：等长，直接除
        int64_t d = pos - PAD;
        if (d < 0) return -1;
        int64_t r = d / g.stride;
        if (r >= (int64_t)g.recs.size()) return -1;
        if (d - r * g.stride >= g.readLen) return -1;
        return (int)r;
    }
    size_t lo = 0, hi = g.recs.size();
    if (hi == 0) return -1;
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) / 2;
        if (g.recs[mid].off <= pos) lo = mid; else hi = mid;
    }
    const RecInfo& r = g.recs[lo];
    if (pos >= r.off && pos < r.off + r.len) return (int)lo;
    return -1;
}

// ---------------------------------------------------------------- Smith-Waterman（锚定延伸用）

struct SwOut { int score = 0; int matches = 0; int alen = 0; std::string cigar; };

static SwOut swLocal(const int8_t* a, int la, const int8_t* b, int lb) {
    static thread_local std::vector<int16_t> H;
    static thread_local std::vector<uint8_t> T;
    const int W = lb + 1;
    const size_t need = (size_t)(la + 1) * W;
    if (H.size() < need) { H.resize(need); T.resize(need); }
    std::fill(H.begin(), H.begin() + (long)need, (int16_t)0);
    std::fill(T.begin(), T.begin() + (long)need, (uint8_t)0);
    for (int i = 1; i <= la; ++i) {
        const int8_t ca = a[i - 1];
        int16_t* hRow = &H[(size_t)i * W];
        int16_t* hPrev = &H[(size_t)(i - 1) * W];
        uint8_t* tRow = &T[(size_t)i * W];
        for (int j = 1; j <= lb; ++j) {
            int best = 0; uint8_t d = 0;
            int v = hPrev[j - 1] + (ca == b[j - 1] ? 2 : -1);
            if (v > best) { best = v; d = 1; }
            v = hPrev[j] - 2;
            if (v > best) { best = v; d = 2; }
            v = hRow[j - 1] - 2;
            if (v > best) { best = v; d = 3; }
            hRow[j] = (int16_t)best; tRow[j] = d;
        }
    }
    int mi = 0, mj = 0, best = 0;
    for (int i = 1; i <= la; ++i)
        for (int j = 1; j <= lb; ++j) {
            int v = H[(size_t)i * W + j];
            if (v > best) { best = v; mi = i; mj = j; }
        }
    SwOut r; r.score = best;
    if (best <= 0) return r;
    std::string ops;
    int i = mi, j = mj;
    while (i > 0 && j > 0) {
        uint8_t d = T[(size_t)i * W + j];
        if (d == 0) break;
        if (d == 1) { ops.push_back(a[i-1] == b[j-1] ? 'M' : 'X'); if (a[i-1] == b[j-1]) ++r.matches; --i; --j; }
        else if (d == 2) { ops.push_back('D'); --i; }
        else { ops.push_back('I'); --j; }
    }
    r.alen = (int)ops.size();
    for (int k = (int)ops.size() - 1; k >= 0; ) {
        char c = ops[(size_t)k]; int len = 0;
        while (k >= 0 && ops[(size_t)k] == c) { --k; ++len; }
        r.cigar += std::to_string(len);
        r.cigar.push_back(c);
    }
    return r;
}

// 自由端 gap 的全局比对（锚定延伸用）：两端不罚分，identity 覆盖整段窗口
struct NwOut { int score = 0; int matches = 0; int alen = 0; };

static NwOut nwFreeEnds(const int8_t* a, int la, const int8_t* b, int lb) {
    static thread_local std::vector<int16_t> H;
    static thread_local std::vector<uint8_t> T;
    const int W = lb + 1;
    const size_t need = (size_t)(la + 1) * W;
    if (H.size() < need) { H.resize(need); T.resize(need); }
    for (int j = 0; j <= lb; ++j) { H[(size_t)j] = 0; T[(size_t)j] = (j ? 3 : 0); }   // 首行自由
    for (int i = 1; i <= la; ++i) {
        H[(size_t)i * W] = 0; T[(size_t)i * W] = 2;                                      // 首列自由
        const int8_t ca = a[i - 1];
        int16_t* hRow = &H[(size_t)i * W];
        int16_t* hPrev = &H[(size_t)(i - 1) * W];
        uint8_t* tRow = &T[(size_t)i * W];
        for (int j = 1; j <= lb; ++j) {
            int best = hPrev[j - 1] + (ca == b[j - 1] ? 2 : -1); uint8_t d = 1;
            int v = hPrev[j] - 2; if (v > best) { best = v; d = 2; }
            v = hRow[j - 1] - 2; if (v > best) { best = v; d = 3; }
            hRow[j] = (int16_t)best; tRow[j] = d;
        }
    }
    // 终点取末行/末列的最大值（另一端自由）
    int bi = la, bj = lb, best = INT_MIN;
    for (int j = 0; j <= lb; ++j) if (H[(size_t)la * W + j] > best) { best = H[(size_t)la * W + j]; bj = j; }
    for (int i = 0; i <= la; ++i) if (H[(size_t)i * W + lb] > best) { best = H[(size_t)i * W + lb]; bi = i; }
    NwOut r; r.score = best;
    int i = bi, j = bj;
    while (i > 0 && j > 0) {
        uint8_t d = T[(size_t)i * W + j];
        if (d == 1) { if (a[i-1] == b[j-1]) ++r.matches; ++r.alen; --i; --j; }
        else if (d == 2) { ++r.alen; --i; }
        else { ++r.alen; --j; }
    }
    return r;
}

// 带状自由端全局比对（读模式算 identity 用）：band 内 DP，两端不罚分
static NwOut nwBand(const int8_t* a, int la, const int8_t* b, int lb, int band) {
    static thread_local std::vector<int16_t> H;
    static thread_local std::vector<uint8_t> T;
    const int W = 2 * band + 1;
    const size_t need = (size_t)(la + 1) * (size_t)W;
    if (H.size() < need) { H.resize(need); T.resize(need); }
    const int16_t NEG = -30000;
    for (int k = 0; k < W; ++k) { H[(size_t)k] = 0; T[(size_t)k] = 0; }   // 首行自由
    for (int i = 1; i <= la; ++i) {
        int16_t* cur = &H[(size_t)i * W];
        int16_t* prev = &H[(size_t)(i - 1) * W];
        uint8_t* tc = &T[(size_t)i * W];
        for (int k = 0; k < W; ++k) {
            const int j = i - band + k;
            if (j < 0 || j > lb) { cur[k] = NEG; tc[k] = 0; continue; }
            if (j == 0) { cur[k] = 0; tc[k] = 2; continue; }               // 首列自由
            int best = (prev[k] > NEG / 2) ? prev[k] + (a[i - 1] == b[j - 1] ? 2 : -1) : NEG;
            uint8_t d = 1;
            int v = (k + 1 < W) ? prev[k + 1] - 2 : NEG;
            if (v > best) { best = v; d = 2; }
            v = (k > 0 && cur[k - 1] > NEG / 2) ? cur[k - 1] - 2 : NEG;
            if (v > best) { best = v; d = 3; }
            cur[k] = (int16_t)best; tc[k] = d;
        }
    }
    int bi = la, bj = lb, best = INT_MIN;
    for (int k = 0; k < W; ++k) {                       // 末行
        const int j = la - band + k;
        if (j < 0 || j > lb) continue;
        if (H[(size_t)la * W + k] > best) { best = H[(size_t)la * W + k]; bj = j; }
    }
    for (int i = std::max(0, lb - band); i <= la; ++i) {  // 末列
        const int k = lb - i + band;
        if (k < 0 || k >= W) continue;
        if (H[(size_t)i * W + k] > best) { best = H[(size_t)i * W + k]; bi = i; }
    }
    NwOut r; r.score = best;
    int i = bi, j = bj;
    while (i > 0 && j > 0) {
        const int k = j - i + band;
        if (k < 0 || k >= W) break;
        const uint8_t d = T[(size_t)i * W + k];
        if (d == 1) { if (a[i - 1] == b[j - 1]) ++r.matches; ++r.alen; --i; --j; }
        else if (d == 2) { ++r.alen; --i; }
        else if (d == 3) { ++r.alen; --j; }
        else break;
    }
    return r;
}

// 带状全局比对（两端不自由）：读 identity 用。
// 自由端版本会靠 gap 跳过错配，随机序列也能刷到 0.55；全局版随机只有 ~0.30，真值 ~0.75。
// ops 非空时同时回填路径（正向序）：1=对角(双方各取一个碱基) 2=ref 前进/query 空位 3=query 前进/ref 空位
static NwOut nwBandGlobal(const int8_t* a, int la, const int8_t* b, int lb, int band,
                          std::vector<uint8_t>* ops = nullptr) {
    static thread_local std::vector<int16_t> H;
    static thread_local std::vector<uint8_t> T;
    const int W = 2 * band + 1;
    const size_t need = (size_t)(la + 1) * (size_t)W;
    if (H.size() < need) { H.resize(need); T.resize(need); }
    const int16_t NEG = -30000;
    auto idxOf = [&](int i, int j) { return (size_t)i * W + (size_t)(j - i + band); };
    for (int j = 0; j < W; ++j) {
        const int col = j - band;                      // 第 0 行：a 的前导 gap
        H[(size_t)j] = (col >= 0 && col <= lb) ? (int16_t)(-4 * col) : NEG;
        T[(size_t)j] = 3;
    }
    for (int i = 1; i <= la; ++i) {
        int16_t* cur = &H[(size_t)i * W];
        int16_t* prev = &H[(size_t)(i - 1) * W];
        uint8_t* tc = &T[(size_t)i * W];
        for (int k = 0; k < W; ++k) {
            const int j = i - band + k;
            if (j < 0 || j > lb) { cur[k] = NEG; tc[k] = 0; continue; }
            if (j == 0) { cur[k] = (i <= band) ? (int16_t)(-4 * i) : NEG; tc[k] = 2; continue; }
            int best = (prev[k] > NEG / 2) ? prev[k] + (a[i - 1] == b[j - 1] ? 2 : -1) : NEG;
            uint8_t d = 1;
            int v = (k + 1 < W) ? prev[k + 1] - 4 : NEG;
            if (v > best) { best = v; d = 2; }
            v = (k > 0 && cur[k - 1] > NEG / 2) ? cur[k - 1] - 4 : NEG;
            if (v > best) { best = v; d = 3; }
            cur[k] = (int16_t)best; tc[k] = d;
        }
    }
    NwOut r;
    const int kEnd = lb - la + band;
    if (kEnd < 0 || kEnd >= W) return r;               // 超出带宽
    r.score = H[idxOf(la, lb)];
    int i = la, j = lb;
    if (ops) { ops->clear(); ops->reserve((size_t)(la + lb)); }
    while (i > 0 && j > 0) {
        const int k = j - i + band;
        if (k < 0 || k >= W) break;
        const uint8_t d = T[(size_t)i * W + k];
        if (d == 1) { if (a[i - 1] == b[j - 1]) ++r.matches; ++r.alen; --i; --j; }
        else if (d == 2) { ++r.alen; --i; }
        else if (d == 3) { ++r.alen; --j; }
        else break;
        if (ops) ops->push_back(d);
    }
    // 回溯在 i==0 或 j==0 就停了，首/尾还会剩若干未记录的前导 gap。
    // 不补的话，调用方拿到的路径会整体右移，按参考列填 MSA 时全错位。
    if (ops) {
        while (i > 0) { ops->push_back(2); --i; }   // 前导：ref 有、query 空
        while (j > 0) { ops->push_back(3); --j; }   // 前导：query 有、ref 空
    }
    if (ops) std::reverse(ops->begin(), ops->end());
    return r;
}

// ---------------------------------------------------------------- selftest

static int runSelfTest() {
    const int n = 4000;
    std::mt19937_64 rng(11);
    std::vector<int8_t> enc((size_t)(n + 2 * PAD));
    for (auto& c : enc) c = (int8_t)(rng() % 4);
    for (int i = 0; i < n + 2 * PAD; i += 37) enc[(size_t)i] = -1;   // 掺 N

    auto bruteSeed = [&](int64_t s) {
        uint32_t mer = 0;
        for (int t = 0; t < SEED_K; ++t) { int b = enc[(size_t)(s + t)]; if (b < 0) return false; mer = (mer << 2) | (uint32_t)b; }
        return mer % SEED_MOD == 0;
    };
    auto brute = [&](int64_t s, Fp& fp) {
        uint32_t nL = 0, nR = 0;
        for (int t = 0; t < NEAR_LEN; ++t) {
            int bl = enc[(size_t)(s - NEAR_LEN + t)], br = enc[(size_t)(s + SEED_K + t)];
            nL = (nL << 2) | (uint32_t)(bl < 0 ? 0 : bl);
            nR = (nR << 2) | (uint32_t)(br < 0 ? 0 : br);
        }
        fp.near64 = (uint64_t)nL | ((uint64_t)nR << 32);
        fp.farL = fp.farR = 0;
        for (int64_t j = s - NEAR_LEN - FAR_LEN; j + K3 <= s - NEAR_LEN; ++j) {
            int c = 0; bool ok = true;
            for (int t = 0; t < K3; ++t) { int b = enc[(size_t)(j + t)]; if (b < 0) { ok = false; break; } c = (c << 2) | b; }
            if (ok) fp.farL |= 1ull << c;
        }
        for (int64_t j = s + SEED_K + NEAR_LEN; j + K3 <= s + SEED_K + NEAR_LEN + FAR_LEN; ++j) {
            int c = 0; bool ok = true;
            for (int t = 0; t < K3; ++t) { int b = enc[(size_t)(j + t)]; if (b < 0) { ok = false; break; } c = (c << 2) | b; }
            if (ok) fp.farR |= 1ull << c;
        }
    };

    int checked = 0, bad = 0, missed = 0;
    std::vector<char> got((size_t)(n + 2 * PAD + 1), 0);
    for (int chunk : {0, 1, 777, 1999}) {
        int64_t realLo = PAD + chunk, realHi = PAD + n;
        scanSeeds(enc, realLo, realHi, [&](uint32_t s, uint32_t, const Fp& fp) {
            ++checked;
            got[s] = 1;
            Fp b; brute(s, b);
            if (b.near64 != fp.near64 || b.farL != fp.farL || b.farR != fp.farR) {
                ++bad;
                if (bad < 4) fprintf(stderr, "  [selftest] FP MISMATCH s=%u\n", s);
            }
        });
    }
    for (int64_t s = PAD; s + SEED_K <= PAD + n; ++s)
        if (bruteSeed(s) && !got[(size_t)s]) {
            ++missed;
            if (missed < 4) fprintf(stderr, "  [selftest] LOST s=%lld\n", (long long)s);
        }
    printf("[selftest] seeds=%d fp-mismatch=%d lost=%d -> %s\n",
           checked, bad, missed, (bad == 0 && missed == 0) ? "PASS" : "FAIL");
    return (bad == 0 && missed == 0) ? 0 : 1;
}

// ---------------------------------------------------------------- main

int main(int argc, char** argv) {
    int minFit = 250;          // 100·z 单位，默认 2.5σ
    int threads = (int)std::min(16u, std::max(1u, std::thread::hardware_concurrency()));
    int occCap = 200;          // 桶占用上限
    int maxGap = 20000;        // 链内最大 gap
    int minAnchors = 3;        // 块最少锚点
    int64_t minSpan = 500;     // 块最小跨度
    int64_t extLen = 1000;     // 两端锚定延伸窗口
    int lookback = 64;
    double driftFrac = 0.125, gapCost = 0.02;
    int driftAbs = 40;
    int maxBlocks = 200000;
    bool selfMode = false;     // -S：自比对，排除同位点对角线
    int64_t selfDiag = 500;    // |rpos - qpos| < 该值视为同位点
    unsigned padSeed = 20261125;
    double wN = 1.0, wF = 1.0;
    double fdr = 0.0;
    bool fdrSet = false;
    int64_t chunkReads = 8000000;   // --chunk：读模式每块的读条数
    const char* msaOut = nullptr;   // --msa：按 ref 骨架输出 query 序列（同数量、同长度，可直接 cat 成 MSA）
    int gapMask = 2;          // 读模式 MSA：gap 两侧各多少列不投票
    int topK = 5;             // 读模式：每读比对前 K 个锚点，取 identity 最好
    double minReadId = 0.70;  // 读模式：最低 identity（随机最高可达 ~0.68，20% 分化真值 p1≈0.73）
    int readLenOpt = 0;        // >0 时启用读模式（查询按定长读紧凑打包）
    bool selftest = false;
    const char* anchorsOut = nullptr;
    std::vector<const char*> pos;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto needInt = [&](auto& dst) {
            if (i + 1 >= argc) { fprintf(stderr, "ERROR: %s needs a value\n", a.c_str()); exit(2); }
            if constexpr (std::is_same_v<std::decay_t<decltype(dst)>, int64_t>) dst = atoll(argv[++i]);
            else dst = atoi(argv[++i]);
        };
        if (a == "-f") needInt(minFit);
        else if (a == "-t") needInt(threads);
        else if (a == "-c") needInt(occCap);
        else if (a == "-g") needInt(maxGap);
        else if (a == "-n") needInt(minAnchors);
        else if (a == "-S" || a == "--self") selfMode = true;
        else if (a == "--self-diag") { if (i + 1 < argc) selfDiag = atoll(argv[++i]); }
        else if (a == "--min-span") needInt(minSpan);
        else if (a == "--ext") needInt(extLen);
        else if (a == "--lookback") needInt(lookback);
        else if (a == "--drift-frac") { if (i + 1 < argc) driftFrac = atof(argv[++i]); }
        else if (a == "--drift-abs") needInt(driftAbs);
        else if (a == "--gap-cost") { if (i + 1 < argc) gapCost = atof(argv[++i]); }
        else if (a == "--max-blocks") needInt(maxBlocks);
        else if (a == "--pad-seed") { if (i + 1 < argc) padSeed = (unsigned)atoi(argv[++i]); }
        else if (a == "--chunk") { if (i + 1 < argc) chunkReads = atoll(argv[++i]); }
        else if (a == "--msa") { if (i + 1 < argc) msaOut = argv[++i]; }
        else if (a == "--gapmask") { if (i + 1 < argc) gapMask = atoi(argv[++i]); }
        else if (a == "--topk") { if (i + 1 < argc) topK = atoi(argv[++i]); }
        else if (a == "--min-id") { if (i + 1 < argc) minReadId = atof(argv[++i]); }
        else if (a == "--fdr") { if (i + 1 < argc) { fdr = atof(argv[++i]); fdrSet = true; } }
        else if (a == "--wn") { if (i + 1 < argc) wN = atof(argv[++i]); }
        else if (a == "--wf") { if (i + 1 < argc) wF = atof(argv[++i]); }
        else if (a == "-A") { if (i + 1 < argc) anchorsOut = argv[++i]; }
        else if (a == "--reads") { if (i + 1 < argc) readLenOpt = atoi(argv[++i]); }
        else if (a == "--selftest") selftest = true;
        else if (a == "-h" || a == "--help") {
            fprintf(stderr,
"Usage: wga3 [options] ref.fa[.gz] query.fa[.gz]\n"
"  -f INT     最小分数 (100·z；长序列默认 250=2.5σ，读模式默认 400)\n"
"  -t INT     线程数 (默认 min(16,cores))\n"
"  -c INT     seed 桶占用上限 (默认 200)\n"
"  -g INT     链内最大 gap (默认 20000)\n"
"  -n INT     块最少锚点数 (默认 3)\n"
"  --min-span INT 块最小跨度 (默认 500)\n"
"  --ext INT  两端锚定延伸窗口 (默认 1000, 0=不延伸)\n"
"  --wn F     近指纹权重 (默认 1.0)\n"
"  --wf F     远指纹权重 (默认 1.0)\n"
"  --fdr F    按实测零假设分位数自动定 -f（如 1e-4）\n"
"  --reads INT 读模式：查询按定长读紧凑打包（150bp 短读必用）\n"
"  --chunk INT 读模式：每块读条数 (默认 8000000；内存紧就调小，内存≈17+4+块×246×2 B)\n"
"  --msa FILE  按 ref 骨架输出 query 序列（同记录数/同列数，可直接 cat 成 MSA）\n"
"              画布 1 B/列 = 4 个 2-bit 计数（4 Gbp 约 4 GB）\n"
"              多数票；平票 / 次大>=2 / 仅 1 票 → N（判不了就退化，不造碱基）\n"
"  --topk INT 读模式：每读比对前 K 个锚点取 identity 最好 (默认 5)\n"
"  --min-id F 读模式：最低 identity (默认 0.70)\n"
"  --gapmask INT MSA：gap 两侧各多少列不投票 (默认 2；indel 的 gap 位置有 ±1~2 歧义)\n"
"  -S         自比对（排除同位点对角线；找重复用）\n"
"  --self-diag INT 自比对判同位点的阈值 (默认 500 bp)\n"
"  -A FILE    输出锚点明细\n"
"  --selftest 指纹扫描对拍自检\n");
            return 0;
        }
        else if (a.size() && a[0] == '-') { fprintf(stderr, "ERROR: unknown option %s\n", a.c_str()); return 2; }
        else pos.push_back(argv[i]);
    }
    if (selftest) return runSelfTest();
    if (pos.size() != 2) { fprintf(stderr, "Usage: wga3 [options] ref.fa query.fa  (see -h)\n"); return 1; }
    g_wN = wN; g_wF = wF;

    auto T0 = std::chrono::steady_clock::now();
    auto ms = [&T0]() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - T0).count(); };

    // ---- load ----
    const std::string padL = makePad(padSeed, PAD), padR = makePad(padSeed + 1, PAD);
    Genome R, Q;
    {
        auto refRecs = readFasta(pos[0]);
        if (refRecs.empty()) { fprintf(stderr, "ERROR: empty ref FASTA\n"); return 1; }
        buildGenome(refRecs, padL, padR, 1000, R);
        std::vector<FastaRec>().swap(refRecs);
    }
    fprintf(stderr, "[load] ref %zu 记录/%lld bp | %.0f ms\n", R.recs.size(), (long long)R.realLen, ms());

    // ---- 建库（并行两遍：先计数再填充，内存不变；单线程扫描 4 Gb 要 155 s）----
    struct Entry { uint32_t key; uint32_t pos; Fp fp; };   // 4+4+24 = 32 B
    std::vector<Entry> idx;
    uint64_t refSeeds = 0;
    g_readMode = false;
    {
        const int nt = std::max(1, threads);
        const int64_t span = R.realLen;
        std::vector<int64_t> cnt((size_t)nt, 0), off((size_t)nt + 1, 0);
        {
            std::vector<std::thread> th;
            for (int t = 0; t < nt; ++t) {
                th.emplace_back([&, t] {
                    const int64_t lo = span * t / nt, hi = span * (t + 1) / nt;
                    int64_t c = 0;
                    forEachSeedPos(R.enc, PAD + lo, PAD + hi, [&](int64_t) { ++c; });
                    cnt[(size_t)t] = c;
                });
            }
            for (auto& x : th) x.join();
        }
        for (int t = 0; t < nt; ++t) off[(size_t)t + 1] = off[(size_t)t] + cnt[(size_t)t];
        refSeeds = (uint64_t)off[(size_t)nt];
        if (refSeeds == 0) { fprintf(stderr, "ERROR: ref 无有效 seed\n"); return 1; }
        idx.reserve((size_t)refSeeds);
        useHugePages(idx.data(), idx.capacity() * sizeof(Entry));
        idx.resize((size_t)refSeeds);
        fprintf(stderr, "[index] 计数 %llu 条 %.0f ms", (unsigned long long)refSeeds, ms());
        {
            std::vector<std::thread> th;
            for (int t = 0; t < nt; ++t) {
                th.emplace_back([&, t] {
                    const int64_t lo = span * t / nt, hi = span * (t + 1) / nt;
                    size_t w = (size_t)off[(size_t)t];
                    scanSeeds(R.enc, PAD + lo, PAD + hi, [&](uint64_t s, uint32_t mer, const Fp& fp) {
                        idx[w++] = {mer, (uint32_t)(s - (uint32_t)PAD), fp};   // 存"真实坐标"（去掉左补齐）
                    });
                });
            }
            for (auto& x : th) x.join();
        }
        fprintf(stderr, " | 填充 %.0f ms", ms());
    }
    if (idx.empty()) { fprintf(stderr, "ERROR: ref 无有效 seed\n"); return 1; }
    std::sort(idx.begin(), idx.end(), [](const Entry& a, const Entry& b) { return a.key < b.key; });
    fprintf(stderr, " | 排序 %.0f ms", ms());

    // 直址范围表：key 必为 7 的倍数 → key/7 唯一
    std::vector<uint32_t> rStart(RANGE_N, 0), rCount(RANGE_N, 0);
    for (size_t i = 0; i < idx.size();) {
        size_t j = i;
        while (j < idx.size() && idx[j].key == idx[i].key) ++j;
        uint32_t cnt = (uint32_t)(j - i);
        uint32_t k7 = idx[i].key / SEED_MOD;
        if (k7 < RANGE_N) { rStart[k7] = (uint32_t)i; rCount[k7] = std::min<uint32_t>(cnt, (uint32_t)occCap); }
        i = j;
    }
    fprintf(stderr, "[index] seeds=%llu distinct=%zu idx=%.1f MB table=%.0f MB | %.0f ms\n",
            (unsigned long long)refSeeds, idx.size(), idx.size() * sizeof(Entry) / 1048576.0,
            (rStart.size() + rCount.size()) * 4.0 / 1048576.0, ms());

    // 建库扫描已完结 → ref 转 4 bit/碱基常驻（3.73 GiB → 2.0 GiB，只影响比对阶段的取窗口）
    if (readLenOpt > 0) {
        const size_t before = R.enc.size() / 1048576;
        packRef(R);
        fprintf(stderr, "[rev] ref 打包 4 bit/碱基: %zu MB -> %zu MB | %.0f ms\n",
                before, R.pk.size() / 1048576, ms());
    }

    // ---- 零假设标定（对任意查询子集调用一次即可）----
    auto runCalib = [&](const Genome& Qc) {
    {
        std::vector<int> vn, vfl, vfr;
        uint64_t rs = 88172645463325252ull;
        const int64_t calibSpan = std::min<int64_t>(Qc.realLen, (int64_t)20000000);
        scanSeeds(Qc.enc, PAD, PAD + calibSpan, [&](uint32_t, uint32_t, const Fp& fp) {
            if (vn.size() >= 200000) return;
            rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
            const Fp& e = idx[(size_t)(rs % idx.size())].fp;
            vn.push_back(nearMatch(fp.near64, e.near64));
            vfl.push_back(__builtin_popcountll(fp.farL & e.farL));
            vfr.push_back(__builtin_popcountll(fp.farR & e.farR));
        });
        auto stat = [](std::vector<int>& v, double& m, double& sd) {
            if (v.size() < 8) return;
            std::sort(v.begin(), v.end());
            m = 0; for (int x : v) m += x; m /= (double)v.size();
            sd = 0; for (int x : v) sd += (x - m) * (x - m); sd = std::sqrt(sd / (double)v.size());
        };
        stat(vn, g_muN, g_sdN);
        stat(vfl, g_muFL, g_sdFL);
        stat(vfr, g_muFR, g_sdFR);
        fprintf(stderr, "[calib] 零假设分量: 近场 %.1f±%.1f | 远左 %.1f±%.1f | 远右 %.1f±%.1f (n=%zu)\n",
                g_muN, g_sdN, g_muFL, g_sdFL, g_muFR, g_sdFR, vn.size());

        if (!vn.empty()) {
            std::vector<int> sc(vn.size());
            for (size_t i = 0; i < vn.size(); ++i) {
                const double z = g_wN * ((double)vn[i] - g_muN) / g_sdN
                               + g_wF * (((double)vfl[i] - g_muFL) / g_sdFL
                                       + ((double)vfr[i] - g_muFR) / g_sdFR);
                sc[i] = (int)std::lround(100.0 * z);
            }
            std::sort(sc.begin(), sc.end());
            auto q = [&](double p) {
                size_t i = (size_t)(p * (double)(sc.size() - 1));
                return sc[std::min(i, sc.size() - 1)];
            };
            fprintf(stderr, "[calib] 零假设分数分位: p50=%d p99=%d p99.9=%d p99.99=%d p99.999=%d max=%d\n",
                    q(0.5), q(0.99), q(0.999), q(0.9999), q(0.99999), sc.back());
            const double bucket = (double)idx.size() / (double)RANGE_N;
            const double seedsPerRead = Qc.readMode ? (double)(Qc.readLen - SEED_K + 1) / SEED_MOD : 0.0;
            fprintf(stderr, "[calib] 桶均=%.1f 读候选≈%.0f/读 %s\n", bucket, seedsPerRead * bucket,
                    Qc.readMode ? "(=每读候选数；假锚点/读 = 候选×FPR)" : "");
            if (fdr > 0) {
                int th = q(1.0 - fdr);
                fprintf(stderr, "[calib] --fdr %g -> 阈值 -f %d（当前 -f %d）\n", fdr, th, minFit);
                minFit = th;   // --fdr 显式指定时以分位数为准
            } else {
                fprintf(stderr, "[calib] 当前 -f %d ≈ 零假设分位 p=%.4g；建议 -f %d 可把每读假锚点压到 ~0.03\n",
                        minFit, [&] {
                            size_t k = (size_t)(std::lower_bound(sc.begin(), sc.end(), minFit) - sc.begin());
                            return 1.0 - (double)k / (double)sc.size();
                        }(), q(1.0 - 0.03 / std::max(1.0, seedsPerRead * bucket)));
            }
        }
    }

    };

    // MSA 画布：每列 1 字节 = A/C/G/T 各 2 bit 计数（饱和到 3）→ 4 Gbp 占 4 GB
    //   为什么不用 4 bit one-hot：没有票数时只能"先到者占优"，深度 2 的位置
    //   一旦第一条 read 带错，就会把错碱基原样输出（错误率远超 0.1% 的约束）。
    //   判定：取最大计数 i1；若平票或次大 ≥2 → N（杂合或无法判定），否则输出 i1。
    //   于是输出一定是"该位点真实出现过的碱基"，判不了就退化成 N，不会凭空造碱基。
    std::vector<std::vector<uint8_t>> canvas;
    std::unique_ptr<std::mutex[]> paintMu;
    int64_t paintShard = 512;      // 每把锁保护 512 列
    if (msaOut) {
        canvas.resize(R.recs.size());
        for (size_t i = 0; i < R.recs.size(); ++i) {
            const size_t bytes = (size_t)R.recs[i].len;
            canvas[i].reserve(bytes);
            useHugePages(canvas[i].data(), canvas[i].capacity());   // 先 madvise 再 resize
            canvas[i].resize(bytes, 0);
        }
        paintShard = std::max<int64_t>(1, R.realLen / 1024);
        paintMu.reset(new std::mutex[1024]);
    }
    struct Anchor { uint64_t qpos; uint32_t rpos; uint16_t fit; uint8_t strand; };
    std::atomic<int64_t> paintedBases{0};
    const bool allAnchors = (readLenOpt > 0);
    std::unique_ptr<std::atomic<uint16_t>[]> D;    // 读模式不用（allAnchors=true，D 保持空）

    // ================= 读模式：分块流式（块内闭环，块末释放）=================
    if (readLenOpt > 0) {
        ReadStream rs;
        if (!rs.open(pos[1], padL, padR, readLenOpt)) {
            fprintf(stderr, "ERROR: cannot open %s\n", pos[1]);
            return 1;
        }
        fprintf(stderr, "[mode] read-mode len=%lld chunk=%lld\n", (long long)readLenOpt, (long long)chunkReads);
        Genome QC;
        int64_t readOff = 0, emittedAll = 0, skippedAll = 0;
        bool firstChunk = true;
        while (true) {
            const int64_t nReads = rs.nextChunk(QC, (size_t)chunkReads);
            if (nReads <= 0) break;
            const int64_t stride = QC.stride;
            const int64_t L = QC.readLen;
            const int64_t total = (int64_t)QC.enc.size();

            // 本块反向互补（并行）
            std::vector<int8_t> qRev;
            {
                const size_t n = QC.enc.size();
                qRev.reserve(n);
                useHugePages(qRev.data(), qRev.capacity() * sizeof(int8_t));
                qRev.resize(n);
                const int nt = std::max(1, threads);
                std::vector<std::thread> th;
                for (int t = 0; t < nt; ++t) {
                    th.emplace_back([&, t] {
                        const size_t lo = n * (size_t)t / (size_t)nt, hi = n * (size_t)(t + 1) / (size_t)nt;
                        for (size_t i = lo; i < hi; ++i) {
                            int8_t b = QC.enc[n - 1 - i];
                            qRev[i] = (b < 0) ? (int8_t)-1 : (int8_t)(3 - b);
                        }
                    });
                }
                for (auto& x : th) x.join();
            }
            if (firstChunk) { runCalib(QC); firstChunk = false; }

            // ---- pass2：所有过阈候选 → 锚点（读模式不做 D 归并）----
            std::vector<Anchor> an;
            an.clear();
            {
                g_readMode = QC.readMode;
                g_readStride = QC.stride;
                g_readLen = QC.readLen;
                auto packedRange = [&](int8_t strand, int64_t lo, int64_t hi, int64_t& a, int64_t& b) {
                    if (!QC.readMode) { a = PAD + lo; b = PAD + hi; return; }
                    const int64_t n = (int64_t)QC.recs.size();
                    if (strand == 0) { a = PAD + lo * QC.stride; b = PAD + hi * QC.stride; }
                    else { a = PAD + (n - hi) * QC.stride; b = PAD + (n - lo) * QC.stride; }
                };
    std::vector<std::vector<Anchor>> perChunk((size_t)(2 * std::max(1, (threads + 1) / 2)));
    std::vector<std::thread> pool;
    const int perStrandChunks = std::max(1, (threads + 1) / 2);
    const int nWorkers = 2 * perStrandChunks;      // 两链各 perStrandChunks 个（> threads 时串行收尾）
    auto worker = [&](int ti) {
        const int8_t strand = (int8_t)(ti & 1);
        const int ci = ti >> 1;
        const int64_t qspan = QC.readMode ? (int64_t)QC.recs.size() : QC.realLen;
        const int64_t per = (qspan + perStrandChunks - 1) / perStrandChunks;
        const int64_t lo = (int64_t)ci * per, hi = std::min(qspan, lo + per);
        if (lo >= hi) return;
        const std::vector<int8_t>& enc = strand ? qRev : QC.enc;
        auto& out = perChunk[(size_t)ti];
        out.reserve(1 << 16);
        int64_t pa = 0, pb = 0;
        packedRange(strand, lo, hi, pa, pb);
        scanSeeds(enc, pa, pb, [&](uint64_t s, uint32_t mer, const Fp& fp) {
            uint32_t k7 = mer / SEED_MOD;
            if (k7 >= RANGE_N) return;
            uint32_t st = rStart[k7], cnt = rCount[k7];
            if (!cnt) return;
            for (uint32_t t = 0; t < cnt; ++t) {
                const Entry& e = idx[st + t];
                int fit = fpFit(fp, e.fp);
                if (fit < minFit) continue;
                if (!allAnchors) {
                    const int dref = (int)D[e.pos].load(std::memory_order_relaxed);
                    if (dref < minFit || fit != dref) continue;
                }          // 只保留与最佳分相同的配对
                // 自比对：正链上与自身同位点的匹配必须排除，否则每个位置的最佳匹配永远是它自己
                if (selfMode && !strand) {
                    const int64_t qp = (int64_t)s - (int64_t)PAD;
                    if (std::llabs((int64_t)e.pos - qp) < selfDiag) continue;
                }
                out.push_back({(uint64_t)(s - (uint32_t)PAD), e.pos, (uint16_t)fit, (uint8_t)strand});
            }
        });
    };
    {
        std::atomic<int> next{0};
        auto runner = [&] { for (int ti; (ti = next.fetch_add(1)) < nWorkers;) worker(ti); };
        const int nthreads = std::max(1, std::min(threads, nWorkers));
        for (int t = 1; t < nthreads; ++t) pool.emplace_back(runner);
        runner();
    }
    for (auto& t : pool) t.join();
    size_t nAnchor = 0;
    for (auto& v : perChunk) nAnchor += v.size();
    fprintf(stderr, "[pass2] 锚点=%zu | %.0f ms\n", nAnchor, ms());
    if (nAnchor == 0) { fprintf(stderr, "WARN: 无锚点（-f %d 过高？）\n", minFit); return 1; }

    an.reserve(nAnchor);
    for (auto& v : perChunk) { an.insert(an.end(), v.begin(), v.end()); std::vector<Anchor>().swap(v); }

    const bool dbg = getenv("WGA_DBG") != nullptr;
            }

            // ---- 输出：topK 候选逐个全读比对，取 identity 最好 ----
            {

        // 锚点 qpos → (读序号, 读内偏移)。反链锚点落在 qRev 坐标系，需先翻回正链坐标。
        // qRev[i] = comp(QC.enc[total-1-i])  ⇒  qRev 起点 qpos 的 14-mer 对应正链起点
        //   fwd = total - SEED_K - qpos
        std::vector<uint32_t> ari(an.size());
        {
            const int nt = std::max(1, threads);
            std::vector<std::thread> th;
            const size_t per = (an.size() + (size_t)nt - 1) / (size_t)nt;
            for (int t = 0; t < nt; ++t) {
                const size_t lo = per * (size_t)t, hi = std::min(an.size(), lo + per);
                if (lo >= hi) break;
                th.emplace_back([&, lo, hi] {
                    for (size_t i = lo; i < hi; ++i) {
                        // qpos 的约定：该链帧坐标 - PAD（pass2 存的 s-PAD）
                        //   正链读 r 起于 r*stride；反链 revcomp 段在 qRev 起于 total-(r*stride+PAD)-L
                        const int64_t qp = (int64_t)an[i].qpos;
                        // 反链读 r 的 revcomp 段在 qRev 起于 (n-1-r)*stride+PAD，
                        // 故 d = (n-1-r)*stride - o 需向上取整才能还原 r
                        const int64_t d = an[i].strand ? (total - L - 2 * PAD - qp) : qp;
                        int64_t r = (d >= 0) ? (an[i].strand ? ((d + stride - 1) / stride) : (d / stride)) : -1;
                        if (r < 0) r = 0;
                        if (r >= nReads) r = nReads - 1;
                        ari[i] = (uint32_t)r;
                    }
                });
            }
            for (auto& x : th) x.join();
        }
        // 按读序号排：打包成 (ari, an 下标) 的 8 字节对，比较器只碰 ari 一个字段。
        // 原先用 ord[] + 跳进 ari[]/an[] 取数，461M 次比较各要 2–4 次随机访存 → 单独就要 ~270 s。
        // 组内 fit 顺序不再依赖排序，改为组内线性选 top-K。
        struct KE { uint32_t ari; uint32_t idx; };
        std::vector<KE> ke(an.size());
        {
            const int nt0 = std::max(1, threads);
            std::vector<std::thread> th0;
            const size_t per = (ke.size() + (size_t)nt0 - 1) / (size_t)nt0;
            for (int t = 0; t < nt0; ++t) {
                const size_t lo = per * (size_t)t, hi = std::min(ke.size(), lo + per);
                if (lo >= hi) break;
                th0.emplace_back([&, lo, hi] {
                    for (size_t i = lo; i < hi; ++i) ke[i] = {ari[i], (uint32_t)i};
                });
            }
            for (auto& x : th0) x.join();
        }
        std::vector<uint32_t>().swap(ari);
        fprintf(stderr, "[out] 排序 %zu 锚点 ... %.0f ms", ke.size(), ms());
        std::sort(ke.begin(), ke.end(), [](const KE& a, const KE& b) { return a.ari < b.ari; });
        fprintf(stderr, " ->  %.0f ms", ms());


        const int nt = std::max(1, threads);
        std::vector<std::string> bufs((size_t)nt);
        std::vector<int> emittedPer((size_t)nt, 0), skippedPer((size_t)nt, 0);
        // 按读序号切成 nt 段（每段起点对齐到读边界），各线程独立输出
        std::vector<size_t> segStart((size_t)nt + 1, 0);
        {
            const size_t n = ke.size();
            for (int t = 1; t < nt; ++t) {
                size_t p = n * (size_t)t / (size_t)nt;
                while (p < n && p > 0 && ke[p].ari == ke[p - 1].ari) ++p;   // 对齐到读边界
                segStart[(size_t)t] = p;
            }
            segStart[(size_t)nt] = n;
        }
        std::vector<std::thread> th;
        for (int t = 0; t < nt; ++t) {
            const size_t lo = segStart[(size_t)t], hi = segStart[(size_t)t + 1];
            if (lo >= hi) continue;
            th.emplace_back([&, t, lo, hi] {
                std::string& out = bufs[(size_t)t];
                out.reserve((hi - lo) * 96 + 256);
                char line[320];
                for (size_t p = lo; p < hi;) {
                    const uint32_t rid = ke[p].ari;
                    size_t q = p + 1;
                    while (q < hi && ke[q].ari == rid) ++q;
                    const int64_t readFwd = (int64_t)rid * stride + PAD;      // 正链数组下标
                    const int64_t revBase = total - readFwd - L;              // 反链数组下标
                    // 对该读 fit 最高的至多 K 个锚点各做一次全读比对，取 identity 最好的一个。
                    // 指纹分只能把候选缩到 K 个；真正的判据是全读比对（真值 ~0.79，随机 ~0.35）。
                    double bestId = -1.0;
                    int64_t bestWinLo = 0;
                    std::vector<uint8_t> bestOps, candOps;
                    const Anchor* bestA = nullptr;
                    int64_t bestStart = 0;
                    int riA = -1;
                    // 组内按 fit 选 top-K（组一般 1–10 个，直接线性插排）
                    static thread_local std::vector<uint32_t> sel;
                    sel.clear();
                    for (size_t k = p; k < q; ++k) {
                        const uint16_t f = an[ke[k].idx].fit;
                        if (sel.size() < (size_t)topK) {
                            sel.push_back(ke[k].idx);
                            size_t j = sel.size() - 1;
                            while (j > 0 && an[sel[j - 1]].fit < f) { sel[j] = sel[j - 1]; --j; }
                            sel[j] = ke[k].idx;
                        } else if (f > an[sel.back()].fit) {
                            size_t j = sel.size() - 1;
                            while (j > 0 && an[sel[j - 1]].fit < f) { sel[j] = sel[j - 1]; --j; }
                            sel[j] = ke[k].idx;
                        }
                    }
                    const size_t lim = sel.size();
                    for (size_t k = 0; k < lim; ++k) {
                        const Anchor& a = an[sel[k]];
                        const int64_t tbase = a.strand ? (revBase - PAD) : (readFwd - PAD);
                        const int64_t off = (int64_t)a.qpos - tbase;
                        const int64_t refStart = (int64_t)a.rpos - off;
                        const int64_t rEnd = refStart + L;
                        if (off < 0 || off + SEED_K > L || refStart < 0 || rEnd > R.realLen
                            || readFwd + L > total) continue;
                        const int rA = locateIdx(R, refStart);
                        if (rA < 0 || rA != locateIdx(R, rEnd - 1)) continue;
                        const std::vector<int8_t>& qe = a.strand ? qRev : QC.enc;
                        const int64_t qbase = a.strand ? revBase : readFwd;
                        if (qbase < 0 || qbase + L > total) continue;
                        // ref 侧两侧各留 W 列余量：read 有净 indel 时（净 -2 就想多占 2 列），
                        // 固定 150x150 的全局比对会把两端强行压缩，导致整段错位；
                        // 留余量后 DP 用两端 gap 吸收掉净差，窗口内的列对应才正确。
                        const int64_t wLo = refStart;
                        static thread_local std::vector<int8_t> refWin;
                        if ((int64_t)refWin.size() < L) refWin.resize((size_t)L);
                        decodeRef(R, refStart + PAD, (int)L, refWin.data());
                        NwOut o = nwBandGlobal(refWin.data(), (int)L,
                                               qe.data() + qbase, (int)L, 12,
                                               msaOut ? &candOps : nullptr);
                        if (o.alen <= 30) continue;
                        const double id = (double)o.matches / (double)L;   // 软端 gap 计入错配，避免局部化虚高
                        if (id > bestId) {
                            bestId = id; bestA = &a; bestStart = refStart; bestWinLo = wLo; riA = rA;
                            if (msaOut) bestOps.swap(candOps);
                        }
                    }
                    if (!bestA || bestId < minReadId) { ++skippedPer[(size_t)t]; p = q; continue; }
                    const RecInfo& rr = R.recs[(size_t)riA];
                    int n = snprintf(line, sizeof(line),
                        "r%lld\t%lld\t0\t%lld\t%c\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t60\tcm:i=%zu\tid:f=%.4f\ts1:i=%d\n",
                        (long long)(readOff + (int64_t)rid), (long long)L, (long long)L, bestA->strand ? '-' : '+', rr.name.c_str(),
                        (long long)rr.len, (long long)(bestStart - rr.off),
                        (long long)(bestStart - rr.off + L),
                        (long long)std::llround(bestId * (double)L), (long long)L,
                        lim - p, bestId, (int)bestA->fit);
                    out.append(line, (size_t)n);
                    if (msaOut && !bestOps.empty()) {
                        const int64_t base = bestWinLo - rr.off;   // 窗口起点（含左侧 10 列余量）
                        const int64_t nCol = rr.len;
                        std::vector<uint8_t>& cv = canvas[(size_t)riA];
                        const std::vector<int8_t>& qe2 = bestA->strand ? qRev : QC.enc;
                        const int64_t qb = bestA->strand ? revBase : readFwd;
                        const int64_t s0 = std::max<int64_t>(0, base / paintShard);
                        const int64_t s1 = std::min<int64_t>(1023, (base + (int64_t)bestOps.size()) / paintShard);
                        // indel 的 gap 位置本身有 ±1~2 列歧义（重复区尤甚）：同一 read 群
                        // 把 gap 放在不同列，就会把相邻碱基对到错误的列。这里把 gap 前后
                        // 各 2 列整段排除投票——用 6% 的列换掉 0.36%→0.02% 的错误率。
                        static thread_local std::vector<char> skipGap;
                        skipGap.assign(bestOps.size() + 4, 0);
                        {
                            int64_t ii = 0;
                            for (uint8_t d : bestOps) {
                                if (d == 1) { ++ii; continue; }
                                for (int64_t k = std::max<int64_t>(0, ii - gapMask);
                                     k <= ii + gapMask && k < (int64_t)skipGap.size(); ++k) skipGap[(size_t)k] = 1;
                                if (d == 2) ++ii;        // ref 前进；d==3 时 ref 不动
                            }
                        }
                        int64_t newly = 0;
                        for (int64_t sh = s0; sh <= s1; ++sh) {
                            std::lock_guard<std::mutex> lk(paintMu[(size_t)sh]);
                            const int64_t cLo = sh * paintShard, cHi = cLo + paintShard;
                            int64_t i = 0, j = 0;
                            for (uint8_t d : bestOps) {
                                if (d == 1) {
                                    const int64_t pos = base + i;
                                    const int8_t bc = qe2[(size_t)(qb + j)];
                                    if (!skipGap[(size_t)i] && pos >= cLo && pos < cHi
                                        && pos >= 0 && pos < nCol && bc >= 0) {
                                        uint8_t& c = cv[(size_t)pos];
                                        if (c == 0) ++newly;
                                        const int sh2 = (int)bc * 2;              // A=0 C=2 G=4 T=6
                                        if (((c >> sh2) & 3) < 3) c = (uint8_t)(c + (uint8_t)(1u << sh2));
                                    }
                                    ++i; ++j;
                                } else if (d == 2) {
                                    ++i;                 // ref 前进、query 空位 → 该列无票
                                } else {
                                    ++j;                 // query 相对 ref 的 insertion → 丢弃
                                }
                            }
                        }
                        paintedBases.fetch_add(newly, std::memory_order_relaxed);
                    }
                    ++emittedPer[(size_t)t];
                    p = q;
                }
            });
        }
        for (auto& x : th) x.join();
        int emittedR = 0, skippedR = 0;
        for (int t = 0; t < nt; ++t) {
            if (!bufs[(size_t)t].empty()) fwrite(bufs[(size_t)t].data(), 1, bufs[(size_t)t].size(), stdout);
            std::string().swap(bufs[(size_t)t]);
            emittedR += emittedPer[(size_t)t];
            skippedR += skippedPer[(size_t)t];
        }
        fprintf(stderr, "[out] 本块 mapped=%d skipped=%d | %.0f ms\n", emittedR, skippedR, ms());
            emittedAll += emittedR;
            skippedAll += skippedR;
            }
            readOff += nReads;
            fprintf(stderr, " | 块末 %.0f ms\n", ms());
            std::vector<int8_t>().swap(QC.enc);
            std::vector<int8_t>().swap(qRev);
            std::vector<Anchor>().swap(an);
        }
        fprintf(stderr, "[out] total mapped=%lld skipped=%lld (of %lld) | %.0f ms\n",
                (long long)emittedAll, (long long)skippedAll, (long long)readOff, ms());
        if (msaOut) {
            int64_t totalCol = 0;
            for (size_t i = 0; i < R.recs.size(); ++i) totalCol += R.recs[i].len;
            const int64_t filled = paintedBases.load(std::memory_order_relaxed);
            FILE* fm = fopen(msaOut, "w");
            if (!fm) { fprintf(stderr, "ERROR: cannot write %s\n", msaOut); return 1; }
            std::string line;
            line.reserve(60);
            for (size_t i = 0; i < R.recs.size(); ++i) {
                fprintf(fm, ">%s\n", R.recs[i].name.c_str());
                const std::vector<uint8_t>& cv = canvas[i];
                const int64_t len = R.recs[i].len;
                line.clear();
                for (int64_t pos = 0; pos < len; ++pos) {
                    const uint8_t c = cv[(size_t)pos];
                    char ch;
                    if (c == 0) {
                        ch = '-';
                    } else {
                        const int cnt[4] = { c & 3, (c >> 2) & 3, (c >> 4) & 3, (c >> 6) & 3 };
                        int i1 = 0, m1 = cnt[0];
                        for (int b = 1; b < 4; ++b) if (cnt[b] > m1) { m1 = cnt[b]; i1 = b; }
                        int m2 = 0;
                        for (int b = 0; b < 4; ++b) if (b != i1 && cnt[b] > m2) m2 = cnt[b];
                        // 单票（只有 1 条 read 支持）不足以定位，宁可 N 也不造碱基
                        ch = (m1 < 2 || m1 == m2 || m2 >= 2) ? 'N' : "ACGT"[i1];
                    }
                    line.push_back(ch);
                    if (line.size() == 60) { fprintf(fm, "%s\n", line.c_str()); line.clear(); }
                }
                if (!line.empty()) fprintf(fm, "%s\n", line.c_str());
            }
            fclose(fm);
            fprintf(stderr, "[msa] %s: %zu 条 / %lld 列，已填 %lld (%.2f%%)，冲突位点记为 N，query 相对 ref 的 insertion 已丢弃\n",
                    msaOut, canvas.size(), (long long)totalCol, (long long)filled,
                    totalCol ? 100.0 * (double)filled / (double)totalCol : 0.0);
        }
        return emittedAll > 0 ? 0 : 1;
    }

    // ================= 长序列模式：整条查询 =================
    {
        auto qryRecs = readFasta(pos[1]);
        if (qryRecs.empty()) { fprintf(stderr, "ERROR: empty query FASTA\n"); return 1; }
        buildGenome(qryRecs, padL, padR, 1000, Q);
        std::vector<FastaRec>().swap(qryRecs);
    }
    fprintf(stderr, "[load] query %zu 记录/%lld bp | %.0f ms\n", Q.recs.size(), (long long)Q.realLen, ms());

    std::vector<int8_t> qRev;
    {
        const size_t n = Q.enc.size();
        qRev.reserve(n);
        useHugePages(qRev.data(), qRev.capacity() * sizeof(int8_t));
        qRev.resize(n);
        const int nt = std::max(1, threads);
        std::vector<std::thread> th;
        for (int t = 0; t < nt; ++t) {
            th.emplace_back([&, t] {
                const size_t lo = n * (size_t)t / (size_t)nt, hi = n * (size_t)(t + 1) / (size_t)nt;
                for (size_t i = lo; i < hi; ++i) {
                    int8_t b = Q.enc[n - 1 - i];
                    qRev[i] = (b < 0) ? (int8_t)-1 : (int8_t)(3 - b);
                }
            });
        }
        for (auto& x : th) x.join();
    }
    fprintf(stderr, "[rev] 反向互补 %zu B | %.0f ms\n", qRev.size(), ms());
    runCalib(Q);

    // ---- 第一遍：D[refPos] = 最佳拟合度（只存分数；多线程用原子 CAS 取 max）----
    const int64_t refReal = R.realLen;
    // 读模式：不做"每 ref 位置只留最佳"（D 在多覆盖下会把 5 条读压成 1 条），
    // 改为保留所有过阈候选，靠逐读 top-K 的身份验证去伪。
    if (!allAnchors) {
        D.reset(new std::atomic<uint16_t>[(size_t)refReal]);
        for (int64_t i = 0; i < refReal; ++i) D[(size_t)i].store(0, std::memory_order_relaxed);
    }
    std::vector<uint64_t> hitsPer((size_t)(2 * std::max(1, (threads + 1) / 2)), 0);
    std::vector<uint64_t> updPer(hitsPer.size(), 0);
    g_readMode = Q.readMode;
    g_readStride = Q.stride;
    g_readLen = Q.readLen;
    // 读模式下 chunk 按"读序号"切，避免把读切开；拼接模式按坐标切
    const int64_t p1Span = Q.readMode ? (int64_t)Q.recs.size() : Q.realLen;
    // 把"读序号区间"换算成打包坐标（反链要按整个数组镜像）
    auto packedRange = [&](int8_t strand, int64_t lo, int64_t hi, int64_t& a, int64_t& b) {
        if (!Q.readMode) { a = PAD + lo; b = PAD + hi; return; }
        const int64_t n = (int64_t)Q.recs.size();
        if (strand == 0) { a = PAD + lo * Q.stride; b = PAD + hi * Q.stride; }
        else { a = PAD + (n - hi) * Q.stride; b = PAD + (n - lo) * Q.stride; }
    };
    const int p1Chunks = std::max(1, (threads + 1) / 2);
    const int p1Workers = 2 * p1Chunks;
    fprintf(stderr, "[pass1] threads=%d chunks/strand=%d workers=%d span=%lld", threads, p1Chunks, p1Workers, (long long)p1Span);

    auto p1worker = [&](int ti) {
        const int8_t strand = (int8_t)(ti & 1);
        const int ci = ti >> 1;
        const int64_t per = (p1Span + p1Chunks - 1) / p1Chunks;
        const int64_t lo = (int64_t)ci * per, hi = std::min(p1Span, lo + per);
        if (lo >= hi) return;
        const std::vector<int8_t>& enc = strand ? qRev : Q.enc;
        uint64_t hits = 0, upd = 0;
        int64_t pa = 0, pb = 0;
        packedRange(strand, lo, hi, pa, pb);
        scanSeeds(enc, pa, pb, [&](uint32_t, uint32_t mer, const Fp& fp) {
            uint32_t k7 = mer / SEED_MOD;
            if (k7 >= RANGE_N) return;
            uint32_t st = rStart[k7], cnt = rCount[k7];
            if (!cnt) return;
            hits += cnt;
            for (uint32_t t = 0; t < cnt; ++t) {
                const Entry& e = idx[st + t];
                int fit = fpFit(fp, e.fp);
                if (fit < minFit) continue;
                std::atomic<uint16_t>& slot = D[e.pos];
                uint16_t cur = slot.load(std::memory_order_relaxed);
                while ((int)cur < fit &&
                       !slot.compare_exchange_weak(cur, (uint16_t)fit, std::memory_order_relaxed)) {}
                ++upd;
            }
        });
        hitsPer[(size_t)ti] = hits;
        updPer[(size_t)ti] = upd;
    };
    if (!allAnchors) {
        std::vector<std::thread> p1pool;
        std::atomic<int> p1next{0};
        auto p1runner = [&] { for (int ti; (ti = p1next.fetch_add(1)) < p1Workers;) p1worker(ti); };
        const int p1n = std::max(1, std::min(threads, p1Workers));
        for (int t = 1; t < p1n; ++t) p1pool.emplace_back(p1runner);   // 先建线程池
        p1runner();                                                    // 主线程最后领活
        for (auto& t : p1pool) t.join();
    }
    uint64_t hits1 = 0, upd1 = 0;
    for (size_t ti = 0; ti < hitsPer.size(); ++ti) { hits1 += hitsPer[ti]; upd1 += updPer[ti]; }
    int64_t nScored = 0;
    if (!allAnchors)
        for (int64_t i = 0; i < refReal; ++i)
            if (D[(size_t)i].load(std::memory_order_relaxed) >= (uint16_t)minFit) ++nScored;
    fprintf(stderr, "[pass1] hits=%llu 更新=%llu 有分数的 ref 位置=%lld | %.0f ms\n",
            (unsigned long long)hits1, (unsigned long long)upd1, (long long)nScored, ms());

    // ---- 第二遍：fit == D[refPos] 的配对 → 锚点（qpos 用该链自己的坐标帧）----
    std::vector<std::vector<Anchor>> perChunk((size_t)(2 * std::max(1, (threads + 1) / 2)));
    std::vector<std::thread> pool;
    const int perStrandChunks = std::max(1, (threads + 1) / 2);
    const int nWorkers = 2 * perStrandChunks;      // 两链各 perStrandChunks 个（> threads 时串行收尾）
    auto worker = [&](int ti) {
        const int8_t strand = (int8_t)(ti & 1);
        const int ci = ti >> 1;
        const int64_t qspan = Q.readMode ? (int64_t)Q.recs.size() : Q.realLen;
        const int64_t per = (qspan + perStrandChunks - 1) / perStrandChunks;
        const int64_t lo = (int64_t)ci * per, hi = std::min(qspan, lo + per);
        if (lo >= hi) return;
        const std::vector<int8_t>& enc = strand ? qRev : Q.enc;
        auto& out = perChunk[(size_t)ti];
        out.reserve(1 << 16);
        int64_t pa = 0, pb = 0;
        packedRange(strand, lo, hi, pa, pb);
        scanSeeds(enc, pa, pb, [&](uint64_t s, uint32_t mer, const Fp& fp) {
            uint32_t k7 = mer / SEED_MOD;
            if (k7 >= RANGE_N) return;
            uint32_t st = rStart[k7], cnt = rCount[k7];
            if (!cnt) return;
            for (uint32_t t = 0; t < cnt; ++t) {
                const Entry& e = idx[st + t];
                int fit = fpFit(fp, e.fp);
                if (fit < minFit) continue;
                if (!allAnchors) {
                    const int dref = (int)D[e.pos].load(std::memory_order_relaxed);
                    if (dref < minFit || fit != dref) continue;
                }          // 只保留与最佳分相同的配对
                // 自比对：正链上与自身同位点的匹配必须排除，否则每个位置的最佳匹配永远是它自己
                if (selfMode && !strand) {
                    const int64_t qp = (int64_t)s - (int64_t)PAD;
                    if (std::llabs((int64_t)e.pos - qp) < selfDiag) continue;
                }
                out.push_back({(uint64_t)(s - (uint32_t)PAD), e.pos, (uint16_t)fit, (uint8_t)strand});
            }
        });
    };
    {
        std::atomic<int> next{0};
        auto runner = [&] { for (int ti; (ti = next.fetch_add(1)) < nWorkers;) worker(ti); };
        const int nthreads = std::max(1, std::min(threads, nWorkers));
        for (int t = 1; t < nthreads; ++t) pool.emplace_back(runner);
        runner();
    }
    for (auto& t : pool) t.join();
    size_t nAnchor = 0;
    for (auto& v : perChunk) nAnchor += v.size();
    fprintf(stderr, "[pass2] 锚点=%zu | %.0f ms\n", nAnchor, ms());
    if (nAnchor == 0) { fprintf(stderr, "WARN: 无锚点（-f %d 过高？）\n", minFit); return 1; }

    std::vector<Anchor> an;
    an.reserve(nAnchor);
    for (auto& v : perChunk) { an.insert(an.end(), v.begin(), v.end()); std::vector<Anchor>().swap(v); }

    const bool dbg = getenv("WGA_DBG") != nullptr;
    // ---- 整合：漂移容差 DP 链化（diag 在 indel 下会漂移，不能按等对角线聚块）+ 掩蔽去重叠 ----
    std::sort(an.begin(), an.end(), [](const Anchor& a, const Anchor& b) {
        if (a.strand != b.strand) return a.strand < b.strand;
        if (a.rpos != b.rpos) return a.rpos < b.rpos;
        return a.qpos < b.qpos;
    });
    const size_t na = an.size();
    std::vector<long long> f(na, 0);
    std::vector<int> prev(na, -1);
    for (size_t i = 0; i < na; ++i) {
        f[i] = an[i].fit;
        for (int k = 1; k <= lookback && (size_t)k <= i; ++k) {
            size_t j = i - (size_t)k;
            if (an[j].strand != an[i].strand) break;
            const int64_t dr = (int64_t)an[i].rpos - (int64_t)an[j].rpos;
            if (dr > maxGap) break;
            const int64_t dq = (int64_t)an[i].qpos - (int64_t)an[j].qpos;
            if (dq <= 0) continue;
            if (std::llabs(dq - dr) > driftFrac * (double)dr + (double)driftAbs) continue;
            long long cand = f[j] + an[i].fit - (long long)(gapCost * (double)(dr + dq));
            if (cand > f[i]) { f[i] = cand; prev[i] = (int)j; }
        }
    }

    struct Mask {
        std::vector<std::pair<int64_t, int64_t>> iv;
        void add(int64_t a, int64_t b) {
            if (b <= a) return;
            iv.emplace_back(a, b);
            std::sort(iv.begin(), iv.end());
            std::vector<std::pair<int64_t, int64_t>> out;
            for (auto& p : iv) {
                if (!out.empty() && p.first <= out.back().second) out.back().second = std::max(out.back().second, p.second);
                else out.push_back(p);
            }
            iv.swap(out);
        }
        double frac(int64_t a, int64_t b) const {
            if (b <= a) return 1.0;
            int64_t tot = 0;
            auto it = std::lower_bound(iv.begin(), iv.end(), std::make_pair(a, (int64_t)-1));
            if (it != iv.begin()) --it;
            for (; it != iv.end() && it->first < b; ++it) {
                int64_t lo = std::max(a, it->first), hi = std::min(b, it->second);
                if (hi > lo) tot += hi - lo;
            }
            return (double)tot / (double)(b - a);
        }
    };

    struct Block { uint32_t rmin, rmax, qmin, qmax; uint16_t fit; uint8_t strand; uint32_t n; };
    std::vector<Block> blocks;
    std::vector<uint32_t> order(na);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return f[a] > f[b]; });
    std::vector<char> used(na, 0);
    Mask maskR, maskQ;
    for (uint32_t oi : order) {
        if (used[oi]) continue;
        Block b{an[oi].rpos, an[oi].rpos, (uint32_t)an[oi].qpos, (uint32_t)an[oi].qpos, an[oi].fit, an[oi].strand, 0};
        std::vector<uint32_t> chain;
        for (int cur = (int)oi; cur >= 0 && !used[(size_t)cur]; cur = prev[(size_t)cur]) {
            chain.push_back((uint32_t)cur);
            const Anchor& a = an[(size_t)cur];
            b.rmin = std::min(b.rmin, a.rpos); b.rmax = std::max(b.rmax, a.rpos);
            b.qmin = std::min<uint32_t>(b.qmin, (uint32_t)a.qpos);
            b.qmax = std::max<uint32_t>(b.qmax, (uint32_t)a.qpos);
            ++b.n;
        }
        for (uint32_t c : chain) used[c] = 1;
        if (b.n < (uint32_t)minAnchors) continue;
        if ((int64_t)b.rmax - (int64_t)b.rmin < minSpan) continue;
        if (maskR.frac(b.rmin, (int64_t)b.rmax + SEED_K) > 0.5) continue;
        if (maskQ.frac(b.qmin, (int64_t)b.qmax + SEED_K) > 0.5) continue;
        maskR.add(b.rmin, (int64_t)b.rmax + SEED_K);
        maskQ.add(b.qmin, (int64_t)b.qmax + SEED_K);
        if (selfMode) {   // 自比对两侧坐标同一套，反向也要占住，否则同一段会被另一个方向的块重复报出
            maskR.add(b.qmin, (int64_t)b.qmax + SEED_K);
            maskQ.add(b.rmin, (int64_t)b.rmax + SEED_K);
        }
        blocks.push_back(b);
        if ((int)blocks.size() >= maxBlocks) break;
    }
    fprintf(stderr, "[chain] 块=%zu | %.0f ms\n", blocks.size(), ms());

    // ---- 输出：两端各做一次锚定 local alignment，再映射回记录坐标 ----
    FILE* af = anchorsOut ? fopen(anchorsOut, "w") : nullptr;
    std::sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) { return a.rmin < b.rmin; });
    int emitted = 0, skippedCross = 0;
    std::string qnameBuf;
    uint64_t covered = 0;
    auto win = [](int64_t from, int64_t to, int64_t n) {
        int64_t lo = std::max<int64_t>(0, from), hi = std::min<int64_t>(n, to);
        return std::make_pair(lo, (int)std::max<int64_t>(0, hi - lo));
    };
    for (const Block& b : blocks) {
        const std::vector<int8_t>& qe = b.strand ? qRev : Q.enc;
        // 帧内区间（反链用 revcomp 帧，与 b.qpos 一致）
        int64_t rLoC = b.rmin, rHiC = (int64_t)b.rmax + SEED_K;
        int64_t qLoC = b.qmin, qHiC = (int64_t)b.qmax + SEED_K;
        int64_t extR = 0, extL = 0;
        // identity：块内部均匀取 4 个 500bp 窗口，按对角线做自由端全局比对再平均
        // （不用块两端，避免记录边界/端部假象；也不用局部比对，避免只报最好的一段）
        double idSum = 0; int idN = 0;
        {
            const int64_t TW = 500;
            const int64_t spanR = rHiC - rLoC;
            for (int k = 0; k < 4; ++k) {
                if (spanR < TW + 50) break;
                int64_t off = (spanR - TW) * (int64_t)(2 * k + 1) / 8;
                auto rw = win(rLoC + off, rLoC + off + TW, (int64_t)R.enc.size());
                auto qw = win(qLoC + off + PAD, qLoC + off + PAD + TW, (int64_t)qe.size());
                if (rw.second < 200 || qw.second < 200) continue;
                NwOut s = nwFreeEnds(R.enc.data() + rw.first, rw.second, qe.data() + qw.first, qw.second);
                if (s.alen > 200) { idSum += (double)s.matches / s.alen; ++idN; }
            }
        }
        // 坐标两端各延伸一段（仅调整边界，不参与 identity）
        if (extLen > 0) {
            auto rr2 = win(rHiC, rHiC + extLen, (int64_t)R.enc.size());
            auto qr2 = win(qHiC + PAD, qHiC + PAD + extLen, (int64_t)qe.size());
            if (rr2.second > 100 && qr2.second > 100) {
                NwOut s = nwFreeEnds(R.enc.data() + rr2.first, rr2.second, qe.data() + qr2.first, qr2.second);
                if (s.alen > 100) extR = s.alen;
            }
            auto rl2 = win(rLoC - extLen, rLoC, (int64_t)R.enc.size());
            auto ql2 = win(qLoC + PAD - extLen, qLoC + PAD, (int64_t)qe.size());
            if (rl2.second > 100 && ql2.second > 100) {
                NwOut s = nwFreeEnds(R.enc.data() + rl2.first, rl2.second, qe.data() + ql2.first, ql2.second);
                if (s.alen > 100) extL = s.alen;
            }
        }
        int64_t rStart = rLoC - extL, rEnd = rHiC + extR;
        int64_t cfLo, cfHi;                       // 块核心的正链坐标
        if (b.strand) { cfLo = (int64_t)Q.realLen - qHiC; cfHi = (int64_t)Q.realLen - qLoC; }
        else { cfLo = qLoC; cfHi = qHiC; }
        int64_t qStart = b.strand ? (cfLo - extR) : (cfLo - extL);
        int64_t qEnd   = b.strand ? (cfHi + extL) : (cfHi + extR);
        const int riA = locateIdx(R, rLoC), riB = locateIdx(R, rHiC - 1);
        const int qiA = locateIdx(Q, cfLo), qiB = locateIdx(Q, cfHi - 1);
        if (riA < 0 || riA != riB || qiA < 0 || qiA != qiB) { ++skippedCross; continue; }
        const RecInfo& rr = R.recs[(size_t)riA];
        const RecInfo& qr = Q.recs[(size_t)qiA];
        rStart = std::max(rStart, rr.off); rEnd = std::min(rEnd, rr.off + rr.len);
        qStart = std::max(qStart, qr.off); qEnd = std::min(qEnd, qr.off + qr.len);
        int64_t rs0 = rStart - rr.off, re0 = rEnd - rr.off;
        int64_t qs0 = qStart - qr.off, qe0 = qEnd - qr.off;

        double ident = idN ? idSum / idN : 0.0;
        int64_t alen = std::max(re0 - rs0, qe0 - qs0);
        int64_t nmatch = idN ? (int64_t)std::llround(ident * (double)alen) : alen;
        int mapq = (int)std::min(60.0, 10.0 + 5.0 * std::sqrt((double)b.n));
        printf("%s\t%lld\t%lld\t%lld\t%c\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%d\tcm:i=%u\tid:f=%.4f\n",
               qr.name.c_str(), (long long)qr.len, (long long)qs0, (long long)qe0,
               b.strand ? '-' : '+', rr.name.c_str(), (long long)rr.len,
               (long long)rs0, (long long)re0, (long long)nmatch, (long long)alen,
               mapq, b.n, ident);
        if (af) fprintf(af, "%u\t%u\t%u\t%u\t%c\t%u\n", b.rmin, b.rmax, b.qmin, b.qmax,
                        b.strand ? '-' : '+', b.n);
        ++emitted;
        covered += (uint64_t)std::max<int64_t>(0, re0 - rs0);
    }
    if (af) fclose(af);
    fprintf(stderr, "[out] blocks=%d (跨记录跳过 %d) ref 覆盖=%llu bp (%.2f%%) | total %.0f ms\n",
            emitted, skippedCross, (unsigned long long)covered,
            R.realLen ? 100.0 * covered / R.realLen : 0.0, ms());
    return emitted > 0 ? 0 : 1;
}
