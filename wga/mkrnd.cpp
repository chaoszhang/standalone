// mkrnd.cpp — 流式生成随机 FASTA（不占内存）
// 用法: mkrnd <bp> <seed> [lineLen=60]
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <bp> <seed> [lineLen]\n", argv[0]); return 1; }
    const long long bp = atoll(argv[1]);
    std::mt19937_64 rng((unsigned long long)atoll(argv[2]));
    const int lineLen = argc > 3 ? atoi(argv[3]) : 60;
    const char* B = "ACGT";
    fputs(">ref\n", stdout);
    std::string line;
    line.reserve((size_t)lineLen);
    for (long long i = 0; i < bp; ++i) {
        line.push_back(B[rng() & 3]);
        if ((int)line.size() == lineLen) { fputs(line.c_str(), stdout); fputc('\n', stdout); line.clear(); }
    }
    if (!line.empty()) { fputs(line.c_str(), stdout); fputc('\n', stdout); }
    return 0;
}
