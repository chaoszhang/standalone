// NJ —— 经典邻居连接（Neighbor-Joining），朴素扫描版
// 编译: g++ -O2 -std=c++20 -o NJ NJ.cpp
#include "limbNJ.hpp"

int main(int argc, char** argv) {
    return limbNJ::run(argc, argv, limbNJ::Mode::NJ, "NJ");
}
