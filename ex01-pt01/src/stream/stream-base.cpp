#include <chrono>

#include "../util.h"
#include "stream-util.h"

inline void stream(size_t nx, const double *__restrict__ src, double *__restrict__ dest) {
    for (int i = 0; i < nx; ++i)
        dest[i] = src[i] + 1;
}


int main(int argc, char *argv[]) {
    size_t nx, nItWarmUp, nIt;
    parseCLA_1d(argc, argv, nx, nItWarmUp, nIt);

    auto src = new double[nx];
    auto dest = new double[nx];

    // init
    initStream(src, nx);

    // warm-up
    for (int i = 0; i < nItWarmUp; ++i) {
        stream(nx, src, dest);
        std::swap(src, dest);
    }

    // measurement
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < nIt; ++i) {
        stream(nx, src, dest);
        std::swap(src, dest);
    }

    auto end = std::chrono::steady_clock::now();

    printStats(end - start, nx, nIt, streamNumReads, streamNumWrites);

    // check solution
    checkSolutionStream(src, nx, nIt + nItWarmUp);

    delete[] src;
    delete[] dest;

    return 0;
}
