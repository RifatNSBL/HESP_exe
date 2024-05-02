#pragma once

#include <cstdlib>
#include <iostream>
#include <chrono>

void parseCLA_1d(int argc, char *const *argv, size_t &nx, size_t &nItWarmUp, size_t &nIt) {
    // default values
    nx = 1024 * 1024;
    nItWarmUp = 2;
    nIt = 10;

    // override with command line arguments
    int i = 1;
    if (argc > i) nx = atoi(argv[i]);
    ++i;
    if (argc > i) nItWarmUp = atoi(argv[i]);
    ++i;
    if (argc > i) nIt = atoi(argv[i]);
    ++i;
}


void printStats(const std::chrono::duration<double> elapsedSeconds, size_t nCells, size_t nIt, unsigned int numReads, unsigned int numWrites) {
    std::cout << "  #cells / #it:  " << nCells << " / " << nIt << "\n";
    std::cout << "  elapsed time:  " << 1e3 * elapsedSeconds.count() << " ms\n";
    std::cout << "  per iteration: " << 1e3 * elapsedSeconds.count() / nIt << " ms\n";
    std::cout << "  MLUP/s:        " << 1e-6 * nCells * nIt / elapsedSeconds.count() << "\n";
    std::cout << "  bandwidth:     " << 1e-9 * (numReads + numWrites) * sizeof(double) * nCells * nIt / elapsedSeconds.count() << " GB/s\n";
}
