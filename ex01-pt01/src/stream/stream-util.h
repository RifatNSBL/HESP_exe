#pragma once

#include <iostream>

constexpr unsigned int streamNumWrites = 1;
constexpr unsigned int streamNumReads = 1 + 0 * streamNumWrites; // assume non-temporal stores

void initStream(double *vec, size_t nx) {
    for (size_t i = 0; i < nx; ++i)
        vec[i] = (double) i;
}

void checkSolutionStream(const double *const vec, size_t nx, size_t nIt) {
    for (size_t i = 0; i < nx; ++i)
        if ((double) (i + nIt) != vec[i]) {
            std::cerr << "Stream check failed for element " << i << " (expected " << i + nIt << " but got " << vec[i] << ")" << std::endl;
            return;
        }

    std::cout << "  Passed result check" << std::endl;
}
