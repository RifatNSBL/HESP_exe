#include <chrono>

#include "../util.h"
#include "stream-util.h"

__global__ void stream(size_t nx, const double *__restrict__ src, double *__restrict__ dest) {
    int i;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    for (i = index; i < nx; i += stride)
    {
        dest[i] = src[i] + 1;
    }
}


int main(int argc, char *argv[]) {
    size_t nx, nItWarmUp, nIt;
    parseCLA_1d(argc, argv, nx, nItWarmUp, nIt);

    double* src;
    double* dest;

    size_t size = nx * sizeof(double);
    cudaMallocManaged(&src, size);
    cudaMallocManaged(&dest, size);
    
    cudaMemPrefetchAsync(src, size, 0);
    cudaMemPrefetchAsync(dest, size, 0);

    int deviceId;
    cudaGetDevice(&deviceId);
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, deviceId);
    int multiProcessorCount = props.multiProcessorCount;
    int warpSize = props.warpSize;

    size_t threads_per_block = 256;
    size_t number_of_blocks = warpSize * multiProcessorCount / 2;

    // init
    initStream(src, nx);

    // warm-up
    for (int i = 0; i < nItWarmUp; ++i) {
        stream<<<number_of_blocks, threads_per_block>>>(nx, src, dest);
        cudaDeviceSynchronize();
        std::swap(src, dest);
    }


    // measurement
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < nIt; ++i) {
        stream<<<number_of_blocks, threads_per_block>>>(nx, src, dest);
        cudaDeviceSynchronize();
        std::swap(src, dest);
    }

    
    auto end = std::chrono::steady_clock::now();

    printStats(end - start, nx, nIt, streamNumReads, streamNumWrites);

    // check solution
    checkSolutionStream(src, nx, nIt + nItWarmUp);

    cudaFree(src);
    cudaFree(dest);
    return 0;
}