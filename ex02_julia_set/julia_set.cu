
#include "CImg.h"
using namespace cimg_library;
#include <iostream>
#include <cuda/std/complex>
#include <complex>
#include <vector>
#include "lodepng.h"

using my_Complex = cuda::std::complex<double>; 

void encodeOneStep(const char* filename, const unsigned char *image, unsigned int width, unsigned int height) {

    std::vector<unsigned char> png;
    lodepng::State state;
        
    unsigned error = lodepng::encode(png, image, width, height, state);
    if(!error) lodepng::save_file(png, filename);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}


__global__ void initWindow(my_Complex *window, int *iter_matrix, int length,
                           int width, double min_value, double factor_x, double factor_y){

    int row_index = blockIdx.x * blockDim.x + threadIdx.x;
    int col_index = blockIdx.y * blockDim.y + threadIdx.y;

    // size_t strideX = gridDim.x * blockDim.x;
    // size_t strideY = gridDim.y * blockDim.y;
    // for(int i = row_index; i < length; i += strideX)
    //     for(int j = col_index; j < width; j += strideY){
    //         if (row_index < length && col_index < width){
    //             window[i * width + j] = my_Complex(i * factor_x, j * factor_y);
    //             iter_matrix[i * width + j] = 0;
    //         }
    //     }

    // size_t strideX = gridDim.x * blockDim.x;
    // size_t strideY = gridDim.y * blockDim.y;

    if (row_index < length && col_index < width){
        //for(int i = col_index; i < length; i++){
            window[row_index * width + col_index] = my_Complex( min_value + row_index * factor_x, min_value + col_index * factor_y);
            iter_matrix[row_index * width + col_index] = 0;
        //}
    }
}


__global__ void applyIteration(my_Complex *window, int *iter_matrix, 
                               int length, int width, my_Complex constant, double MAX_VAL, int MAX_ITER_NUM){

    int row_index = blockIdx.x * blockDim.x + threadIdx.x;
    int col_index = blockIdx.y * blockDim.y + threadIdx.y;

    // size_t strideX = gridDim.x * blockDim.x;
    // size_t strideY = gridDim.y * blockDim.y;
    // for(int i = row_index; i < length; i += strideX)
    //     for(int j = col_index; j < width; j += strideY){
            if (cuda::std::abs(window[row_index * width + col_index]) < MAX_VAL  && (row_index < length && col_index < width)){
                    window[row_index * width + col_index] = cuda::std::pow(window[row_index * width + col_index], 2) + constant;
                    iter_matrix[row_index * width + col_index] += 1;
                }
    }

int main(int argc, char *argv[]) {
    // default values
    int length = 100;
    int width = 1000;

    double max_val = 2.0;
    double min_val = -2.0;
    int MAX_ITER_NUM = 128;
    double MAX_ABS_VAL = 10.0;
    double const_real = -0.8;
    double const_imag = 0.156;

    int size = length * width * sizeof (my_Complex);
    
    double factor_x = (max_val - min_val) / (length - 1);
    double factor_y = (max_val - min_val) / (width - 1); 

    my_Complex constant(const_real, const_imag); 
    my_Complex *window;

    int *iter_matrix;
    int size_2 = length * width * sizeof(int);

    cudaMallocManaged (&window, size);
    cudaMemPrefetchAsync(window, size, 0);
    

    cudaMallocManaged (&iter_matrix, size_2);
    cudaMemPrefetchAsync(iter_matrix, size_2, 0);


    dim3 threads_per_block (8, 40, 1);

    dim3 number_of_blocks (25, 25, 1);

    initWindow<<<number_of_blocks, threads_per_block>>>(window, iter_matrix, length, width, min_val, factor_x, factor_y);
 

    for(int i = 0; i < MAX_ITER_NUM; i++)
        applyIteration<<<number_of_blocks, threads_per_block>>>(window, iter_matrix, length, width, constant, MAX_ABS_VAL, MAX_ITER_NUM);

    cudaDeviceSynchronize();


    // for(int i = 0; i < length; i++){
    //     for(int j = 0; j < width; j++)
    //         std::cout << int(iter_matrix[i * width + j]) << " ";
    //     std::cout << "\n";
    //     }

    const char* name= "image.png";
    
    const unsigned char *to_print = (unsigned char *)iter_matrix;

    encodeOneStep(name, to_print, width, length);

    
    // const CImg<unsigned char> image = CImg<>(name).normalize(0,255);

    //CImgDisplay main_disp(image,"Julia set image",0);

    // while (!main_disp.is_closed() && !main_disp.is_keyESC()
    //        && !main_disp.is_keyQ() ) {
    //     main_disp.display(image);
    //     cimg::wait(20);
    // }

    cudaFree(window);
    cudaFree(iter_matrix);

}
