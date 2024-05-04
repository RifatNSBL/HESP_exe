#include <iostream>
#include <cuda/std/complex>
#include <complex>
#include <vector>
#include <fstream>
#include "lodepng.h"

using my_Complex = cuda::std::complex<double>; 


// should be in separate file, but no time to implement/debug it
void makeImagePPM(int *iter_matrix, int height, int width, const char* name);
void encodeOneStep(const char* filename, std::vector<unsigned char> &image, unsigned int width, unsigned int height);


__global__ void initWindow(my_Complex *window, int *iter_matrix, int height, int width, 
                           double min_value, double factor_x, double factor_y){

    int row_index = blockIdx.x * blockDim.x + threadIdx.x;
    int col_index = blockIdx.y * blockDim.y + threadIdx.y;

    if (row_index < height && col_index < width){
        window[row_index * width + col_index] = my_Complex( min_value + row_index * factor_x, min_value + col_index * factor_y);
        iter_matrix[row_index * width + col_index] = 0;
    }
}


__global__ void applyIteration(my_Complex *window, int *iter_matrix, int height, int width, 
                               my_Complex constant, double MAX_VAL, int MAX_ITER_NUM){

    int row_index = blockIdx.x * blockDim.x + threadIdx.x;
    int col_index = blockIdx.y * blockDim.y + threadIdx.y;

    if (cuda::std::abs(window[row_index * width + col_index]) < MAX_VAL  && (row_index < height && col_index < width)){
            window[row_index * width + col_index] = cuda::std::pow(window[row_index * width + col_index], 2) + constant;
            iter_matrix[row_index * width + col_index] += 1;
    }
}


__global__ void makeImagePNG(int* iter_matrix, int *image, int height, int width){

    int row_index = blockIdx.x * blockDim.x + threadIdx.x;
    int col_index = blockIdx.y * blockDim.y + threadIdx.y;

    if (row_index < height && col_index < width){
        image[4 * width * row_index + 4 * col_index + 0] = iter_matrix[row_index * width + col_index];
        image[4 * width * row_index + 4 * col_index + 1] = 0;
        image[4 * width * row_index + 4 * col_index + 2] = 20;
        image[4 * width * row_index + 4 * col_index + 3] = 255;
    }
}


int main(int argc, char *argv[]) {
    // default values
    int height = 4800;
    int width = 3600;
    double max_val = 2.0;
    double min_val = -2.0;
    int MAX_ITER_NUM = 128;
    double MAX_ABS_VAL = 10.0;
    double const_real = -0.8;
    double const_imag = 0.178;


    my_Complex constant(const_real, const_imag);
    double factor_x = (max_val - min_val) / (height - 1);
    double factor_y = (max_val - min_val) / (width - 1); 

    my_Complex *window;
    int size = height * width * sizeof (my_Complex); 
    cudaMallocManaged (&window, size);
    cudaMemPrefetchAsync(window, size, 0);
    
    int *iter_matrix;
    int size_2 = height * width * sizeof(int);
    cudaMallocManaged (&iter_matrix, size_2);
    cudaMemPrefetchAsync(iter_matrix, size_2, 0);

    int *image;
    int size_img = 4 * height * width * sizeof (int);
    cudaMallocManaged (&image, size_img);
    cudaMemPrefetchAsync(image, size_img, 0);

    // x*y*z <= 1024
    dim3 threads_per_block (32, 32, 1);
    dim3 number_of_blocks ((height / (double)threads_per_block.x), (width / (double)threads_per_block.y), 1);


    initWindow<<<number_of_blocks, threads_per_block>>>(window, iter_matrix, height, width, min_val, factor_x, factor_y);
 
    for(int i = 0; i < MAX_ITER_NUM; i++)
        applyIteration<<<number_of_blocks, threads_per_block>>>(window, iter_matrix, height, width, constant, MAX_ABS_VAL, MAX_ITER_NUM);
    
    makeImagePNG<<<number_of_blocks, threads_per_block>>>(iter_matrix, image, height, width);

    // all matrices in GPU, so no need to synchronize every time
    cudaDeviceSynchronize();


    std::vector<unsigned char> image_png(image, image + height * width * 4);

    // .png image is ~70-100 times cheaper in space than .ppl
    // but might encounter artifacts
    const char* name = "build/image.png";
    encodeOneStep(name, image_png, width, height);


    // use this to make .ppl image - 100% correct image
    // const char *name_ppm = "build/image.ppm";
    // makeImagePPM(iter_matrix, height, width, name_ppm);

    cudaFree(window);
    cudaFree(iter_matrix);
    cudaFree(image);

}


void makeImagePPM(int *iter_matrix, int height, int width, const char* name){
    std::ofstream image_ppm(name);
    if (!image_ppm) {
        std::cerr << "Error: Unable to create image file" << std::endl;
        return;
    }
    image_ppm << "P3\n" << width << " " << height << "\n255\n";
    for (int x = 0; x < height; ++x) {
        for (int y = 0; y < width; ++y) {
            int red = (iter_matrix[x * width + y]);
            int green = 0;
            int blue = 20;
            image_ppm << red << " " << green << " " << blue << " ";

        }
        image_ppm << std::endl;
    }
    image_ppm.close();
}

void encodeOneStep(const char* filename, std::vector<unsigned char> &image, unsigned int width, unsigned int height) {
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    if(!error) lodepng::save_file(png, filename);
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}