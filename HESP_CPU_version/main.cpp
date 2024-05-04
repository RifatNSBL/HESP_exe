#include <iostream>
#include <fstream>
#include <complex>

// Define constants
const std::complex<double> c(-0.8, 0.2); // Constant c value
const double threshold = 10;              // Absolute value threshold
const int max_iterations = 1000;          // Maximum number of iterations
const int width = 800;                    // Image width
const int height = 800;                   // Image height


int main() {
    // Create output file
    std::ofstream image("julia_set.ppm");
    if (!image) {
        std::cerr << "Error: Unable to create image file" << std::endl;
        return 1;
    }

    // Write PPM header
    image << "P3\n" << width << " " << height << "\n255\n";

    // Iterate over each pixel
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Scale pixel positions to a suitable range
            double real = (4.0 * x - width) / width;
            double imag = (4.0 * y - height) / height;

            // Initialize starting value z_0
            std::complex<double> z(real, imag);
            int iterations = 0;

            // Apply iteration rule
            while (std::abs(z) < threshold && iterations < max_iterations) {
                z = z * z + c;
                iterations++;
            }

             // Color mapping
            int red = (iterations * 7) % 256;
            int green = (iterations * 13) % 256;
            int blue = (iterations * 23) % 256;

            image << red << " " << green << " " << blue << " ";

        }
        image << std::endl;
    }

    // Close the file
    image.close();

    std::cout << "Image generation completed." << std::endl;

    return 0;
}
