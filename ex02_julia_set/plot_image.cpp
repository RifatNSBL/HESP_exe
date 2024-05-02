#include <complex>
#include <vector>
#include "CImg.h"
using namespace cimg_library;

int main(int argc, char *argv[]) {
    CImg<unsigned char> image;
    image.load("image.png");
    //image.resize(1920, 1080);
    //image.display();

    CImgDisplay main_disp(image, "The image");
    while (!main_disp.is_closed())
    {
        main_disp.wait();
    }

    return 0;
}


// #include "CImg.h"
// using namespace cimg_library;
// int main() {
//   CImg<unsigned char> img(640,400,1,3);  // Define a 640x400 color image with 8 bits per color component.
//   img.fill(0);                           // Set pixel values to 0 (color : black)
//   unsigned char purple[] = { 255,0,255 };        // Define a purple color
//   img.draw_text(100,100,"Hello World",purple); // Draw a purple "Hello world" at coordinates (100,100).
//   img.display("My first CImg code");             // Display the image in a display window.
//   return 0;
// }