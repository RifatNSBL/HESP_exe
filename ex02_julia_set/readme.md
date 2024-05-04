Notes for Exercise 2: Julia Set

1. Was tested on Ubuntu 22.04 with Kernel version 6.5.0-28-generic, cuda version 12.2, g++ version 11.4.0, driver version 535.171.04;
2. Hardware: NVIDIA RTX 4060 mobile (https://www.techpowerup.com/gpu-specs/geforce-rtx-4060-mobile.c3946)
3. To run it in Huber CIP computers, need try to add -ccbin clang-14 -allow-unsupported-compiler -lstdc++ into line 23 in Makefile,
and remove -arch=sm_89 from it;
4. To build just run make in terminal;
5. To change initial parametes, go to julia_set.cu -> lines 58-65
6. If successfully built, the png image is located in build folder.
7. Sometimes, the image.png looks weird, to get correct image - uncomment lines 112-113 in julia_set.cu
