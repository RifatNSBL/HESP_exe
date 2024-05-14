Notes for Exercise 3: Molecular Dynamics

1. Was tested on Ubuntu 22.04 with Kernel version 6.5.0-28-generic, cuda version 12.2, g++ version 11.4.0, driver version 535 171.04;
2. Hardware: NVIDIA RTX 4060 mobile (https://www.techpowerup.com/gpu-specs/geforce-rtx-4060-mobile.c3946)
3. To run it on Huber CIP computers, need to add -ccbin clang-14 -allow-unsupported-compiler -lstdc++ -ln into line 23 in Makefile, and remove -arch=sm_89 from it;
4. To build just run make -B in terminal;
5. To change initial parametes (sigma, eps, step size and number of steps), go to Makefile -> lines 9-13
6. If successfully built, the executable will start automatically, every 100 steps(molecular_dynamics -> line 145) total energy will be printed in the console.
7. Simulation files can be found in visualization folder, to see the particle animation use Paraview.
8. All particles initial values are written in points.txt with the following structure:
8.1. First line always represent number of particles.
8.2. Every next line provide particle values: mass, position x, position y, position z, velocity x, velocity y, velocity z, acceleration x, acceleration y, acceleration z.
