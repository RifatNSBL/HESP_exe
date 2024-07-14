#ifndef GRID_H
#define GRID_H

#include <vector>

struct Grid
{
    size_t box_size;
    double cell_length;
    int cells_per_side;
    int* cells;
    int* particle;
    int* cells_h;
    int* particles_h;
    size_t number_of_cells;
    int number_of_particles;

    Grid(size_t box_size, double sigma, int number_of_particles, double max_particle_size, double cell_length_mulltiplier)
        : box_size(box_size), cell_length(cell_length_mulltiplier * max_particle_size) {
        this->cells_per_side = static_cast<int>(box_size / cell_length);
        this->number_of_cells = cells_per_side * cells_per_side * cells_per_side;
        this->number_of_particles = number_of_particles;
    }

    void allocateDeviceMemory() { 
        cudaMallocManaged(&this->cells, this->number_of_cells * sizeof(int));
        for(int i = 0; i < this->number_of_cells; i++)
            this->cells[i] = -1;
        cudaMallocManaged(&this->particle, this->number_of_particles * sizeof(int));
        for(int i = 0; i < this->number_of_particles; i++)
            this->particle[i] = i;
    }

    ~Grid(){
        cudaFreeHost(cells_h);
        cudaFreeHost(particle);
        cudaFree(cells);
        cudaFree(particle);
    }

};

#endif
