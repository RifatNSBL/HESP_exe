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

    Grid(size_t box_size, double sigma, int number_of_particles)
        : box_size(box_size), cell_length(2.5 * sigma) { // *2.5
        this->cells_per_side = static_cast<int>(box_size / cell_length);
        this->number_of_cells = cells_per_side * cells_per_side * cells_per_side;
        this->number_of_particles = number_of_particles;
    }

    void allocateDeviceMemory() { 
        // cudaMallocManaged(&this->cells, this->number_of_cells * sizeof(int));
        // for(int i = 0; i < this->number_of_cells; i++)
        //     this->cells[i] = -1;
        // cudaMallocManaged(&this->particles, this->number_of_particles * sizeof(int));
        // for(int i = 0; i < this->number_of_particles; i++)
        //     this->particles[i] = i;
        cudaMalloc(&this->cells, this->number_of_cells * sizeof(int));
        cudaMemcpy(this->cells, this->cells_h, this->number_of_cells * sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc(&this->particle, this->number_of_particles * sizeof(int));
        cudaMemcpy(this->particle, this->particles_h, this->number_of_particles * sizeof(int), cudaMemcpyHostToDevice);

    }

    void construct_first_iteration(Molecule* particles){
        cudaMallocHost(&this->cells_h, this->number_of_cells * sizeof(int));
        cudaMallocHost(&this->particles_h, this->number_of_particles * sizeof(int));

        for(int i = 0; i < this->number_of_particles; i++)
            this->particles_h[i] = i;
        for(int i = 0; i < this->number_of_cells; i++)
            this->cells_h[i] = -1;

        for(int i = 0; i < this->number_of_particles; i++){

            if (this->cells_h[particles[i].cell_id] == -1){
                this->particles_h[i] = -1;
                this->cells_h[particles[i].cell_id] = i;
            }

            else{
                int value_in_cell = this->cells_h[particles[i].cell_id];
                this->cells_h[particles[i].cell_id] = i;
                this->particles_h[i] = value_in_cell;
            }

        }
        printf("Filling particles finished\n");
        for(int i = 0; i < this->number_of_cells; i++)
            printf("%d ", this->cells_h[i]);
        printf("\n\n\n");
        for(int i = 0; i < this->number_of_particles; i++)
            printf("%d ", this->particles_h[i]);
        printf("\n\n\n");
    }

    ~Grid(){
        cudaFreeHost(cells_h);
        cudaFreeHost(particle);
        cudaFree(cells);
        cudaFree(particle);
    }

};

#endif
