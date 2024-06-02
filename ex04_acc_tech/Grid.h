#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cmath>
#include <iostream>
#include <cuda_runtime.h>

struct Grid
{
    size_t box_size;
    double cell_length;
    int cells_per_side;
    size_t grid_size;
    std::vector<std::vector<int>> grid;
    int* d_grid;
    int* d_cell_start;
    int* d_cell_end;

    Grid(size_t box_size, double sigma)
        : box_size(box_size), cell_length(2.5 * sigma), d_grid(nullptr), d_cell_start(nullptr), d_cell_end(nullptr) {
        cells_per_side = static_cast<int>(box_size / cell_length);
        grid_size = cells_per_side * cells_per_side * cells_per_side;
        grid.resize(grid_size);
    }

    ~Grid() {
        if (d_grid) cudaFree(d_grid);
        if (d_cell_start) cudaFree(d_cell_start);
        if (d_cell_end) cudaFree(d_cell_end);
    }

    int getCellIndex(int x, int y, int z) const {
        return x * cells_per_side * cells_per_side + y * cells_per_side + z;
    }

    void addParticle(double px, double py, double pz, int particle_index, Molecule* particles) {
        int cell_x = static_cast<int>(px / cell_length);
        int cell_y = static_cast<int>(py / cell_length);
        int cell_z = static_cast<int>(pz / cell_length);
        int cell_index = getCellIndex(cell_x, cell_y, cell_z);

        particles[particle_index].cell_id = cell_index; // Assign cell_id to the particle
        grid[cell_index].push_back(particle_index);
    }

    const std::vector<int>& getParticlesInCell(int cell_index) const {
        return grid[cell_index];
    }

    void printGridInfo() const {
        std::cout << "Grid dimensions: " << cells_per_side << " x " << cells_per_side << " x " << cells_per_side << "\n";
        std::cout << "Total cells: " << grid_size << "\n";
    }

    void printGridContents( Molecule* particles) {
        for (int x = 0; x < cells_per_side; ++x) {
            for (int y = 0; y < cells_per_side; ++y) {
                for (int z = 0; z < cells_per_side; ++z) {
                    int cell_index = getCellIndex(x, y, z);
                    std::vector<int>& cell_particles = grid[cell_index];
                    std::cout << "Cell (" << x << ", " << y << ", " << z << "): ";
                    for (int particle_index : cell_particles) {
                        particles[particle_index].cell_id = cell_index;
                        std::cout << particles[particle_index].unique << " " << particles[particle_index].cell_id << " ";
                    
                    }
                    std::cout << "\n";
                }
            }
        }
    }

    void allocateDeviceMemory() {
        size_t total_elements = 0;
        for (const auto& cell : grid) {
            total_elements += cell.size();
        }
        
        cudaMalloc(&d_grid, total_elements * sizeof(int));
        cudaMalloc(&d_cell_start, grid_size * sizeof(int));
        cudaMalloc(&d_cell_end, grid_size * sizeof(int));
    }

    void copyToDevice() {
        std::vector<int> h_grid;
        std::vector<int> h_cell_start(grid_size, 0);
        std::vector<int> h_cell_end(grid_size, 0);
        
        int offset = 0;
        for (size_t i = 0; i < grid.size(); ++i) {
            h_cell_start[i] = offset;
            h_grid.insert(h_grid.end(), grid[i].begin(), grid[i].end());
            offset += grid[i].size();
            h_cell_end[i] = offset;
        }

        cudaMemcpy(d_grid, h_grid.data(), h_grid.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cell_start, h_cell_start.data(), grid_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cell_end, h_cell_end.data(), grid_size * sizeof(int), cudaMemcpyHostToDevice);
    }
};

#endif
