#include <chrono>
#include <iostream>
#include <cmath>

#include "Molecule.h"
#include "Force.h"
#include "util.h"
#include "Grid.h"

__device__ void update_acc(Molecule &particle, Force &force){
    Molecule _particle = particle;
    Force    _force = force;

    _particle.xa = _force.x / _particle.mass;
    _particle.ya = _force.y / _particle.mass;
    _particle.za = _force.z / _particle.mass;

    particle = _particle;
    force    = _force;
}

__device__ double min_dist(double particle_i, double particle_j, size_t box_size){
    if      (abs(particle_i - particle_j) < ((double)box_size / 2))
            return particle_i - particle_j;

    else if (particle_i < particle_j)
            return particle_i + (box_size - particle_j);

    else    return particle_i - (box_size + particle_j);
}

__global__ void calculate_forces(Molecule *particles, size_t num_particles,
                                 double sigma, double eps, size_t box_size,
                                 double cutoff_dist, int cells_per_side,
                                 int* particle_list, int* cells_list) {

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    
    if (index < num_particles) {
        Molecule _particle = particles[index];

        Force _force_i;
        _force_i.x = 0;
        _force_i.y = 0;
        _force_i.z = 0;

        double sigma_sqr = sigma * sigma;
        double cutoff_sqr = cutoff_dist * cutoff_dist;
        int cell_x = static_cast<int>(_particle.x / cutoff_dist);
        int cell_y = static_cast<int>(_particle.y / cutoff_dist);
        int cell_z = static_cast<int>(_particle.z / cutoff_dist);
        
        //Iterate over the 3x3x3 grid around the current cell
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    int neighbor_cell_x = (cell_x + dx + cells_per_side) % cells_per_side;
                    int neighbor_cell_y = (cell_y + dy + cells_per_side) % cells_per_side;
                    int neighbor_cell_z = (cell_z + dz + cells_per_side) % cells_per_side;

                    int neighbor_cell_index = neighbor_cell_x * cells_per_side * cells_per_side
                                              + neighbor_cell_y * cells_per_side
                                              + neighbor_cell_z;
                                              
                    int first_particle = cells_list[neighbor_cell_index];
                    if (first_particle != -1 && first_particle != index){
                        Molecule _particle_i = particles[first_particle];

                        double x_proj = min_dist(_particle.x, _particle_i.x, box_size);
                        double y_proj = min_dist(_particle.y, _particle_i.y, box_size);
                        double z_proj = min_dist(_particle.z, _particle_i.z, box_size);

                        double dist_sqr = x_proj * x_proj + y_proj * y_proj + z_proj * z_proj;
                        if (dist_sqr < cutoff_sqr){

                            double factor_sqr = sigma_sqr / dist_sqr;
                            double factor_hex = factor_sqr * factor_sqr * factor_sqr;

                            double res = 24 * eps * factor_hex * (2 * factor_hex - 1) / dist_sqr;

                            _force_i.x += res * x_proj;
                            _force_i.y += res * y_proj;
                            _force_i.z += res * z_proj;
                        }
                    }
                    int sliding_index = particle_list[first_particle];
                    while (sliding_index != -1) {
                        //int neighbor_index = grid[i];
                        if (sliding_index == index) sliding_index = particle_list[sliding_index];
                        else {
                            Molecule _particle_i = particles[sliding_index];
                            double x_proj = min_dist(_particle.x, _particle_i.x, box_size);
                            double y_proj = min_dist(_particle.y, _particle_i.y, box_size);
                            double z_proj = min_dist(_particle.z, _particle_i.z, box_size);

                            double dist_sqr = x_proj * x_proj + y_proj * y_proj + z_proj * z_proj;
                            if (dist_sqr < cutoff_sqr){

                                double factor_sqr = sigma_sqr / dist_sqr;
                                double factor_hex = factor_sqr * factor_sqr * factor_sqr;

                                double res = 24 * eps * factor_hex * (2 * factor_hex - 1) / dist_sqr;

                                _force_i.x += res * x_proj;
                                _force_i.y += res * y_proj;
                                _force_i.z += res * z_proj;

                            }
                            sliding_index = particle_list[sliding_index];
                        }
                    }
                }
            }
        }
        update_acc(particles[index], _force_i); // a(t + 1/2 dt)
    }
}


__global__ void integration_step_begin(Molecule *particles, int num_particles,
                                       double time_step, size_t box_size){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles){
        double time_step_sqr = time_step * time_step;
        Molecule _particle = particles[index];


        // Update particle position
        _particle.x += _particle.xv * time_step + (_particle.xa * time_step_sqr) / 2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        _particle.y += _particle.yv * time_step + (_particle.ya * time_step_sqr) / 2;
        _particle.z += _particle.zv * time_step + (_particle.za * time_step_sqr) / 2;

        // Handle boundary conditions
        if (_particle.x < 0) _particle.x += box_size;
        if (_particle.y < 0) _particle.y += box_size;
        if (_particle.z < 0) _particle.z += box_size;
        if (_particle.x > box_size) _particle.x -= box_size;
        if (_particle.y > box_size) _particle.y -= box_size;
        if (_particle.z > box_size) _particle.z -= box_size;

        // Update particle velocity
        _particle.xv += _particle.xa * time_step / 2; // v(t + 1/2 dt) = v(t) + a(t) * dt / 2
        _particle.yv += _particle.ya * time_step / 2;
        _particle.zv += _particle.za * time_step / 2;
        // Write updated particle back to global memory
        particles[index] = _particle;
    }
}

__global__ void update_cells(Molecule *particles, int num_particles,int cells_per_side, double cell_length,
                                       int* particle_list, int* cells_list){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles){
        Molecule _particle = particles[index];

        // Calculate new cell indices based on the updated position
        int cell_x = static_cast<int>(_particle.x / cell_length);
        int cell_y = static_cast<int>(_particle.y / cell_length);
        int cell_z = static_cast<int>(_particle.z / cell_length);
        int new_cell_id = cell_x * cells_per_side * cells_per_side + cell_y * cells_per_side + cell_z;

        particle_list[index] = atomicExch(&cells_list[new_cell_id], index);
        _particle.cell_id = new_cell_id;
        
    }
}

__global__ void show_list(int cells_per_side, int* particle_list, int* cells_list, size_t num_particles, int iter_num){
    printf("ITER NUM %d\n\n\n", iter_num);
    int cell_triple = cells_per_side * cells_per_side * cells_per_side;
    int min_ones = 0;
    for(int i = 0; i < cell_triple; i++){
        printf("%d ", cells_list[i]);
        if (cells_list[i] == -1) min_ones +=1;
    }
    printf("TOTAL NUMBER OF CELLS %d\n\n\n", cell_triple);
    for(int i = 0; i  < num_particles; i++){
        printf("%d ", particle_list[i]);
        if (particle_list[i] == -1) min_ones +=1;
    }
    printf("\n\nNUMBER OF -1's %d\n", min_ones);
}



__global__ void regenerate_lists(int* cells_list, size_t number_of_cells){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < number_of_cells){
        cells_list[index] = -1;
    }
}


__global__ void integration_step_end(Molecule *particle, int num_particles,
                                     double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles){
        Molecule _particle = particle[index];

        _particle.xv += _particle.xa * time_step / 2; // v(t + dt) = v(t + 1/2 dt) + a(t + dt) * dt / 2
        _particle.yv += _particle.ya * time_step / 2;
        _particle.zv += _particle.za * time_step / 2;

        particle[index] = _particle;
    }
}

__global__ void make_me_updated(){
    printf("Cells updated\n");

}

int main(int argc, char *argv[]) {
    // Parse command line arguments
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " <time_step> <num_steps> <sigma> <eps> <particle_datafile> <cutoff_dist> <box_size>" << std::endl;
        return -1;
    }

    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    double sigma = atof(argv[3]);
    double eps = atof(argv[4]);
    std::string particle_datafile = argv[5];
    double cutoff_dist = atof(argv[6]);
    size_t box_size = atoi(argv[7]);
    double cell_length = 2.5 * sigma;

    size_t num_particles = 0;

    std::ifstream file(particle_datafile);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << particle_datafile << std::endl;
        return -1;
    }

    std::string line;
    std::getline(file, line);
    num_particles = std::stoi(line);

    size_t size_molecule = num_particles * sizeof(Molecule);

    Molecule *particles;
    cudaMallocManaged(&particles, size_molecule);
    if (particles == nullptr) {
        std::cerr << "Failed to allocate managed memory for particles!" << std::endl;
        return -1;
    }

    fill_particles(particles, num_particles, cell_length, box_size / cell_length, file);

    Grid grid(box_size, sigma, num_particles);
    grid.construct_first_iteration(particles);
    grid.allocateDeviceMemory();

    int NUM_THREAD = 256;
    int NUM_BLOCK = (num_particles + NUM_THREAD - 1) / NUM_THREAD;

    auto start_global = std::chrono::steady_clock::now();

    for (int i = 0; i < num_steps; i++) {
        auto start = std::chrono::steady_clock::now();

        integration_step_begin<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, box_size);

        regenerate_lists<<<NUM_BLOCK, NUM_THREAD>>>(grid.cells, grid.number_of_cells);

        update_cells<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, grid.cells_per_side, grid.cell_length, grid.particle, grid.cells);

        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, sigma, eps, box_size, cutoff_dist, grid.cells_per_side, grid.particle, grid.cells);

        integration_step_end<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, sigma, eps);

        if (i % 10 == 0) {
            auto end = std::chrono::steady_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Iteration took " << time.count() << " milliseconds\n";
            writeVTK(i, num_particles, particles);
        }
    }

    auto end_global = std::chrono::steady_clock::now();
    auto time_global = std::chrono::duration_cast<std::chrono::seconds>(end_global - start_global);
    std::cout << "Total iteration time " << time_global.count() << " seconds\n";

    cudaFree(particles);
}
