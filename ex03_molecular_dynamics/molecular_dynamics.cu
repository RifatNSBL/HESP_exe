#include <chrono>

#include "Molecule.h"
#include "Force.h"
#include "util.h"

__device__ void update_acc(Molecule &particle, Force &force){
    
    particle.xa = force.x / particle.mass;
    particle.ya = force.y / particle.mass;
    particle.za = force.z / particle.mass;
}

__global__ void calculate_forces(Molecule *particle, size_t num_molecules,
                                 double sigma, double eps){
        int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < num_molecules){
            Molecule _particle = particle[index];

            Force _force_i;
            
            _force_i.x = 0;
            _force_i.y = 0;
            _force_i.z = 0;

            double sigma_sqr  = sigma * sigma;
            for(auto i = 0; i < num_molecules; i++){
                if (i == index) continue;
                Molecule _particle_i = particle[i];
                
                double x_proj = _particle.x - _particle_i.x;
                double y_proj = _particle.y - _particle_i.y;
                double z_proj = _particle.z - _particle_i.z;

                double dist_sqr = x_proj * x_proj  + y_proj * y_proj + z_proj * z_proj;

                double factor_sqr = sigma_sqr / dist_sqr;
                double factor_hex = factor_sqr * factor_sqr * factor_sqr;

                double res = 24 * eps * factor_hex * (2 * factor_hex - 1) / dist_sqr;

                _force_i.x +=  res * x_proj;
                _force_i.y +=  res * y_proj;
                _force_i.z +=  res * z_proj;
            }

            update_acc(particle[index], _force_i); // a(t + 1/2 dt)
        }
}

__global__ void integration_step_begin(Molecule *particle, int num_molecules, 
                                 double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        double time_step_sqr = time_step * time_step;
        particle[index].x += particle[index].xv * time_step + ( particle[index].xa * time_step_sqr)/2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        particle[index].y += particle[index].yv * time_step + ( particle[index].ya * time_step_sqr)/2; 
        particle[index].z += particle[index].zv * time_step + ( particle[index].za * time_step_sqr)/2;

        particle[index].xv += particle[index].xa * time_step / 2; // v(t + 1/2 dt) = v(t) + a(t) * dt / 2
        particle[index].yv += particle[index].ya * time_step / 2;
        particle[index].zv += particle[index].za * time_step / 2;

        }
}

__global__ void integration_step_end(Molecule *particle, int num_molecules, 
                                 double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){

        particle[index].xv += particle[index].xa * time_step / 2; // v(t + dt) = v(t + 1/2 dt) + a(t + dt) * dt / 2
        particle[index].yv += particle[index].ya * time_step / 2;
        particle[index].zv += particle[index].za * time_step / 2;
    }
}

int main(int argc, char *argv[]){
    
    double time_step = 0.003;
    int num_steps = 5;
    double sigma = 1; // stable distance = abs_dist * pow(0.5, 1/6)
    double eps = 2;

    size_t num_molecules = 1000;


    size_t size_molecule = num_molecules * sizeof(Molecule);

    Molecule *grid;

    cudaMallocManaged(&grid, size_molecule);

    fill_particles_on(grid, num_molecules);

    cudaMemPrefetchAsync(grid, size_molecule, 0);

    int NUM_THREAD = 512;
    int NUM_BLOCK = (num_molecules + NUM_THREAD - 1) / NUM_THREAD;

    for(int i = 0; i < num_steps; i++){
        auto start = std::chrono::steady_clock::now();
        // most stuff should be done here
        integration_step_begin<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, time_step, sigma, eps);
        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, sigma, eps);
        integration_step_end<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, time_step, sigma, eps);

        if (i % 10 == 0){
            //cudaDeviceSynchronize();
            auto end = std::chrono::steady_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // don't know if it's correct

            std::cout << "Iteration took " << time.count() << " milliseconds\n";
            //calculate_energy(grid, num_molecules, sigma, eps);
            //writeVTK(i, num_molecules, grid);
        }
    }

    cudaFree(grid);
}
