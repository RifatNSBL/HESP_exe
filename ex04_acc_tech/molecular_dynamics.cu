#include <chrono>

#include "Molecule.h"
#include "Force.h"
#include "util.h"

__device__ void update_acc(Molecule &particle, Force &force){
    
    particle.xa = force.x / particle.mass;
    particle.ya = force.y / particle.mass;
    particle.za = force.z / particle.mass;
}

__global__ void calculate_forces(Molecule *particle, Force *force, size_t num_molecules,
                                 double sigma, double eps){
        int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < num_molecules){
            
            double particle_pos_x = particle[index].x;
            double particle_pos_y = particle[index].y;
            double particle_pos_z = particle[index].z;
            
            double f_x = 0;
            double f_y = 0;
            double f_z = 0;
            double sigma_sqr = sigma * sigma;
            double cutoff_sqr = 2.5 * 2.5 * sigma_sqr;
            for(auto i = 0; i < num_molecules; i++){
                if (i == index) continue;
                
                double x_proj = particle_pos_x - particle[i].x;
                double y_proj = particle_pos_y - particle[i].y;
                double z_proj = particle_pos_z - particle[i].z;
                
                double dist_sqr = x_proj * x_proj  + y_proj * y_proj + z_proj * z_proj;
                if (dist_sqr < cutoff_sqr){
                    double factor_sqr = sigma_sqr / dist_sqr;
                    double factor_hex = factor_sqr * factor_sqr * factor_sqr;

                    double res = 24 * eps * factor_hex * (2 * factor_hex - 1) / dist_sqr;
                    f_x +=  res * x_proj;
                    f_y +=  res * y_proj;
                    f_z +=  res * z_proj;
                }
            }

            force[index].x = f_x;
            force[index].y = f_y;
            force[index].z = f_z;

            update_acc(particle[index], force[index]); // a(t + 1/2 dt)
        }
}

__global__ void integration_step_begin(Molecule *particle, int num_molecules, 
                                       double time_step, double sigma, double eps, 
                                       double max_pos_x, double max_pos_y, double max_pos_z){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        double time_step_sqr = time_step * time_step;
        particle[index].x += particle[index].xv * time_step + (particle[index].xa * time_step_sqr)/2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        particle[index].y += particle[index].yv * time_step + (particle[index].ya * time_step_sqr)/2; 
        particle[index].z += particle[index].zv * time_step + (particle[index].za * time_step_sqr)/2;

        if (particle[index].x < 0) particle[index].x +=  max_pos_x;
        if (particle[index].y < 0) particle[index].y +=  max_pos_y;
        if (particle[index].z < 0) particle[index].z +=  max_pos_z;

        if (particle[index].x > max_pos_x) particle[index].x -=  max_pos_x;
        if (particle[index].y > max_pos_y) particle[index].y -=  max_pos_y;
        if (particle[index].z > max_pos_z) particle[index].z -=  max_pos_z;

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
    
    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    double sigma = atof(argv[3]); // stable distance = abs_dist * pow(0.5, 1/6)
    double eps = atof(argv[4]);
    std::string patricle_datafile = argv[5];
    size_t box_size = atoi(argv[6]);

    size_t num_molecules = 0;
    double max_x = 0;
    double max_y = 0;
    double max_z = 0;

    std::ifstream file;
    file.open(patricle_datafile);

    if (file.is_open()){
        std::string line;
        std::getline(file, line);
        num_molecules = std::stoi(line);    
    }

    size_t size_molecule = num_molecules * sizeof(Molecule);
    size_t size_forces =  num_molecules * sizeof(Force);

    Molecule *grid_host, *grid_device;
    Force *forces_device;

    cudaMallocHost(&grid_host, size_molecule);
    cudaMalloc(&grid_device, size_molecule);

    cudaMalloc(&forces_device, size_forces);

    fill_particles(grid_host, num_molecules, file, max_x, max_y, max_z);

    size_t num_blocks_x = ceil(max_x / box_size);
    size_t num_blocks_y = ceil(max_y / box_size);
    size_t num_blocks_z = ceil(max_z / box_size);

    double max_pos_x = num_blocks_x * box_size;
    double max_pos_y = num_blocks_y * box_size;
    double max_pos_z = num_blocks_z * box_size;

    cudaMemcpy(grid_device, grid_host, size_molecule, cudaMemcpyHostToDevice);

    int NUM_THREAD = 512;
    int NUM_BLOCK = (num_molecules + NUM_THREAD - 1) / NUM_THREAD;

    for(int i = 0; i < num_steps; i++){
        auto start = std::chrono::steady_clock::now();
        // most stuff should be done here
        integration_step_begin<<<NUM_BLOCK, NUM_THREAD>>>(grid_device, num_molecules, time_step, 
                                                          sigma, eps, max_pos_x, max_pos_y, max_pos_z);
        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(grid_device, forces_device, num_molecules, sigma, eps);
        integration_step_end<<<NUM_BLOCK, NUM_THREAD>>>(grid_device, num_molecules, time_step, sigma, eps);

        if (i % 10 == 0){
            cudaDeviceSynchronize();
            auto end = std::chrono::steady_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // don't know if it's correct

            std::cout << "Iteration took " << time.count() << " milliseconds\n";
            cudaMemcpy(grid_host, grid_device, size_molecule, cudaMemcpyDeviceToHost);
            // calculate_energy(grid_host, num_molecules, sigma, eps);
            writeVTK(i, num_molecules, grid_host);
        }
    }

    cudaFree(grid_device);
    cudaFree(forces_device);
    cudaFreeHost(grid_host);
}
