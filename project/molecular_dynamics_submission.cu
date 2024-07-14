#include <chrono>
#include <iostream>
#include <cmath>
#include <curand_kernel.h>

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

__device__ int sign(double value){ return value > 0 ? 1 : -1; }

__global__ void calculate_forces(Molecule *particles, size_t num_particles, double time_step,
                                 size_t box_size, double eff_young_modulus, double eff_shear_modulus){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles) {
        Molecule _particle = particles[index];
        if (_particle.flag != 0)
        {
            Force _force_i;
            _force_i.x = 0;
            _force_i.y = 0;
            _force_i.z = -4.5;

            for(int i = 0; i < num_particles; i++) {
                Molecule _particle_i = particles[i];
                if (i == index) continue; 
                double x_proj = _particle.x - _particle_i.x;
                double y_proj = _particle.y - _particle_i.y;
                double z_proj = _particle.z - _particle_i.z;

                double dist_sqr = x_proj * x_proj + y_proj * y_proj + z_proj * z_proj;
                double dist = sqrt(dist_sqr);

                double overlap_distance = (_particle.diameter + _particle_i.diameter) / 2.0;
                if (dist < overlap_distance) { // Check if particles are overlapping

                    double overlap = overlap_distance - dist;
                    double sqrt_overlap = sqrt(overlap);
                    //double cube_sqrt_overlap = sqrt_overlap * sqrt_overlap * sqrt_overlap;

                    double normalized_x_proj = x_proj / dist;
                    double normalized_y_proj = y_proj / dist;
                    double normalized_z_proj = z_proj / dist;

                    double eff_radius = (_particle.diameter * _particle_i.diameter) / (2 * (_particle.diameter + _particle_i.diameter));
                    double sqrt_eff_radius = sqrt(eff_radius);

                    //find projections of overlap distance
                    double overlap_normalized_x = overlap * normalized_x_proj;
                    double overlap_normalized_y = overlap * normalized_y_proj;
                    double overlap_normalized_z = overlap * normalized_z_proj;

                    // Calculate damping force components
                    double relative_velocity_x = _particle.xv - _particle_i.xv;
                    double relative_velocity_y = _particle.yv - _particle_i.yv;
                    double relative_velocity_z = _particle.zv - _particle_i.zv;
                    
                    double normal_relative_velocity_x = relative_velocity_x * normalized_x_proj * normalized_x_proj;
                    double normal_relative_velocity_y = relative_velocity_y * normalized_y_proj * normalized_y_proj;
                    double normal_relative_velocity_z = relative_velocity_z * normalized_z_proj * normalized_z_proj;

                    double normal_force_abs = 1.333 * eff_young_modulus * sqrt_eff_radius;

                    double normal_force_x = normal_force_abs * sqrt(overlap_normalized_x * overlap_normalized_x) * sign(x_proj);
                    double normal_force_y = normal_force_abs * sqrt(overlap_normalized_y * overlap_normalized_y) * sign(y_proj);
                    double normal_force_z = normal_force_abs * sqrt(overlap_normalized_z * overlap_normalized_z) * sign(z_proj);
                    
                    double tangential_relative_velocity_x = relative_velocity_x - normal_relative_velocity_x;
                    double tangential_relative_velocity_y = relative_velocity_y - normal_relative_velocity_y;
                    double tangential_relative_velocity_z = relative_velocity_z - normal_relative_velocity_z;

                    double tangential_displacement_x = tangential_relative_velocity_x * time_step;
                    double tangential_displacement_y = tangential_relative_velocity_y * time_step;
                    double tangential_displacement_z = tangential_relative_velocity_z * time_step;

                    double factor_tangential_damping = -8 * eff_shear_modulus * sqrt_eff_radius;
                    double sqrt_factor_tangential_damping = sqrt(factor_tangential_damping);

                    double tangential_force_x = factor_tangential_damping * sqrt(abs(tangential_displacement_x)) * sign(tangential_displacement_x);
                    double tangential_force_y = factor_tangential_damping * sqrt(abs(tangential_displacement_y)) * sign(tangential_displacement_y);
                    double tangential_force_z = factor_tangential_damping * sqrt(abs(tangential_displacement_z)) * sign(tangential_displacement_z);


                    // add damping forces
                    double restitution = 0.8;
                    double log_rest = log(restitution);
                    double sqr_log_rest = log_rest * log_rest;
                    double sum = sqr_log_rest + 9.8696;
                    double beta = - log_rest / sum; 
                    double eff_mass = (_particle.mass * _particle_i.mass) / (_particle.mass + _particle_i.mass);
                    double sqrt_eff_mass = sqrt(eff_mass);
                    double factor_normal_damping = 2 * eff_young_modulus * sqrt_eff_radius;

                    double S_normal_damping_x = factor_normal_damping * sqrt(abs(overlap_normalized_x));
                    double S_normal_damping_y = factor_normal_damping * sqrt(abs(overlap_normalized_y));
                    double S_normal_damping_z = factor_normal_damping * sqrt(abs(overlap_normalized_z));

                    double factor_normal_damping_force = -1.825741858 * beta * sqrt_eff_mass;
                    double normal_damping_force_x = factor_normal_damping_force * sqrt(S_normal_damping_x)  * relative_velocity_x;
                    double normal_damping_force_y = factor_normal_damping_force * sqrt(S_normal_damping_y)  * relative_velocity_y;
                    double normal_damping_force_z = factor_normal_damping_force * sqrt(S_normal_damping_z)  * relative_velocity_z;

                    double factor_tangential_damping_force = -1.825741858 * beta * sqrt(eff_mass);
                    double tangential_damping_force_x = factor_tangential_damping_force * sqrt(abs(factor_tangential_damping) * sqrt(abs(overlap_normalized_x))) * tangential_relative_velocity_x;
                    double tangential_damping_force_y = factor_tangential_damping_force * sqrt(abs(factor_tangential_damping) * sqrt(abs(overlap_normalized_y))) * tangential_relative_velocity_y;
                    double tangential_damping_force_z = factor_tangential_damping_force * sqrt(abs(factor_tangential_damping) * sqrt(abs(overlap_normalized_z))) * tangential_relative_velocity_z;

                    double friction_restriction_x = 1.15 * abs(normal_force_x + normal_damping_force_x);
                    double friction_restriction_y = 1.15 * abs(normal_force_y + normal_damping_force_y);
                    double friction_restriction_z = 1.15 * abs(normal_force_z + normal_damping_force_z);

                    // double total_damping_force_x = tangential_damping_force_x + tangential_force_x < friction_restriction_x ? tangential_damping_force_x + tangential_force_x : friction_restriction_x;
                    // double total_damping_force_y = tangential_damping_force_y + tangential_force_y < friction_restriction_y ? tangential_damping_force_y + tangential_force_y : friction_restriction_y;
                    // double total_damping_force_z = tangential_damping_force_z + tangential_force_z < friction_restriction_z ? tangential_damping_force_z + tangential_force_z : friction_restriction_z;

                    double total_tangential_force_x = abs(tangential_force_x + tangential_damping_force_x) < friction_restriction_x ? tangential_force_x + tangential_damping_force_x : friction_restriction_x * sign(tangential_force_x + tangential_damping_force_x);
                    double total_tangential_force_y = abs(tangential_force_y + tangential_damping_force_x) < friction_restriction_y ? tangential_force_y + tangential_damping_force_y : friction_restriction_y * sign(tangential_force_y + tangential_damping_force_y);
                    double total_tangential_force_z = abs(tangential_force_z + tangential_damping_force_z) < friction_restriction_z ? tangential_force_z + tangential_damping_force_z : friction_restriction_z * sign(tangential_force_z + tangential_damping_force_z);
                    
                    // Update forces on particle
                    _force_i.x += normal_force_x + total_tangential_force_x + normal_damping_force_x;
                    _force_i.y += normal_force_y + total_tangential_force_y + normal_damping_force_y;
                    _force_i.z += normal_force_z + total_tangential_force_z + normal_damping_force_z;
                }
            }
            if(_particle.flag == 1)
            {
            update_acc(particles[index], _force_i); // a(t + 1/2 dt)
            }
        }
    }
}

__global__ void position_update(Molecule *particles, int num_particles, double time_step, size_t box_size) { // spring constant and damping coefficient

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles) {
        Molecule _particle = particles[index];
        if (_particle.flag != 0) {
            double time_step_sqr = time_step * time_step;
            double particle_diameter = _particle.diameter;

            // Update particle position
            _particle.x += _particle.xv * time_step + (_particle.xa * time_step_sqr) / 2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
            _particle.y += _particle.yv * time_step + (_particle.ya * time_step_sqr) / 2;
            _particle.z += _particle.zv * time_step + (_particle.za * time_step_sqr) / 2;

            double k = 10.0;
            double b = 0.95;
            // Wall collisions
            // X direction
            if (_particle.x < particle_diameter || _particle.x > box_size - particle_diameter) {
                double overlap = _particle.x < particle_diameter ? particle_diameter - _particle.x : _particle.x - (box_size - particle_diameter);
                double relative_velocity = _particle.xv;
                double spring_force = k * overlap;
                double damping_force = b * relative_velocity;
                double total_force = spring_force - damping_force;
                double acceleration = total_force / _particle.mass; 

                // Update velocity due to collision
                _particle.xv += (_particle.x < particle_diameter ? acceleration : -acceleration) * time_step;
                
                // Reflect position to simulate bounce
                if (_particle.x < 0) {
                    _particle.x = 0;
                } 
                if (_particle.x > box_size)
                {
                    _particle.x = box_size;
                }
            }

            // Y direction
            if (_particle.y < particle_diameter || _particle.y > box_size - particle_diameter) {
                double overlap = _particle.y < particle_diameter ? particle_diameter - _particle.y : _particle.y - (box_size - particle_diameter);
                double relative_velocity = _particle.yv;
                double spring_force = k * overlap;
                double damping_force = b * relative_velocity;
                double total_force = spring_force - damping_force;
                double acceleration = total_force / _particle.mass;

                // Update velocity due to collision
                _particle.yv += (_particle.y < particle_diameter ? acceleration : -acceleration) * time_step;
                
                // Reflect position to simulate bounce
                if (_particle.y < 0) {
                    _particle.y = 0;
                } 
                if (_particle.y > box_size)
                {
                    _particle.y = box_size;
                }
            }

            // Z direction
            if (_particle.z < particle_diameter || _particle.z > box_size - particle_diameter) {
                double overlap = _particle.z < particle_diameter ? particle_diameter - _particle.z : _particle.z - (box_size - particle_diameter);
                double relative_velocity = _particle.zv;
                double spring_force = k * overlap;
                double damping_force = b*relative_velocity;
                double total_force = spring_force - damping_force;
                double acceleration = total_force / _particle.mass;

                // Update velocity due to collision
                _particle.zv += (_particle.z < particle_diameter ? acceleration : -acceleration) * time_step;
                
                // Reflect position to simulate bounce
                if (_particle.z < 0) {
                    _particle.z = 0;
                } 
                if (_particle.z > box_size)
                {
                    _particle.z = box_size;
                }
            }

            // Write updated particle back to global memory
            particles[index] = _particle;
        }
    }
}

__global__ void velocity_update(Molecule *particle, int num_particles,
                                     double time_step, curandState *my_curandstate){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles){
        Molecule _particle = particle[index];
        if (_particle.flag != 0){
            _particle.xv += _particle.xa * time_step;
            if (abs(_particle.xv) < 0.01) _particle.xv +=  (curand_uniform(my_curandstate+index)) - 0.5;
            _particle.yv += _particle.ya * time_step;
            _particle.zv += _particle.za * time_step;

            particle[index] = _particle;
        }
    }
}

__global__ void setup_kernel(curandState *state){

  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  curand_init(1234, idx, 0, &state[idx]);
}

int main(int argc, char *argv[]) {
    // if (argc < 8) {
    //     std::cerr << "Usage: " << argv[0] << " <time_step> <num_steps> <k> <b> <particle_datafile> <cutoff_dist> <box_size>" << std::endl;
    //     return -1;
    // }

    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    std::string particle_datafile = argv[3];
    size_t box_size = atoi(argv[4]);
    double poisson_ratio = atof(argv[5]);
    double young_modulus = atof(argv[6]) * pow(10, 4);

    double cell_length_mulltiplier = 3.0;

    double eff_young_modulus = young_modulus / (2 - 2 * poisson_ratio * poisson_ratio);
    double eff_shear_modulus = young_modulus / (4 + 4 * poisson_ratio);
    std::cout << eff_young_modulus << "  " << eff_shear_modulus << "\n";

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

    double max_particle_size = fill_particles(particles, num_particles, cell_length_mulltiplier, box_size, file);

    writeBoxVTK(box_size);  // Write the box boundary to a VTK file

    int NUM_THREAD = 256;
    int NUM_BLOCK = (num_particles + NUM_THREAD - 1) / NUM_THREAD;

    auto start_global = std::chrono::steady_clock::now();
    cudaError_t syncErr, asyncErr;

    curandState *d_state;
    cudaMalloc(&d_state, sizeof(curandState));
    setup_kernel<<<NUM_BLOCK, NUM_THREAD>>>(d_state);

    for (int i = 0; i < num_steps; i++) {
        //auto start = std::chrono::steady_clock::now();

        position_update<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, box_size);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, box_size, eff_young_modulus, eff_shear_modulus);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

        velocity_update<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, d_state);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

        if (i % 10 == 0) {
            // auto end = std::chrono::steady_clock::now();
            // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            // std::cout << "Iteration took " << time.count() << " milliseconds\n";
            writeVTK(i, num_particles, particles);
        }
    }

    auto end_global = std::chrono::steady_clock::now();
    auto time_global = std::chrono::duration_cast<std::chrono::seconds>(end_global - start_global);
    std::cout << "Total iteration time " << time_global.count() << " seconds\n";

    cudaFree(particles);
}



