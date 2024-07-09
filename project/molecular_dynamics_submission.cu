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
    _particle.alpha = _force.torque / _particle.inertia;

    particle = _particle;
    force    = _force;
}

__global__ void calculate_forces(Molecule *particles, size_t num_particles,
                                 double k, double b, double mu, size_t box_size){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles) {
        Molecule _particle = particles[index];

        Force _force_i;
        _force_i.x = 0;
        _force_i.y = -9.81;
        _force_i.z = 0;

        //Rotational
        _force_i.torque = 0.0;

        for(int i = 0; i < num_particles; i++) {
            if (i == index) continue; 
            Molecule _particle_i = particles[i];
            double x_proj = _particle.x - _particle_i.x;
            double y_proj = _particle.y - _particle_i.y;
            double z_proj = _particle.z - _particle_i.z;

            double dist_sqr = x_proj * x_proj + y_proj * y_proj + z_proj * z_proj;
            double dist = sqrt(dist_sqr);

            double combined_radius = (_particle.diameter + _particle_i.diameter) / 2;
            if (dist < combined_radius) { // Check if particles are overlapping
                double overlap = combined_radius - dist;

                // Normalize the projection to get the normal direction
                double normal_x = x_proj / dist;
                double normal_y = y_proj / dist;
                double normal_z = z_proj / dist;

                double spring_force_magnitude = k * overlap;
                double spring_force_x = spring_force_magnitude * normal_x;
                double spring_force_y = spring_force_magnitude * normal_y;
                double spring_force_z = spring_force_magnitude * normal_z;

                                    // Calculate damping force components
                double relative_velocity_x = _particle.xv - _particle_i.xv;
                double relative_velocity_y = _particle.yv - _particle_i.yv;
                double relative_velocity_z = _particle.zv - _particle_i.zv;
                double relative_velocity_dot_proj = (relative_velocity_x * normal_x + relative_velocity_y * normal_y + relative_velocity_z * normal_z);

                double damping_force_magnitude = b * relative_velocity_dot_proj;
                double damping_force_x = damping_force_magnitude * normal_x;
                double damping_force_y = damping_force_magnitude * normal_y;
                double damping_force_z = damping_force_magnitude * normal_z;

                                    // Calculate total force components
                double total_force_x = spring_force_x - damping_force_x;
                double total_force_y = spring_force_y - damping_force_y;
                double total_force_z = spring_force_z - damping_force_z;

                //Calculate tangential force and torque
                double tangential_velocity_x = relative_velocity_x - relative_velocity_dot_proj * normal_x;
                double tangential_velocity_y = relative_velocity_y - relative_velocity_dot_proj * normal_y;
               
                double tangential_velocity_magnitude = sqrt(tangential_velocity_x * tangential_velocity_x +
                                                            tangential_velocity_y * tangential_velocity_y);

                double tangential_force_magnitude = mu * sqrt(total_force_x * total_force_x + total_force_y * total_force_y);
                double tangential_force_x = tangential_force_magnitude * (tangential_velocity_magnitude > 0 ? tangential_velocity_x / tangential_velocity_magnitude : 0);
                double tangential_force_y = tangential_force_magnitude * (tangential_velocity_magnitude > 0 ? tangential_velocity_y / tangential_velocity_magnitude : 0);

                // Compute the torque due to the tangential force
                double torque = (_particle.diameter / 2) * (tangential_force_y * normal_x - tangential_force_x * normal_y);

                // Update torque on particle
                _force_i.torque += torque;

                                    // Update forces on particle
                _force_i.x += total_force_x;
                _force_i.y += total_force_y;
                _force_i.z += total_force_z;

                //Update torque on particle
                _force_i.torque+= torque;


            }
        }
        update_acc(particles[index], _force_i); // a(t + 1/2 dt)
    }
}

__global__ void position_update(Molecule *particles, int num_particles, double time_step, size_t box_size, double k, double b, double mu) { // spring constant and damping coefficient

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles) {
        double time_step_sqr = time_step * time_step;
        Molecule _particle = particles[index];
        double particle_diameter = _particle.diameter;

        // Update particle position
        _particle.x += _particle.xv * time_step + (_particle.xa * time_step_sqr) / 2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        _particle.y += _particle.yv * time_step + (_particle.ya * time_step_sqr) / 2;
        _particle.z += _particle.zv * time_step + (_particle.za * time_step_sqr) / 2;

        // Update particle angular_position
        _particle.theta = _particle.omega * time_step + (_particle.alpha* time_step_sqr) / 2;

        //Check for angle wrap
        _particle.theta = fmod(_particle.theta, 360.0);  // Wrap-around to [0, 360)
        if (_particle.theta < 0) {
        _particle.theta += 360.0;  // Adjust negative values to be within [0, 360)
        }

        // Wall collisions
        // X direction
        if (_particle.x < particle_diameter || _particle.x > box_size - particle_diameter) {
            double overlap = _particle.x < particle_diameter ? particle_diameter - _particle.x : _particle.x - (box_size - particle_diameter);
            double relative_velocity = _particle.xv;
            double spring_force = k * overlap;
            double damping_force = b * relative_velocity;
            double total_force = spring_force - damping_force;
            double acceleration = total_force / _particle.mass; 

            //Angular Components
            double frictional_force = total_force * mu;
            double torque = frictional_force * particle_diameter/2 * (_particle.xv >= 0 ? -1 : 1);  // Negative for positive y-velocity, positive for negative y-velocity
            double alpha = torque /  _particle.inertia;

            // Update velocity due to collision
            _particle.xv += (_particle.x < particle_diameter ? acceleration : -acceleration) * time_step;

             // Update angular velocity due to collision
            _particle.omega += alpha * time_step;
            
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

             //Angular Components
            double frictional_force = total_force * mu;
            double torque = frictional_force * particle_diameter/2 * (_particle.yv >= 0 ? -1 : 1);  // Negative for positive y-velocity, positive for negative y-velocity
            double alpha = torque /  _particle.inertia;

            // Update velocity due to collision
            _particle.yv += (_particle.y < particle_diameter ? acceleration : -acceleration) * time_step;

             // Update angular velocity due to collision
            _particle.omega += alpha * time_step;
            
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

__global__ void velocity_update(Molecule *particle, int num_particles,
                                     double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_particles){
        Molecule _particle = particle[index];

        _particle.xv += _particle.xa * time_step; 
        _particle.yv += _particle.ya * time_step;
        _particle.zv += _particle.za * time_step;

        particle[index] = _particle;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " <time_step> <num_steps> <k> <b> <particle_datafile> <cutoff_dist> <box_size>" << std::endl;
        return -1;
    }

    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    double k = atof(argv[3]);
    double b = atof(argv[4]);

    //Add friction-coefficient
    double mu = atof(argv[5]);

    std::string particle_datafile = argv[6];
    double cutoff_dist = atof(argv[7]);
    size_t box_size = atoi(argv[8]);
    double cell_length_mulltiplier = 3.0;

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

    for (int i = 0; i < num_steps; i++) {
        auto start = std::chrono::steady_clock::now();

        position_update<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, box_size, k, b, mu);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, k, b, mu, box_size);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

        velocity_update<<<NUM_BLOCK, NUM_THREAD>>>(particles, num_particles, time_step, k, b);
        syncErr = cudaGetLastError();
        asyncErr = cudaDeviceSynchronize();
        if (syncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(syncErr));
        if (asyncErr != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(asyncErr));

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



