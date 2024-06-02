#include <chrono>
#include <list>
#include <vector>
#include "Molecule.h"
#include "Force.h"
#include "util.h"
#include "neighbourList.cuh"

__device__ void update_acc(Molecule &particle, Force &force){
    Molecule _particle = particle;
    Force    _force = force;

    _particle.xa = _force.x / _particle.mass;
    _particle.ya = _force.y / _particle.mass;
    _particle.za = _force.z / _particle.mass;

    particle = _particle;
}

__device__ double min_dist(double particle_i, double particle_j, size_t box_size){
    if      (abs(particle_i - particle_j) < ((double)box_size / 2)) 
            return particle_i - particle_j;

    else if (particle_i < particle_j) 
            return particle_i + (box_size - particle_j);

    else    return particle_i - (box_size + particle_j);
}

__global__ void calculate_forces(Molecule *particle, size_t num_molecules,neighbourList* neighbour_list,
                                 double sigma, double eps, size_t box_size){
        int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < num_molecules){
            Molecule _particle = particle[index];

            Force _force_i;
            
            _force_i.x = 0;
            _force_i.y = 0;
            _force_i.z = 0;

            double sigma_sqr  = sigma * sigma;
            neighbourList _neighbour_list = neighbour_list[index];
            Node_list* temp = _neighbour_list.head;

            while (temp != NULL){
                Molecule _particle_i = particle[temp->id];
                
                double x_proj = min_dist(_particle.x , _particle_i.x, box_size);
                double y_proj = min_dist(_particle.y , _particle_i.y, box_size);
                double z_proj = min_dist(_particle.z , _particle_i.z, box_size);

                double dist_sqr = x_proj * x_proj  + y_proj * y_proj + z_proj * z_proj;

                double factor_sqr = sigma_sqr / dist_sqr;
                double factor_hex = factor_sqr * factor_sqr * factor_sqr;

                double res = 24 * eps * factor_hex * (2 * factor_hex - 1) / dist_sqr;

                _force_i.x +=  res * x_proj;
                _force_i.y +=  res * y_proj;
                _force_i.z +=  res * z_proj;

                temp = temp->next;
            }

            update_acc(particle[index], _force_i); // a(t + 1/2 dt)
        }
}

__global__ void  generate_neighbour_list(Molecule *particle, neighbourList* neighbour_list, 
                                   const size_t num_molecules, double sigma, size_t box_size){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        neighbourList _neighbour_list = neighbour_list[index];
        if (!isEmpty(&_neighbour_list)){
            freeList(&_neighbour_list);
            initNode_list(&_neighbour_list);
        }
        Molecule _particle = particle[index];

        double sigma_sqr  = sigma * sigma;
        double cutoff_sqr = 2.5 * 2.5 * sigma_sqr;

        for(auto i = 0; i < num_molecules; i++){
            if (i == index) continue;
            
            Molecule _particle_i = particle[i];
                    
            double x_proj = min_dist(_particle.x , _particle_i.x, box_size);
            double y_proj = min_dist(_particle.y , _particle_i.y, box_size);
            double z_proj = min_dist(_particle.z , _particle_i.z, box_size);

            double dist_sqr = x_proj * x_proj  + y_proj * y_proj + z_proj * z_proj;

            if (dist_sqr < cutoff_sqr){
                appendNode_list(&_neighbour_list, i);
                }
            }
        neighbour_list[index] = _neighbour_list;
        if (index == 999) printList(&neighbour_list[index]);
    }
}

__global__ void integration_step_begin(Molecule *particle, int num_molecules, 
                                       double time_step, double sigma, double eps, 
                                       size_t box_size){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        double time_step_sqr = time_step * time_step;
        Molecule _particle = particle[index];

        _particle.x += _particle.xv * time_step + (_particle.xa * time_step_sqr)/2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        _particle.y += _particle.yv * time_step + (_particle.ya * time_step_sqr)/2; 
        _particle.z += _particle.zv * time_step + (_particle.za * time_step_sqr)/2;

        // if particle flies out from left -> move it to right
        if (_particle.x < 0) _particle.x +=  box_size; // implicit type cast
        if (_particle.y < 0) _particle.y +=  box_size;
        if (_particle.z < 0) _particle.z +=  box_size;
        // if particle flies out from right -> move it to left
        if (_particle.x > box_size) _particle.x -=  box_size;
        if (_particle.y > box_size) _particle.y -=  box_size;
        if (_particle.z > box_size) _particle.z -=  box_size;

        _particle.xv += _particle.xa * time_step / 2; // v(t + 1/2 dt) = v(t) + a(t) * dt / 2
        _particle.yv += _particle.ya * time_step / 2;
        _particle.zv += _particle.za * time_step / 2;

        particle[index] = _particle;
        }
}

__global__ void integration_step_end(Molecule *particle, int num_molecules, 
                                     double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        Molecule _particle = particle[index];
        
        _particle.xv += _particle.xa * time_step / 2; // v(t + dt) = v(t + 1/2 dt) + a(t + dt) * dt / 2
        _particle.yv += _particle.ya * time_step / 2;
        _particle.zv += _particle.za * time_step / 2;

        particle[index] = _particle;
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

    std::ifstream file;
    file.open(patricle_datafile);

    if (file.is_open()){
        std::string line;
        std::getline(file, line);
        num_molecules = std::stoi(line);    
    }

    size_t size_molecule = num_molecules * sizeof(Molecule);

    Molecule *grid;
    size_t heapsize = sizeof(neighbourList) * num_molecules * num_molecules * sizeof(Node_list);
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapsize);

    cudaMallocManaged(&grid, size_molecule);

    fill_particles(grid, num_molecules, file);
    
    Node_list *last_nodes;
    size_t size_of_nodes = num_molecules * sizeof(Node_list);
    cudaMallocManaged(&last_nodes, size_of_nodes);

    neighbourList* neighbour_list = initializeArrayOfLists(num_molecules);
    size_t size_neighbour_list = sizeof(neighbour_list);
    cudaMallocManaged(&neighbour_list, size_neighbour_list);

    int NUM_THREAD = 512;
    int NUM_BLOCK = (num_molecules + NUM_THREAD - 1) / NUM_THREAD;
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < num_steps; i++){
        
        // most stuff should be done here
        integration_step_begin<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, time_step, 
                                                         sigma, eps, box_size);
        if (i % 10 == 0){ 
            generate_neighbour_list<<<NUM_BLOCK, NUM_THREAD>>>(grid, neighbour_list, num_molecules, 
                                                               sigma, box_size);
        }
        calculate_forces<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, neighbour_list, 
                                                   sigma, eps, box_size);
        integration_step_end<<<NUM_BLOCK, NUM_THREAD>>>(grid, num_molecules, time_step, sigma, eps);

        if (i % 10 == 0){
            cudaDeviceSynchronize();
            // calculate_energy(grid, num_molecules, sigma, eps);
            writeVTK(i, num_molecules, grid);
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::seconds>(end - start); // don't know if it's correct

    std::cout << "Iteration took " << time.count() << " seconds\n";
    freeArrayOfLists<<<NUM_BLOCK, NUM_THREAD>>>(neighbour_list, num_molecules);
    cudaFree(grid);
}
