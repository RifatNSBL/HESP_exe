#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

struct Molecule
{
    double mass;
    double x, y, z;
    double xv, yv, zv;
    double xa, ya, za;
};

struct Force
{
    double x;
    double y;
    double z;
};

void fill_particles(Molecule*, size_t, std::ifstream &);
void get_values(std::vector<double> &, std::string);
void writeVTK(int, size_t, Molecule *);

__device__ void calculate_forces(Molecule *particle, Force *force, size_t index, size_t num_molecules,
                                 double sigma, double eps){

        for(auto i = 0; i < num_molecules; i++){
            if (i == index) continue;
            double x_proj = particle[index].x - particle[i].x;
            double y_proj = particle[index].y - particle[i].y;
            double z_proj = particle[index].z - particle[i].z;
            
            double abs_dist = sqrt(pow(x_proj, 2) + pow(y_proj, 2) + pow(z_proj, 2));
            double factor = pow(sigma / abs_dist, 6);
            double res = 24 * eps * factor * (2 * factor - 1) / pow(abs_dist, 2);
            
            // x
            force[index].x += res * x_proj;
            //y
            force[index].y += res * y_proj;
            //z
            force[index].z += res * z_proj;
        }
}

__device__ void update_acc(Molecule &particle, Force &force){
    particle.xa = force.x / particle.mass;
    particle.ya = force.y / particle.mass;
    particle.za = force.z / particle.mass;
}

__device__ void reset_forces(Force &force){
    force.x = 0;
    force.y = 0;
    force.z = 0;
}

__global__ void integration_step(Molecule *particle, Force *forces, int num_molecules, 
                                 double time_step, double sigma, double eps){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        //calculate_forces(particle, forces, num_molecules, sigma, eps, index);
        particle[index].x += particle[index].xv * time_step + ( particle[index].xa * pow(time_step, 2) )/2; // x(t + dt) = x(t) + v(t) * dt + a(t) * dt² / 2
        particle[index].y += particle[index].yv * time_step + ( particle[index].ya * pow(time_step, 2) )/2; 
        particle[index].z += particle[index].zv * time_step + ( particle[index].za * pow(time_step, 2) )/2;

        particle[index].xv += particle[index].xa * time_step / 2; // v(t + 1/2 dt) = v(t) + a(t) * dt / 2
        particle[index].yv += particle[index].ya * time_step / 2;
        particle[index].zv += particle[index].za * time_step / 2;
        calculate_forces(particle, forces, index, num_molecules, sigma, eps);

        // calculation of forces  
        update_acc(particle[index], forces[index]); // a(t + 1/2 dt)

        particle[index].xv += particle[index].xa * time_step / 2; // v(t + dt) = v(t + 1/2 dt) + a(t + dt) * dt / 2
        particle[index].yv += particle[index].ya * time_step / 2;
        particle[index].zv += particle[index].za * time_step / 2;
        reset_forces(forces[index]);
        }
}


void calculate_energy(Molecule *points, size_t num_points, double sigma, double eps){
    double energy = 0.0;
    for(int i = 0; i < num_points; i++){
        energy += 0.5 * (points[i].mass * (pow(points[i].xv, 2) + pow(points[i].yv, 2) + pow(points[i].zv, 2)));
        for(int j = i + 1; j < num_points; j++){
            double x_proj = points[i].x - points[j].x;
            double y_proj = points[i].y - points[j].y;
            double z_proj = points[i].z - points[j].z;
            double abs_dist = sqrt(pow(x_proj, 2) + pow(y_proj, 2) + pow(z_proj, 2));
            double factor = pow(sigma / abs_dist, 6);
            double V = 4 * eps * ( pow(factor, 2) - factor );
            energy += V;
        }
    }
    std::cout << energy << "\n";
}

int main(int argc, char *argv[]){
    
    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    double sigma = atof(argv[3]); // stable distance = abs_dist * pow(0.5, 1/6)
    double eps = atof(argv[4]);
    std::string name = argv[5];

    size_t num_molecules = 0;
    
    std::ifstream file;
    file.open(name);

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

    fill_particles(grid_host, num_molecules, file);

    cudaMemcpy(grid_device, grid_host, size_molecule, cudaMemcpyHostToDevice);

    int NUM_THREAD = num_molecules;
    int NUM_BLOCK = (num_molecules + NUM_THREAD - 1) / NUM_THREAD;

    for(int i = 0; i < num_steps; i++){
        integration_step<<<NUM_BLOCK, NUM_THREAD>>>(grid_device, forces_device, num_molecules, time_step, sigma, eps);
        cudaDeviceSynchronize();

        if (i % 100 == 0){
            cudaMemcpy(grid_host, grid_device, size_molecule, cudaMemcpyDeviceToHost);
            calculate_energy(grid_host, num_molecules, sigma, eps);
            writeVTK(i, num_molecules, grid_host);
        }
    }

    cudaFree(grid_device);
    cudaFree(forces_device);
    cudaFreeHost(grid_host);
}

void fill_particles(Molecule *grid, size_t num_molecules, std::ifstream &file){

    for(auto i = 0; i < num_molecules; i++){
        if(file.good()){
            std::string line;
            std::getline(file, line);
            
            std::vector<double> values;
            get_values(values, line);
            
            grid[i].mass = values[0];
            grid[i].x = values[1];
            grid[i].y = values[2];
            grid[i].z = values[3];
            grid[i].xv = values[4];
            grid[i].yv = values[5];
            grid[i].zv = values[6];
            grid[i].xa = values[7];
            grid[i].ya = values[8];
            grid[i].za = values[9];
        }
        else{
            std::cout << "Error: file has reached the end\n";
            break; 
        }
    }
    file.close();
}

void get_values(std::vector<double> &values, std::string line){
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        values.push_back(val);
    }

}

void writeVTK(int index, size_t num_molecules, Molecule *points){
    std::string name = "visualization/output_" + std::to_string(index) + ".vtk";
    std::ofstream vtk_output(name);
    vtk_output.precision(6);
        if (!vtk_output) {
        std::cerr << "Error: Unable to create .vtk" << std::endl;
        return;
    }

    vtk_output << "# vtk DataFile Version 4.2\n";
    vtk_output << "hesp visualization file\n";
    vtk_output << "ASCII\n";
    vtk_output << "DATASET UNSTRUCTURED_GRID\n";
    vtk_output << "POINTS " <<  num_molecules << " double\n";
    for(int i = 0; i < num_molecules; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].x
                   << " " << points[i].y << " " << points[i].z << "\n";
    }
    vtk_output << "CELLS 0 0\n";
    vtk_output << "CELL_TYPES 0\n";
    vtk_output << "POINT_DATA " << num_molecules << "\n";
    vtk_output << "SCALARS m double\n";
    vtk_output << "LOOKUP_TABLE default\n";
    for(int i = 0; i < num_molecules; i++)
        vtk_output << std::fixed << std::setprecision(6) << 1 << "\n";
        
    vtk_output << "VECTORS v double\n";

    for(int i = 0; i < num_molecules; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].xv
                   << " " << points[i].yv << " " << points[i].zv << "\n";
    }
    vtk_output.close();
}