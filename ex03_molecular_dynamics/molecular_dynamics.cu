#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

struct Molecule
{
    double mass;
    double x;
    double y;
    double z;
    double xv;
    double yv;
    double zv;
    double xa = 0;
    double ya = 0; // -9.8?
    double za = 0;
};

struct Force
{
    double x = 0;
    double y = 0;
    double z = 0;
};

void fill_points(Molecule*, int, std::ifstream &);
void get_values(std::vector<double> &, std::string);
void writeVTK(int index, int num_molecules, Molecule *points);

__global__ void change_data(Molecule *data, int num_molecules, double time_step){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        data[index].x += time_step;
        data[index].y += time_step;
        data[index].z += time_step;
        }
}


int main(int argc, char *argv[]){
    
    double time_step = atof(argv[1]);
    int num_steps = atoi(argv[2]);
    double delta = atof(argv[3]);
    double eps = atof(argv[4]);

     int num_molecules = 0;
    
    std::ifstream file;
    file.open("points.txt");

    if (file.is_open()){
        std::string line;
        std::getline(file, line);
        num_molecules = std::stoi(line);    
    }

    int size_molecule = num_molecules * sizeof(Molecule);
    int size_forces =  num_molecules * sizeof(Force);

    Molecule *grid_host, *grid_device;
    Force *forces_host, *forces_device;

    cudaMallocHost(&grid_host, size_molecule);
    cudaMalloc(&grid_device, size_molecule);

    cudaMallocHost(&forces_host, size_forces);
    cudaMalloc(&forces_device, size_forces);

    fill_points(grid_host, num_molecules, file);
    file.close();

    cudaMemcpy(grid_device, grid_host, size_molecule, cudaMemcpyHostToDevice);
    cudaMemcpy(forces_device, forces_host, size_forces, cudaMemcpyHostToDevice);



    int NUM_THREAD = 5;
    int NUM_BLOCK = 1;
    for(int i = 0; i < num_steps; i++){
        change_data<<<NUM_BLOCK, NUM_THREAD>>>(grid_device, num_molecules, time_step);
        cudaDeviceSynchronize();

        cudaMemcpy(grid_host, grid_device, size_molecule, cudaMemcpyDeviceToHost);
        cudaMemcpy(forces_host, forces_device, size_forces, cudaMemcpyDeviceToHost);
        if (i % 10000 == 0)
            writeVTK(i / 10000, num_molecules, grid_host);
    }


    cudaFree(grid_device);
    cudaFree(forces_device);
}

void fill_points(Molecule *grid, int num_molecules, std::ifstream &file){

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
        }
        else{
            std::cout << "Error: file has reached the end\n";
            break; 
        }
    }
}

void get_values(std::vector<double> &values, std::string line){
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        values.push_back(val);
    }
    std::cout << std::endl;

}

void writeVTK(int index, int num_molecules, Molecule *points){
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
        vtk_output << std::fixed << std::setprecision(6) << 1.0 << "\n";
        
    vtk_output << "VECTORS v double\n";

    for(int i = 0; i < num_molecules; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].xv
                   << " " << points[i].yv << " " << points[i].zv << "\n";
    }
    vtk_output.close();
}