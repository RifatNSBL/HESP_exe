#ifndef UTIL_H
#define UTIL_H
#include "Molecule.h"
#include "Force.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

void get_values(std::vector<double> &values, std::string line){
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        values.push_back(val);
    }

}

void find_max(std::vector<double> values, double &max_x, double &max_y, double &max_z){
    max_x = values[1] > max_x ? values[1] : max_x;
    max_y = values[2] > max_y ? values[2] : max_y;
    max_z = values[3] > max_z ? values[3] : max_z;
}

void fill_particles(Molecule *grid, size_t num_molecules, std::ifstream &file, 
                    double &max_x, double &max_y, double &max_z){

    for(auto i = 0; i < num_molecules; i++){
        if(file.good()){
            std::string line;
            std::getline(file, line);
            
            std::vector<double> values;
            get_values(values, line);
            find_max(values, max_x, max_y, max_z);
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
        vtk_output << std::fixed << std::setprecision(6) << points[i].mass  << "\n";
        
    vtk_output << "VECTORS v double\n";

    for(int i = 0; i < num_molecules; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].xv
                   << " " << points[i].yv << " " << points[i].zv << "\n";
    }
    vtk_output.close();
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
#endif