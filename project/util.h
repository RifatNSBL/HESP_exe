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
#include <math.h>

void get_values(std::vector<double> &values, std::string line){
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        values.push_back(val);
    }
}

double fill_particles(Molecule *grid, size_t num_molecules, double cell_length_multiplier, double box_size, std::ifstream &file){
    double max_size = 0;
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
            grid[i].diameter = values[10];  // Add size
            grid[i].orientation = Quaternion(0.0, 1.0, 0.0, 0.0);
            //grid[i].orientation = Quaternion(0.0, 0.0, 90.0 * 2 * M_PI / 360);
            grid[i].inertia = 0.4 * grid[i].mass * grid[i].diameter * grid[i].diameter / 4;
            if (values[10] > max_size) max_size = values[0];

        }
        else{
            std::cout << "Error: file has reached the end\n";
            break; 
        }
    }
    file.close();
    int cells_per_side = box_size / (cell_length_multiplier * max_size);

    for(auto i = 0; i < num_molecules; i++){
            int cell_x = static_cast<int>(grid[i].x / cells_per_side);
            int cell_y = static_cast<int>(grid[i].y / cells_per_side);
            int cell_z = static_cast<int>(grid[i].z / cells_per_side);
            int flatten_cel_id = cell_x * cells_per_side * cells_per_side + cell_y * cells_per_side + cell_z;
            grid[i].cell_id = flatten_cel_id;
            // Debugging output
            printf("Particle %d: x=%f, y=%f, z=%f, cell_id=%d\n", i, grid[i].x, grid[i].y, grid[i].z, flatten_cel_id);
    }
    return max_size;
}



void writeBoxVTK(size_t box_size) {
    std::string name = "visualization/box.vtk";
    std::ofstream vtk_output(name);
    if (!vtk_output) {
        std::cerr << "Error: Unable to create box.vtk" << std::endl;
        return;
    }

    vtk_output << "# vtk DataFile Version 4.2\n";
    vtk_output << "Box boundary\n";
    vtk_output << "ASCII\n";
    vtk_output << "DATASET POLYDATA\n";

    // Define the 8 vertices of the box
    vtk_output << "POINTS 8 float\n";
    vtk_output << "0 0 0\n";
    vtk_output << box_size << " 0 0\n";
    vtk_output << box_size << " " << box_size << " 0\n";
    vtk_output << "0 " << box_size << " 0\n";
    vtk_output << "0 0 " << box_size << "\n";
    vtk_output << box_size << " 0 " << box_size << "\n";
    vtk_output << box_size << " " << box_size << " " << box_size << "\n";
    vtk_output << "0 " << box_size << " " << box_size << "\n";

    // Define the 12 edges of the box
    vtk_output << "LINES 12 36\n";
    vtk_output << "2 0 1\n";
    vtk_output << "2 1 2\n";
    vtk_output << "2 2 3\n";
    vtk_output << "2 3 0\n";
    vtk_output << "2 4 5\n";
    vtk_output << "2 5 6\n";
    vtk_output << "2 6 7\n";
    vtk_output << "2 7 4\n";
    vtk_output << "2 0 4\n";
    vtk_output << "2 1 5\n";
    vtk_output << "2 2 6\n";
    vtk_output << "2 3 7\n";

    vtk_output.close();
}


void writeVTK(int index, size_t num_molecules, Molecule *points) {
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
    vtk_output << "POINTS " << num_molecules << " double\n";
    for (size_t i = 0; i < num_molecules; i++) {
        vtk_output << std::fixed << std::setprecision(6) << points[i].x
                   << " " << points[i].y << " " << points[i].z << "\n";
    }
    vtk_output << "CELLS 0 0\n";
    vtk_output << "CELL_TYPES 0\n";
    vtk_output << "POINT_DATA " << num_molecules << "\n";

    vtk_output << "SCALARS size double\n";
    vtk_output << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < num_molecules; i++) {
        vtk_output << std::fixed << std::setprecision(6) << points[i].diameter << "\n";
    }

    vtk_output << "SCALARS mass double\n";
    vtk_output << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < num_molecules; i++) {
        vtk_output << std::fixed << std::setprecision(6) << points[i].mass << "\n";
    }
    vtk_output << "COLOR_SCALARS color 3\n";
    for (size_t i = 0; i < num_molecules; i++) {
        vtk_output << std::fixed << std::setprecision(6) << double(i) / (num_molecules + 1) + 0.000001
                   << " " << double(i) / (num_molecules + 1) + 0.000001 << " " << double(i) / (num_molecules + 1) + 0.000001 << "\n";
    }

    vtk_output << "VECTORS orientation double\n";
    // // for (size_t i = 0; i < num_molecules; i++) {
    // //     double x = sin(2 * M_PI * double(index) / 3000);
    // //     double y = cos(2 * M_PI * double(index) / 3000);
    // //     vtk_output << std::fixed << std::setprecision(6) << x
    // //         << " " << y << " " << 0.0 << "\n";
    // // }
    // double alpha = 45.0 * 2 * M_PI / 360.0; // angle from X in XY plane
    // double polar = 45.0 * 2 * M_PI / 360.0; // angle from XY plane

    // double y = sin(alpha) * cos(polar);
    // double x = cos(alpha) * cos(polar);
    // double z = sin(polar);
    
    vtk_output << std::fixed << std::setprecision(6) << points[0].orientation.q1
        << " " << points[0].orientation.q2 << " " << points[0].orientation.q3 << "\n";

    double rotation[3] = {0.0, 0.0, 1.0};

    points[1].orientation = rotate( Quaternion(90.0 * 2 * M_PI / 360, rotation), points[0].orientation );
    vtk_output << std::fixed << std::setprecision(6) << points[1].orientation.q1
        << " " << points[1].orientation.q2 << " " << points[1].orientation.q3 << "\n";

    vtk_output << "VECTORS velocity double\n";
    for (size_t i = 0; i < num_molecules; i++) {
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
