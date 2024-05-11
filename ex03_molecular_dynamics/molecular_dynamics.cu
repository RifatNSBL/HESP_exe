#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

struct molecule
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
    double x;
    double y;
    double z;
};

std::vector<molecule> fill_points(int, std::ifstream &);
void get_values(std::vector<double> &, std::string);
void writeVTK(const char* name, int num_points, std::vector<molecule> points);


int main(int argc, char *argv[]){
    std::ifstream file;
    file.open("points.txt");

    double time_step = atoi(argv[1]);
    int num_steps = atoi(argv[2]);
    double delta = atoi(argv[3]);
    double eps = atoi(argv[4]);

    int num_points = 0;

    if (file.is_open()){
        std::string line;
        std::getline(file, line);
        num_points = std::stoi(line);    
    }
    std::vector<molecule> grid = fill_points(num_points, file);   
    file.close();

    // here do steps



    writeVTK("output.vtk", num_points, grid);
}

std::vector<molecule> fill_points(int num_points, std::ifstream &file){
    std::vector<molecule> grid(num_points);
    for(auto i = 0; i < num_points; i++){
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
    return grid;
}

void get_values(std::vector<double> &values, std::string line){
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        values.push_back(val);
    }
    std::cout << std::endl;

}

void writeVTK(const char* name, int num_points, std::vector<molecule> points){
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
    vtk_output << "POINTS " <<  num_points << " double\n";
    for(int i = 0; i < num_points; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].x
                   << " " << points[i].y << " " << points[i].z << "\n";
    }
    vtk_output << "CELLS 0 0\n";
    vtk_output << "CELL_TYPES 0\n";
    vtk_output << "POINT_DATA " << num_points << "\n";
    vtk_output << "SCALARS m double\n";
    vtk_output << "LOOKUP_TABLE default\n";
    for(int i = 0; i < num_points; i++)
        vtk_output << std::fixed << std::setprecision(6) << 1.0 << "\n";
        
    vtk_output << "VECTORS v double\n";

    for(int i = 0; i < num_points; i++){
        vtk_output << std::fixed << std::setprecision(6) << points[i].xv
                   << " " << points[i].yv << " " << points[i].zv << "\n";
    }
    vtk_output.close();
}