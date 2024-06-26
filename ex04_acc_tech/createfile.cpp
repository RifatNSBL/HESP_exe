#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip> // For std::fixed and std::setprecision
#include <random>
#include <ctime>

double generateRandomNumber(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

void writeLineToFile(std::ofstream& file, const std::vector<double>& line) {

    // Write each double in the vector to the file, separated by spaces
    for (const double& value : line) {
        file << std::fixed << std::setprecision(6) << value << " ";
    }

    // Write a newline at the end of the line
    file << std::endl;

}

int main(int argc, char *argv[]) {
    // Example usage
    size_t cube_size = atoi(argv[1]) / 2;
    std::string filePath = "points.txt";
    std::ofstream file(filePath);
    double min = -1.0;
    double max = 1.0;
    file << pow(cube_size, 3) << "\n";
    for(auto x = 0; x < cube_size; x++){
        for(auto y = 0; y < cube_size; y++){
            for(auto z = 0; z < cube_size; z++){
                std::vector<double> line = {1, 2* static_cast<double>(x) + 0.5, 2* static_cast<double>(y) + 0.5, 
                                            2* static_cast<double>(z) + 0.5, generateRandomNumber(min, max), 
                                            generateRandomNumber(min, max), generateRandomNumber(min, max), 0, 0, 0};
                writeLineToFile(file, line);
                }
        }
    }

    // for(auto x = 15; x < 25; x++){
    //     for(auto y = 0; y < 10; y++){
    //         for(auto z = 0; z < 10; z++){
    //             std::vector<double> line = {3, static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5, 
    //                                         static_cast<double>(z) + 0.5, 1, 0, 0, 0, 0, 0};
    //             writeLineToFile(file, line);
    //             }
    //     }
    // }
    file.close();
    return 0;
}