#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip> // For std::fixed and std::setprecision

void writeLineToFile(std::ofstream& file, const std::vector<double>& line) {

    // Write each double in the vector to the file, separated by spaces
    for (const double& value : line) {
        file << std::fixed << std::setprecision(6) << value << " ";
    }

    // Write a newline at the end of the line
    file << std::endl;

}

int main() {
    // Example usage
    std::string filePath = "points.txt";
    std::ofstream file(filePath);
    file << 2000 << "\n"; // change this if use different number of particles
    for(auto x = 0; x < 10; x++){
        for(auto y = 0; y < 10; y++){
            for(auto z = 0; z < 10; z++){
                std::vector<double> line = {1, static_cast<double>(x), static_cast<double>(y), 
                                            static_cast<double>(z), -1, 0, 0, 0, 0, 0};
                writeLineToFile(file, line);
                }
        }
    }

    for(auto x = -15; x < -5; x++){
        for(auto y = 0; y < 10; y++){
            for(auto z = 0; z < 10; z++){
                std::vector<double> line = {3, static_cast<double>(x), static_cast<double>(y), 
                                            static_cast<double>(z), 1, 0, 0, 0, 0, 0};
                writeLineToFile(file, line);
                }
        }
    }
    file.close();
    return 0;
}