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
    // Ensure the correct number of arguments are provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <cube_size> <ball_size>\n";
        return 1;
    }

    // Example usage
    double cube_size = std::atof(argv[1]) / 25.0; // Use atof to convert to double
    double ball_size = std::atof(argv[2]); // Use atof to convert to double

    double min = -1.0;
    double max = 1.0;

    // Debug prints
    std::cout << "Parsed cube_size: " << cube_size << std::endl;
    std::cout << "Parsed ball_size: " << ball_size << std::endl;

    std::string filePath = "points.txt";
    std::ofstream file(filePath);

    // Ensure the file is open
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filePath << "\n";
        return 1;
    }

    file << 3177 << "\n";

    for (double x = 0; x < 20 ; x++) {                                                           //pegs
        for (double y = 0; y < 47 ; y++) {
            std::vector<double> line = {1, 3.0 + 2.0 * y + std::fmod(x, 2.0), 12.0, 
                                        -std::pow(3.0, 0.5) * x + 50.0, 0, 0, 0, 0, 0, 0, ball_size, 0};
            writeLineToFile(file, line);
        }
    }   

    for (double x = 0; x < 14 ; x++) {                                                           //collecting wall
        for (double y = 0; y < 24 ; y++) {
            std::vector<double> line = {1, 2.0 + 4.0 * y , 12.0, 
                                        x + 2.05, 0, 0, 0, 0, 0, 0, ball_size, 0};
            writeLineToFile(file, line);
        }
    }


    for (double x = 0; x < 1 ; x++) {                                                           //collecting wall
        for (double y = 0; y < 336 ; y++) {
            std::vector<double> line = {1, 2.75 + 0.275*y, 12.0, 
                                        1.05, 0, 0, 0, 0, 0, 0, ball_size, 0};
            writeLineToFile(file, line);
        }
    }


    for (double x = 0; x < 170 ; x++) {                                                          //left handle
        for (double y = 0; y < 1 ; y++) {
            std::vector<double> line = {1, 45 -  0.275*x , 12.0, 
                                        55.0 + 0.275*x, 0, 0, 0, 0, 0, 0, ball_size, 0};
            writeLineToFile(file, line);
        }
    }

    for (double x = 0; x < 170 ; x++) {                                                          //right handle
        for (double y = 0; y < 1 ; y++) {
            std::vector<double> line = {1, 55 +  0.275*x , 12.0, 
                                        55.0 + 0.275*x, 0, 0, 0, 0, 0, 0, ball_size, 0};
            writeLineToFile(file, line);
        }
    }


/*
    for (double x = 0; x < 1 ; x++) {
        for (double y = 0; y < 10 ; y++) {
            double random_val = generateRandomValue();
            std::vector<double> line = {1, 11 + random_val, 12.0, 
                                        24 + 20 * x + 1.5 * y, 0, 0, 0, 0, 0, 0, ball_size, 1};
            writeLineToFile(file, line);
        }
    }
*/

    // for (double x = 0; x < 15 ; x++) {                                                          //left handle
    //     for (double y = 0; y < 1 ; y++) {
    //         std::vector<double> line = {1, 48 -  0.88*x , 12.0, 
    //                                     55.0 + 0.68*x, 0, 0, 0, 0, 0, 0, ball_size, 0};
    //         writeLineToFile(file, line);
    //     }
    // }

    // for (double x = 0; x < 15 ; x++) {                                                          //right handle
    //     for (double y = 0; y < 1 ; y++) {
    //         std::vector<double> line = {1, 51 +  0.88*x , 12.0, 
    //                                     55.0 + 0.68*x, 0, 0, 0, 0, 0, 0, ball_size, 0};
    //         writeLineToFile(file, line);
    //    }
    //}
    // before - 140


 /*   for (double x = 0; x < 30 ; x++) {                                                           //dropping balls
        for (double y = 0; y < 15 ; y++) {
            //double random_val = generateRandomNumber(min, max);
            std::vector<double> line = {1, 27.5 + x*1.5, 12.0, 
                                        75 + 1.1 * y , 0, 0, 0, 0, 0, 0, 0.85*ball_size, 1};
            writeLineToFile(file, line);
        }
    }
*/
    for (double x = 0; x < 35 ; x++) {                                                           //dropping balls
        for (double y = 0; y < 2*x + 1 ; y++) {
            //double random_val = generateRandomNumber(min, max);
            std::vector<double> line = {1, 50 - x*1.1 + y*1.1, 12.0, 
                                        56 + 1.1 * x , 0, 0, 0, 0, 0, 0, 0.85*ball_size, 1};
            writeLineToFile(file, line);
        }
    }


    file.close();
    return 0;
}

// int main(int argc, char *argv[]) {
//     // Example usage
//     size_t cube_size = atoi(argv[1]) / 6;
//     double ball_size = atoi(argv[2]);
//     std::string filePath = "points.txt";
//     std::ofstream file(filePath);
//     double min = -1.0;
//     double max = 1.0;
//     file << pow(atoi(argv[1]) / 6, 3) << "\n";
//     for(auto x = 0; x < cube_size; x++){
//         for(auto y = 0; y < cube_size; y++){
//             for(auto z = 0; z < cube_size; z++){
            
//                 std::vector<double> line = {1,  
//                                             3.0*static_cast<double>(x) + 3.0,  // x
//                                             3.0*static_cast<double>(y) + 3.0,  // y
//                                             3.0*static_cast<double>(z) + 3.0,  // z
//                                             0 ,                                // v_x
//                                             -0.5*y,                            // v_y
//                                             0,                                 // v_z
//                                             0,                                 // a_x
//                                             0,                                 // a_x
//                                             0, ball_size};
//                 writeLineToFile(file, line);
                
//         }
//         }
//     }

//     // for(auto x = 15; x < 25; x++){
//     //     for(auto y = 0; y < 10; y++){
//     //         for(auto z = 0; z < 10; z++){
//     //             std::vector<double> line = {3, static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5, 
//     //                                         static_cast<double>(z) + 0.5, 1, 0, 0, 0, 0, 0};
//     //             writeLineToFile(file, line);
//     //             }
//     //     }
//     // }
//     file.close();
//     return 0;
// }