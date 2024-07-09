#ifndef MOLECULE_H
#define MOLECULE_H
#include "Quaternion.h"
struct Molecule
{
    double mass;
    double x, y, z;
    double xv, yv, zv;
    double xa, ya, za;
    double diameter;  // Add size attribute
    int cell_id;

    // Quaternion orientation;
    
    double theta; 
    double omega;
    double alpha;
    double inertia;
};
#endif
