#ifndef MOLECULE_H
#define MOLECULE_H
struct Molecule
{
    double mass;
    double x, y, z;
    double xv, yv, zv;
    double xa, ya, za;
    int cell_id;
};

#endif