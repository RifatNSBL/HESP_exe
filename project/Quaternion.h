#ifndef QUATERNOIN_H
#define QUATERNOIN_H

#include <math.h> 
#include <iostream>
struct Quaternion{
    double q0, q1, q2, q3;

    Quaternion(double alpha, double axis[3]){ // Constructor given basic definition
        q0 = cos(alpha / 2.0);
        q1 = axis[0] * sin(alpha / 2);
        q2 = axis[1] * sin(alpha / 2);
        q3 = axis[2] * sin(alpha / 2);
    }

    Quaternion(double phi, double theta, double psi){ // constructor for giver Euler angles
        q0 = cos(theta / 2.0) * cos( (phi + psi) / 2);
        q1 = sin(theta / 2.0) * sin( (phi - psi) / 2);
        q2 = sin(theta / 2.0) * cos( (phi - psi) / 2);
        q3 = cos(theta / 2.0) * sin( (phi + psi) / 2);
    }

    Quaternion(double q0, double q1, double q2, double q3): // copy costructor for multiplication
        q0(q0), q1(q1), q2(q2), q3(q3){}

    Quaternion(): // default quaternion
        q0(1), q1(0), q2(0), q3(0){}

    bool norm_is_unit() {
        double norm = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        return norm == 1;
    }
};

Quaternion conjugate(Quaternion q){
    return Quaternion(q.q0, -q.q1, -q.q2, -q.q3);
} 

Quaternion multiply(Quaternion l, Quaternion r){ // operation is not commutative (q1*q2 != q2*q1)
    double q0 = l.q0 * r.q0 - l.q1 * r.q1 - l.q2 * r.q2 - l.q3 * r.q3;
    double q1 = l.q0 * r.q1 + l.q1 * r.q0 + l.q2 * r.q3 - l.q3 * r.q2;
    double q2 = l.q0 * r.q2 - l.q1 * r.q3 + l.q2 * r.q0 + l.q3 * r.q1;
    double q3 = l.q0 * r.q3 + l.q1 * r.q2 - l.q2 * r.q1 + l.q3 * r.q0;
    return Quaternion(q0, q1, q2, q3);
}

Quaternion rotate(Quaternion q, Quaternion curr_orient){
    Quaternion conj_q = conjugate(q);
    Quaternion temp = multiply(q, curr_orient);
    std::cout << temp.q0 << " " << temp.q1 << " " << temp.q2 << " " << temp.q3 << "\n";
    Quaternion result = multiply(temp, conj_q);
    std::cout << result.q0 << " " << result.q1 << " " << result.q2 << " " << result.q3 << "\n";
    return Quaternion(result.q0, result.q1, result.q2, result.q3);
}

#endif