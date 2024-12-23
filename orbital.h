#ifndef ORBITAL_H
#define ORBITAL_H

#include <vector>
#include <cmath>
#include <iostream>

class Orbital {
public:
    Orbital(double exponent, const std::vector<double>& coefficients);
    double evaluate(double x, double y, double z) const;
    void print() const;

private:
    double exponent;
    std::vector<double> coefficients;
};

#endif // ORBITAL_H