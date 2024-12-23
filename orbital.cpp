#include "Orbital.h"

Orbital::Orbital(double exponent, const std::vector<double>& coefficients)
    : exponent(exponent), coefficients(coefficients) {}

double Orbital::evaluate(double x, double y, double z) const {
    double r2 = x*x + y*y + z*z;
    return exp(-exponent * r2) * (coefficients[0] + coefficients[1] * x + coefficients[2] * y + coefficients[3] * z);
}

void Orbital::print() const {
    std::cout << "Orbital with exponent: " << exponent << " and coefficients: ";
    for (double coef : coefficients) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;
}