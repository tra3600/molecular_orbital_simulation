#ifndef MOLECULARORBITAL_H
#define MOLECULARORBITAL_H

#include "Orbital.h"
#include <vector>

class MolecularOrbital {
public:
    MolecularOrbital(const Orbital& ao1, const Orbital& ao2);
    double evaluate(double x, double y, double z) const;
    void print() const;

private:
    Orbital ao1;
    Orbital ao2;
};

#endif // MOLECULARORBITAL_H