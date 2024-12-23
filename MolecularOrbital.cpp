#include "MolecularOrbital.h"

MolecularOrbital::MolecularOrbital(const Orbital& ao1, const Orbital& ao2)
    : ao1(ao1), ao2(ao2) {}

double MolecularOrbital::evaluate(double x, double y, double z) const {
    return ao1.evaluate(x, y, z) + ao2.evaluate(x, y, z);
}

void MolecularOrbital::print() const {
    std::cout << "Molecular Orbital formed by: " << std::endl;
    ao1.print();
    ao2.print();
}
