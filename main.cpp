#include <iostream>
#include "Orbital.h"
#include "MolecularOrbital.h"

int main() {
    // Définir deux orbitales atomiques
    Orbital ao1(1.0, {1.0, 0.0, 0.0, 0.0});  // Orbital 1s
    Orbital ao2(1.0, {0.0, 1.0, 0.0, 0.0});  // Orbital 2px

    // Créer une orbitale moléculaire en combinant les deux orbitales atomiques
    MolecularOrbital mo(ao1, ao2);

    // Afficher les orbitales atomiques et l'orbitale moléculaire
    ao1.print();
    ao2.print();
    mo.print();

    // Évaluer l'orbitale moléculaire en un point donné
    double x = 1.0, y = 0.0, z = 0.0;
    std::cout << "Molecular Orbital value at (" << x << ", " << y << ", " << z << "): " 
              << mo.evaluate(x, y, z) << std::endl;

    return 0;
}