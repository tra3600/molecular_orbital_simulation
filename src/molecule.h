// Atomes, molécules, bibliothèque de géométries et lecture de fichiers XYZ.
#pragma once
#include <array>
#include <string>
#include <vector>

namespace mo {

constexpr double ANGSTROM_TO_BOHR = 1.8897259886;
constexpr double HARTREE_TO_EV = 27.211386246;
constexpr double AU_TO_DEBYE = 2.541746473;

struct Atom {
    int Z = 0;
    double x = 0, y = 0, z = 0;  // bohr
};

struct Molecule {
    std::string name;
    std::string description;
    std::vector<Atom> atoms;
    int charge = 0;
    int multiplicity = 1;

    int nElectrons() const;
    double nuclearRepulsion() const;
    std::string formula() const;
    double distance(int a, int b) const;  // bohr
};

std::string elementSymbol(int Z);
int atomicNumber(const std::string& symbol);  // -1 si inconnu
double covalentRadiusAngstrom(int Z);

// Bibliothèque de molécules prédéfinies.
std::vector<Molecule> builtinMolecules();
bool findBuiltin(const std::string& name, Molecule& out);

// Lecture d'un fichier XYZ (Å). Les lignes de commentaire peuvent contenir
// "charge=<n>" et "mult=<n>".
Molecule readXYZ(const std::string& path);
// Analyse d'une description "O 0 0 0; H 0.75 0.58 0" (Å).
Molecule parseAtoms(const std::string& spec);

}  // namespace mo
