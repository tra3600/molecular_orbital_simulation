// Base d'orbitales atomiques gaussiennes contractées (STO-3G, H–Ne).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "molecule.h"

namespace mo {

struct Primitive {
    double alpha;  // exposant
    double coef;   // coefficient incluant les normalisations
};

struct BasisFunction {
    int atom = 0;
    std::array<double, 3> R{};  // centre (bohr)
    std::array<int, 3> l{};     // (lx, ly, lz)
    std::vector<Primitive> prims;
    std::string label;  // ex. "O1 2px"
};

struct Basis {
    std::vector<BasisFunction> fns;
    int size() const { return static_cast<int>(fns.size()); }
};

Basis buildBasis(const Molecule& mol, const std::string& name = "sto-3g");
double evalBasisFunction(const BasisFunction& f, double x, double y, double z);

}  // namespace mo
