// Évaluation des orbitales sur des grilles 3D / plans 2D, export .cube.
#pragma once
#include <array>
#include <string>
#include <vector>

#include "basis.h"
#include "linalg.h"
#include "molecule.h"

namespace mo {

struct Grid {
    std::array<double, 3> origin{};
    double step = 0.25;
    std::array<int, 3> n{};
    size_t size() const { return size_t(n[0]) * n[1] * n[2]; }
    std::array<double, 3> point(int i, int j, int k) const {
        return {origin[0] + i * step, origin[1] + j * step, origin[2] + k * step};
    }
};

Grid makeGrid(const Molecule& mol, double margin, double spacing, size_t maxPoints = 2500000);

// Valeurs de chaque colonne de C (orbitales) + densité (matrice P) aux points donnés.
struct FieldSet {
    std::vector<std::vector<float>> orbitals;  // [orbitale][point]
    std::vector<float> density;                // [point] (vide si P non fourni)
};
FieldSet evaluateFields(const Basis& basis, const std::vector<const Mat*>& coeffs, const std::vector<int>& orbIdx,
                        const Mat* density, const std::vector<std::array<double, 3>>& pts);
FieldSet evaluateOnGrid(const Basis& basis, const Mat& C, const std::vector<int>& orbIdx, const Mat* density,
                        const Grid& g);

void writeCube(const std::string& path, const Molecule& mol, const Grid& g, const std::vector<float>& values,
               const std::string& comment);

}  // namespace mo
