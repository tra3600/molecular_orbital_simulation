// Intégrales moléculaires sur gaussiennes cartésiennes (s et p) par la méthode
// de McMurchie–Davidson.
#pragma once
#include <array>
#include <vector>

#include "basis.h"
#include "linalg.h"
#include "molecule.h"

namespace mo {

double overlapPair(const BasisFunction& a, const BasisFunction& b);

Mat overlapMatrix(const Basis& b);
Mat kineticMatrix(const Basis& b);
Mat nuclearMatrix(const Basis& b, const Molecule& mol);
std::array<Mat, 3> dipoleMatrices(const Basis& b);  // <μ|r|ν>, origine (0,0,0)

// Intégrales biélectroniques (μν|λσ), stockage compact à 8 symétries.
class ERITensor {
public:
    explicit ERITensor(int n = 0) : n_(n) {
        size_t np = static_cast<size_t>(n) * (n + 1) / 2;
        v_.assign(np * (np + 1) / 2, 0.0);
    }
    static size_t pair(int i, int j) { return i >= j ? size_t(i) * (i + 1) / 2 + j : size_t(j) * (j + 1) / 2 + i; }
    static size_t quad(size_t ij, size_t kl) { return ij >= kl ? ij * (ij + 1) / 2 + kl : kl * (kl + 1) / 2 + ij; }
    double operator()(int i, int j, int k, int l) const { return v_[quad(pair(i, j), pair(k, l))]; }
    double& at(size_t ij, size_t kl) { return v_[quad(ij, kl)]; }
    int n() const { return n_; }

private:
    int n_;
    std::vector<double> v_;
};

ERITensor computeERI(const Basis& b);

// Fonction de Boys F_n(x) pour n = 0..nmax.
void boysFunction(int nmax, double x, double* F);

}  // namespace mo
