// Champ auto-cohérent Hartree–Fock (RHF / UHF) avec accélération DIIS.
#pragma once
#include <array>
#include <string>
#include <vector>

#include "basis.h"
#include "integrals.h"
#include "linalg.h"
#include "molecule.h"

namespace mo {

struct SCFOptions {
    bool uhf = false;
    int maxIter = 200;
    double energyTol = 1e-10;
    double errorTol = 1e-7;
    bool breakSymmetry = true;  // UHF singulet : mélange HOMO/LUMO bêta pour casser la symétrie
    bool verbose = false;
};

struct SCFResult {
    bool converged = false;
    bool uhf = false;
    int iterations = 0;
    int nAlpha = 0, nBeta = 0;
    double energy = 0, eNuc = 0, eElec = 0;
    std::vector<double> epsA, epsB;  // énergies orbitalaires (hartree)
    Mat CA, CB;                      // coefficients (colonnes = orbitales)
    Mat PA, PB;                      // matrices densité de spin
    Mat S, H;
    std::vector<double> charges;       // charges de Mulliken par atome
    std::vector<double> spinPop;       // populations de spin (UHF)
    Mat bondOrders;                    // indices de liaison de Mayer
    std::array<double, 3> dipole{};    // bohr·e → u.a.
    double s2 = 0;                     // <S²> (UHF)
    std::vector<double> energyHistory;
};

SCFResult runSCF(const Molecule& mol, const Basis& basis, const SCFOptions& opt);

}  // namespace mo
