#include "basis.h"

#include <cmath>
#include <stdexcept>

#include "integrals.h"

namespace mo {

namespace {
// Développement universel STO-3G (ζ = 1) des orbitales de Slater en 3 gaussiennes.
const double A1S[3] = {2.227660584, 0.4057711562, 0.1098175104};
const double C1S[3] = {0.1543289673, 0.5353281423, 0.4446345422};
const double A2SP[3] = {0.9942027296, 0.2310313333, 0.0751386610};
const double C2S[3] = {-0.09996722919, 0.3995128261, 0.7001154689};
const double C2P[3] = {0.1559162750, 0.6076837186, 0.3919573931};

// Facteurs d'échelle ζ de Slater (Hehre, Stewart, Pople) : {ζ(1s), ζ(2sp)}.
const double ZETA[11][2] = {{0, 0},       {1.24, 0},     {1.69, 0},     {2.69, 0.80}, {3.68, 1.15}, {4.68, 1.50},
                            {5.67, 1.72}, {6.67, 1.95}, {7.66, 2.25}, {8.65, 2.55}, {9.64, 2.88}};

double dfact(int n) {  // (n)!! avec (-1)!! = 1
    double r = 1;
    for (int k = n; k > 1; k -= 2) r *= k;
    return r;
}

double primNorm(double a, int lx, int ly, int lz) {
    int L = lx + ly + lz;
    return std::pow(2 * a / M_PI, 0.75) * std::pow(4 * a, L / 2.0) /
           std::sqrt(dfact(2 * lx - 1) * dfact(2 * ly - 1) * dfact(2 * lz - 1));
}

void addFunction(Basis& b, const Molecule& mol, int atom, std::array<int, 3> l, const double* alpha0,
                 const double* c, double zeta, const char* tag) {
    BasisFunction f;
    f.atom = atom;
    f.R = {mol.atoms[atom].x, mol.atoms[atom].y, mol.atoms[atom].z};
    f.l = l;
    for (int k = 0; k < 3; ++k) {
        double a = alpha0[k] * zeta * zeta;
        f.prims.push_back({a, c[k] * primNorm(a, l[0], l[1], l[2])});
    }
    std::string suffix;
    if (l[0] + l[1] + l[2] == 1) suffix = l[0] ? "x" : (l[1] ? "y" : "z");
    f.label = elementSymbol(mol.atoms[atom].Z) + std::to_string(atom + 1) + " " + tag + suffix;
    // normalisation de la fonction contractée
    double s = overlapPair(f, f);
    double n = 1.0 / std::sqrt(s);
    for (Primitive& p : f.prims) p.coef *= n;
    b.fns.push_back(f);
}
}  // namespace

Basis buildBasis(const Molecule& mol, const std::string& name) {
    if (name != "sto-3g" && name != "STO-3G") throw std::runtime_error("Base non supportée : " + name + " (seule STO-3G est disponible)");
    Basis b;
    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        int Z = mol.atoms[i].Z;
        if (Z < 1 || Z > 10) throw std::runtime_error("Élément non supporté (H–Ne uniquement)");
        int a = static_cast<int>(i);
        addFunction(b, mol, a, {0, 0, 0}, A1S, C1S, ZETA[Z][0], "1s");
        if (Z >= 3) {
            addFunction(b, mol, a, {0, 0, 0}, A2SP, C2S, ZETA[Z][1], "2s");
            addFunction(b, mol, a, {1, 0, 0}, A2SP, C2P, ZETA[Z][1], "2p");
            addFunction(b, mol, a, {0, 1, 0}, A2SP, C2P, ZETA[Z][1], "2p");
            addFunction(b, mol, a, {0, 0, 1}, A2SP, C2P, ZETA[Z][1], "2p");
        }
    }
    return b;
}

double evalBasisFunction(const BasisFunction& f, double x, double y, double z) {
    double dx = x - f.R[0], dy = y - f.R[1], dz = z - f.R[2];
    double r2 = dx * dx + dy * dy + dz * dz;
    double ang = 1.0;
    for (int k = 0; k < f.l[0]; ++k) ang *= dx;
    for (int k = 0; k < f.l[1]; ++k) ang *= dy;
    for (int k = 0; k < f.l[2]; ++k) ang *= dz;
    double s = 0;
    for (const Primitive& p : f.prims) s += p.coef * std::exp(-p.alpha * r2);
    return ang * s;
}

}  // namespace mo
