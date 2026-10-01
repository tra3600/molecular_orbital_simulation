#include "grid.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

namespace mo {

Grid makeGrid(const Molecule& mol, double margin, double spacing, size_t maxPoints) {
    Grid g;
    double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9};
    for (const Atom& a : mol.atoms) {
        double c[3] = {a.x, a.y, a.z};
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], c[d]);
            hi[d] = std::max(hi[d], c[d]);
        }
    }
    for (;;) {
        for (int d = 0; d < 3; ++d) {
            g.origin[d] = lo[d] - margin;
            g.n[d] = std::max(8, int(std::ceil((hi[d] - lo[d] + 2 * margin) / spacing)) + 1);
        }
        g.step = spacing;
        if (g.size() <= maxPoints) break;
        spacing *= 1.1;
    }
    return g;
}

FieldSet evaluateFields(const Basis& basis, const std::vector<const Mat*>& coeffs, const std::vector<int>& orbIdx,
                        const Mat* density, const std::vector<std::array<double, 3>>& pts) {
    // « coeffs » : une matrice C par orbitale demandée (permet de mélanger α et β)
    const int n = basis.size();
    FieldSet fs;
    fs.orbitals.assign(orbIdx.size(), std::vector<float>(pts.size()));
    if (density) fs.density.assign(pts.size(), 0.f);
#pragma omp parallel
    {
        std::vector<double> phi(n);
#pragma omp for schedule(static)
        for (long long p = 0; p < static_cast<long long>(pts.size()); ++p) {
            for (int m = 0; m < n; ++m) phi[m] = evalBasisFunction(basis.fns[m], pts[p][0], pts[p][1], pts[p][2]);
            for (size_t o = 0; o < orbIdx.size(); ++o) {
                double s = 0;
                const Mat& C = *coeffs[o];
                for (int m = 0; m < n; ++m) s += C(m, orbIdx[o]) * phi[m];
                fs.orbitals[o][p] = static_cast<float>(s);
            }
            if (density) {
                double s = 0;
                for (int i = 0; i < n; ++i) {
                    double r = 0;
                    for (int j = 0; j < n; ++j) r += (*density)(i, j) * phi[j];
                    s += phi[i] * r;
                }
                fs.density[p] = static_cast<float>(s);
            }
        }
    }
    return fs;
}

FieldSet evaluateOnGrid(const Basis& basis, const Mat& C, const std::vector<int>& orbIdx, const Mat* density,
                        const Grid& g) {
    std::vector<std::array<double, 3>> pts;
    pts.reserve(g.size());
    for (int i = 0; i < g.n[0]; ++i)
        for (int j = 0; j < g.n[1]; ++j)
            for (int k = 0; k < g.n[2]; ++k) pts.push_back(g.point(i, j, k));
    std::vector<const Mat*> cs(orbIdx.size(), &C);
    return evaluateFields(basis, cs, orbIdx, density, pts);
}

void writeCube(const std::string& path, const Molecule& mol, const Grid& g, const std::vector<float>& v,
               const std::string& comment) {
    std::ofstream f(path);
    f << "molecular_orbital_simulation\n" << comment << "\n";
    char buf[128];
    std::snprintf(buf, sizeof buf, "%5d %12.6f %12.6f %12.6f\n", int(mol.atoms.size()), g.origin[0], g.origin[1],
                  g.origin[2]);
    f << buf;
    for (int d = 0; d < 3; ++d) {
        double vec[3] = {0, 0, 0};
        vec[d] = g.step;
        std::snprintf(buf, sizeof buf, "%5d %12.6f %12.6f %12.6f\n", g.n[d], vec[0], vec[1], vec[2]);
        f << buf;
    }
    for (const Atom& a : mol.atoms) {
        std::snprintf(buf, sizeof buf, "%5d %12.6f %12.6f %12.6f %12.6f\n", a.Z, double(a.Z), a.x, a.y, a.z);
        f << buf;
    }
    size_t idx = 0;
    for (int i = 0; i < g.n[0]; ++i)
        for (int j = 0; j < g.n[1]; ++j) {
            for (int k = 0; k < g.n[2]; ++k) {
                std::snprintf(buf, sizeof buf, " %12.5E", v[idx++]);
                f << buf;
                if (k % 6 == 5) f << "\n";
            }
            if (g.n[2] % 6 != 0) f << "\n";
        }
}

}  // namespace mo
