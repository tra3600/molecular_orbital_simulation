#include "scf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <stdexcept>

namespace mo {

namespace {

// J[P] et K[P] à partir du tenseur ERI compact.
void buildJK(const ERITensor& eri, const Mat& P, Mat& J, Mat& K) {
    int n = P.r;
    J = Mat(n, n);
    K = Mat(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double j_ = 0, k_ = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l) {
                    j_ += P(k, l) * eri(i, j, k, l);
                    k_ += P(k, l) * eri(i, k, j, l);
                }
            J(i, j) = J(j, i) = j_;
            K(i, j) = K(j, i) = k_;
        }
}

Mat densityFromC(const Mat& C, int nocc) {
    int n = C.r;
    Mat P(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < nocc; ++k) s += C(i, k) * C(j, k);
            P(i, j) = s;
        }
    return P;
}

// Orthogonalisation canonique : X = U s^{-1/2} (colonnes de faible valeur propre supprimées)
Mat orthogonalizer(const Mat& S) {
    std::vector<double> w;
    Mat U;
    jacobiEigen(S, w, U);
    int n = S.r, m = 0;
    for (double v : w)
        if (v > 1e-7) ++m;
    Mat X(n, m);
    int col = 0;
    for (int k = 0; k < n; ++k) {
        if (w[k] <= 1e-7) continue;
        double f = 1.0 / std::sqrt(w[k]);
        for (int i = 0; i < n; ++i) X(i, col) = U(i, k) * f;
        ++col;
    }
    return X;
}

void diagonalize(const Mat& F, const Mat& X, std::vector<double>& eps, Mat& C) {
    Mat Fp = transpose(X) * F * X;
    Mat Cp;
    jacobiEigen(Fp, eps, Cp);
    C = X * Cp;
}

struct Diis {
    std::deque<std::vector<double>> focks, errs;
    size_t maxVec = 8;
    void add(const std::vector<double>& f, const std::vector<double>& e) {
        focks.push_back(f);
        errs.push_back(e);
        if (focks.size() > maxVec) {
            focks.pop_front();
            errs.pop_front();
        }
    }
    bool extrapolate(std::vector<double>& out) const {
        int m = static_cast<int>(focks.size());
        if (m < 2) return false;
        std::vector<std::vector<double>> B(m + 1, std::vector<double>(m + 1, 0.0));
        std::vector<double> rhs(m + 1, 0.0);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j) {
                double s = 0;
                for (size_t k = 0; k < errs[i].size(); ++k) s += errs[i][k] * errs[j][k];
                B[i][j] = s;
            }
        for (int i = 0; i < m; ++i) B[i][m] = B[m][i] = -1.0;
        rhs[m] = -1.0;
        std::vector<double> c;
        if (!solveLinear(B, rhs, c)) return false;
        out.assign(focks[0].size(), 0.0);
        for (int i = 0; i < m; ++i)
            for (size_t k = 0; k < out.size(); ++k) out[k] += c[i] * focks[i][k];
        return true;
    }
};

}  // namespace

SCFResult runSCF(const Molecule& mol, const Basis& basis, const SCFOptions& opt) {
    SCFResult res;
    const int n = basis.size();
    const int nel = mol.nElectrons();
    if (nel < 0) throw std::runtime_error("Nombre d'électrons négatif");
    int unpaired = mol.multiplicity - 1;
    if ((nel - unpaired) % 2 != 0 || unpaired > nel)
        throw std::runtime_error("Charge " + std::to_string(mol.charge) + " et multiplicité " +
                                 std::to_string(mol.multiplicity) + " incompatibles avec " + std::to_string(nel) +
                                 " électrons");
    res.nBeta = (nel - unpaired) / 2;
    res.nAlpha = res.nBeta + unpaired;
    if (res.nAlpha > n) throw std::runtime_error("Base trop petite pour ce nombre d'électrons");
    const bool uhf = opt.uhf || unpaired > 0;
    res.uhf = uhf;

    res.S = overlapMatrix(basis);
    Mat T = kineticMatrix(basis);
    Mat V = nuclearMatrix(basis, mol);
    res.H = T + V;
    ERITensor eri = computeERI(basis);
    res.eNuc = mol.nuclearRepulsion();
    Mat X = orthogonalizer(res.S);

    // Estimation initiale : Hamiltonien de Wolfsberg–Helmholz généralisé
    Mat F0(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            F0(i, j) = (i == j) ? res.H(i, i) : 0.875 * res.S(i, j) * (res.H(i, i) + res.H(j, j));
    diagonalize(F0, X, res.epsA, res.CA);
    res.CB = res.CA;
    res.epsB = res.epsA;
    if (uhf && opt.breakSymmetry && unpaired == 0 && res.nBeta > 0 && res.nBeta < static_cast<int>(res.epsA.size())) {
        int h = res.nBeta - 1, l = res.nBeta;
        double c = std::cos(M_PI / 8), s = std::sin(M_PI / 8);  // mélange partiel
        for (int i = 0; i < n; ++i) {
            double a = res.CB(i, h), b = res.CB(i, l);
            res.CB(i, h) = c * a + s * b;
            res.CB(i, l) = -s * a + c * b;
        }
    }
    res.PA = densityFromC(res.CA, res.nAlpha);
    res.PB = densityFromC(res.CB, res.nBeta);

    Diis diis;
    double eOld = 0;
    Mat FA, FB;
    for (int it = 1; it <= opt.maxIter; ++it) {
        Mat PT = res.PA + res.PB;
        Mat J, KA, KB, tmp;
        buildJK(eri, PT.r ? PT : PT, J, tmp);  // J[Ptotal]
        if (uhf) {
            Mat J2;
            buildJK(eri, res.PA, J2, KA);
            buildJK(eri, res.PB, J2, KB);
            FA = res.H + J - KA;
            FB = res.H + J - KB;
        } else {
            Mat J2, K;
            buildJK(eri, res.PA, J2, K);  // PA = P/2
            FA = res.H + J - K;
            FB = FA;
        }
        double eel = 0.5 * (trace_prod(PT, res.H) + trace_prod(res.PA, FA) + trace_prod(res.PB, FB));
        double etot = eel + res.eNuc;

        // erreur DIIS : X^T (F P S − S P F) X
        std::vector<double> fv, ev;
        auto errOf = [&](const Mat& F, const Mat& P) {
            Mat FPS = F * P * res.S;
            Mat e = transpose(X) * (FPS - transpose(FPS)) * X;
            return e;
        };
        Mat eA = errOf(FA, res.PA);
        fv.insert(fv.end(), FA.d.begin(), FA.d.end());
        ev.insert(ev.end(), eA.d.begin(), eA.d.end());
        if (uhf) {
            Mat eB = errOf(FB, res.PB);
            fv.insert(fv.end(), FB.d.begin(), FB.d.end());
            ev.insert(ev.end(), eB.d.begin(), eB.d.end());
        }
        double rms = 0;
        for (double v : ev) rms += v * v;
        rms = std::sqrt(rms / ev.size());
        res.energyHistory.push_back(etot);
        if (opt.verbose) std::printf("  iter %3d  E = %.10f  ΔE = %+.2e  err = %.2e\n", it, etot, etot - eOld, rms);

        res.iterations = it;
        res.eElec = eel;
        res.energy = etot;
        if (it > 1 && std::fabs(etot - eOld) < opt.energyTol && rms < opt.errorTol) {
            res.converged = true;
            break;
        }
        eOld = etot;

        diis.add(fv, ev);
        std::vector<double> ex;
        if (diis.extrapolate(ex)) {
            std::copy(ex.begin(), ex.begin() + FA.d.size(), FA.d.begin());
            if (uhf) std::copy(ex.begin() + FA.d.size(), ex.end(), FB.d.begin());
            else FB = FA;
        }
        diagonalize(FA, X, res.epsA, res.CA);
        if (uhf) diagonalize(FB, X, res.epsB, res.CB);
        else {
            res.CB = res.CA;
            res.epsB = res.epsA;
        }
        res.PA = densityFromC(res.CA, res.nAlpha);
        res.PB = densityFromC(res.CB, res.nBeta);
    }
    // dernière diagonalisation propre (sans extrapolation) pour des énergies orbitalaires cohérentes
    diagonalize(FA, X, res.epsA, res.CA);
    if (uhf) diagonalize(FB, X, res.epsB, res.CB);
    else {
        res.CB = res.CA;
        res.epsB = res.epsA;
    }
    res.PA = densityFromC(res.CA, res.nAlpha);
    res.PB = densityFromC(res.CB, res.nBeta);

    // ---- propriétés
    const int na = static_cast<int>(mol.atoms.size());
    Mat PT = res.PA + res.PB;
    Mat PS = PT * res.S;
    res.charges.assign(na, 0.0);
    res.spinPop.assign(na, 0.0);
    Mat PSs = (res.PA - res.PB) * res.S;
    for (int a = 0; a < na; ++a) res.charges[a] = mol.atoms[a].Z;
    for (int mu = 0; mu < n; ++mu) {
        res.charges[basis.fns[mu].atom] -= PS(mu, mu);
        res.spinPop[basis.fns[mu].atom] += PSs(mu, mu);
    }
    // indices de liaison de Mayer
    res.bondOrders = Mat(na, na);
    Mat PAS = res.PA * res.S, PBS = res.PB * res.S;
    for (int mu = 0; mu < n; ++mu)
        for (int nu = 0; nu < n; ++nu) {
            int A = basis.fns[mu].atom, B = basis.fns[nu].atom;
            if (A == B) continue;
            res.bondOrders(A, B) += 2.0 * (PAS(mu, nu) * PAS(nu, mu) + PBS(mu, nu) * PBS(nu, mu));
        }
    // moment dipolaire
    auto D = dipoleMatrices(basis);
    for (int d = 0; d < 3; ++d) {
        double el = -trace_prod(PT, D[d]);
        double nuc = 0;
        for (const Atom& a : mol.atoms) nuc += a.Z * (d == 0 ? a.x : d == 1 ? a.y : a.z);
        res.dipole[d] = el + nuc;
    }
    // <S²> pour UHF
    if (uhf) {
        double sz = 0.5 * (res.nAlpha - res.nBeta);
        double ov = 0;
        Mat M = transpose(res.CA) * res.S * res.CB;
        for (int i = 0; i < res.nAlpha; ++i)
            for (int j = 0; j < res.nBeta; ++j) ov += M(i, j) * M(i, j);
        res.s2 = sz * (sz + 1) + res.nBeta - ov;
    }
    return res;
}

}  // namespace mo
