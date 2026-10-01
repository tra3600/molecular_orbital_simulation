#include "integrals.h"

#include <algorithm>
#include <cmath>

namespace mo {

// ---------------------------------------------------------------- Boys
void boysFunction(int nmax, double x, double* F) {
    if (x < 1e-13) {
        for (int n = 0; n <= nmax; ++n) F[n] = 1.0 / (2 * n + 1);
        return;
    }
    if (x > 35.0) {
        // développement asymptotique : F_n(x) ≈ (2n-1)!!/2^{n+1} sqrt(π/x^{2n+1})
        F[0] = 0.5 * std::sqrt(M_PI / x);
        for (int n = 1; n <= nmax; ++n) F[n] = F[n - 1] * (2 * n - 1) / (2 * x);
        return;
    }
    // série pour F_nmax puis récurrence descendante (stable)
    double ex = std::exp(-x);
    double term = 1.0 / (2 * nmax + 1), sum = term;
    for (int k = 1; k < 300; ++k) {
        term *= 2 * x / (2 * nmax + 2 * k + 1);
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    F[nmax] = ex * sum;
    for (int n = nmax; n > 0; --n) F[n - 1] = (2 * x * F[n] + ex) / (2 * n - 1);
}

// ---------------------------------------------------- Coefficients de Hermite
static double Ecoef(int i, int j, int t, double Q, double a, double b) {
    if (t < 0 || t > i + j) return 0.0;
    double p = a + b, q = a * b / p;
    if (i == 0 && j == 0 && t == 0) return std::exp(-q * Q * Q);
    if (j == 0)
        return Ecoef(i - 1, j, t - 1, Q, a, b) / (2 * p) - (q * Q / a) * Ecoef(i - 1, j, t, Q, a, b) +
               (t + 1) * Ecoef(i - 1, j, t + 1, Q, a, b);
    return Ecoef(i, j - 1, t - 1, Q, a, b) / (2 * p) + (q * Q / b) * Ecoef(i, j - 1, t, Q, a, b) +
           (t + 1) * Ecoef(i, j - 1, t + 1, Q, a, b);
}

// Recouvrement 1D (sans facteur sqrt(π/p) → inclus ici)
static double S1(int i, int j, double Q, double a, double b) {
    if (j < 0 || i < 0) return 0.0;
    return Ecoef(i, j, 0, Q, a, b) * std::sqrt(M_PI / (a + b));
}

double overlapPair(const BasisFunction& A, const BasisFunction& B) {
    double s = 0;
    for (const Primitive& pa : A.prims)
        for (const Primitive& pb : B.prims) {
            double v = 1.0;
            for (int d = 0; d < 3; ++d) v *= S1(A.l[d], B.l[d], A.R[d] - B.R[d], pa.alpha, pb.alpha);
            s += pa.coef * pb.coef * v;
        }
    return s;
}

static double kineticPair(const BasisFunction& A, const BasisFunction& B) {
    double s = 0;
    for (const Primitive& pa : A.prims)
        for (const Primitive& pb : B.prims) {
            double a = pa.alpha, b = pb.alpha;
            double S[3], T[3];
            for (int d = 0; d < 3; ++d) {
                int i = A.l[d], j = B.l[d];
                double Q = A.R[d] - B.R[d];
                S[d] = S1(i, j, Q, a, b);
                T[d] = -0.5 * (j * (j - 1) * S1(i, j - 2, Q, a, b) - 2 * b * (2 * j + 1) * S[d] +
                               4 * b * b * S1(i, j + 2, Q, a, b));
            }
            s += pa.coef * pb.coef * (T[0] * S[1] * S[2] + S[0] * T[1] * S[2] + S[0] * S[1] * T[2]);
        }
    return s;
}

// ----------------------------------------------- Paires de primitives (ERI, V)
namespace {
constexpr int NR = 5;  // ordre de Hermite max pour s,p : 4
struct PrimPair {
    double p, cc;
    double P[3];
    int lmax[3];            // la+lb par dimension
    double E[3][3];         // E[dim][t], t<=2
};

struct ContractedPair {
    std::vector<PrimPair> pp;
    double Qschwarz = 0;
};

ContractedPair makePair(const BasisFunction& A, const BasisFunction& B) {
    ContractedPair cp;
    for (const Primitive& pa : A.prims)
        for (const Primitive& pb : B.prims) {
            PrimPair q;
            double a = pa.alpha, b = pb.alpha;
            q.p = a + b;
            q.cc = pa.coef * pb.coef;
            for (int d = 0; d < 3; ++d) {
                q.P[d] = (a * A.R[d] + b * B.R[d]) / q.p;
                q.lmax[d] = A.l[d] + B.l[d];
                for (int t = 0; t < 3; ++t) q.E[d][t] = Ecoef(A.l[d], B.l[d], t, A.R[d] - B.R[d], a, b);
            }
            cp.pp.push_back(q);
        }
    return cp;
}

// R^0_{tuv}(α, X) pour t+u+v <= L ; stocké dans R[t][u][v]
void hermiteR(int L, double alpha, double X, double Y, double Z, double R[NR][NR][NR]) {
    static thread_local double Rn[NR][NR][NR][NR];
    double F[NR + 1];
    boysFunction(L, alpha * (X * X + Y * Y + Z * Z), F);
    double pw = 1.0;
    for (int n = 0; n <= L; ++n) {
        Rn[n][0][0][0] = pw * F[n];
        pw *= -2 * alpha;
    }
    for (int T = 1; T <= L; ++T)
        for (int n = 0; n <= L - T; ++n)
            for (int t = 0; t <= T; ++t)
                for (int u = 0; u <= T - t; ++u) {
                    int v = T - t - u;
                    double val;
                    if (t > 0)
                        val = (t > 1 ? (t - 1) * Rn[n + 1][t - 2][u][v] : 0.0) + X * Rn[n + 1][t - 1][u][v];
                    else if (u > 0)
                        val = (u > 1 ? (u - 1) * Rn[n + 1][t][u - 2][v] : 0.0) + Y * Rn[n + 1][t][u - 1][v];
                    else
                        val = (v > 1 ? (v - 1) * Rn[n + 1][t][u][v - 2] : 0.0) + Z * Rn[n + 1][t][u][v - 1];
                    Rn[n][t][u][v] = val;
                }
    for (int t = 0; t <= L; ++t)
        for (int u = 0; t + u <= L; ++u)
            for (int v = 0; t + u + v <= L; ++v) R[t][u][v] = Rn[0][t][u][v];
}

double eriPrim(const PrimPair& a, const PrimPair& b) {
    double alpha = a.p * b.p / (a.p + b.p);
    double R[NR][NR][NR];
    int L = a.lmax[0] + a.lmax[1] + a.lmax[2] + b.lmax[0] + b.lmax[1] + b.lmax[2];
    hermiteR(L, alpha, a.P[0] - b.P[0], a.P[1] - b.P[1], a.P[2] - b.P[2], R);
    double s = 0;
    for (int t = 0; t <= a.lmax[0]; ++t)
        for (int u = 0; u <= a.lmax[1]; ++u)
            for (int v = 0; v <= a.lmax[2]; ++v) {
                double Eb = a.E[0][t] * a.E[1][u] * a.E[2][v];
                if (Eb == 0.0) continue;
                for (int tt = 0; tt <= b.lmax[0]; ++tt)
                    for (int uu = 0; uu <= b.lmax[1]; ++uu)
                        for (int vv = 0; vv <= b.lmax[2]; ++vv) {
                            double Ek = b.E[0][tt] * b.E[1][uu] * b.E[2][vv];
                            double sign = ((tt + uu + vv) & 1) ? -1.0 : 1.0;
                            s += Eb * Ek * sign * R[t + tt][u + uu][v + vv];
                        }
            }
    return a.cc * b.cc * s * 2.0 * std::pow(M_PI, 2.5) / (a.p * b.p * std::sqrt(a.p + b.p));
}

double eriContracted(const ContractedPair& a, const ContractedPair& b) {
    double s = 0;
    for (const PrimPair& x : a.pp)
        for (const PrimPair& y : b.pp) s += eriPrim(x, y);
    return s;
}
}  // namespace

// ------------------------------------------------------------- Matrices
Mat overlapMatrix(const Basis& b) {
    int n = b.size();
    Mat S(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) S(i, j) = S(j, i) = overlapPair(b.fns[i], b.fns[j]);
    return S;
}

Mat kineticMatrix(const Basis& b) {
    int n = b.size();
    Mat T(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) T(i, j) = T(j, i) = kineticPair(b.fns[i], b.fns[j]);
    return T;
}

Mat nuclearMatrix(const Basis& b, const Molecule& mol) {
    int n = b.size();
    Mat V(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            ContractedPair cp = makePair(b.fns[i], b.fns[j]);
            double s = 0;
            for (const PrimPair& q : cp.pp)
                for (const Atom& at : mol.atoms) {
                    double R[NR][NR][NR];
                    int L = q.lmax[0] + q.lmax[1] + q.lmax[2];
                    hermiteR(L, q.p, q.P[0] - at.x, q.P[1] - at.y, q.P[2] - at.z, R);
                    double v = 0;
                    for (int t = 0; t <= q.lmax[0]; ++t)
                        for (int u = 0; u <= q.lmax[1]; ++u)
                            for (int w = 0; w <= q.lmax[2]; ++w) v += q.E[0][t] * q.E[1][u] * q.E[2][w] * R[t][u][w];
                    s -= at.Z * q.cc * 2 * M_PI / q.p * v;
                }
            V(i, j) = V(j, i) = s;
        }
    return V;
}

std::array<Mat, 3> dipoleMatrices(const Basis& b) {
    int n = b.size();
    std::array<Mat, 3> M{Mat(n, n), Mat(n, n), Mat(n, n)};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            const BasisFunction &A = b.fns[i], &B = b.fns[j];
            double acc[3] = {0, 0, 0};
            for (const Primitive& pa : A.prims)
                for (const Primitive& pb : B.prims) {
                    double a = pa.alpha, bb = pb.alpha, p = a + bb;
                    double S[3], D[3];
                    for (int d = 0; d < 3; ++d) {
                        double Q = A.R[d] - B.R[d], P = (a * A.R[d] + bb * B.R[d]) / p;
                        double sq = std::sqrt(M_PI / p);
                        S[d] = Ecoef(A.l[d], B.l[d], 0, Q, a, bb) * sq;
                        D[d] = (Ecoef(A.l[d], B.l[d], 1, Q, a, bb) + P * Ecoef(A.l[d], B.l[d], 0, Q, a, bb)) * sq;
                    }
                    double c = pa.coef * pb.coef;
                    acc[0] += c * D[0] * S[1] * S[2];
                    acc[1] += c * S[0] * D[1] * S[2];
                    acc[2] += c * S[0] * S[1] * D[2];
                }
            for (int d = 0; d < 3; ++d) M[d](i, j) = M[d](j, i) = acc[d];
        }
    return M;
}

ERITensor computeERI(const Basis& b) {
    const int n = b.size();
    ERITensor T(n);
    const size_t np = static_cast<size_t>(n) * (n + 1) / 2;
    std::vector<ContractedPair> pairs(np);
    std::vector<std::pair<int, int>> ij(np);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            size_t k = ERITensor::pair(i, j);
            pairs[k] = makePair(b.fns[i], b.fns[j]);
            ij[k] = {i, j};
        }
    // Schwarz : Q_ij = sqrt((ij|ij))
    for (size_t k = 0; k < np; ++k) pairs[k].Qschwarz = std::sqrt(std::fabs(eriContracted(pairs[k], pairs[k])));
#pragma omp parallel for schedule(dynamic)
    for (long long k = 0; k < static_cast<long long>(np); ++k)
        for (size_t l = 0; l <= static_cast<size_t>(k); ++l) {
            if (pairs[k].Qschwarz * pairs[l].Qschwarz < 1e-11) continue;
            T.at(k, l) = eriContracted(pairs[k], pairs[l]);
        }
    return T;
}

}  // namespace mo
