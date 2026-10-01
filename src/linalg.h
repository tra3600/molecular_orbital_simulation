// Algèbre linéaire minimale (matrices denses, diagonalisation de Jacobi).
#pragma once
#include <cmath>
#include <vector>

namespace mo {

struct Mat {
    int r = 0, c = 0;
    std::vector<double> d;
    Mat() = default;
    Mat(int rows, int cols, double v = 0.0) : r(rows), c(cols), d(static_cast<size_t>(rows) * cols, v) {}
    double& operator()(int i, int j) { return d[static_cast<size_t>(i) * c + j]; }
    double operator()(int i, int j) const { return d[static_cast<size_t>(i) * c + j]; }
    static Mat identity(int n) {
        Mat m(n, n);
        for (int i = 0; i < n; ++i) m(i, i) = 1.0;
        return m;
    }
};

Mat operator*(const Mat& a, const Mat& b);
Mat operator+(const Mat& a, const Mat& b);
Mat operator-(const Mat& a, const Mat& b);
Mat operator*(double s, const Mat& a);
Mat transpose(const Mat& a);
double trace_prod(const Mat& a, const Mat& b);  // Tr(A^T B) = somme_ij A_ij B_ij

// Diagonalisation d'une matrice symétrique : valeurs propres croissantes,
// vecteurs propres en colonnes de V.
void jacobiEigen(const Mat& A, std::vector<double>& w, Mat& V);

// Résout A x = b (élimination de Gauss avec pivot). Retourne false si singulier.
bool solveLinear(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& x);

}  // namespace mo
