#include "linalg.h"

#include <algorithm>
#include <numeric>

namespace mo {

Mat operator*(const Mat& a, const Mat& b) {
    Mat m(a.r, b.c);
    for (int i = 0; i < a.r; ++i)
        for (int k = 0; k < a.c; ++k) {
            double aik = a(i, k);
            if (aik == 0.0) continue;
            for (int j = 0; j < b.c; ++j) m(i, j) += aik * b(k, j);
        }
    return m;
}
Mat operator+(const Mat& a, const Mat& b) {
    Mat m = a;
    for (size_t i = 0; i < m.d.size(); ++i) m.d[i] += b.d[i];
    return m;
}
Mat operator-(const Mat& a, const Mat& b) {
    Mat m = a;
    for (size_t i = 0; i < m.d.size(); ++i) m.d[i] -= b.d[i];
    return m;
}
Mat operator*(double s, const Mat& a) {
    Mat m = a;
    for (double& v : m.d) v *= s;
    return m;
}
Mat transpose(const Mat& a) {
    Mat m(a.c, a.r);
    for (int i = 0; i < a.r; ++i)
        for (int j = 0; j < a.c; ++j) m(j, i) = a(i, j);
    return m;
}
double trace_prod(const Mat& a, const Mat& b) {
    double s = 0;
    for (size_t i = 0; i < a.d.size(); ++i) s += a.d[i] * b.d[i];
    return s;
}

void jacobiEigen(const Mat& A0, std::vector<double>& w, Mat& V) {
    const int n = A0.r;
    Mat A = A0;
    V = Mat::identity(n);
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) off += A(i, j) * A(i, j);
        if (off < 1e-24) break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                double apq = A(p, q);
                if (std::fabs(apq) < 1e-300) continue;
                double theta = (A(q, q) - A(p, p)) / (2.0 * apq);
                double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; ++k) {
                    double akp = A(k, p), akq = A(k, q);
                    A(k, p) = c * akp - s * akq;
                    A(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    double apk = A(p, k), aqk = A(q, k);
                    A(p, k) = c * apk - s * aqk;
                    A(q, k) = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    double vkp = V(k, p), vkq = V(k, q);
                    V(k, p) = c * vkp - s * vkq;
                    V(k, q) = s * vkp + c * vkq;
                }
            }
    }
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return A(a, a) < A(b, b); });
    w.assign(n, 0.0);
    Mat Vs(n, n);
    for (int j = 0; j < n; ++j) {
        w[j] = A(idx[j], idx[j]);
        for (int i = 0; i < n; ++i) Vs(i, j) = V(i, idx[j]);
    }
    V = Vs;
}

bool solveLinear(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col])) piv = r;
        if (std::fabs(A[piv][col]) < 1e-14) return false;
        std::swap(A[piv], A[col]);
        std::swap(b[piv], b[col]);
        for (int r = col + 1; r < n; ++r) {
            double f = A[r][col] / A[col][col];
            for (int c = col; c < n; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return true;
}

}  // namespace mo
