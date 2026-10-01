#include "report.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "grid.h"
#include "png.h"
#include "report_template.h"

namespace mo {

namespace fs = std::filesystem;

std::string OrbitalRef::name() const {
    std::string s = std::to_string(index + 1);
    if (spin == 1) s += "α";
    if (spin == 2) s += "β";
    return s;
}

static std::string fmt(const char* f, double v) {
    char b[64];
    std::snprintf(b, sizeof b, f, v);
    return b;
}

std::vector<OrbitalRef> listOrbitals(const Basis& basis, const SCFResult& r) {
    std::vector<OrbitalRef> out;
    const int n = basis.size();
    int nm = static_cast<int>(r.epsA.size());
    auto make = [&](int spin, int i, const Mat& C, double eps, double occ, int nocc) {
        OrbitalRef o;
        o.spin = spin;
        o.index = i;
        o.eps = eps;
        o.occ = occ;
        if (i == nocc - 1) o.tag = "HOMO";
        if (i == nocc) o.tag = "LUMO";
        std::vector<std::pair<std::string, double>> comp;
        for (int m = 0; m < n; ++m)
            if (std::fabs(C(m, i)) > 0.12) comp.push_back({basis.fns[m].label, C(m, i)});
        std::sort(comp.begin(), comp.end(), [](auto& a, auto& b) { return std::fabs(a.second) > std::fabs(b.second); });
        if (comp.size() > 8) comp.resize(8);
        o.comp = comp;
        return o;
    };
    if (!r.uhf) {
        for (int i = 0; i < nm; ++i) out.push_back(make(0, i, r.CA, r.epsA[i], i < r.nAlpha ? 2 : 0, r.nAlpha));
    } else {
        for (int i = 0; i < nm; ++i) out.push_back(make(1, i, r.CA, r.epsA[i], i < r.nAlpha ? 1 : 0, r.nAlpha));
        for (int i = 0; i < nm; ++i) out.push_back(make(2, i, r.CB, r.epsB[i], i < r.nBeta ? 1 : 0, r.nBeta));
    }
    return out;
}

// ------------------------------------------------------------ console
void printSummary(std::ostream& os, const Molecule& mol, const Basis& basis, const SCFResult& r) {
    auto orbs = listOrbitals(basis, r);
    os << "\n=== " << mol.name << "  (" << mol.formula() << ") ===\n";
    if (!mol.description.empty()) os << mol.description << "\n";
    os << "Méthode : " << (r.uhf ? "UHF" : "RHF") << "/STO-3G   charge " << mol.charge << ", multiplicité "
       << mol.multiplicity << "   " << basis.size() << " fonctions de base\n";
    os << std::fixed << std::setprecision(8);
    os << "SCF " << (r.converged ? "convergé" : "NON CONVERGÉ") << " en " << r.iterations << " itérations\n";
    os << "  Énergie nucléaire    : " << std::setw(16) << r.eNuc << " Eh\n";
    os << "  Énergie électronique : " << std::setw(16) << r.eElec << " Eh\n";
    os << "  Énergie totale       : " << std::setw(16) << r.energy << " Eh   (" << std::setprecision(4)
       << r.energy * HARTREE_TO_EV << " eV)\n";
    if (r.uhf) os << std::setprecision(4) << "  <S²> = " << r.s2 << "  (théorique " << 0.25 * (r.nAlpha - r.nBeta) * (r.nAlpha - r.nBeta + 2) << ")\n";

    os << "\nOrbitales moléculaires (LCAO) :\n";
    os << "   MO   énergie (Eh)   énergie (eV)  occ  composition dominante\n";
    for (const OrbitalRef& o : orbs) {
        os << std::setw(5) << o.name() << std::setw(15) << std::setprecision(5) << o.eps << std::setw(15)
           << std::setprecision(3) << o.eps * HARTREE_TO_EV << std::setw(5) << std::setprecision(0) << o.occ << "  ";
        int shown = 0;
        for (auto& c : o.comp) {
            if (shown++ >= 4) break;
            os << (c.second >= 0 ? "+" : "-") << std::setprecision(2) << std::fabs(c.second) << "·" << c.first << " ";
        }
        if (!o.tag.empty()) os << " <" << o.tag << ">";
        os << "\n";
    }
    // diagramme ASCII (valence)
    os << "\nDiagramme d'énergie (orbitales de valence) :\n";
    {
        std::vector<OrbitalRef> v;
        for (auto& o : orbs)
            if (o.spin != 2 && o.eps > -2.0 - 1e-9) v.push_back(o);
        std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a.eps > b.eps; });
        for (auto& o : v) {
            std::string occ = o.occ >= 2 ? "↑↓" : (o.occ == 1 ? "↑ " : "  ");
            os << "  " << std::setw(8) << std::setprecision(3) << o.eps << " Eh  ──[" << occ << "]──  MO " << o.name();
            if (!o.tag.empty()) os << "  (" << o.tag << ")";
            os << "\n";
        }
    }
    int h = r.nAlpha - 1;
    if (h >= 0 && h + 1 < static_cast<int>(r.epsA.size()))
        os << std::setprecision(3) << "\nGap HOMO–LUMO : " << (r.epsA[h + 1] - r.epsA[h]) * HARTREE_TO_EV << " eV\n";

    os << "\nCharges de Mulliken" << (r.uhf ? " et populations de spin" : "") << " :\n";
    for (size_t a = 0; a < mol.atoms.size(); ++a) {
        os << "  " << std::setw(2) << elementSymbol(mol.atoms[a].Z) << (a + 1) << std::setw(10) << std::setprecision(4)
           << std::showpos << r.charges[a] << std::noshowpos;
        if (r.uhf) os << std::setw(10) << std::showpos << r.spinPop[a] << std::noshowpos;
        os << "\n";
    }
    bool any = false;
    for (size_t a = 0; a < mol.atoms.size(); ++a)
        for (size_t b = a + 1; b < mol.atoms.size(); ++b)
            if (r.bondOrders(int(a), int(b)) > 0.15) {
                if (!any) os << "\nIndices de liaison de Mayer :\n";
                any = true;
                os << "  " << elementSymbol(mol.atoms[a].Z) << a + 1 << "–" << elementSymbol(mol.atoms[b].Z) << b + 1
                   << "  (" << std::setprecision(3) << mol.distance(int(a), int(b)) / ANGSTROM_TO_BOHR << " Å)  "
                   << std::setprecision(3) << r.bondOrders(int(a), int(b)) << "\n";
            }
    double dn = std::sqrt(r.dipole[0] * r.dipole[0] + r.dipole[1] * r.dipole[1] + r.dipole[2] * r.dipole[2]);
    if (mol.charge == 0)
        os << "\nMoment dipolaire : " << std::setprecision(3) << dn * AU_TO_DEBYE << " D\n";
    os << std::defaultfloat;
}

// ------------------------------------------------------------ SVG
static const char* SVG_FONT = "font-family=\"DejaVu Sans,Arial,sans-serif\"";

std::string energyDiagramSVG(const Molecule& mol, const std::vector<OrbitalRef>& orbs, bool uhf) {
    auto f = [](double e) { return e >= -1.5 ? e : -1.5 - 0.17 * std::log(1 + (-1.5 - e)); };
    double lo = 1e9, hi = -1e9;
    for (auto& o : orbs) {
        lo = std::min(lo, f(o.eps));
        hi = std::max(hi, f(o.eps));
    }
    if (hi - lo < 1e-6) hi = lo + 1;
    const int W = uhf ? 560 : 380, top = 50, bot = 40;
    const double GAP = 17;
    struct Group { std::vector<const OrbitalRef*> m; double y; };
    std::vector<std::vector<Group>> cols;
    double ymax = 0;
    for (int spin : uhf ? std::vector<int>{1, 2} : std::vector<int>{0}) {
        std::vector<const OrbitalRef*> v;
        for (auto& o : orbs) if (o.spin == spin) v.push_back(&o);
        std::sort(v.begin(), v.end(), [](auto* a, auto* b) { return a->eps > b->eps; });
        std::vector<Group> gs;
        for (auto* o : v) {
            if (!gs.empty() && std::fabs(gs.back().m[0]->eps - o->eps) < 2e-3) gs.back().m.push_back(o);
            else gs.push_back({{o}, 0});
        }
        double prev = -1e9;
        for (auto& g : gs) {
            g.y = std::max(top + (hi - f(g.m[0]->eps)) / (hi - lo) * 480, prev + GAP);
            prev = g.y;
        }
        ymax = std::max(ymax, prev);
        cols.push_back(gs);
    }
    const int H = int(std::max(520.0, ymax + bot));
    std::ostringstream s;
    s << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W
      << " " << H << "\" " << SVG_FONT << " font-size=\"12\">\n<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    s << "<text x=\"" << W / 2 << "\" y=\"26\" text-anchor=\"middle\" font-size=\"15\" font-weight=\"bold\">Niveaux d'énergie de "
      << mol.formula() << "</text>\n";
    for (size_t c = 0; c < cols.size(); ++c) {
        double cx = uhf ? (c == 0 ? W * 0.27 : W * 0.73) : W * 0.45;
        if (uhf) s << "<text x=\"" << cx << "\" y=\"" << top - 14 << "\" text-anchor=\"middle\" fill=\"#555\">" << (c == 0 ? "spin α ↑" : "spin β ↓") << "</text>\n";
        for (auto& g : cols[c]) {
            double tot = double(g.m.size()), sw = std::min(46.0, (W - 120) / tot), xs = cx - tot * sw / 2 + 3;
            std::string names, tag;
            for (size_t k = 0; k < g.m.size(); ++k) {
                const OrbitalRef* o = g.m[k];
                double x0 = xs + sw * k;
                const char* col = o->occ > 0 ? "#c0392b" : "#2c6fbb";
                s << "<line x1=\"" << x0 << "\" x2=\"" << x0 + sw - 8 << "\" y1=\"" << g.y << "\" y2=\"" << g.y << "\" stroke=\"" << col << "\" stroke-width=\"3\" stroke-linecap=\"round\"/>\n";
                std::string ar = o->occ >= 2 ? "↑↓" : (o->occ == 1 ? (o->spin == 2 ? "↓" : "↑") : "");
                if (!ar.empty()) s << "<text x=\"" << x0 + (sw - 8) / 2 << "\" y=\"" << g.y - 4 << "\" text-anchor=\"middle\" fill=\"#222\">" << ar << "</text>\n";
                if (tot <= 3) names += (k ? " " : "") + o->name();
                if (!o->tag.empty()) tag = " · " + o->tag;
            }
            if (tot > 3) names = g.m.front()->name() + "–" + g.m.back()->name();
            s << "<text x=\"" << xs - 6 << "\" y=\"" << g.y + 4 << "\" text-anchor=\"end\" font-size=\"11\" fill=\"#555\">" << fmt("%.3f", g.m[0]->eps) << "</text>\n";
            s << "<text x=\"" << xs + tot * sw << "\" y=\"" << g.y + 4 << "\" font-size=\"11\" fill=\"#777\">" << names << tag << "</text>\n";
        }
    }
    s << "<text x=\"10\" y=\"" << H - 12 << "\" font-size=\"10\" fill=\"#888\">énergies en hartree (axe comprimé sous −1,5 Eh)</text>\n</svg>\n";
    return s.str();
}

std::string svgPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel,
                    const std::vector<Series>& series) {
    double x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18;
    for (auto& s : series)
        for (size_t i = 0; i < s.x.size(); ++i) {
            x0 = std::min(x0, s.x[i]); x1 = std::max(x1, s.x[i]);
            y0 = std::min(y0, s.y[i]); y1 = std::max(y1, s.y[i]);
        }
    if (y1 - y0 < 1e-9) y1 = y0 + 1;
    if (x1 - x0 < 1e-9) x1 = x0 + 1;
    double pad = (y1 - y0) * 0.06;
    y0 -= pad; y1 += pad;
    const int W = 760, H = 480, L = 80, R = 20, T = 50, B = 60;
    auto X = [&](double x) { return L + (x - x0) / (x1 - x0) * (W - L - R); };
    auto Y = [&](double y) { return T + (y1 - y) / (y1 - y0) * (H - T - B); };
    std::ostringstream s;
    s << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W << " " << H << "\" " << SVG_FONT << " font-size=\"12\">\n<rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>\n";
    s << "<text x=\"" << W / 2 << "\" y=\"28\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">" << title << "</text>\n";
    for (int i = 0; i <= 5; ++i) {
        double xv = x0 + (x1 - x0) * i / 5, yv = y0 + (y1 - y0) * i / 5;
        s << "<line x1=\"" << X(xv) << "\" x2=\"" << X(xv) << "\" y1=\"" << T << "\" y2=\"" << H - B << "\" stroke=\"#e4e4e4\"/>\n";
        s << "<line x1=\"" << L << "\" x2=\"" << W - R << "\" y1=\"" << Y(yv) << "\" y2=\"" << Y(yv) << "\" stroke=\"#e4e4e4\"/>\n";
        s << "<text x=\"" << X(xv) << "\" y=\"" << H - B + 18 << "\" text-anchor=\"middle\">" << fmt("%.2f", xv) << "</text>\n";
        s << "<text x=\"" << L - 8 << "\" y=\"" << Y(yv) + 4 << "\" text-anchor=\"end\">" << fmt("%.3f", yv) << "</text>\n";
    }
    s << "<rect x=\"" << L << "\" y=\"" << T << "\" width=\"" << W - L - R << "\" height=\"" << H - T - B << "\" fill=\"none\" stroke=\"#888\"/>\n";
    s << "<text x=\"" << (L + W - R) / 2 << "\" y=\"" << H - 15 << "\" text-anchor=\"middle\">" << xlabel << "</text>\n";
    s << "<text transform=\"translate(18," << (T + H - B) / 2 << ") rotate(-90)\" text-anchor=\"middle\">" << ylabel << "</text>\n";
    int ly = T + 16;
    for (auto& se : series) {
        std::ostringstream pts;
        for (size_t i = 0; i < se.x.size(); ++i) pts << X(se.x[i]) << "," << Y(se.y[i]) << " ";
        s << "<polyline fill=\"none\" stroke=\"" << se.color << "\" stroke-width=\"2\" points=\"" << pts.str() << "\"/>\n";
        if (se.markers)
            for (size_t i = 0; i < se.x.size(); ++i) s << "<circle cx=\"" << X(se.x[i]) << "\" cy=\"" << Y(se.y[i]) << "\" r=\"2.6\" fill=\"" << se.color << "\"/>\n";
        s << "<line x1=\"" << W - R - 170 << "\" x2=\"" << W - R - 145 << "\" y1=\"" << ly << "\" y2=\"" << ly << "\" stroke=\"" << se.color << "\" stroke-width=\"3\"/><text x=\"" << W - R - 138 << "\" y=\"" << ly + 4 << "\">" << se.name << "</text>\n";
        ly += 18;
    }
    s << "</svg>\n";
    return s.str();
}

// ------------------------------------------------------------ coupes 2D
namespace {
struct Plane {
    std::array<double, 3> c, u, v;
    double umin, umax, vmin, vmax;
    int w, h;
};

Plane makePlane(const Molecule& mol, int which, int px) {
    Plane p;
    int na = static_cast<int>(mol.atoms.size());
    std::array<double, 3> c{0, 0, 0};
    for (auto& a : mol.atoms) { c[0] += a.x / na; c[1] += a.y / na; c[2] += a.z / na; }
    Mat cov(3, 3);
    for (auto& a : mol.atoms) {
        double d[3] = {a.x - c[0], a.y - c[1], a.z - c[2]};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) cov(i, j) += d[i] * d[j];
    }
    std::vector<double> w;
    Mat V;
    jacobiEigen(cov, w, V);  // croissant
    auto col = [&](int k) { return std::array<double, 3>{V(0, k), V(1, k), V(2, k)}; };
    p.c = c;
    p.u = col(2);
    p.v = which == 0 ? col(1) : col(0);
    // moléculaire plane/linéaire : le second axe est bien défini par l'orthonormalité de Jacobi
    double lo_u = 1e9, hi_u = -1e9, lo_v = 1e9, hi_v = -1e9;
    for (auto& a : mol.atoms) {
        double d[3] = {a.x - c[0], a.y - c[1], a.z - c[2]};
        double uu = d[0] * p.u[0] + d[1] * p.u[1] + d[2] * p.u[2];
        double vv = d[0] * p.v[0] + d[1] * p.v[1] + d[2] * p.v[2];
        lo_u = std::min(lo_u, uu); hi_u = std::max(hi_u, uu);
        lo_v = std::min(lo_v, vv); hi_v = std::max(hi_v, vv);
    }
    const double m = 3.6;
    p.umin = lo_u - m; p.umax = hi_u + m; p.vmin = lo_v - m; p.vmax = hi_v + m;
    // le plan n'est pas toujours carré : on garde une échelle uniforme
    double s = px / std::max(p.umax - p.umin, p.vmax - p.vmin);
    p.w = std::max(40, int(std::round((p.umax - p.umin) * s)));
    p.h = std::max(40, int(std::round((p.vmax - p.vmin) * s)));
    return p;
}

std::array<double, 3> planePoint(const Plane& p, int ix, int iy) {
    double u = p.umin + (ix + 0.5) / p.w * (p.umax - p.umin);
    double v = p.vmax - (iy + 0.5) / p.h * (p.vmax - p.vmin);
    return {p.c[0] + u * p.u[0] + v * p.v[0], p.c[1] + u * p.u[1] + v * p.v[1], p.c[2] + u * p.u[2] + v * p.v[2]};
}

void shade(std::vector<uint8_t>& img, int w, int h, const std::vector<float>& val, bool signedField) {
    double vmax = 1e-12;
    for (float v : val) vmax = std::max(vmax, double(std::fabs(v)));
    static const double lev[] = {0.03, 0.1, 0.25, 0.5, 0.8};
    std::vector<int> L(val.size());
    std::vector<double> t(val.size());
    for (size_t i = 0; i < val.size(); ++i) {
        double a = std::fabs(val[i]) / vmax;
        int k = 0;
        for (double l : lev) if (a >= l) ++k;
        L[i] = (val[i] < 0 && signedField) ? -k : k;
        t[i] = a;
    }
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            size_t i = size_t(y) * w + x;
            double r, g, b;
            if (signedField) {
                double a = std::pow(t[i], 0.55);
                if (val[i] >= 0) { r = 255 - a * (255 - 214); g = 255 - a * (255 - 52); b = 255 - a * (255 - 40); }
                else { r = 255 - a * (255 - 38); g = 255 - a * (255 - 96); b = 255 - a * (255 - 206); }
            } else {
                double a = std::pow(t[i], 0.35);
                static const double st[5][3] = {{255, 255, 255}, {255, 232, 150}, {238, 120, 55}, {130, 30, 100}, {25, 10, 60}};
                double q = a * 4;
                int k = std::min(3, int(q));
                double f = q - k;
                r = st[k][0] + f * (st[k + 1][0] - st[k][0]);
                g = st[k][1] + f * (st[k + 1][1] - st[k][1]);
                b = st[k][2] + f * (st[k + 1][2] - st[k][2]);
            }
            bool edge = (x + 1 < w && L[i] != L[i + 1]) || (y + 1 < h && L[i] != L[i + w]);
            if (edge) { r *= 0.45; g *= 0.45; b *= 0.45; }
            img[i * 3] = uint8_t(r); img[i * 3 + 1] = uint8_t(g); img[i * 3 + 2] = uint8_t(b);
        }
}

void drawAtoms(std::vector<uint8_t>& img, const Plane& p, const Molecule& mol) {
    for (auto& a : mol.atoms) {
        double d[3] = {a.x - p.c[0], a.y - p.c[1], a.z - p.c[2]};
        double u = d[0] * p.u[0] + d[1] * p.u[1] + d[2] * p.u[2];
        double v = d[0] * p.v[0] + d[1] * p.v[1] + d[2] * p.v[2];
        double nrm[3] = {p.u[1] * p.v[2] - p.u[2] * p.v[1], p.u[2] * p.v[0] - p.u[0] * p.v[2], p.u[0] * p.v[1] - p.u[1] * p.v[0]};
        double off = std::fabs(d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2]);
        double cx = (u - p.umin) / (p.umax - p.umin) * p.w, cy = (p.vmax - v) / (p.vmax - p.vmin) * p.h;
        bool inPlane = off < 0.4;
        double rad = inPlane ? 5.5 : 3.5;
        for (int y = int(cy - rad - 1); y <= int(cy + rad + 1); ++y)
            for (int x = int(cx - rad - 1); x <= int(cx + rad + 1); ++x) {
                if (x < 0 || y < 0 || x >= p.w || y >= p.h) continue;
                double r = std::hypot(x + 0.5 - cx, y + 0.5 - cy);
                size_t i = (size_t(y) * p.w + x) * 3;
                if (r <= rad && inPlane) { img[i] = img[i + 1] = img[i + 2] = a.Z == 1 ? 255 : 30; }
                if (r <= rad && r > rad - 1.6) { img[i] = img[i + 1] = img[i + 2] = 20; }
            }
    }
}
}  // namespace

// ------------------------------------------------------------ JSON
static std::string jstr(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (c == '\n') o += "\\n";
        else o += c;
    }
    return o + "\"";
}
static std::string jnum(double v) {
    if (!std::isfinite(v)) return "0";
    char b[64];
    std::snprintf(b, sizeof b, "%.8g", v);
    return b;
}

std::string writeOutputs(const Molecule& mol, const Basis& basis, const SCFResult& r, const OutputOptions& opt) {
    fs::create_directories(opt.outDir);
    auto orbs = listOrbitals(basis, r);
    std::string svg = energyDiagramSVG(mol, orbs, r.uhf);
    { std::ofstream(opt.outDir + "/energy_levels.svg") << svg; }

    Grid grid = makeGrid(mol, opt.margin, opt.spacing, 1500000);
    Plane planes[2] = {makePlane(mol, 0, opt.slicePx), makePlane(mol, 1, opt.slicePx)};
    Mat PT = r.PA + r.PB;

    // champs sur les plans (orbitales + densité)
    std::vector<std::string> sliceB64[2];
    std::string densB64[2];
    std::vector<int> idxs;
    std::vector<const Mat*> cs;
    for (auto& o : orbs) { idxs.push_back(o.index); cs.push_back(o.spin == 2 ? &r.CB : &r.CA); }
    for (int pl = 0; pl < 2; ++pl) {
        const Plane& P = planes[pl];
        std::vector<std::array<double, 3>> pts;
        for (int y = 0; y < P.h; ++y) for (int x = 0; x < P.w; ++x) pts.push_back(planePoint(P, x, y));
        FieldSet fsx = evaluateFields(basis, cs, idxs, &PT, pts);
        std::string stem = opt.outDir + "/slice" + std::string(pl == 0 ? "A" : "B");
        for (size_t k = 0; k < orbs.size(); ++k) {
            std::vector<uint8_t> img(size_t(P.w) * P.h * 3);
            shade(img, P.w, P.h, fsx.orbitals[k], true);
            drawAtoms(img, P, mol);
            sliceB64[pl].push_back(base64(encodePNG(P.w, P.h, img)));
            writePNG(stem + "_mo" + orbs[k].name() + ".png", P.w, P.h, img);
        }
        std::vector<uint8_t> img(size_t(P.w) * P.h * 3);
        shade(img, P.w, P.h, fsx.density, false);
        drawAtoms(img, P, mol);
        densB64[pl] = base64(encodePNG(P.w, P.h, img));
        writePNG(stem + "_density.png", P.w, P.h, img);
    }

    // grille 3D : cubes + échantillonnage du nuage de points
    std::mt19937 rng(12345);
    std::vector<std::vector<int>> clouds(orbs.size());
    const size_t batch = 6;
    for (size_t b0 = 0; b0 < orbs.size() + 1; b0 += batch) {
        // le dernier lot contient aussi la densité
        std::vector<int> ix;
        std::vector<const Mat*> cc;
        for (size_t k = b0; k < std::min(orbs.size(), b0 + batch); ++k) { ix.push_back(idxs[k]); cc.push_back(cs[k]); }
        bool lastBatch = b0 + batch >= orbs.size();
        if (ix.empty() && !lastBatch) break;
        std::vector<std::array<double, 3>> pts;
        pts.reserve(grid.size());
        for (int i = 0; i < grid.n[0]; ++i) for (int j = 0; j < grid.n[1]; ++j) for (int k = 0; k < grid.n[2]; ++k) pts.push_back(grid.point(i, j, k));
        FieldSet f = evaluateFields(basis, cc, ix, lastBatch ? &PT : nullptr, pts);
        for (size_t k = 0; k < ix.size(); ++k) {
            size_t o = b0 + k;
            if (opt.cubes)
                writeCube(opt.outDir + "/mo_" + orbs[o].name() + ".cube", mol, grid, f.orbitals[k],
                          "Orbitale moleculaire " + orbs[o].name() + " E=" + fmt("%.6f", orbs[o].eps) + " Eh");
            std::vector<double> w(f.orbitals[k].size());
            for (size_t i = 0; i < w.size(); ++i) w[i] = double(f.orbitals[k][i]) * f.orbitals[k][i];
            std::discrete_distribution<size_t> dist(w.begin(), w.end());
            std::uniform_real_distribution<double> jit(-0.5, 0.5);
            for (int s = 0; s < opt.samples; ++s) {
                size_t id = dist(rng);
                int kk = int(id % grid.n[2]), jj = int((id / grid.n[2]) % grid.n[1]), ii = int(id / (size_t(grid.n[2]) * grid.n[1]));
                auto p = grid.point(ii, jj, kk);
                for (int d = 0; d < 3; ++d) clouds[o].push_back(int(std::lround((p[d] + jit(rng) * grid.step) * 100)));
                clouds[o].push_back(f.orbitals[k][id] >= 0 ? 1 : -1);
            }
        }
        if (lastBatch && opt.cubes) writeCube(opt.outDir + "/density.cube", mol, grid, f.density, "Densite electronique totale");
        if (lastBatch) break;
    }

    // résumé texte
    {
        std::ofstream t(opt.outDir + "/resume.txt");
        printSummary(t, mol, basis, r);
    }
    if (!opt.html) return "";

    // JSON pour la page
    std::ostringstream j;
    double dn = std::sqrt(r.dipole[0] * r.dipole[0] + r.dipole[1] * r.dipole[1] + r.dipole[2] * r.dipole[2]);
    int h = r.nAlpha - 1;
    double gap = (h >= 0 && h + 1 < int(r.epsA.size())) ? (r.epsA[h + 1] - r.epsA[h]) * HARTREE_TO_EV : 0;
    j << "{\"title\":" << jstr(mol.name) << ",\"desc\":" << jstr(mol.description) << ",\"formula\":" << jstr(mol.formula())
      << ",\"method\":" << jstr(r.uhf ? "UHF/STO-3G" : "RHF/STO-3G") << ",\"charge\":" << mol.charge << ",\"mult\":" << mol.multiplicity
      << ",\"nbf\":" << basis.size() << ",\"uhf\":" << (r.uhf ? "true" : "false") << ",\"converged\":" << (r.converged ? "true" : "false")
      << ",\"iterations\":" << r.iterations << ",\"energy\":" << jnum(r.energy) << ",\"eNuc\":" << jnum(r.eNuc) << ",\"eElec\":" << jnum(r.eElec)
      << ",\"gap\":" << jnum(gap) << ",\"dipole\":" << jnum(dn * AU_TO_DEBYE) << ",\"s2\":" << jnum(r.s2) << ",\"history\":[";
    for (size_t i = 0; i < r.energyHistory.size(); ++i) j << (i ? "," : "") << jnum(r.energyHistory[i]);
    j << "],\"atoms\":[";
    for (size_t a = 0; a < mol.atoms.size(); ++a) {
        const Atom& at = mol.atoms[a];
        j << (a ? "," : "") << "{\"s\":" << jstr(elementSymbol(at.Z)) << ",\"Z\":" << at.Z << ",\"x\":" << jnum(at.x) << ",\"y\":" << jnum(at.y)
          << ",\"z\":" << jnum(at.z) << ",\"q\":" << jnum(r.charges[a]) << ",\"spin\":" << jnum(r.spinPop[a]) << ",\"r\":" << jnum(covalentRadiusAngstrom(at.Z) * ANGSTROM_TO_BOHR) << "}";
    }
    j << "],\"bonds\":[";
    bool first = true;
    for (size_t a = 0; a < mol.atoms.size(); ++a)
        for (size_t b = a + 1; b < mol.atoms.size(); ++b) {
            double order = r.bondOrders(int(a), int(b));
            if (order > 0.15) { j << (first ? "" : ",") << "[" << a << "," << b << "," << jnum(order) << "," << jnum(mol.distance(int(a), int(b)) / ANGSTROM_TO_BOHR) << "]"; first = false; }
        }
    j << "],\"orbitals\":[";
    for (size_t k = 0; k < orbs.size(); ++k) {
        const OrbitalRef& o = orbs[k];
        j << (k ? "," : "") << "{\"name\":" << jstr(o.name()) << ",\"spin\":" << o.spin << ",\"eps\":" << jnum(o.eps) << ",\"occ\":" << o.occ
          << ",\"tag\":" << jstr(o.tag) << ",\"comp\":[";
        for (size_t c = 0; c < o.comp.size(); ++c) j << (c ? "," : "") << "[" << jstr(o.comp[c].first) << "," << jnum(o.comp[c].second) << "]";
        j << "],\"sa\":" << jstr(sliceB64[0][k]) << ",\"sb\":" << jstr(sliceB64[1][k]) << ",\"pts\":[";
        for (size_t c = 0; c < clouds[k].size(); ++c) j << (c ? "," : "") << clouds[k][c];
        j << "]}";
    }
    j << "],\"density\":[" << jstr(densB64[0]) << "," << jstr(densB64[1]) << "]}";

    std::string html = HTML_REPORT_TEMPLATE;
    size_t pos = html.find("__DATA__");
    html.replace(pos, 8, j.str());
    std::string path = opt.outDir + "/rapport.html";
    std::ofstream(path) << html;
    return path;
}

}  // namespace mo
