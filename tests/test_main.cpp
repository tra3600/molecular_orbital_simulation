// Tests de non-régression : valeurs de référence de la littérature (HF/STO-3G).
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include "huckel.h"
#include "scf.h"

using namespace mo;

static int failures = 0, checks = 0;
static void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::printf("  ECHEC : %s\n", what.c_str()); }
}
static void near(double a, double b, double tol, const std::string& what) {
    char buf[200];
    std::snprintf(buf, sizeof buf, "%s : %.8f attendu %.8f (tol %.1e)", what.c_str(), a, b, tol);
    check(std::fabs(a - b) <= tol, buf);
}

static SCFResult scf(const Molecule& m, bool uhf = false) {
    SCFOptions o;
    o.uhf = uhf;
    return runSCF(m, buildBasis(m), o);
}
static Molecule lib(const std::string& n) {
    Molecule m;
    if (!findBuiltin(n, m)) { std::printf("molécule absente : %s\n", n.c_str()); std::exit(2); }
    return m;
}

int main() {
    // Fonction de Boys
    { double F[3]; boysFunction(2, 1.0, F);
      near(F[0], 0.7468241328, 1e-9, "F0(1)"); near(F[1], 0.1894723458, 1e-9, "F1(1)");
      boysFunction(2, 0.0, F); near(F[2], 0.2, 1e-12, "F2(0)");
      boysFunction(1, 50.0, F); near(F[0], 0.5 * std::sqrt(M_PI / 50.0), 1e-12, "F0(50)"); }

    // H2 à 1,4 bohr (Szabo & Ostlund) : S12 = 0.6593, E = −1.1167
    { Molecule h2 = parseAtoms("H 0 0 0; H 0 0 " + std::to_string(1.4 / ANGSTROM_TO_BOHR));
      Basis b = buildBasis(h2);
      Mat S = overlapMatrix(b);
      near(S(0, 0), 1.0, 1e-10, "S11"); near(S(0, 1), 0.6593, 5e-4, "S12 H2");
      SCFResult r = scf(h2);
      check(r.converged, "H2 converge");
      near(r.energy, -1.1167, 5e-4, "E(H2)");
      near(r.epsA[0], -0.5782, 1e-3, "eps(σg) H2");
      near(r.bondOrders(0, 1), 1.0, 1e-6, "ordre de liaison H2"); }

    // Énergies totales de référence (géométries de la bibliothèque)
    near(scf(lib("LiH")).energy, -7.8620, 1e-3, "E(LiH)");
    near(scf(lib("CH4")).energy, -39.7267, 1e-3, "E(CH4)");
    near(scf(lib("NH3")).energy, -55.4540, 1e-3, "E(NH3)");
    near(scf(lib("H2O")).energy, -74.9629, 1e-3, "E(H2O)");
    near(scf(lib("N2")).energy, -107.4960, 1e-3, "E(N2)");
    near(scf(lib("C2H4")).energy, -77.0721, 1e-3, "E(C2H4)");

    // Invariance par rotation + translation (teste les fonctions p)
    { Molecule a = lib("H2O"), b = a;
      double c = std::cos(0.7), s = std::sin(0.7);
      for (Atom& t : b.atoms) {
          double x = t.x, y = t.y, z = t.z;
          t.x = c * x - s * z + 1.3; t.z = s * x + c * z - 0.4; t.y = y + 0.2;
      }
      near(scf(a).energy, scf(b).energy, 1e-8, "invariance rotationnelle H2O");
      SCFResult ra = scf(a), rb = scf(b);
      double da = std::hypot(ra.dipole[0], ra.dipole[1], ra.dipole[2]), db = std::hypot(rb.dipole[0], rb.dipole[1], rb.dipole[2]);
      near(da, db, 1e-7, "norme du dipôle invariante");
      near(da * AU_TO_DEBYE, 1.73, 0.05, "dipôle H2O (D)"); }

    // Charges de Mulliken : somme = charge totale
    { Molecule m = lib("HeH+"); SCFResult r = scf(m);
      double q = 0; for (double v : r.charges) q += v;
      near(q, 1.0, 1e-8, "somme des charges HeH+"); }
    { SCFResult r = scf(lib("H2O")); double q = 0; for (double v : r.charges) q += v;
      near(q, 0.0, 1e-8, "somme des charges H2O"); check(r.charges[0] < 0, "O négatif dans H2O"); }

    // UHF
    { Molecule h2 = lib("H2");
      near(scf(h2, true).energy, scf(h2, false).energy, 1e-7, "UHF = RHF pour H2 à l'équilibre");
      Molecule far = parseAtoms("H 0 0 0; H 0 0 3.0");
      double eu = scf(far, true).energy, er = scf(far, false).energy;
      check(eu < er - 0.1, "UHF sous le RHF à grande distance");
      near(eu, -0.9332, 2e-3, "limite de dissociation (2 H)");
      SCFResult t = scf(lib("O2"));
      check(t.uhf, "O2 triplet utilise UHF"); near(t.s2, 2.0, 0.02, "<S²> O2 triplet");
      near(t.energy, -147.634, 1e-3, "E(O2 triplet)"); }

    // Hückel
    { HuckelSystem s; findHuckel("benzene", s); HuckelResult r = solveHuckel(s);
      near(r.x[0], 2, 1e-9, "benzène x1"); near(r.x[5], -2, 1e-9, "benzène x6");
      near(r.piEnergy, 8.0, 1e-9, "énergie π benzène"); near(r.delocalization, 2.0, 1e-9, "délocalisation benzène");
      findHuckel("butadiene", s); r = solveHuckel(s);
      near(r.x[0], 1.6180339887, 1e-8, "butadiène x1"); near(r.piEnergy, 4.4721359550, 1e-8, "énergie π butadiène");
      near(r.bondOrders(0, 1), 0.8944271910, 1e-8, "ordre π C1–C2 butadiène"); }

    std::printf("%d vérifications, %d échec(s)\n", checks, failures);
    return failures ? 1 : 0;
}
