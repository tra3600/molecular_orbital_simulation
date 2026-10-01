// molorb — simulation d'orbitales moléculaires (LCAO) : Hartree–Fock/STO-3G et Hückel.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include "huckel.h"
#include "report.h"
#include "scf.h"

using namespace mo;
namespace fs = std::filesystem;

static void usage() {
    std::cout << R"(molorb — orbitales moléculaires par LCAO (Hartree–Fock / STO-3G, Hückel)

Usage :
  molorb list                          liste les molécules et systèmes π disponibles
  molorb scf <molécule> [options]      calcul SCF + rapport HTML/cube/PNG/SVG
  molorb scan <diatomique> [options]   courbe d'énergie potentielle (RHF vs UHF)
  molorb huckel <système> [--charge n] méthode de Hückel (π) + rapport HTML
  molorb gallery [--fast]              calcule toute la bibliothèque et crée out/index.html
  molorb <molécule>                    raccourci pour « scf »

<molécule> : nom de la bibliothèque (H2O, CH4, N2…), fichier .xyz (Å),
             ou description inline "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"

Options scf :
  --uhf              Hartree–Fock non restreint (obligatoire si multiplicité > 1, automatique)
  --charge Q         charge totale            --mult M   multiplicité de spin (2S+1)
  --out DIR          dossier de sortie (défaut : out/<nom>)
  --grid H           pas de grille en bohr (défaut 0.25 ; plus grand = plus rapide)
  --samples N        points du nuage 3D par orbitale (défaut 4000)
  --no-cubes         ne pas écrire les fichiers .cube        --no-html   pas de rapport HTML
  -v, --verbose      affiche les itérations SCF

Options scan :
  --from R0 --to R1 --steps N   distances interatomiques en Å (défaut 0.4 → 3.5, 32 points)
)";
}

struct Args {
    std::vector<std::string> pos;
    std::map<std::string, std::string> opt;
    bool flag(const std::string& k) const { return opt.count(k) > 0; }
    std::string get(const std::string& k, const std::string& d = "") const { auto it = opt.find(k); return it == opt.end() ? d : it->second; }
    double num(const std::string& k, double d) const { return opt.count(k) ? std::atof(opt.at(k).c_str()) : d; }
};

static Args parseArgs(int argc, char** argv, int from) {
    static const std::vector<std::string> withValue = {"--charge", "--mult", "--out", "--grid", "--samples", "--from", "--to", "--steps"};
    Args a;
    for (int i = from; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--", 0) == 0 || s == "-v") {
            if (s == "-v") s = "--verbose";
            if (std::find(withValue.begin(), withValue.end(), s) != withValue.end()) {
                if (i + 1 >= argc) throw std::runtime_error("Valeur manquante pour " + s);
                a.opt[s] = argv[++i];
            } else a.opt[s] = "1";
        } else a.pos.push_back(s);
    }
    return a;
}

static Molecule loadMolecule(const std::string& spec, const Args& a) {
    Molecule m;
    if (findBuiltin(spec, m)) {}
    else if (fs::exists(spec)) m = readXYZ(spec);
    else {
        m = parseAtoms(spec);
        if (m.atoms.empty()) throw std::runtime_error("Molécule inconnue : « " + spec + " » (voir « molorb list »)");
        m.name = m.formula();
    }
    if (a.flag("--charge")) m.charge = int(a.num("--charge", 0));
    if (a.flag("--mult")) m.multiplicity = int(a.num("--mult", 1));
    else if (a.flag("--charge") && m.nElectrons() % 2 != 0) m.multiplicity = 2;
    else if (a.flag("--charge") && m.nElectrons() % 2 == 0 && m.multiplicity % 2 == 0) m.multiplicity = 1;
    return m;
}

static std::string safeName(std::string s) {
    for (char& c : s) if (!std::isalnum(static_cast<unsigned char>(c)) && c != '+' && c != '-') c = '_';
    return s;
}

static SCFResult runAndReport(const Molecule& mol, const Args& a, bool quiet, std::string* htmlOut = nullptr) {
    Basis basis = buildBasis(mol);
    SCFOptions so;
    so.uhf = a.flag("--uhf");
    so.verbose = a.flag("--verbose");
    SCFResult r = runSCF(mol, basis, so);
    if (!quiet) printSummary(std::cout, mol, basis, r);
    OutputOptions oo;
    oo.outDir = a.get("--out", "out/" + safeName(mol.name));
    oo.spacing = a.num("--grid", 0.25);
    oo.samples = int(a.num("--samples", 4000));
    oo.cubes = !a.flag("--no-cubes");
    oo.html = !a.flag("--no-html");
    std::string html = writeOutputs(mol, basis, r, oo);
    if (!quiet) {
        std::cout << "\nFichiers écrits dans " << oo.outDir << "/ :\n";
        if (!html.empty()) std::cout << "  rapport.html          rapport interactif (nuage 3D rotatif, niveaux cliquables)\n";
        std::cout << "  energy_levels.svg     diagramme d'énergie\n  slice*_mo*.png        coupes 2D de chaque orbitale\n  resume.txt            ce résumé\n";
        if (oo.cubes) std::cout << "  mo_*.cube, density.cube   orbitales et densité 3D (Avogadro, VMD, Jmol…)\n";
    }
    if (htmlOut) *htmlOut = html;
    return r;
}

static int cmdScan(const Args& a) {
    if (a.pos.empty()) throw std::runtime_error("scan : indiquer une molécule diatomique");
    Molecule base = loadMolecule(a.pos[0], a);
    if (base.atoms.size() != 2) throw std::runtime_error("scan : la molécule doit être diatomique");
    double r0 = a.num("--from", 0.4), r1 = a.num("--to", 3.5);
    int steps = std::max(3, int(a.num("--steps", 32)));
    std::vector<double> R, Er, Eu;
    bool singlet = base.multiplicity == 1;
    std::cout << "Courbe d'énergie potentielle de " << base.name << " (STO-3G)\n";
    std::cout << "   R (Å)        E(RHF)           E(UHF)\n";
    std::ofstream csv;
    std::string dir = a.get("--out", "out/scan_" + safeName(base.name));
    fs::create_directories(dir);
    csv.open(dir + "/scan.csv");
    csv << "R_angstrom,E_RHF,E_UHF\n";
    for (int i = 0; i < steps; ++i) {
        double R_ = r0 + (r1 - r0) * i / (steps - 1);
        Molecule m = base;
        m.atoms[0].x = m.atoms[0].y = m.atoms[0].z = 0;
        m.atoms[1].x = m.atoms[1].y = 0;
        m.atoms[1].z = R_ * ANGSTROM_TO_BOHR;
        Basis b = buildBasis(m);
        SCFOptions so;
        double er = std::nan(""), eu;
        if (singlet) { so.uhf = false; SCFResult rr = runSCF(m, b, so); if (rr.converged) er = rr.energy; }
        so.uhf = true;
        SCFResult ru = runSCF(m, b, so);
        eu = ru.converged ? ru.energy : std::nan("");
        if (!singlet) er = eu;
        R.push_back(R_); Er.push_back(er); Eu.push_back(eu);
        std::cout << std::fixed << std::setprecision(4) << std::setw(8) << R_ << std::setprecision(8) << std::setw(18) << er << std::setw(17) << eu << "\n";
        csv << R_ << "," << er << "," << eu << "\n";
    }
    int im = 0;
    for (size_t i = 0; i < R.size(); ++i) if (!std::isnan(Eu[i]) && (Eu[im] > Eu[i] || std::isnan(Eu[im]))) im = int(i);
    std::cout << std::setprecision(4) << "\nMinimum (UHF) vers R = " << R[im] << " Å, E = " << std::setprecision(6) << Eu[im] << " Eh\n";
    Series sr{"RHF", "#c0392b", {}, {}, true}, su{"UHF", "#2c6fbb", {}, {}, true};
    for (size_t i = 0; i < R.size(); ++i) {
        if (!std::isnan(Er[i])) { sr.x.push_back(R[i]); sr.y.push_back(Er[i]); }
        if (!std::isnan(Eu[i])) { su.x.push_back(R[i]); su.y.push_back(Eu[i]); }
    }
    std::vector<Series> ss = singlet ? std::vector<Series>{sr, su} : std::vector<Series>{su};
    std::ofstream(dir + "/scan.svg") << svgPlot("Courbe d'énergie potentielle — " + base.name + " (STO-3G)", "distance interatomique (Å)", "énergie (hartree)", ss);
    std::cout << "Écrits : " << dir << "/scan.csv, " << dir << "/scan.svg\n";
    if (singlet) std::cout << "Note : au-delà du point de Coulson–Fischer l'UHF (symétrie brisée) passe sous le RHF, qui dissocie mal en fragments ioniques.\n";
    return 0;
}

static int cmdHuckel(const Args& a) {
    if (a.pos.empty()) throw std::runtime_error("huckel : indiquer un système (voir « molorb list »)");
    HuckelSystem s;
    if (!findHuckel(a.pos[0], s)) throw std::runtime_error("Système de Hückel inconnu : " + a.pos[0]);
    s.charge = int(a.num("--charge", 0));
    HuckelResult r = solveHuckel(s);
    printHuckel(std::cout, s, r);
    std::string dir = a.get("--out", "out/huckel_" + s.name);
    fs::create_directories(dir);
    std::ofstream(dir + "/huckel.html") << huckelHTML(s, r);
    std::cout << "\nRapport : " << dir << "/huckel.html\n";
    return 0;
}

static int cmdList() {
    std::cout << "Molécules (Hartree–Fock / STO-3G) :\n";
    for (const Molecule& m : builtinMolecules())
        std::cout << "  " << std::left << std::setw(8) << m.name << m.description << "\n";
    std::cout << "\nSystèmes π (Hückel) :\n";
    for (auto& s : huckelLibrary()) std::cout << "  " << std::left << std::setw(18) << s.name << s.description << "\n";
    return 0;
}

static int cmdGallery(const Args& a) {
    std::string root = a.get("--out", "out");
    fs::create_directories(root);
    std::ostringstream cards;
    for (const Molecule& m : builtinMolecules()) {
        if (a.flag("--fast") && m.name == "C6H6") continue;
        std::cout << "→ " << m.name << " ..." << std::flush;
        Args b = a;
        b.opt["--out"] = root + "/" + safeName(m.name);
        b.opt["--no-cubes"] = "1";
        SCFResult r = runAndReport(m, b, true);
        std::cout << " E = " << std::fixed << std::setprecision(6) << r.energy << " Eh\n";
        cards << "<a class=c href=\"" << safeName(m.name) << "/rapport.html\"><b>" << m.formula() << "</b><span>" << m.description
              << "</span><i>E = " << std::fixed << std::setprecision(5) << r.energy << " Eh · " << (r.uhf ? "UHF" : "RHF") << "</i></a>\n";
    }
    std::ostringstream hu;
    for (auto& s : huckelLibrary()) {
        HuckelResult r = solveHuckel(s);
        std::string d = root + "/huckel_" + s.name;
        fs::create_directories(d);
        std::ofstream(d + "/huckel.html") << huckelHTML(s, r);
        hu << "<a class=c href=\"huckel_" << s.name << "/huckel.html\"><b>" << s.name << "</b><span>" << s.description << "</span></a>\n";
    }
    std::ofstream(root + "/index.html") << R"(<!doctype html><html lang="fr"><meta charset="utf-8"><title>Galerie d'orbitales moléculaires</title>
<style>body{background:#0f1420;color:#e6ebf5;font:15px system-ui,sans-serif;margin:0;padding:24px}h1{margin:0 0 4px}h2{color:#7cc4ff;margin-top:26px}
p{color:#8b97b3}.g{display:grid;grid-template-columns:repeat(auto-fill,minmax(250px,1fr));gap:12px}
.c{background:#171e2e;border:1px solid #2a3550;border-radius:12px;padding:14px;color:inherit;text-decoration:none;display:flex;flex-direction:column;gap:4px}
.c:hover{border-color:#7cc4ff}.c b{font-size:19px}.c span{color:#8b97b3;font-size:13px}.c i{font-size:12px;color:#7cc4ff;font-style:normal}</style>
<h1>Galerie d'orbitales moléculaires</h1><p>Hartree–Fock / STO-3G et méthode de Hückel — cliquez sur une molécule pour explorer ses orbitales en 3D.</p>
<h2>Hartree–Fock / STO-3G</h2><div class=g>)" << cards.str() << "</div><h2>Méthode de Hückel (π)</h2><div class=g>" << hu.str() << "</div></html>\n";
    std::cout << "\nGalerie : " << root << "/index.html\n";
    return 0;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) { usage(); return 0; }
        std::string cmd = argv[1];
        if (cmd == "-h" || cmd == "--help" || cmd == "help") { usage(); return 0; }
        if (cmd == "list") return cmdList();
        if (cmd == "scan") return cmdScan(parseArgs(argc, argv, 2));
        if (cmd == "huckel") return cmdHuckel(parseArgs(argc, argv, 2));
        if (cmd == "gallery") return cmdGallery(parseArgs(argc, argv, 2));
        Args a = parseArgs(argc, argv, cmd == "scf" ? 2 : 1);
        if (a.pos.empty()) { usage(); return 1; }
        Molecule m = loadMolecule(a.pos[0], a);
        runAndReport(m, a, false);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << "\n";
        return 1;
    }
}
