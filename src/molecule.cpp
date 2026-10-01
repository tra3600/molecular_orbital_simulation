#include "molecule.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace mo {

static const char* SYMBOLS[] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"};
static const double RADII[] = {0, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58};

std::string elementSymbol(int Z) { return (Z >= 0 && Z <= 10) ? SYMBOLS[Z] : "?"; }
double covalentRadiusAngstrom(int Z) { return (Z >= 1 && Z <= 10) ? RADII[Z] : 0.8; }

int atomicNumber(const std::string& s) {
    std::string t = s;
    if (!t.empty()) t[0] = static_cast<char>(std::toupper(t[0]));
    for (size_t i = 1; i < t.size(); ++i) t[i] = static_cast<char>(std::tolower(t[i]));
    for (int z = 1; z <= 10; ++z)
        if (t == SYMBOLS[z]) return z;
    if (!t.empty() && std::isdigit(static_cast<unsigned char>(t[0]))) {
        int z = std::atoi(t.c_str());
        if (z >= 1 && z <= 10) return z;
    }
    return -1;
}

int Molecule::nElectrons() const {
    int n = -charge;
    for (const Atom& a : atoms) n += a.Z;
    return n;
}

double Molecule::distance(int a, int b) const {
    const Atom &p = atoms[a], &q = atoms[b];
    return std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
}

double Molecule::nuclearRepulsion() const {
    double e = 0;
    for (size_t i = 0; i < atoms.size(); ++i)
        for (size_t j = i + 1; j < atoms.size(); ++j)
            e += double(atoms[i].Z) * atoms[j].Z / distance(int(i), int(j));
    return e;
}

std::string Molecule::formula() const {
    std::map<int, int> count;
    for (const Atom& a : atoms) count[a.Z]++;
    std::string s;
    auto emit = [&](int z) {
        s += elementSymbol(z);
        if (count[z] > 1) s += std::to_string(count[z]);
        count.erase(z);
    };
    // Convention de Hill : C, H, puis ordre alphabétique.
    if (count.count(6)) {
        emit(6);
        if (count.count(1)) emit(1);
    }
    std::vector<std::pair<std::string, int>> rest;
    for (auto& kv : count) rest.push_back({elementSymbol(kv.first), kv.first});
    std::sort(rest.begin(), rest.end());
    for (auto& r : rest) emit(r.second);
    if (charge > 0) s += "^" + (charge > 1 ? std::to_string(charge) : std::string()) + "+";
    if (charge < 0) s += "^" + (charge < -1 ? std::to_string(-charge) : std::string()) + "-";
    return s;
}

Molecule parseAtoms(const std::string& spec) {
    Molecule m;
    std::string s = spec;
    for (char& c : s)
        if (c == ';' || c == ',') c = '\n';
    std::istringstream in(s);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string sym;
        double x, y, z;
        if (!(ls >> sym >> x >> y >> z)) continue;
        int Z = atomicNumber(sym);
        if (Z < 0) throw std::runtime_error("Élément non supporté (H–Ne uniquement) : " + sym);
        m.atoms.push_back({Z, x * ANGSTROM_TO_BOHR, y * ANGSTROM_TO_BOHR, z * ANGSTROM_TO_BOHR});
    }
    return m;
}

Molecule readXYZ(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Impossible d'ouvrir " + path);
    std::string l1, l2;
    std::getline(f, l1);
    std::getline(f, l2);
    int n = std::atoi(l1.c_str());
    std::string rest, line;
    while (std::getline(f, line)) rest += line + "\n";
    // les lignes sont séparées par '\n' ; parseAtoms accepte aussi ';'
    Molecule m = parseAtoms(rest);
    if (n > 0 && int(m.atoms.size()) != n) throw std::runtime_error("Fichier XYZ incohérent : " + path);
    auto grab = [&](const char* key, int& out) {
        size_t p = l2.find(key);
        if (p != std::string::npos) out = std::atoi(l2.c_str() + p + std::string(key).size());
    };
    grab("charge=", m.charge);
    grab("mult=", m.multiplicity);
    size_t slash = path.find_last_of('/');
    std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    m.name = dot == std::string::npos ? base : base.substr(0, dot);
    return m;
}

namespace {
struct Entry {
    const char *name, *desc;
    int charge, mult;
    const char* atoms;
};
const Entry LIB[] = {
    {"H2", "Dihydrogène (liaison σ)", 0, 1, "H 0 0 0; H 0 0 0.7414"},
    {"HeH+", "Hydrure d'hélium (cation)", 1, 1, "He 0 0 0; H 0 0 0.7743"},
    {"H3+", "Cation trihydrogène (triangle équilatéral)", 1, 1,
     "H 0.4500 0 0; H -0.2250 0.3897 0; H -0.2250 -0.3897 0"},
    {"LiH", "Hydrure de lithium (liaison très polaire)", 0, 1, "Li 0 0 0; H 0 0 1.595"},
    {"Li2", "Dilithium", 0, 1, "Li 0 0 0; Li 0 0 2.673"},
    {"HF", "Fluorure d'hydrogène", 0, 1, "F 0 0 0; H 0 0 0.917"},
    {"N2", "Diazote (triple liaison)", 0, 1, "N 0 0 0; N 0 0 1.098"},
    {"O2", "Dioxygène (état triplet, radical diatomique)", 0, 3, "O 0 0 0; O 0 0 1.208"},
    {"F2", "Difluor", 0, 1, "F 0 0 0; F 0 0 1.412"},
    {"CO", "Monoxyde de carbone", 0, 1, "C 0 0 0; O 0 0 1.128"},
    {"H2O", "Eau", 0, 1, "O 0 0 0; H 0.7570 0.5859 0; H -0.7570 0.5859 0"},
    {"NH3", "Ammoniac", 0, 1,
     "N 0 0 0; H 0.9375 0 -0.3810; H -0.4688 0.8119 -0.3810; H -0.4688 -0.8119 -0.3810"},
    {"CH4", "Méthane (tétraèdre)", 0, 1,
     "C 0 0 0; H 0.6287 0.6287 0.6287; H 0.6287 -0.6287 -0.6287; H -0.6287 0.6287 -0.6287; H -0.6287 -0.6287 0.6287"},
    {"BH3", "Borane (trigonal plan, lacune électronique)", 0, 1,
     "B 0 0 0; H 1.19 0 0; H -0.5950 1.0306 0; H -0.5950 -1.0306 0"},
    {"CO2", "Dioxyde de carbone (linéaire)", 0, 1, "C 0 0 0; O 0 0 1.162; O 0 0 -1.162"},
    {"HCN", "Cyanure d'hydrogène", 0, 1, "H 0 0 -1.064; C 0 0 0; N 0 0 1.156"},
    {"C2H2", "Éthyne (acétylène)", 0, 1, "C 0.6015 0 0; C -0.6015 0 0; H 1.6645 0 0; H -1.6645 0 0"},
    {"C2H4", "Éthène (éthylène, liaison π)", 0, 1,
     "C 0.6695 0 0; C -0.6695 0 0; H 1.2342 0.9288 0; H 1.2342 -0.9288 0; H -1.2342 0.9288 0; H -1.2342 -0.9288 0"},
    {"H2CO", "Méthanal (formaldéhyde)", 0, 1, "C 0 0 0; O 0 0 1.205; H 0.9427 0 -0.5879; H -0.9427 0 -0.5879"},
    {"C6H6", "Benzène (36 fonctions de base, calcul plus long)", 0, 1,
     "C 1.3970 0 0; C 0.6985 1.2098 0; C -0.6985 1.2098 0; C -1.3970 0 0; C -0.6985 -1.2098 0; C 0.6985 -1.2098 0;"
     "H 2.4810 0 0; H 1.2405 2.1486 0; H -1.2405 2.1486 0; H -2.4810 0 0; H -1.2405 -2.1486 0; H 1.2405 -2.1486 0"},
};
}  // namespace

std::vector<Molecule> builtinMolecules() {
    std::vector<Molecule> v;
    for (const Entry& e : LIB) {
        Molecule m = parseAtoms(e.atoms);
        m.name = e.name;
        m.description = e.desc;
        m.charge = e.charge;
        m.multiplicity = e.mult;
        v.push_back(m);
    }
    return v;
}

bool findBuiltin(const std::string& name, Molecule& out) {
    auto up = [](std::string s) {
        for (char& c : s) c = static_cast<char>(std::toupper(c));
        return s;
    };
    for (const Molecule& m : builtinMolecules())
        if (up(m.name) == up(name)) {
            out = m;
            return true;
        }
    return false;
}

}  // namespace mo
