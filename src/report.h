// Sorties : résumé console, fichiers .cube, coupes PNG, diagrammes SVG,
// rapport HTML interactif (nuage 3D rotatif, diagramme d'énergie cliquable).
#pragma once
#include <iosfwd>
#include <string>
#include <vector>

#include "basis.h"
#include "molecule.h"
#include "scf.h"

namespace mo {

struct OutputOptions {
    std::string outDir = "out";
    double spacing = 0.25;  // pas de grille (bohr)
    double margin = 4.0;    // marge autour de la molécule (bohr)
    int samples = 4000;     // points du nuage 3D par orbitale
    int slicePx = 320;      // taille des coupes 2D
    bool cubes = true;
    bool html = true;
};

struct OrbitalRef {
    int spin;   // 0 = RHF (α=β), 1 = α, 2 = β
    int index;  // numéro (0-based) dans la liste de ce spin
    double eps;
    double occ;
    std::string tag;  // HOMO, LUMO, ...
    std::vector<std::pair<std::string, double>> comp;  // contributions AO dominantes
    std::string name() const;
};

std::vector<OrbitalRef> listOrbitals(const Basis& basis, const SCFResult& r);

void printSummary(std::ostream& os, const Molecule& mol, const Basis& basis, const SCFResult& r);
std::string energyDiagramSVG(const Molecule& mol, const std::vector<OrbitalRef>& orbs, bool uhf);

// Génère tous les fichiers dans outDir ; retourne le chemin du rapport HTML (ou "").
std::string writeOutputs(const Molecule& mol, const Basis& basis, const SCFResult& r, const OutputOptions& opt);

struct Series {
    std::string name, color;
    std::vector<double> x, y;
    bool markers = true;
};
std::string svgPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel,
                    const std::vector<Series>& series);

}  // namespace mo
