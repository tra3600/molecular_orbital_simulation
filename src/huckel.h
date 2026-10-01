// Méthode de Hückel simple (systèmes π) : version pédagogique du LCAO.
#pragma once
#include <string>
#include <vector>

#include "linalg.h"

namespace mo {

struct HuckelAtom {
    std::string sym = "C";
    double x = 0, y = 0;  // coordonnées 2D de dessin
    double h = 0;         // α_i = α + h β
    int electrons = 1;    // électrons π apportés
};
struct HuckelBond { int a, b; double k = 1.0; };  // β_ab = k β

struct HuckelSystem {
    std::string name, description;
    std::vector<HuckelAtom> atoms;
    std::vector<HuckelBond> bonds;
    int charge = 0;
};

struct HuckelResult {
    std::vector<double> x;        // énergies E = α + xβ (x décroissant = énergie croissante)
    Mat C;                        // coefficients (colonne = MO)
    std::vector<double> occ;
    double piEnergy = 0;          // en unités de β
    double delocalization = 0;    // vs. liaisons π localisées (système tout carbone)
    bool hasDelocalization = false;
    std::vector<double> charges;  // charge π nette
    Mat bondOrders;
};

std::vector<HuckelSystem> huckelLibrary();
bool findHuckel(const std::string& name, HuckelSystem& out);
HuckelResult solveHuckel(const HuckelSystem& s);
void printHuckel(std::ostream& os, const HuckelSystem& s, const HuckelResult& r);
std::string huckelHTML(const HuckelSystem& s, const HuckelResult& r);

}  // namespace mo
