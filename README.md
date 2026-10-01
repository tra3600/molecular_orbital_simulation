# molecular_orbital_simulation

Simulateur d'**orbitales moléculaires** en C++17 (sans dépendance externe) : les orbitales moléculaires (MO) sont construites par **combinaison linéaire d'orbitales atomiques (LCAO)** et calculées par un vrai **champ auto-cohérent Hartree–Fock** en base **STO-3G** (H → Ne), plus un mode **Hückel** pédagogique pour les systèmes π.

Chaque calcul produit un **rapport HTML interactif** : nuage 3D rotatif des orbitales, diagramme d'énergie cliquable, coupes 2D, charges, ordres de liaison, convergence du SCF.

## Compilation

```sh
make            # ou : cmake -S . -B build && cmake --build build
make tests      # 35 vérifications contre des valeurs de la littérature
```

Nécessite uniquement un compilateur C++17 (OpenMP utilisé s'il est disponible).

## Utilisation

```sh
./molorb list                      # molécules et systèmes π disponibles
./molorb H2O                       # SCF + rapport complet dans out/H2O/
./molorb scf C2H4 --grid 0.3       # pas de grille plus grossier (plus rapide)
./molorb scf molecule.xyz --charge 1
./molorb "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"   # géométrie en ligne (Å)
./molorb scan H2                   # courbe d'énergie RHF vs UHF (dissociation)
./molorb huckel benzene            # méthode de Hückel
./molorb gallery                   # tout calculer et créer out/index.html
```

Options de `scf` : `--uhf` (automatique si multiplicité > 1), `--charge Q`, `--mult M`, `--out DIR`, `--grid H`, `--samples N`, `--no-cubes`, `--no-html`, `-v`.

Ouvrez `out/<molécule>/rapport.html` dans un navigateur : cliquez sur un niveau d'énergie pour afficher l'orbitale correspondante (lobes + rouge, − bleu), faites-la tourner à la souris, zoomez à la molette, parcourez les orbitales avec ← / →.

### Fichiers produits (`out/<molécule>/`)

| Fichier | Contenu |
|---|---|
| `rapport.html` | rapport interactif autonome (aucune connexion nécessaire) |
| `mo_<n>.cube`, `density.cube` | orbitales et densité sur grille 3D, lisibles par Avogadro, VMD, Jmol… |
| `sliceA_mo<n>.png`, `sliceB_mo<n>.png` | coupes 2D signées avec lignes de contour (deux plans) |
| `energy_levels.svg` | diagramme d'énergie avec occupations |
| `resume.txt` | résumé texte : énergies, composition LCAO, charges, dipôle |

## Ce que calcule le programme

* **Base** : STO-3G (H–Ne), gaussiennes cartésiennes s et p contractées.
* **Intégrales** : recouvrement, cinétique, attraction nucléaire, répulsion électronique et dipôle par la méthode de McMurchie–Davidson (fonction de Boys, criblage de Schwarz, parallélisation OpenMP).
* **SCF** : RHF (couches fermées) et UHF (radicaux, brisure de symétrie), estimation initiale de Wolfsberg–Helmholz, accélération **DIIS**, orthogonalisation canonique.
* **Propriétés** : énergie totale, énergies orbitalaires, composition LCAO de chaque MO, gap HOMO–LUMO, charges et populations de spin de Mulliken, indices de liaison de Mayer, moment dipolaire, ⟨S²⟩.
* **Hückel** : énergies α + xβ, coefficients, énergie π, énergie de délocalisation, charges et ordres de liaison π (éthylène, allyle, butadiène, benzène, naphtalène, pyridine, pyrrole, cyclobutadiène…).

### Molécules disponibles

H₂, HeH⁺, H₃⁺, LiH, Li₂, HF, N₂, O₂ (triplet), F₂, CO, H₂O, NH₃, CH₄, BH₃, CO₂, HCN, C₂H₂, C₂H₄, H₂CO, C₆H₆ — ou n'importe quelle géométrie en fichier XYZ / en ligne (éléments H–Ne).

## Validation

| Système | E(HF/STO-3G) calculée | Référence |
|---|---|---|
| H₂ (1,4 bohr) | −1,11668 Eh | −1,1167 (Szabo & Ostlund) |
| LiH | −7,86202 Eh | −7,8620 |
| CH₄ | −39,72675 Eh | −39,7268 |
| N₂ | −107,49598 Eh | −107,496 |
| O₂ (³Σ, UHF) | −147,63397 Eh, ⟨S²⟩ = 2,003 | −147,634 |

Les tests vérifient aussi l'invariance par rotation/translation (fonctions p), la somme des charges de Mulliken, la limite de dissociation de H₂ en UHF (−0,9332 Eh = 2 atomes H) et les énergies de Hückel exactes (benzène : 2, 1, 1, −1, −1, −2 β).

## Idée physique

Une orbitale moléculaire s'écrit ψᵢ = Σ_μ C_μi φ_μ, où φ_μ sont les orbitales atomiques. Les coefficients C sont obtenus en résolvant les équations de Roothaan **F C = S C ε**, où F dépend elle-même de C (d'où l'itération auto-cohérente). Par exemple pour H₂ on retrouve l'orbitale liante σg = (φ₁ + φ₂)/√(2+2S) et l'antiliante σu* = (φ₁ − φ₂)/√(2−2S).

## Structure du code

`src/` : `molecule` (géométries, XYZ) · `basis` (STO-3G) · `integrals` (McMurchie–Davidson) · `scf` (RHF/UHF/DIIS) · `grid` (évaluation + .cube) · `png` (encodeur PNG/deflate) · `report` + `report_template.h` (HTML, SVG, console) · `huckel` · `main` (CLI). `tests/test_main.cpp` : tests de régression.

## Limites

Base minimale STO-3G uniquement, éléments H–Ne, pas d'optimisation de géométrie ni de corrélation électronique (les énergies HF/STO-3G sont qualitatives).
