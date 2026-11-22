# A polydrafter with isohedral number 9, and seven with isohedral number 6

The *isohedral number* of a shape that tiles the plane is the smallest k such that some tiling by
congruent copies of the shape has exactly k transitivity classes of tiles under the symmetry group
of the tiling (k = 1: the shape tiles isohedrally). J. Myers' tables of
[polyform tilings](https://www.polyomino.org.uk/mathematics/polyform-tiling/) give isohedral
numbers for polyominoes, polyhexes, polyiamonds and polykites; the values 6 and above are rare
(the tables contain every value from 1 to 10 except 7).

This repository reports isohedral numbers for **polydrafters** (unions of 30°–60°–90° triangles of
the kisrhombille grid, the triangles obtained by cutting each equilateral triangle of the triangular grid into six), a grid not
covered by those tables. Among the holeless free polydrafters with at most 22 drafters:

- one **20-drafter has isohedral number 9**;
- one **16-drafter** and six **20-drafters** have isohedral number 6. Two of these 20-drafters are
  unions of ten kites: they are very probably the two 10-kites of isohedral number 6 already in
  Myers' polykite table, seen on the finer grid. The shapes we believe to be new are therefore the
  20-drafter of number 9 and five shapes of number 6 (the 16-drafter and four 20-drafters), none of
  which is a polykite;
- all other polydrafters in that range that tile periodically have isohedral number at most 4;
- no shape of isohedral number 7 was found.

All statements are for tilings whose tiles are aligned with one common drafter grid (the convention
Myers uses for polykites).

> **Status.** Not checked by a third party; novelty not confirmed (unpublished computations may
> exist). The witness tilings (upper bounds) are re-checked by an independent program
> (`verification/`). The lower bounds (no tiling with fewer classes) rely on a single program,
> `search/iso.cpp`, supported by the controls below but not by a second implementation. Only
> grid-aligned tilings are considered; tilings with tiles off the grid are not excluded. The list of
> candidate shapes comes from a Heesch-number census made with
> [heesch-sat](https://github.com/isohedral/heesch-sat) and relies on its classification. Corrections
> and pointers to prior work are welcome: please open an issue.

*Résumé en français : un polydrafter de 20 drafters de nombre isoédrique 9 (tout pavage aligné sur
la grille a au moins 9 classes de tuiles, et un pavage en a exactement 9), et sept de nombre
isoédrique 6 (un de 16 drafters, six de 20, dont deux sont les deux polykites de 10 cerfs-volants
déjà connus de Myers). Les pavages témoins sont revérifiés par un programme indépendant ; les
bornes inférieures reposent sur un seul programme de recherche. Pas encore vérifié par un tiers ;
nouveauté non confirmée.*

Figures (French / English): https://totocaedres.fransys.io/pistes/isoedrique-7/notre-travail

## Files

- `results/witnesses.jsonl` — one JSON object per shape: `id`, `grid`, `n`, `cells` (coordinates
  of the drafters in the drafter grid of heesch-sat), `isohedral_number`, `union_of_kites`,
  `shape_symmetries` (order of the symmetry group of the shape), `search` (nodes explored and
  number of symmetry groups tried for k = 1, 2, …) and `witness`:
  - `lattice` — two vectors generating the translation lattice L of the witness tiling;
  - `group` — generators, modulo L, of a symmetry group G of the tiling;
  - `tiles` — one tile per class: the tile is the image of `cells` under the map.

  Every map is written `[a, b, c, d, e, f]` for `(x, y) ↦ (a x + b y + c, d x + e y + f)` in grid
  coordinates. The tiling is the set of images of the listed tiles under G and L.
- `results/shapes.txt` — the same eight shapes in the input format of `iso` (`drafter x y x y …`).
- `results/drafter_isohedral_numbers.csv` — for each size, the number of candidate shapes and how
  many have each isohedral number; last column: shapes with no periodic grid-aligned tiling with
  at most 12 classes.
- `verification/verify_witnesses.py` — independent check of the witnesses (Python 3, no dependency).
- `search/iso.cpp` — the search program (C++17, no dependency; identifiers and comments in French).
- `search/grilles.txt` — description of ten grids (lattice, symmetries, cells and adjacencies), read
  by both programs.

## Results

| Shape | Drafters | Isohedral number | Witness: tiles per period, symmetry group order mod translations | Polykite |
|---|---|---|---|---|
| `D20-k9` | 20 | **9** | 36, 4 | no |
| `D16-k6` | 16 | **6** | 6, 1 | no |
| `D20-k6-a` … `D20-k6-d` | 20 | **6** | 6, 1 | no |
| `D20-k6-e`, `D20-k6-f` | 20 | 6 | 6, 1 | yes (Myers' two 10-kites) |

Candidate shapes: the 266 holeless free polydrafters with at most 22 drafters that the Heesch census
(see [heesch-polykites](https://github.com/fransys-code/heesch-polykites)) left unclassified, that is,
neither isohedral tilers nor non-tilers with a known finite Heesch number: the 262 shapes counted
as inconclusive in its `census.csv` for these sizes, plus 4 shapes on which the classification
timed out (one 16-drafter and three 20-drafters). Their isohedral numbers:

| Isohedral number | 2 | 3 | 4 | 6 | 9 | none up to 12 |
|---|---|---|---|---|---|---|
| Shapes | 127 | 128 | 1 | 7 | 1 | 2 |

The two shapes with no periodic tiling with at most 12 classes are a 16-drafter and a 20-drafter
that are unions of 8 and 10 kites (the sizes of the aperiodic "hat" and "turtle" polykites).

## Method

A tiling with k classes is periodic and its symmetry group G is a plane crystallographic group;
for a grid-aligned tiling, G is a subgroup of the symmetry group of the grid. If G has point group P
and translation lattice L, the torus R²/L contains N = Σ |P|/sᵢ tiles, where sᵢ is the order of the
stabiliser of a tile of the i-th class. The index of L is therefore determined by k, and only
finitely many groups G have to be tried. For k = 1, 2, 3, … `iso`:

1. enumerates all subgroups G compatible with k, up to conjugacy (conjugacy classes of point
   groups P, P-invariant sublattices in Hermite normal form, cocycles P → Λ/L);
2. for each G, solves an exact cover problem (dancing links, depth k): columns are the G-orbits of
   cells of the torus, rows are the G-orbits of placements of the shape whose members are pairwise
   disjoint.

The first k with a solution is the isohedral number, and the solution is the witness. The search is
exhaustive for each k; there is no bound on the size of the period.

## Verification

```bash
python3 verification/verify_witnesses.py results/witnesses.jsonl --controls
```

The script checks that every map is a symmetry of the grid, that the tiles generated by the group
cover a period exactly once, and recounts the classes under the announced group and under the full
symmetry group of the tiling. With `--controls`, five corrupted copies of each witness (a class
removed, a tile shifted, a tile moved off the grid, a tile duplicated, a wrong number) must be
rejected. Output:

```
OK   D20-k9 20-drafter announced 9 {"valid": true, "tiles_per_cell": 36, "group_order": 4, "classes_under_given_group": 9, "tiles_per_primitive_cell": 36, "full_symmetry_order": 4, "classes_under_full_symmetry": 9}
OK   D16-k6 16-drafter announced 6 {"valid": true, "tiles_per_cell": 6, "group_order": 1, "classes_under_given_group": 6, "tiles_per_primitive_cell": 6, "full_symmetry_order": 1, "classes_under_full_symmetry": 6}
OK   D20-k6-a 20-drafter announced 6 {"valid": true, "tiles_per_cell": 6, "group_order": 1, "classes_under_given_group": 6, "tiles_per_primitive_cell": 6, "full_symmetry_order": 1, "classes_under_full_symmetry": 6}
… (same line for D20-k6-b to D20-k6-f)
8 witnesses checked, 0 failure(s)
40 corrupted witnesses, 40 rejected, 0 wrongly accepted
```

To recompute the isohedral numbers, lower bounds included (well under a second):

```bash
g++ -std=c++17 -O2 -o search/iso search/iso.cpp
search/iso -K 12 resout < results/shapes.txt
```

Each output line gives `k` (the isohedral number), `k_min_prouve` (no tiling with fewer classes)
and a witness (`temoin`). Other modes: `search/iso gen GRID N` lists the free polyforms of size N
of a grid (`omino`, `hex`, `iamond`, `kite`, `drafter`, …).

Controls of the search program against Myers' tables (July 2026 version), all passed:

- the numbers of polykites with each isohedral number ≥ 2, for 4 to 17 kites, among the candidate
  shapes of our polykite census, are identical to his table (largest value 6: two 10-kites and three
  14-kites), and his two aperiodic polykites have no periodic tiling up to 12 classes;
- complete censuses of all free polyforms give the same counts per isohedral number as his tables
  for polyominoes of size 8, polyhexes of sizes 6 and 7, polyiamonds of size 9 and polykites of
  sizes 4, 5, 6, 7 and 10 (for instance `search/iso gen kite 10 | search/iso -K 6 resout`: 12 shapes
  of number 3, 2 of number 6, 9,997 without tiling up to 6 classes; under a minute);
- seven anisohedral shapes taken from his lists are found with his value: a 6-hex (2), a 9-omino (4),
  a 13-omino (5), a 14-omino (6), two 10-hexes (8) and an 11-hex (9).

## Acknowledgements

The search and verification code was written with the help of Claude (Anthropic). The grids come
from C. S. Kaplan's heesch-sat; the reference tables are J. Myers'.

## Licence

Code: MIT (`LICENSE`). Data in `results/`: CC BY 4.0 (`LICENSE-DATA.md`).
