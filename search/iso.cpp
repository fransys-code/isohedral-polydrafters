// Nombre isoédrique d'une polyforme : plus petit k tel que la forme admette un pavage du plan, aligné sur la grille,
// dont les tuiles forment k orbites sous le groupe de symétrie du pavage.
// Méthode : un pavage à k orbites est invariant par un groupe cristallographique plan G, sous-groupe du groupe de
// symétrie de la grille, dont le domaine fondamental a une aire fixée par k. On énumère tous ces G à conjugaison près
// et, pour chacun, on résout une couverture exacte des orbites de cellules par des orbites de placements.
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
typedef pair<int, int> pii;
typedef array<int, 3> cel;

static int mod(int a, int m) { int r = a % m; return r < 0 ? r + m : r; }
static int divf(int a, int m) { return (a - mod(a, m)) / m; }

struct Aff { int a, b, c, d, e, f; };

struct Grille {
	string nom; char abrev;
	int v1x, v1y, v2x, v2y, det, no, c, id;
	vector<Aff> ori;
	vector<pii> reps; vector<int> type; vector<vector<pii>> vois;
	map<pii, int> residu;
	vector<array<int, 4>> ML;
	vector<vector<cel>> comp, imgRep;
	vector<int> inv;

	pii coordsReseau(int x, int y, bool& ok) const {
		int A = x * v2y - y * v2x, B = v1x * y - v1y * x;
		ok = (mod(A, det) == 0 && mod(B, det) == 0);
		return { A / det, B / det };
	}
	bool decompose(int x, int y, cel& out) const {
		int A = x * v2y - y * v2x, B = v1x * y - v1y * x;
		auto it = residu.find({ mod(A, det), mod(B, det) });
		if (it == residu.end()) return false;
		int r = it->second, Ar = reps[r].first * v2y - reps[r].second * v2x, Br = v1x * reps[r].second - v1y * reps[r].first;
		out = { r, (A - Ar) / det, (B - Br) / det };
		return true;
	}
	pii xy(const cel& q) const { return { reps[q[0]].first + q[1] * v1x + q[2] * v2x, reps[q[0]].second + q[1] * v1y + q[2] * v2y }; }
	pii lin(int o, int i, int j) const { return { ML[o][0] * i + ML[o][1] * j, ML[o][2] * i + ML[o][3] * j }; }
	cel image(int o, const cel& q) const {
		const cel& b = imgRep[o][q[0]]; pii l = lin(o, q[1], q[2]);
		return { b[0], b[1] + l.first, b[2] + l.second };
	}
	void prepare() {
		det = v1x * v2y - v1y * v2x;
		if (det <= 0) { cerr << "réseau mal orienté\n"; exit(1); }
		c = reps.size(); no = ori.size();
		for (int r = 0; r < c; ++r) {
			int A = reps[r].first * v2y - reps[r].second * v2x, B = v1x * reps[r].second - v1y * reps[r].first;
			residu[{ mod(A, det), mod(B, det) }] = r;
		}
		bool ok;
		for (auto& T : ori) {
			pii c1 = coordsReseau(T.a * v1x + T.b * v1y, T.d * v1x + T.e * v1y, ok); if (!ok) { cerr << "orientation hors réseau\n"; exit(1); }
			pii c2 = coordsReseau(T.a * v2x + T.b * v2y, T.d * v2x + T.e * v2y, ok); if (!ok) { cerr << "orientation hors réseau\n"; exit(1); }
			ML.push_back({ c1.first, c2.first, c1.second, c2.second });
		}
		comp.assign(no, vector<cel>(no)); inv.assign(no, -1);
		for (int i = 0; i < no; ++i) for (int j = 0; j < no; ++j) {
			const Aff& S = ori[i]; const Aff& T = ori[j];
			int a = S.a * T.a + S.b * T.d, b = S.a * T.b + S.b * T.e, d = S.d * T.a + S.e * T.d, e = S.d * T.b + S.e * T.e;
			int cx = S.a * T.c + S.b * T.f + S.c, cy = S.d * T.c + S.e * T.f + S.f, k = -1;
			for (int l = 0; l < no; ++l) if (ori[l].a == a && ori[l].b == b && ori[l].d == d && ori[l].e == e) { if (k >= 0) { cerr << "orientations en double\n"; exit(1); } k = l; }
			if (k < 0) { cerr << "orientations non closes\n"; exit(1); }
			pii t = coordsReseau(cx - ori[k].c, cy - ori[k].f, ok); if (!ok) { cerr << "composition hors réseau\n"; exit(1); }
			comp[i][j] = { k, t.first, t.second };
		}
		int ident = -1;
		for (int l = 0; l < no; ++l) if (ori[l].a == 1 && ori[l].b == 0 && ori[l].d == 0 && ori[l].e == 1) ident = l;
		if (ident != 0 || ori[0].c != 0 || ori[0].f != 0) { cerr << "l'identité doit être l'orientation 0\n"; exit(1); }
		for (int i = 0; i < no; ++i) for (int j = 0; j < no; ++j) if (comp[i][j][0] == 0) inv[i] = j;
		imgRep.assign(no, vector<cel>(c));
		for (int o = 0; o < no; ++o) for (int r = 0; r < c; ++r) {
			const Aff& T = ori[o];
			if (!decompose(T.a * reps[r].first + T.b * reps[r].second + T.c, T.d * reps[r].first + T.e * reps[r].second + T.f, imgRep[o][r])) { cerr << "orientation hors grille\n"; exit(1); }
		}
	}
};

static vector<Grille> lisGrilles(const string& chemin) {
	ifstream f(chemin);
	if (!f) { cerr << "impossible de lire " << chemin << "\n"; exit(1); }
	vector<Grille> gs; string mot;
	while (f >> mot) {
		if (mot == "grille") { gs.emplace_back(); f >> gs.back().nom >> gs.back().abrev; gs.back().id = gs.size() - 1; }
		else if (mot == "page") { double x; f >> x >> x >> x >> x; }
		else if (mot == "reseau") { Grille& g = gs.back(); f >> g.v1x >> g.v1y >> g.v2x >> g.v2y; }
		else if (mot == "indice_reseau_heesch_sat") { int x; f >> x; }
		else if (mot == "orientations") { int n; f >> n; Grille& g = gs.back(); g.ori.resize(n); for (auto& T : g.ori) f >> T.a >> T.b >> T.c >> T.d >> T.e >> T.f; }
		else if (mot == "cellules") {
			int n; f >> n; Grille& g = gs.back();
			for (int i = 0; i < n; ++i) {
				int x, y, t, nv; f >> x >> y >> t >> nv;
				g.reps.push_back({ x, y }); g.type.push_back(t); g.vois.emplace_back();
				for (int j = 0; j < nv; ++j) { int dx, dy; f >> dx >> dy; g.vois.back().push_back({ dx, dy }); }
			}
		}
	}
	for (auto& g : gs) g.prepare();
	return gs;
}

struct Lat {
	int a, b, d;
	int m() const { return a * d; }
	int idx(int i, int j) const { int q = divf(j, d); j -= q * d; i -= q * b; return mod(i, a) + a * j; }
	bool contient(int i, int j) const { return mod(j, d) == 0 && mod(i - (j / d) * b, a) == 0; }
	bool operator<(const Lat& o) const { return tie(a, b, d) < tie(o.a, o.b, o.d); }
	bool operator==(const Lat& o) const { return a == o.a && b == o.b && d == o.d; }
};

static Lat hnf(int x1, int y1, int x2, int y2) {
	while (y2 != 0) { int q = y1 / y2; x1 -= q * x2; y1 -= q * y2; swap(x1, x2); swap(y1, y2); }
	if (y1 < 0) { x1 = -x1; y1 = -y1; }
	Lat L; L.a = abs(x2); L.d = y1; L.b = mod(x1, L.a);
	return L;
}

struct ClasseP { vector<int> elems, gens, normalisateur; };

struct Groupe { int pc; Lat L; vector<int> lam; };

struct Enumerateur {
	const Grille& g;
	vector<ClasseP> classes;
	map<pii, vector<Groupe>> cache;

	int ferme(const vector<int>& gens) const {
		int masque = 1; bool change = true;
		while (change) {
			change = false;
			for (int x = 0; x < g.no; ++x) if (masque >> x & 1) for (int y : gens) {
				int z = g.comp[y][x][0];
				if (!(masque >> z & 1)) { masque |= 1 << z; change = true; }
			}
		}
		return masque;
	}
	int ordre(int o) const { int n = 1, x = o; while (x != 0) { x = g.comp[o][x][0]; ++n; } return n; }
	bool rotation(int o) const { return g.ori[o].a * g.ori[o].e - g.ori[o].b * g.ori[o].d == 1; }

	Enumerateur(const Grille& g_) : g(g_) {
		set<int> masques;
		for (int x = 0; x < g.no; ++x) for (int y = x; y < g.no; ++y) masques.insert(ferme({ x, y }));
		auto conj = [&](int masque, int q) { int r = 0; for (int x = 0; x < g.no; ++x) if (masque >> x & 1) r |= 1 << g.comp[g.comp[q][x][0]][g.inv[q]][0]; return r; };
		set<int> vus;
		for (int masque : masques) {
			if (vus.count(masque)) continue;
			ClasseP cp;
			for (int q = 0; q < g.no; ++q) { int cm = conj(masque, q); vus.insert(cm); if (cm == masque) cp.normalisateur.push_back(q); }
			for (int x = 0; x < g.no; ++x) if (masque >> x & 1) cp.elems.push_back(x);
			int meilleur = -1;
			for (int x : cp.elems) if (ferme({ x }) == masque) { cp.gens = { x }; break; }
			if (cp.gens.empty()) for (int x : cp.elems) for (int y : cp.elems) if (ferme({ x, y }) == masque) {
				int score = (rotation(x) ? 100 : 0) + ordre(x);
				if (score > meilleur) { meilleur = score; cp.gens = { x, y }; }
			}
			classes.push_back(cp);
		}
	}

	bool cloture(const ClasseP& cp, const Lat& L, const vector<int>& lg, vector<int>& t) const {
		int a = L.a;
		t.assign(g.no, -1); t[0] = 0;
		vector<int> file { 0 };
		for (size_t h = 0; h < file.size(); ++h) {
			int x = file[h], xi = t[x] % a, xj = t[x] / a;
			for (size_t k = 0; k < cp.gens.size(); ++k) {
				int y = cp.gens[k]; const cel& cm = g.comp[y][x]; pii l = g.lin(y, xi, xj);
				int z = cm[0], v = L.idx(cm[1] + l.first + lg[k] % a, cm[2] + l.second + lg[k] / a);
				if (t[z] < 0) { t[z] = v; file.push_back(z); } else if (t[z] != v) return false;
			}
		}
		return true;
	}

	const vector<Groupe>& groupes(int pc, int m) {
		auto it = cache.find({ pc, m });
		if (it != cache.end()) return it->second;
		vector<Groupe>& res = cache[{ pc, m }];
		const ClasseP& cp = classes[pc];
		for (int a = 1; a <= m; ++a) if (m % a == 0) for (int b = 0; b < a; ++b) {
			Lat L { a, b, m / a };
			bool inv = true;
			for (int o : cp.gens) { pii u = g.lin(o, L.a, 0), v = g.lin(o, L.b, L.d); if (!L.contient(u.first, u.second) || !L.contient(v.first, v.second)) { inv = false; break; } }
			if (!inv) continue;
			bool canon = true;
			for (int q : cp.normalisateur) { pii u = g.lin(q, L.a, 0), v = g.lin(q, L.b, L.d); if (hnf(u.first, u.second, v.first, v.second) < L) { canon = false; break; } }
			if (!canon) continue;
			vector<vector<int>> cand(cp.gens.size());
			vector<int> noyau(m); iota(noyau.begin(), noyau.end(), 0);
			for (size_t k = 0; k < cp.gens.size(); ++k) {
				int o = cp.gens[k];
				vector<char> dansW(m, 0); vector<int> W, suivant;
				for (int ti : noyau) {
					int i = ti % a, j = ti / a; pii l = g.lin(o, i, j);
					int w = L.idx(i - l.first, j - l.second);
					if (!dansW[w]) { dansW[w] = 1; W.push_back(w); }
					if (w == 0) suivant.push_back(ti);
				}
				vector<char> marque(m, 0);
				for (int v = 0; v < m; ++v) if (!marque[v]) {
					cand[k].push_back(v);
					for (int w : W) marque[L.idx(v % a + w % a, v / a + w / a)] = 1;
				}
				noyau = suivant;
			}
			vector<int> lg(cp.gens.size()), t;
			if (cp.gens.size() == 1) {
				for (int l1 : cand[0]) { lg[0] = l1; if (cloture(cp, L, lg, t)) res.push_back({ pc, L, t }); }
			} else {
				for (int l1 : cand[0]) for (int l2 : cand[1]) { lg[0] = l1; lg[1] = l2; if (cloture(cp, L, lg, t)) res.push_back({ pc, L, t }); }
			}
		}
		return res;
	}
};

struct DLX {
	vector<int> L, R, U, D, C, S, ligne;
	int ncol, maxLigne, limiteProf;
	long long noeuds, limiteNoeuds;
	bool abandon;
	vector<int> choix;
	void init(int n) {
		ncol = n; maxLigne = 0;
		L.assign(n + 1, 0); R.assign(n + 1, 0); U.assign(n + 1, 0); D.assign(n + 1, 0); C.assign(n + 1, 0); S.assign(n + 1, 0); ligne.assign(n + 1, -1);
		for (int i = 0; i <= n; ++i) { L[i] = i == 0 ? n : i - 1; R[i] = i == n ? 0 : i + 1; U[i] = D[i] = i; C[i] = i; }
	}
	void ajoute(int id, const vector<int>& cols) {
		int premier = -1;
		maxLigne = max(maxLigne, (int)cols.size());
		for (int cidx : cols) {
			int c = cidx + 1, x = L.size();
			L.push_back(x); R.push_back(x); U.push_back(U[c]); D.push_back(c); C.push_back(c); ligne.push_back(id);
			D[U[c]] = x; U[c] = x; ++S[c];
			if (premier < 0) premier = x; else { L[x] = L[premier]; R[x] = premier; R[L[premier]] = x; L[premier] = x; }
		}
	}
	void couvre(int c) {
		R[L[c]] = R[c]; L[R[c]] = L[c];
		for (int i = D[c]; i != c; i = D[i]) for (int j = R[i]; j != i; j = R[j]) { D[U[j]] = D[j]; U[D[j]] = U[j]; --S[C[j]]; }
	}
	void decouvre(int c) {
		for (int i = U[c]; i != c; i = U[i]) for (int j = L[i]; j != i; j = L[j]) { ++S[C[j]]; D[U[j]] = j; U[D[j]] = j; }
		R[L[c]] = c; L[R[c]] = c;
	}
	bool cherche(int prof, int restant) {
		if (R[0] == 0) return true;
		if (prof >= limiteProf || restant > (long long)(limiteProf - prof) * maxLigne) return false;
		if (++noeuds > limiteNoeuds) { abandon = true; return false; }
		int c = R[0];
		for (int j = R[0]; j != 0; j = R[j]) if (S[j] < S[c]) c = j;
		if (S[c] == 0) return false;
		couvre(c);
		for (int r = D[c]; r != c; r = D[r]) {
			int nb = 1;
			for (int j = R[r]; j != r; j = R[j]) { couvre(C[j]); ++nb; }
			choix.push_back(ligne[r]);
			if (cherche(prof + 1, restant - nb)) return true;
			choix.pop_back();
			for (int j = L[r]; j != r; j = L[j]) decouvre(C[j]);
			if (abandon) break;
		}
		decouvre(c);
		return false;
	}
};

struct Forme {
	const Grille& g;
	vector<cel> cellules;
	vector<vector<cel>> orientees;
	vector<int> os, sigma;
	vector<pii> delta;
	vector<vector<cel>> act;
	int n, q, smax;
	bool equilibree;

	static pii normalise(vector<cel>& v) {
		cel anc = v[0];
		for (auto& x : v) if (make_tuple(x[2], x[1], x[0]) < make_tuple(anc[2], anc[1], anc[0])) anc = x;
		for (auto& x : v) { x[1] -= anc[1]; x[2] -= anc[2]; }
		sort(v.begin(), v.end());
		return { anc[1], anc[2] };
	}
	Forme(const Grille& g_, const vector<cel>& cs) : g(g_), cellules(cs) {
		n = cs.size();
		map<vector<cel>, int> vus;
		sigma.resize(g.no); delta.resize(g.no);
		for (int o = 0; o < g.no; ++o) {
			vector<cel> v;
			for (auto& x : cellules) v.push_back(g.image(o, x));
			delta[o] = normalise(v);
			auto it = vus.find(v);
			if (it == vus.end()) { vus[v] = orientees.size(); sigma[o] = orientees.size(); orientees.push_back(v); os.push_back(o); } else sigma[o] = it->second;
		}
		q = orientees.size(); smax = g.no / q;
		act.assign(g.no, vector<cel>(q));
		for (int o = 0; o < g.no; ++o) for (int s = 0; s < q; ++s) {
			const cel& cm = g.comp[o][os[s]]; pii l = g.lin(o, delta[os[s]].first, delta[os[s]].second);
			act[o][s] = { sigma[cm[0]], delta[cm[0]].first + cm[1] - l.first, delta[cm[0]].second + cm[2] - l.second };
			vector<cel> v, w;
			for (auto& x : orientees[s]) v.push_back(g.image(o, x));
			for (auto& x : orientees[act[o][s][0]]) w.push_back({ x[0], x[1] + act[o][s][1], x[2] + act[o][s][2] });
			sort(v.begin(), v.end()); sort(w.begin(), w.end());
			if (v != w) { cerr << "table d'action incohérente\n"; exit(1); }
		}
		vector<int> parRep(g.c, 0);
		for (auto& x : cellules) ++parRep[x[0]];
		vector<int> cl(g.c); iota(cl.begin(), cl.end(), 0);
		for (int o = 0; o < g.no; ++o) for (int r = 0; r < g.c; ++r) { int a = cl[r], b = cl[g.imgRep[o][r][0]]; if (a != b) for (int& x : cl) if (x == b) x = a; }
		vector<int> parClasse(g.c, 0), taille(g.c, 0);
		for (int r = 0; r < g.c; ++r) { parClasse[cl[r]] += parRep[r]; ++taille[cl[r]]; }
		equilibree = true;
		for (int r = 0; r < g.c; ++r) if (taille[r] && parClasse[r] * g.c != n * taille[r]) equilibree = false;
	}
};

struct Temoin { Groupe G; vector<pii> placements; };

struct Solveur {
	const Grille& g; Enumerateur& en;
	DLX dlx;
	vector<int> orbCell, tampon; vector<char> vu;
	int estampille = 0;
	long long groupesTestes = 0;
	Solveur(const Grille& g_, Enumerateur& en_) : g(g_), en(en_) {}

	bool resous(const Forme& F, const Groupe& G, int k, long long limiteNoeuds, Temoin& tem) {
		const ClasseP& cp = en.classes[G.pc]; const Lat& L = G.L;
		int m = L.m(), a = L.a, nc = g.c * m, np = F.q * m;
		orbCell.assign(nc, -1);
		int ni = 0;
		for (int x = 0; x < nc; ++x) if (orbCell[x] < 0) {
			cel q { x / m, (x % m) % a, (x % m) / a };
			for (int o : cp.elems) { cel y = g.image(o, q); orbCell[y[0] * m + L.idx(y[1] + G.lam[o] % a, y[2] + G.lam[o] / a)] = ni; }
			++ni;
		}
		dlx.init(ni);
		vu.assign(np, 0);
		if ((int)tampon.size() < nc) tampon.assign(nc, 0);
		vector<pii> options; vector<int> membres, cols;
		for (int p = 0; p < np; ++p) if (!vu[p]) {
			int s = p / m, ti = (p % m) % a, tj = (p % m) / a;
			membres.clear();
			for (int o : cp.elems) {
				const cel& ac = F.act[o][s]; pii l = g.lin(o, ti, tj);
				int id = ac[0] * m + L.idx(ac[1] + l.first + G.lam[o] % a, ac[2] + l.second + G.lam[o] / a);
				if (!vu[id]) { vu[id] = 1; membres.push_back(id); }
			}
			++estampille;
			bool valide = true;
			for (int id : membres) {
				int s2 = id / m, i2 = (id % m) % a, j2 = (id % m) / a;
				for (auto& x : F.orientees[s2]) {
					int cid = x[0] * m + L.idx(x[1] + i2, x[2] + j2);
					if (tampon[cid] == estampille) { valide = false; break; }
					tampon[cid] = estampille;
				}
				if (!valide) break;
			}
			if (!valide) continue;
			cols.clear();
			for (auto& x : F.orientees[s]) cols.push_back(orbCell[x[0] * m + L.idx(x[1] + ti, x[2] + tj)]);
			sort(cols.begin(), cols.end()); cols.erase(unique(cols.begin(), cols.end()), cols.end());
			dlx.ajoute(options.size(), cols);
			options.push_back({ s, p % m });
		}
		++groupesTestes;
		dlx.limiteProf = k; dlx.limiteNoeuds = limiteNoeuds; dlx.abandon = false; dlx.choix.clear();
		if (!dlx.cherche(0, ni)) return false;
		tem.G = G; tem.placements.clear();
		for (int c : dlx.choix) tem.placements.push_back(options[c]);
		return true;
	}

	// Renvoie k trouvé, ou 0 si aucun pavage à au plus K orbites ; kprouve = plus grand k tel que l'absence de pavage à moins de k orbites est établie.
	int nombreIsoedrique(const Forme& F, int kmin, int K, long long limiteNoeuds, Temoin& tem, int& kprouve, vector<long long>& parK) {
		kprouve = kmin; dlx.noeuds = 0; groupesTestes = 0;
		if (!F.equilibree) { kprouve = 1 << 30; return 0; }
		for (int k = kmin; k <= K; ++k) {
			bool complet = true;
			long long avant = groupesTestes;
			for (size_t pc = 0; pc < en.classes.size(); ++pc) {
				int ordreP = en.classes[pc].elems.size();
				set<int> sommes { 0 };
				for (int i = 0; i < k; ++i) { set<int> nv; for (int x : sommes) for (int s = 1; s <= F.smax; ++s) if (F.smax % s == 0 && ordreP % s == 0) nv.insert(x + ordreP / s); sommes = nv; }
				for (int N : sommes) {
					if ((long long)N * F.n % g.c) continue;
					int m = (long long)N * F.n / g.c;
					for (const Groupe& G : en.groupes(pc, m)) {
						if (resous(F, G, k, limiteNoeuds, tem)) { parK.push_back(groupesTestes - avant); return k; }
						if (dlx.abandon) complet = false;
					}
				}
			}
			parK.push_back(groupesTestes - avant);
			if (!complet) return -1;
			kprouve = k + 1;
		}
		return 0;
	}
};

static bool lisForme(const vector<Grille>& gs, const string& ligne, int& gi, vector<cel>& cs, string& err) {
	istringstream is(ligne); string nom;
	if (!(is >> nom)) return false;
	gi = -1;
	for (auto& g : gs) if (g.nom == nom || (nom.size() == 1 && g.abrev == nom[0])) gi = g.id;
	if (gi < 0) { err = "grille inconnue"; return true; }
	cs.clear(); int x, y; set<cel> vus;
	while (is >> x >> y) { cel q; if (!gs[gi].decompose(x, y, q)) { err = "cellule hors grille"; return true; } if (vus.insert(q).second) cs.push_back(q); }
	if (cs.empty()) err = "forme vide";
	return true;
}

static void ecritAff(ostream& os, int a, int b, int c, int d, int e, int f) { os << '[' << a << ',' << b << ',' << c << ',' << d << ',' << e << ',' << f << ']'; }

static void modeResout(const vector<Grille>& gs, int kmin, int K, long long limiteNoeuds) {
	vector<Enumerateur*> ens(gs.size(), nullptr); vector<Solveur*> sols(gs.size(), nullptr);
	string ligne;
	while (getline(cin, ligne)) {
		int gi; vector<cel> cs; string err;
		if (!lisForme(gs, ligne, gi, cs, err)) continue;
		if (!err.empty()) { cout << "{\"erreur\":\"" << err << "\",\"ligne\":\"" << ligne << "\"}" << endl; continue; }
		const Grille& g = gs[gi];
		if (!ens[gi]) { ens[gi] = new Enumerateur(g); sols[gi] = new Solveur(g, *ens[gi]); }
		Forme F(g, cs); Temoin tem; int kprouve; vector<long long> parK;
		int k = sols[gi]->nombreIsoedrique(F, kmin, K, limiteNoeuds, tem, kprouve, parK);
		cout << "{\"grille\":\"" << g.nom << "\",\"n\":" << F.n << ",\"cellules\":[";
		for (size_t i = 0; i < cs.size(); ++i) { pii p = g.xy(cs[i]); cout << (i ? "," : "") << '[' << p.first << ',' << p.second << ']'; }
		cout << "],\"symetries\":" << F.smax << ",\"statut\":\"" << (k > 0 ? "pave" : k < 0 ? "abandon" : !F.equilibree ? "ne_pave_pas" : "aucun") << "\",\"k\":" << (k > 0 ? k : 0)
			 << ",\"k_min_prouve\":" << (k > 0 ? k : F.equilibree ? kprouve : 0) << ",\"K\":" << K << ",\"noeuds\":" << sols[gi]->dlx.noeuds << ",\"groupes\":[";
		for (size_t i = 0; i < parK.size(); ++i) cout << (i ? "," : "") << parK[i];
		cout << ']';
		if (k > 0) {
			const Lat& L = tem.G.L; const ClasseP& cp = ens[gi]->classes[tem.G.pc];
			cout << ",\"temoin\":{\"reseau\":[[" << L.a * g.v1x << ',' << L.a * g.v1y << "],[" << L.b * g.v1x + L.d * g.v2x << ',' << L.b * g.v1y + L.d * g.v2y << "]],\"groupe\":[";
			bool premier = true;
			for (int o : cp.elems) {
				int i = tem.G.lam[o] % L.a, j = tem.G.lam[o] / L.a; const Aff& T = g.ori[o];
				cout << (premier ? "" : ","); premier = false;
				ecritAff(cout, T.a, T.b, T.c + i * g.v1x + j * g.v2x, T.d, T.e, T.f + i * g.v1y + j * g.v2y);
			}
			cout << "],\"tuiles\":[";
			premier = true;
			for (pii pl : tem.placements) {
				int o = F.os[pl.first], i = pl.second % L.a - F.delta[o].first, j = pl.second / L.a - F.delta[o].second; const Aff& T = g.ori[o];
				cout << (premier ? "" : ","); premier = false;
				ecritAff(cout, T.a, T.b, T.c + i * g.v1x + j * g.v2x, T.d, T.e, T.f + i * g.v1y + j * g.v2y);
			}
			cout << "]}";
		}
		cout << "}" << endl;
	}
}

static void modeGen(const Grille& g, int n, bool sansTrous) {
	auto canon = [&](const vector<cel>& v) {
		vector<cel> meilleur;
		for (int o = 0; o < g.no; ++o) {
			vector<cel> w;
			for (auto& x : v) w.push_back(g.image(o, x));
			Forme::normalise(w);
			if (meilleur.empty() || w < meilleur) meilleur = w;
		}
		return meilleur;
	};
	set<vector<cel>> courant;
	set<int> classesVues;
	for (int r = 0; r < g.c; ++r) courant.insert(canon({ cel { r, 0, 0 } }));
	for (int t = 1; t < n; ++t) {
		set<vector<cel>> suivant;
		for (auto& v : courant) {
			set<cel> dedans(v.begin(), v.end());
			for (auto& x : v) {
				pii p = g.xy(x);
				for (pii d : g.vois[x[0]]) {
					cel y; if (!g.decompose(p.first + d.first, p.second + d.second, y)) { cerr << "voisin hors grille\n"; exit(1); }
					if (dedans.count(y)) continue;
					vector<cel> w = v; w.push_back(y);
					suivant.insert(canon(w));
				}
			}
		}
		courant.swap(suivant);
	}
	for (auto& v : courant) {
		if (sansTrous) {
			set<cel> dedans(v.begin(), v.end()), ext;
			int lo1 = 1 << 30, hi1 = -(1 << 30), lo2 = lo1, hi2 = hi1;
			for (auto& x : v) { lo1 = min(lo1, x[1]); hi1 = max(hi1, x[1]); lo2 = min(lo2, x[2]); hi2 = max(hi2, x[2]); }
			lo1 -= 3; lo2 -= 3; hi1 += 3; hi2 += 3;
			vector<cel> pile { cel { 0, lo1, lo2 } }; ext.insert(pile[0]);
			long long total = (long long)g.c * (hi1 - lo1 + 1) * (hi2 - lo2 + 1);
			while (!pile.empty()) {
				cel x = pile.back(); pile.pop_back(); pii p = g.xy(x);
				for (pii d : g.vois[x[0]]) {
					cel y; g.decompose(p.first + d.first, p.second + d.second, y);
					if (y[1] < lo1 || y[1] > hi1 || y[2] < lo2 || y[2] > hi2 || dedans.count(y) || ext.count(y)) continue;
					ext.insert(y); pile.push_back(y);
				}
			}
			if ((long long)ext.size() + (long long)v.size() != total) continue;
		}
		cout << g.nom;
		for (auto& x : v) { pii p = g.xy(x); cout << ' ' << p.first << ' ' << p.second; }
		cout << '\n';
	}
}

int main(int argc, char** argv) {
	string chemin = string(argv[0]); size_t sl = chemin.rfind('/');
	chemin = (sl == string::npos ? string(".") : chemin.substr(0, sl)) + "/grilles.txt";
	vector<string> args(argv + 1, argv + argc);
	int K = 10, kmin = 1; long long limite = 1LL << 60; bool sansTrous = false;
	vector<string> pos;
	for (size_t i = 0; i < args.size(); ++i) {
		if (args[i] == "-g") chemin = args[++i];
		else if (args[i] == "-K") K = atoi(args[++i].c_str());
		else if (args[i] == "-kmin") kmin = atoi(args[++i].c_str());
		else if (args[i] == "-noeuds") limite = atoll(args[++i].c_str());
		else if (args[i] == "-sans-trous") sansTrous = true;
		else pos.push_back(args[i]);
	}
	vector<Grille> gs = lisGrilles(chemin);
	auto trouve = [&](const string& nom) -> const Grille& { for (auto& g : gs) if (g.nom == nom) return g; cerr << "grille inconnue : " << nom << "\n"; exit(1); };
	if (pos.empty() || pos[0] == "resout") modeResout(gs, kmin, K, limite);
	else if (pos[0] == "gen" && pos.size() == 3) modeGen(trouve(pos[1]), atoi(pos[2].c_str()), sansTrous);
	else if (pos[0] == "groupes" && pos.size() == 3) {
		const Grille& g = trouve(pos[1]); Enumerateur en(g); int mmax = atoi(pos[2].c_str());
		for (size_t pc = 0; pc < en.classes.size(); ++pc) {
			cout << "classe " << pc << " ordre " << en.classes[pc].elems.size() << " gens " << en.classes[pc].gens.size() << " :";
			for (int m = 1; m <= mmax; ++m) cout << ' ' << en.groupes(pc, m).size();
			cout << '\n';
		}
	} else { cerr << "usage : iso [-g grilles.txt] [-K 10] [-kmin 1] [-noeuds N] resout < formes | iso gen GRILLE N [-sans-trous] | iso groupes GRILLE MMAX\n"; return 1; }
	return 0;
}
