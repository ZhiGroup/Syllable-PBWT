#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>

using namespace std;

const unsigned long long MOD = (1ull << 63) - 25;

int M, N, n, L, l, Q; // M = # of haplotypes, N = # of individual sites, n = # of merged sites, L = minimum # of individual sites in a match, l = minimum # of merged sites in a match, Q = # of query haplotypes
// the data structures below use merged sites
vector<vector<unsigned long long>> x, h, z; // x = panel, h = prefix hashes of panel haplotypes, z = query panel
vector<vector<int>> a, d, inv_a; // a = positional prefix arrays, d = divergence arrays, inv_a = inverses of positional prefix arrays
vector<unsigned long long> xp = {1}; // xp[i] = pow(2, i * 64) % MOD

int main(int argc, char** argv) {
	string chr = argv[2];
	string vcf = "/data6/victor/panel" + chr + ".vcf";
	ifstream in(vcf);
	ios_base::sync_with_stdio(0); cin.tie(0);

	// find M, N
	string line;
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	while (getline(in, line)) N++;
	n = (N + 63) / 64;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	M = -9;
	while (getline(ss, line, '\t')) M++;
	M <<= 1;
	L = stoi(argv[1]), assert(L >= 127), l = (L - 63) / 64;

	h.resize(M, vector<unsigned long long>(n + 1)); // pad h with a row of 0s so that h[i][j] = the hash of x[i] over sites [0, j)
	x.resize(M + 1, vector<unsigned long long>(n)); // leave an extra space for the query haplotype
	a.push_back(vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0); // initialize a[0] to ordering in panel
	d.push_back(vector<int>(M));

	N = 0, n = 0; // N = current individual site, n = current merged site
	int rem = -1, np1 = 1; // rem = # of remaining sites to pad onto final merged site, np1 = n + 1
	unsigned __int128 lt, rt; // recurring variables for hash comparison
	int lo, hi, g; // recurring variables for binary search

	while (rem != 0) {
		if (rem == -1) {
			if (!getline(in, line)) { // if finished reading panel, go to else statement to pad up the last site
				rem = (N + 63) / 64 * 64 - N;
				continue;
			}

			if (N % 64 == 0) { // if starting new merged site, copy over last site's prefix hashes and xp
				xp.push_back(xp[n]);
				for (int i = 0; i < M; i++) h[i][np1] = h[i][n];
			}
			// update current merged site with current individual site
			ss = stringstream(line);
			for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
			int i = 0;
			while (getline(ss, line, '\t')) {
				x[i][n] = (x[i][n] << 1) | (line[0] != '0');
				h[i][np1] = (h[i][np1] << 1) | (line[0] != '0');
				if (h[i][np1] >= MOD) h[i][np1] -= MOD;
				i++;
				x[i][n] = (x[i][n] << 1) | (line[2] != '0');
				h[i][np1] = (h[i][np1] << 1) | (line[2] != '0');
				if (h[i][np1] >= MOD) h[i][np1] -= MOD;
				i++;
			}
			xp[np1] <<= 1;
			if (xp[np1] >= MOD) xp[np1] -= MOD;

			if (++N % 64 > 0) continue; // if not at end of merged site, keep reading input
		} else { // pad last site with 0s to reach length 64
			for (int i = 0; i < M; i++) {
				x[i][n] <<= 1;
				h[i][np1] <<= 1;
				if (h[i][np1] >= MOD) h[i][np1] -= MOD;
			}
			xp[np1] <<= 1;
			if (xp[np1] >= MOD) xp[np1] -= MOD;

			if (--rem > 0) continue;
		}

		// let inv_a[k][i] = j such that a[k][j] = i
		inv_a.push_back(vector<int>(M + 1));
		for (int i = 0; i < M; i++) {
			inv_a[n][a[n][i]] = i;
		}

		// order haplotypes in a[n + 1] by value at site n and then by position in a[n]
		a.push_back(a[0]);
		sort(a[np1].begin(), a[np1].end(), [](int i, int j) {
			if (x[i][n] != x[j][n]) return x[i][n] < x[j][n];
			return inv_a[n][i] < inv_a[n][j];
		});

		d.push_back(d[0]);
		d[np1][0] = np1;
		for (int i = 1; i < M; i++) {
			int x_cur = a[np1][i], x_above = a[np1][i - 1], a_prev = inv_a[n][x_cur];

			if (x[x_cur][n] != x[x_above][n]) { // no match
				d[np1][i] = np1;
			} else if (a_prev > 0 && x_above == a[n][a_prev - 1]) { // continued match
				d[np1][i] = d[n][a_prev];
			} else { // different match, binary search for divergence value using hashes
				lo = d[n][a_prev], hi = n;
				rt = h[x_cur][np1] < h[x_above][np1] ? MOD - h[x_above][np1] + h[x_cur][np1] : h[x_cur][np1] - h[x_above][np1];

				while (lo < hi) {
					g = (lo + hi) >> 1;
					lt = h[x_cur][g] < h[x_above][g] ? MOD - h[x_above][g] + h[x_cur][g] : h[x_cur][g] - h[x_above][g];

					if (lt * xp[np1 - g] % MOD == rt) {
						hi = g;
					} else {
						lo = g + 1;
					}
				}

				d[np1][i] = lo;
			}
		}

		n++, np1++;
	}
	in.close();

	cout << "merge " << M << ' ' << N << ' ' << L << endl;
	vector<unsigned long long> hz(n + 1); // prefix hashes of query haplotype
	vector<int> up_end(n), dn_end(n); // up_end[i] = # of matches above the query haplotype that (inclusively) end on site i, dn_end[i] = ... below ...

	clock_t START = clock();
	int matches = 0;
	// load in query haplotypes
	string q_vcf = "/data6/victor/query" + chr + ".vcf";
	ifstream q_in(q_vcf);
	while (getline(q_in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	Q = -9;
	ss = stringstream(line);
	while (getline(ss, line, '\t')) Q++;
	Q <<= 1;
	z.resize(n, vector<unsigned long long>(Q));
	int NN = 0;
	while (getline(q_in, line)) {
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int i = 0, nn = NN / 64;
		while (getline(ss, line, '\t')) {
			z[nn][i] = (z[nn][i] << 1) | (line[0] != '0');
			i++;
			z[nn][i] = (z[nn][i] << 1) | (line[2] != '0');
			i++;
		}
		NN++;
	}
	q_in.close();

	for (int q = 0; q < Q; q++) {
	ofstream out("/home/vwang1/PBWT_compression/query064-" + to_string(q) + ".out");
	z[n - 1][q] <<= n * 64 - N; // pad final merged site with necessary 0s
	for (int i = 0, ip1 = 1; i < n; i++, ip1++) {
		x[M][i] = z[i][q]; // copy query haplotype into the panel
		hz[ip1] = ((unsigned __int128) hz[i] * xp[1] + x[M][i]) % MOD; // build prefix hashes of query haplotype
	}

	memset(&up_end[0], 0, n * sizeof(int));
	memset(&dn_end[0], 0, n * sizeof(int));
	for (int i = 0, ip1 = 1, t = 0, up_ct = 0, dn_ct = 0; i < n; i++, ip1++) { // i = current site, ip1 = i + 1, t = position of query haplotype in positional prefix array, up_ct = # of ongoing matches above query haplotype, dn_ct = ... below ...
		// binary search for the query haplotype's position in a[i + 1], according to the sorting by value at site i and then by position in a[i]
		inv_a[i][M] = t;
		t = lower_bound(a[ip1].begin(), a[ip1].end(), M, [i](int i1, int i2) { if (x[i1][i] != x[i2][i]) return x[i1][i] < x[i2][i]; return inv_a[i][i1] < inv_a[i][i2]; }) - a[ip1].begin();
		if (ip1 < l) continue; // no possible matches

		int s_idx = i - l, min_start = max(0, s_idx * 64 + 1);
		// check for match immediately above query haplotype
		bool touch = false;
		if (up_ct == 0 && t > 0) {
			int above = a[ip1][t - 1];
			if (dn_ct > 0) touch = d[ip1][t] <= ip1 - l; // can use divergence value if there are ongoing matches below the query haplotype
			else if (x[M][i] == x[above][i]) { // otherwise use hashes to check if match extends l sites back
				lt = hz[ip1 - l] < h[above][ip1 - l] ? MOD - h[above][ip1 - l] + hz[ip1 - l] : hz[ip1 - l] - h[above][ip1 - l];
				rt = hz[ip1] < h[above][ip1] ? MOD - h[above][ip1] + hz[ip1] : hz[ip1] - h[above][ip1];
				touch = lt * xp[l] % MOD == rt;
			}
		}

		if (up_ct > 0 || touch) {
			int p = t - up_ct; // pointer for finding matches above the query haplotype
			while (d[ip1][p] <= ip1 - l || touch) { // while there are more matches
				touch = false;
				int above = a[ip1][--p];

				// briefly linear search for end of match
				lo = ip1;
				int pivot = min(n, lo + 10);
				while (lo < pivot && x[M][lo] == x[above][lo]) lo++;
				if (lo == pivot) { // switch to binary search if end was not found with linear search
					hi = n;
					lt = hz[i] < h[above][i] ? MOD - h[above][i] + hz[i] : hz[i] - h[above][i];

					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						rt = hz[g] < h[above][g] ? MOD - h[above][g] + hz[g] : hz[g] - h[above][g];

						if (lt * xp[g - i] % MOD == rt) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				}
				up_end[lo - 1]++; // lo = exclusive end of match

				// binary search for individual site resolution ends of match
				int start = 0, end = N, e_idx = lo;
				if (e_idx < n) {
					lo = 1, hi = 64;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][e_idx] >> g) == (x[above][e_idx] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 64 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 1, hi = 64;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][s_idx] & (-1ull >> g)) == (x[above][s_idx] & (-1ull >> g))) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					start = s_idx * 64 + lo;
				}

				if (end - start >= L) { // report match
					out << start << ' ' << end << ' ' << above << '\n';
					matches++;
				}
			}
			up_ct = t - p - up_end[i]; // update # of ongoing matches
		}

		// analogous process for downward

		touch = false;
		if (dn_ct == 0 && t < M) {
			int below = a[ip1][t];
			if (up_ct > 0) touch = d[ip1][t] <= ip1 - l;
			else if (x[below][i] == x[M][i]) {
				lt = hz[ip1 - l] < h[below][ip1 - l] ? MOD - h[below][ip1 - l] + hz[ip1 - l] : hz[ip1 - l] - h[below][ip1 - l];
				rt = hz[ip1] < h[below][ip1] ? MOD - h[below][ip1] + hz[ip1] : hz[ip1] - h[below][ip1];
				touch = lt * xp[l] % MOD == rt;
			}
		}

		if (dn_ct > 0 || touch) {
			int p = t + dn_ct;
			while ((p < M && d[ip1][p] <= ip1 - l) || touch) {
				touch = false;
				int below = a[ip1][p++];

				lo = ip1;
				int pivot = min(n, lo + 10);
				while (lo < pivot && x[M][lo] == x[below][lo]) lo++;
				if (lo == pivot) {
					hi = n;
					lt = hz[i] < h[below][i] ? MOD - h[below][i] + hz[i] : hz[i] - h[below][i];

					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						rt = hz[g] < h[below][g] ? MOD - h[below][g] + hz[g] : hz[g] - h[below][g];

						if (lt * xp[g - i] % MOD == rt) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				} 
				dn_end[lo - 1]++;

				int start = 0, end = N, e_idx = lo;
				if (e_idx < n) {
					lo = 1, hi = 64;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][e_idx] >> g) == (x[below][e_idx] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 64 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 1, hi = 64;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][s_idx] & (-1ull >> g)) == (x[below][s_idx] & (-1ull >> g))) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					start = s_idx * 64 + lo;
				}

				if (end - start >= L) {
					out << start << ' ' << end << ' ' << below << '\n';
					matches++;
				}
			}
			dn_ct = p - t - dn_end[i];
		}
	}
	}

	cout << clock() - START << " time\n";
	cout << matches << " matches\n\n";

	return 0;
}

