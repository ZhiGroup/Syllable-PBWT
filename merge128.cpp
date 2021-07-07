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

int M, N, n, L, l, Q;
vector<vector<unsigned __int128>> x, z;
vector<vector<unsigned long long>> h;
vector<vector<int>> a, d, inv_a;
vector<unsigned long long> xp = {1};

int main(int argc, char** argv) {
	string chr = argv[2];
	string vcf = "/data6/victor/panel" + chr + ".vcf";
	ifstream in(vcf);
	ios_base::sync_with_stdio(0); cin.tie(0);

	string line;
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	while (getline(in, line)) N++;
	n = (N + 127) / 128;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	M = -9;
	while (getline(ss, line, '\t')) M++;
	M <<= 1;
	L = stoi(argv[1]), assert(L >= 255), l = (L - 127) / 128;

	h.resize(M, vector<unsigned long long>(n + 1));
	x.resize(M + 1, vector<unsigned __int128>(n));
	a.push_back(vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0);
	d.push_back(vector<int>(M));

	N = 0, n = 0;
	int rem = -1, np1 = 1;
	unsigned __int128 lt, rt;
	int lo, hi, g;

	while (rem != 0) {
		if (rem == -1) {
			if (!getline(in, line)) {
				rem = (N + 127) / 128 * 128 - N;
				continue;
			}

			if (N % 128 == 0) {
				xp.push_back(xp[n]);
				for (int i = 0; i < M; i++) h[i][np1] = h[i][n];
			}
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

			if (++N % 128 > 0) continue;
		} else {
			for (int i = 0; i < M; i++) {
				x[i][n] <<= 1;
				h[i][np1] <<= 1;
				if (h[i][np1] >= MOD) h[i][np1] -= MOD;
			}
			xp[np1] <<= 1;
			if (xp[np1] >= MOD) xp[np1] -= MOD;

			if (--rem > 0) continue;
		}

		inv_a.push_back(vector<int>(M + 1));
		for (int i = 0; i < M; i++) {
			inv_a[n][a[n][i]] = i;
		}

		a.push_back(a[0]);
		sort(a[np1].begin(), a[np1].end(), [](int i, int j) {
			if (x[i][n] != x[j][n]) return x[i][n] < x[j][n];
			return inv_a[n][i] < inv_a[n][j];
		});

		d.push_back(d[0]);
		d[np1][0] = np1;
		for (int i = 1; i < M; i++) {
			int x_cur = a[np1][i], x_above = a[np1][i - 1], a_prev = inv_a[n][x_cur];

			if (x[x_cur][n] != x[x_above][n]) {
				d[np1][i] = np1;
			} else if (a_prev > 0 && x_above == a[n][a_prev - 1]) {
				d[np1][i] = d[n][a_prev];
			} else {
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

	cout << "merge128 " << M << ' ' << N << ' ' << L << endl;
	vector<unsigned long long> hz(n + 1);
	vector<int> up_end(n), dn_end(n);
	unsigned __int128 ones = -1;

	clock_t START = clock();
	int matches = 0;
	string q_vcf = "/data6/victor/query" + chr + ".vcf";
	ifstream q_in(q_vcf);
	while (getline(q_in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	Q = -9;
	ss = stringstream(line);
	while (getline(ss, line, '\t')) Q++;
	Q <<= 1;
	z.resize(n, vector<unsigned __int128>(Q));
	int NN = 0;
	while (getline(q_in, line)) {
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int i = 0, nn = NN / 128;
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
	ofstream out("/home/vwang1/PBWT_compression/query0128-" + to_string(q) + ".out");
	z[n - 1][q] <<= n * 128 - N;
	for (int i = 0, ip1 = 1; i < n; i++, ip1++) {
		x[M][i] = z[i][q];
		hz[ip1] = ((unsigned __int128) hz[i] * xp[1] + x[M][i] % MOD) % MOD;
	}

	memset(&up_end[0], 0, n * sizeof(int));
	memset(&dn_end[0], 0, n * sizeof(int));
	for (int i = 0, ip1 = 1, t = 0, up_ct = 0, dn_ct = 0; i < n; i++, ip1++) {
		inv_a[i][M] = t;
		t = lower_bound(a[ip1].begin(), a[ip1].end(), M, [i](int i1, int i2) { if (x[i1][i] != x[i2][i]) return x[i1][i] < x[i2][i]; return inv_a[i][i1] < inv_a[i][i2]; }) - a[ip1].begin();
		if (ip1 < l) continue;

		int s_idx = i - l, min_start = max(0, s_idx * 128 + 1);
		bool touch = false;
		if (up_ct == 0 && t > 0) {
			int above = a[ip1][t - 1];
			if (dn_ct > 0) touch = d[ip1][t] <= ip1 - l;
			else if (x[M][i] == x[above][i]) {
				lt = hz[ip1 - l] < h[above][ip1 - l] ? MOD - h[above][ip1 - l] + hz[ip1 - l] : hz[ip1 - l] - h[above][ip1 - l];
				rt = hz[ip1] < h[above][ip1] ? MOD - h[above][ip1] + hz[ip1] : hz[ip1] - h[above][ip1];
				touch = lt * xp[l] % MOD == rt;
			}
		}

		if (up_ct > 0 || touch) {
			int p = t - up_ct;
			while (d[ip1][p] <= ip1 - l || touch) {
				touch = false;
				int above = a[ip1][--p];

				lo = ip1;
				int pivot = min(n, lo + 10);
				while (lo < pivot && x[M][lo] == x[above][lo]) lo++;
				if (lo == pivot) {
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
				up_end[lo - 1]++;

				int start = 0, end = N, e_idx = lo;
				if (e_idx < n) {
					lo = 1, hi = 128;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][e_idx] >> g) == (x[above][e_idx] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 128 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 1, hi = 128;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][s_idx] & (ones >> g)) == (x[above][s_idx] & (ones >> g))) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					start = s_idx * 128 + lo;
				}

				if (end - start >= L) {
					out << start << ' ' << end << ' ' << above << '\n';
					matches++;
				}
			}
			up_ct = t - p - up_end[i];
		}

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
					lo = 1, hi = 128;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][e_idx] >> g) == (x[below][e_idx] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 128 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 1, hi = 128;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[M][s_idx] & (ones >> g)) == (x[below][s_idx] & (ones >> g))) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					start = s_idx * 128 + lo;
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

