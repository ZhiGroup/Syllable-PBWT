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

int M, N, n, L, l;
vector<vector<unsigned long long>> x, h;
vector<vector<int>> a, d, inv_a;
vector<unsigned long long> xp = {1};
vector<vector<bool>> z;

int main(int argc, char** argv) {
	string vcf = "/data6/ukbiobank/genotype_data/autosomal_chroms_june2018/21981/ukb_hap_GP_removed/ukb_hap_chr21_v2.vcf";
	ifstream in(vcf);
	ios_base::sync_with_stdio(0); cin.tie(0);

	string line;
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	M = -9;
	while (getline(ss, line, '\t')) M++;
	M <<= 1;
	int Q = 10;
	assert(Q % 2 == 0);
	M -= Q;
	N = 0, n = 0;
	L = stoi(argv[1]), assert(L >= 127), l = (L - 63) / 64;

	h.push_back(vector<unsigned long long>(M));
	a.push_back(vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0);
	d.push_back(vector<int>(M));

	int rem = -1, np1 = 1;
	unsigned __int128 lt, rt;
	int lo, hi, g;

	while (rem != 0) {
		if (rem == -1) {
			if (!getline(in, line)) {
				rem = (N + 63) / 64 * 64 - N;
				continue;
			}

			z.push_back(vector<bool>(Q));
			if (N % 64 == 0) {
				x.push_back(vector<unsigned long long>(M + 1));
				h.push_back(h[n]);
				xp.push_back(xp[n]);
			}
			ss = stringstream(line);
			for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
			int i = 0;
			while (i < Q) {
				getline(ss, line, '\t');
				z[N][i++] = line[0] != '0';
				z[N][i++] = line[2] != '0';
			}
			i = 0;
			while (i < M) {
				getline(ss, line, '\t');
				x[n][i] = (x[n][i] << 1) | (line[0] != '0');
				h[np1][i] = (h[np1][i] << 1) | (line[0] != '0');
				if (h[np1][i] >= MOD) h[np1][i] -= MOD;
				i++;
				x[n][i] = (x[n][i] << 1) | (line[2] != '0');
				h[np1][i] = (h[np1][i] << 1) | (line[2] != '0');
				if (h[np1][i] >= MOD) h[np1][i] -= MOD;
				i++;
			}
			xp[np1] <<= 1;
			if (xp[np1] >= MOD) xp[np1] -= MOD;

			if (++N % 64) continue;
		} else {
			for (int i = 0; i < M; i++) {
				x[n][i] <<= 1;
				h[np1][i] <<= 1;
				if (h[np1][i] >= MOD) h[np1][i] -= MOD;
			}
			xp[np1] <<= 1;
			if (xp[np1] >= MOD) xp[np1] -= MOD;

			if (--rem) continue;
		}

		inv_a.push_back(vector<int>(M + 1));
		for (int i = 0; i < M; i++) {
			inv_a[n][a[n][i]] = i;
		}

		a.push_back(a[0]);
		sort(a[np1].begin(), a[np1].end(), [](int i, int j) {
			if (x[n][i] != x[n][j]) return x[n][i] < x[n][j];
			return inv_a[n][i] < inv_a[n][j];
		});

		d.push_back(d[0]);
		d[np1][0] = np1;
		for (int i = 1; i < M; i++) {
			int x_cur = a[np1][i], x_above = a[np1][i - 1], a_prev = inv_a[n][x_cur];

			if (x[n][x_cur] != x[n][x_above]) {
				d[np1][i] = np1;
			} else if (a_prev > 0 && x_above == a[n][a_prev - 1]) {
				d[np1][i] = d[n][a_prev];
			} else {
				lo = d[n][a_prev], hi = n;
				rt = h[np1][x_cur] < h[np1][x_above] ? MOD - h[np1][x_above] + h[np1][x_cur] : h[np1][x_cur] - h[np1][x_above];

				while (lo < hi) {
					g = (lo + hi) >> 1;
					lt = h[g][x_cur] < h[g][x_above] ? MOD - h[g][x_above] + h[g][x_cur] : h[g][x_cur] - h[g][x_above];

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

	// ofstream dout("/home/vwang1/outm.txt");
	cout << "merge " << M << ' ' << N << ' ' << L << endl;
	vector<unsigned long long> hz(n + 1);
	vector<int> up_end(n), dn_end(n);

	for (int q = 0; q < Q; q++) {
	clock_t START = clock();
	ofstream out("/home/vwang1/samq" + to_string(q) + ".txt");
	for (int i = 0, ip1 = 1, idx = 0; i < n; i++, ip1++) {
		x[i][M] = 0;
		hz[ip1] = hz[i];
		for (int b = 0; b < 64; b++) {
			if (idx < N) {
				x[i][M] = (x[i][M] << 1) | z[idx][q];
				hz[ip1] = (hz[ip1] << 1) | z[idx][q];
				idx++;
			} else {
				x[i][M] <<= 1;
				hz[ip1] <<= 1;
			}
			if (hz[ip1] >= MOD) hz[ip1] -= MOD;
		}
	}

	memset(&up_end[0], 0, n * sizeof(int));
	memset(&dn_end[0], 0, n * sizeof(int));
	int matches = 0;
	for (int i = 0, ip1 = 1, t = 0, up_ct = 0, dn_ct = 0; i < n; i++, ip1++) {
		inv_a[i][M] = t;
		t = lower_bound(a[ip1].begin(), a[ip1].end(), M, [i](int i1, int i2) { unsigned long long x1 = x[i][i1], x2 = x[i][i2]; if (x1 != x2) return x1 < x2; return inv_a[i][i1] < inv_a[i][i2]; }) - a[ip1].begin();
		if (ip1 < l) continue;

		int s_idx = i - l, min_start = max(0, s_idx * 64 + 1);
		bool touch = false;
		if (up_ct == 0 && t > 0) {
			int above = a[ip1][t - 1];
			lt = hz[ip1 - l] < h[ip1 - l][above] ? MOD - h[ip1 - l][above] + hz[ip1 - l] : hz[ip1 - l] - h[ip1 - l][above];
			rt = hz[ip1] < h[ip1][above] ? MOD - h[ip1][above] + hz[ip1] : hz[ip1] - h[ip1][above];
			touch = lt * xp[l] % MOD == rt;
		}

		if (up_ct > 0 || touch) {
			int p = t - up_ct;
			while (d[ip1][p] <= ip1 - l || touch) {
				touch = false;
				int above = a[ip1][--p];

				lo = ip1;
				int stop = min(n, lo + 10);
				while (lo < stop && x[lo][M] == x[lo][above]) lo++;
				if (lo == stop) {
					hi = n;
					lt = hz[i] < h[i][above] ? MOD - h[i][above] + hz[i] : hz[i] - h[i][above];

					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						rt = hz[g] < h[g][above] ? MOD - h[g][above] + hz[g] : hz[g] - h[g][above];

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
					lo = 1, hi = 64;
					while (lo < hi) {
						g = (lo + hi) >> 1;
						if ((x[e_idx][M] >> g) == (x[e_idx][above] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 64 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 0, hi = 63;
					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						if ((x[s_idx][M] & ((1ull << g) - 1)) == (x[s_idx][above] & ((1ull << g) - 1))) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
					start = (s_idx + 1) * 64 - lo;
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
			lt = hz[ip1 - l] < h[ip1 - l][below] ? MOD - h[ip1 - l][below] + hz[ip1 - l] : hz[ip1 - l] - h[ip1 - l][below];
			rt = hz[ip1] < h[ip1][below] ? MOD - h[ip1][below] + hz[ip1] : hz[ip1] - h[ip1][below];
			touch = lt * xp[l] % MOD == rt;
		}

		if (dn_ct > 0 || touch) {
			int p = t + dn_ct;
			while ((p < M && d[ip1][p] <= ip1 - l) || touch) {
				touch = false;
				int below = a[ip1][p++];

				lo = ip1;
				int stop = min(n, lo + 10);
				while (lo < stop && x[lo][M] == x[lo][below]) lo++;
				if (lo == stop) {
					hi = n;
					lt = hz[i] < h[i][below] ? MOD - h[i][below] + hz[i] : hz[i] - h[i][below];

					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						rt = hz[g] < h[g][below] ? MOD - h[g][below] + hz[g] : hz[g] - h[g][below];

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
						if ((x[e_idx][M] >> g) == (x[e_idx][below] >> g)) {
							hi = g;
						} else {
							lo = g + 1;
						}
					}
					end = (e_idx + 1) * 64 - lo;
				}
				if (end - min_start < L) continue;
				if (s_idx > -1) {
					lo = 0, hi = 63;
					while (lo < hi) {
						g = (lo + hi + 1) >> 1;
						if ((x[s_idx][M] & ((1ull << g) - 1)) == (x[s_idx][below] & ((1ull << g) - 1))) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
					start = (s_idx + 1) * 64 - lo;
				}

				if (end - start >= L) {
					out << start << ' ' << end << ' ' << below << '\n';
					matches++;
				}
			}
			dn_ct = p - t - dn_end[i];
		}
	}

	cout << clock() - START << " time\n";
	cout << matches << " matches\n\n";

	}

	return 0;
}

