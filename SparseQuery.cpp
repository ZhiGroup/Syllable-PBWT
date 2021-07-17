#include "SparseQuery.h"
#include "SparseTable.h"

template<class T>
int SparseQuery<T>::precompute(const char* input_file) {
	ifstream in(input_file);
	if (in.fail()) return 1;

	// get M = # haplotypes, N = # sites, n = # sparse sites, IDs
	string line;
	while (getline(in, line)) {
		if (line.size() < 2u) return 2;
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	M = -9;
	while (getline(ss, line, '\t')) M++;
	if (M < 1) return 2;
	M <<= 1;
	IDs.resize(M);
	while (getline(in, line)) N++;
	if (N < 1) return 2;
	physLocs.resize(N);
	n = (N + B - 1) / B;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < M; i += 2) {
		getline(ss, IDs[i], '\t');
		IDs[i + 1] = IDs[i] + "-1";
		IDs[i] += "-0";
	}

	// resize data structures for sparse panel
	x.resize(M + 1, vector<T>(n));
	inv_a.resize(n, vector<int>(M + 1));
	a.resize(n + 1, vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0);
	d.resize(n + 1, vector<int>(M));

	{ SparseTable st(M);

	// K = current site, k = current sparse site, kp1 = k+1
	for (int K = 0, k = 0, kp1 = 1; K < N; K++) {
		getline(in, line);
		ss = stringstream(line);
		for (int i = 0; i < 9; i++) { // get site's physical location
			getline(ss, line, '\t');
			if (i == 1) {
				try { physLocs[K] = stoi(line); }
				catch (exception& e) { return 2; }
			}
		}
		// update sparse panel with current site
		int index = 0;
		while (getline(ss, line, '\t')) {
			if (index == M || line.size() < 3u) return 2;
			x[index][k] = (x[index][k] << 1) | (line[0] != '0'), index++;
			x[index][k] = (x[index][k] << 1) | (line[2] != '0'), index++;
		}
		if (index != M) return 2;

		if (K == N - 1) { // if last site, pad last sparse site with 0s
			int pad = n * B - N;
			for (int i = 0; i < M; i++) x[i][k] <<= pad;
			K += pad;
		}
		if ((K + 1) % B > 0) continue; // skip if not end of sparse site

		// min query length = max distance across any 2 adjacent sparse sites
		minPhysL = max(minPhysL, physLocs[min(K, N - 1)] - physLocs[(k - 1) * B]);

		// inv_a[k][i] = j where a[k][j] = i
		for (int i = 0; i < M; i++) {
			inv_a[k][a[k][i]] = i;
		}

		// sort haplotypes in a[k+1] by (value at site k, position in a[k])
		iota(a[kp1].begin(), a[kp1].end(), 0);
		sort(a[kp1].begin(), a[kp1].end(), [&, k](int i, int j) {
			if (x[i][k] != x[j][k]) return x[i][k] < x[j][k];
			return inv_a[k][i] < inv_a[k][j];
		});

		// build sparse table and divergence array
		st.build(d[k]);
		d[kp1][0] = kp1;
		for (int i = 1; i < M; i++) {
			int x_cur = a[kp1][i], x_above = a[kp1][i - 1];
			if (x[x_cur][k] != x[x_above][k]) {
				d[kp1][i] = kp1;
			} else {
				d[kp1][i] = st.query(inv_a[k][x_above] + 1, inv_a[k][x_cur]);
			}
		}
		
		k++, kp1++;
	}} // deconstruct sparse table

	xp.resize(n + 1), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1);
	// xp[i] = pow(3, i * B) % MOD
	xp[0] = xp[1] = 1;
	for (int b = 0; b < B; b++) {
		xp[1] = (unsigned __int128) xp[1] * 3 % MOD;
	}
	for (int k = 0; k < n; k++) {
		xp[k+1] = (unsigned __int128) xp[k] * xp[1] % MOD;
	}

	// build prefix hashes; h[i][k+1] = hash of x[i][0, k]
	h.resize(M, vector<unsigned long long>(n + 1));
	for (int i = 0; i < M; i++) {
		for (int k = 0, kp1 = 1; k < n; k++, kp1++) {
			h[i][kp1] = (unsigned __int128) h[i][k] * xp[1] % MOD;
			unsigned __int128 xp3b = 1;
			T x_copy = x[i][k];
			for (int b = 0; b < B; b++) {
				h[i][kp1] = (h[i][kp1] + ((x_copy & 1) + 1) * xp3b) % MOD;
				xp3b = xp3b * 3 % MOD;
				x_copy >>= 1;
			}
		}
	}

	return 0;
}

template<class T>
int SparseQuery<T>::save(const char* save_file) {
	ofstream out(save_file, ios::binary);
	if (out.fail()) return 1;

	out.write((char*) &M, sizeof(M));
	out.write((char*) &N, sizeof(N));
	out.write((char*) &B, sizeof(B));
	out.write((char*) &minPhysL, sizeof(minPhysL));
	for (int loc : physLocs) out.write((char*) &loc, sizeof(loc));
	for (string ID : IDs) {
		size_t sz = ID.size();
		out.write((char*) &sz, sizeof(sz));
		out.write(ID.c_str(), sz);
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < n; j++) {
			out.write((char*) &x[i][j], sizeof(T));
		}
	}
	for (int i = 0; i < M; i++) {
		for (int j = 1; j <= n; j++) {
			out.write((char*) &h[i][j], sizeof(unsigned long long));
		}
	}
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			out.write((char*) &a[i][j], sizeof(int));
		}
	}
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			out.write((char*) &d[i][j], sizeof(int));
		}
	}

	return 0;
}

template<class T>
int SparseQuery<T>::load(const char* load_file) {
	ifstream in(load_file, ios::binary);
	if (in.fail()) return 1;

	if (!in.read((char*) &M, sizeof(M))) return 2;
	if (!in.read((char*) &N, sizeof(N))) return 2;
	if (M < 1 || N < 1) return 2;
	if (!in.read((char*) &B, sizeof(B))) return 2;
	if (B != sizeof(T) * 8) return 3;
	if (!in.read((char*) &minPhysL, sizeof(minPhysL))) return 2;
	physLocs.resize(N);
	for (int i = 0; i < N; i++) if (!in.read((char*) &physLocs[i], sizeof(physLocs[i]))) return 2;
	IDs.resize(M);
	for (int i = 0; i < M; i++) {
		size_t sz;
		if (!in.read((char*) &sz, sizeof(sz)) || sz < 1) return 2;
		IDs[i].resize(sz);
		if (!in.read(&IDs[i][0], sz)) return 2;
	}
	n = (N + B - 1) / B, minSiteL = B * 2 - 1;
	xp.resize(n + 1), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1);
	xp[0] = xp[1] = 1;
	for (int b = 0; b < B; b++) {
		xp[1] = (unsigned __int128) xp[1] * 3 % MOD;
	}
	for (int i = 1; i < n; i++) xp[i + 1] = (unsigned __int128) xp[i] * xp[1] % MOD;
	x.resize(M + 1, vector<T>(n));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < n; j++) {
			if (!in.read((char*) &x[i][j], sizeof(T))) return 2;
		}
	}
	h.resize(M, vector<unsigned long long>(n + 1));
	for (int i = 0; i < M; i++) {
		for (int j = 1; j <= n; j++) {
			if (!in.read((char*) &h[i][j], sizeof(unsigned long long))) return 2;
		}
	}
	a.resize(n + 1, vector<int>(M));
	inv_a.resize(n, vector<int>(M + 1));
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			if (!in.read((char*) &a[i][j], sizeof(int)) || !(-1 < a[i][j] && a[i][j] < M)) return 2;
		}
		if (i == n) break;
		for (int j = 0; j < M; j++) {
			inv_a[i][a[i][j]] = j;
		}
	}
	d.resize(n + 1, vector<int>(M));
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			if (!in.read((char*) &d[i][j], sizeof(int))) return 2;
		}
	}

	return 0;
}

template<class T>
void SparseQuery<T>::show_attributes() {
	cout << "\nData attributes:\n" <<
		"\tM = " << M << " haplotypes, N = " << N << " sites, B = " << B << '\n';
        if (minGeneL == -1) {
		cout << "\tYour minimum allowable query lengths are " << minSiteL << " sites or " << minPhysL << " bps.\n" <<
                        "\tProvide a genetic map to query in units of cM." << endl;
	} else {
                cout << "\tYour minimum allowable query lengths are " << minGeneL << " cM, " << minSiteL << " sites, or " << minPhysL << " bps." << endl;
        }
}

template<class T>
int SparseQuery<T>::set_gen_map(const char* map_file) {
	ifstream in(map_file);
	if (in.fail()) return 1;

	if ((int) geneLocs.size() != N) geneLocs.resize(N);
	minGeneL = -1;
	string line;
	for (int i = 0; i < N; i++) {
		if (!getline(in, line)) return 2;
		stringstream ss(line);
		while (getline(ss, line, '\t')) {}
		try { geneLocs[i] = stod(line); }
		catch (out_of_range &e) {}
		catch (exception &e) { return 3; }
	}
	// min query length = max distance across any 2 adjacent sparse sites
	for (int i = 1; i < n; i++) {
		minGeneL = max(minGeneL, geneLocs[min((i + 1) * B, N) - 1] - geneLocs[(i - 1) * B]);
	}

	return 0;
}

template<class T>
template<class U>
int SparseQuery<T>::query(const char* query_file, const char* output_dir, const U L, const vector<U> &locs) { // variably distributed site locations
	if (sizeof(U) == sizeof(double)) {
		if (minGeneL == -1) return 3;
		if (L < minGeneL) return 4;
	} else {
		if (L < minPhysL) return 4;
	}
	if (access(output_dir, W_OK) != 0) return 2;
	ifstream in(query_file);
	if (in.fail()) return 1;
 
	string line;
	while (getline(in, line)) {
		if (line.size() < 2u) return 5;
		if (line[0] != '#' || line[1] != '#') break;
	}
	int Q = -9;
	stringstream ss(line);
	while (getline(ss, line, '\t')) Q++;
	if (Q < 1) return 5;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	vector<string> qIDs(Q);
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < Q; i++) getline(ss, qIDs[i], '\t');
	Q <<= 1;
	vector<vector<T>> z(n, vector<T>(Q));

	// read query panel; "site" will now refer to "sparse site"
	for (int K = 0; K < N; K++) {
		if (!getline(in, line)) return 5;
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int k = K / B, i = 0;
		while (getline(ss, line, '\t')) {
			if (i == Q || line.size() < 3u) return 5;
			z[k][i] = (z[k][i] << 1) | (line[0] != '0'), i++;
			z[k][i] = (z[k][i] << 1) | (line[2] != '0'), i++;
		}
	}
	in.close();

	unsigned __int128 lt, rt;

	for (int q = 0; q < Q; q++) {
		ofstream out(string(output_dir) + "/" + qIDs[q >> 1] + "-" + to_string(q & 1) + ".txt");

		z[n - 1][q] <<= n * B - N; // pad last site with 0s
		for (int k = 0, kp1 = 1; k < n; k++, kp1++) {
			x[M][k] = z[k][q]; // copy query into panel
			hz[kp1] = ((unsigned __int128) hz[k] * xp[1]) % MOD;
			unsigned __int128 xp3b = 1;
			T x_copy = x[M][k];
			for (int b = 0; b < B; b++) {
				hz[kp1] = (hz[kp1] + ((x_copy & 1) + 1) * xp3b) % MOD;
				xp3b = xp3b * 3 % MOD;
				x_copy >>= 1;
			}
		}
		// up_end[k] = # matches (above the query) that end at site k, dn_end[k] = ... below ... (above/below refer to positions in a[k])
		memset(&up_end[1], 0, n * sizeof(int));
		memset(&dn_end[1], 0, n * sizeof(int));

		/* k = current site, kp1 = k+1, t = virtual location of z in a[k+1], req_idx = rightmost site that must fully match,
		   up_beg = non-inclusive start site of above match, dn_beg = ... below ..., up_on = # ongoing above matches, dn_on = ... below ... */
		for (int k = 0, kp1 = 1, t = 0, req_idx = 0, up_beg = -1, dn_beg = -1, up_on = 0, dn_on = 0; kp1 < n; k++, kp1++) {
			// binary search for query's position in a[k+1], according to sorting by (value at site k, position in a[k])
			inv_a[k][M] = t;
			t = lower_bound(a[kp1].begin(), a[kp1].end(), M, [&, k](int i, int j) {
				if (x[i][k] != x[j][k]) return x[i][k] < x[j][k];
				return inv_a[k][i] < inv_a[k][j];
			}) - a[kp1].begin();
			// increment req_idx to the rightmost index such that matching over [req_idx, k+1] is not long enough
			/* the last individual site can be omitted from site k+1 here since a full match of site k+1
			   means the match will have a chance to be detected next time, unless this is the last site */
			U loc = locs[kp1 == n - 1 ? N - 1 : (kp1 + 1) * B - 2];
			while (req_idx < k && loc - locs[req_idx * B] >= L) req_idx++;
			if (req_idx == 0) continue; // no possible matches

			// subtract terminated matches from ongoing match counts
			up_on -= up_end[k];
			dn_on -= dn_end[k];
			// check for match immediately above query
			bool touch = false;
			if (up_on == 0 && t > 0) {
				int above = a[kp1][t - 1];
				if (dn_on > 0) touch = d[kp1][t] <= req_idx; // can use divergence value if there are ongoing matches below query
				else if (x[M][req_idx] == x[above][req_idx]) { // otherwise use hashes to check for match over sites [req_idx, k]
					lt = hz[req_idx] < h[above][req_idx] ? MOD - h[above][req_idx] + hz[req_idx] : hz[req_idx] - h[above][req_idx];
					rt = hz[kp1] < h[above][kp1] ? MOD - h[above][kp1] + hz[kp1] : hz[kp1] - h[above][kp1];
					touch = lt * xp[kp1 - req_idx] % MOD == rt;
				}

				if (touch) {
					// find beginning of match; total # iterations bounded by n for a given query
					int up_beg_ = req_idx;
					while (--up_beg_ > up_beg && x[M][up_beg_] == x[above][up_beg_]) {}
					up_beg = up_beg_;
				}
			}

			if (up_on > 0 || touch) {
				int p = t - up_on; // a[k+1] pointer for finding matches above query
				while (d[kp1][p] <= req_idx || touch) { // while there's another match
					if (touch) touch = false;
					else up_beg = max(up_beg, d[kp1][p] - 1);
					int above = a[kp1][--p];

					// briefly linear search for end of match
					int lo = kp1, pivot = min(n, lo + 10);
					while (lo < pivot && x[M][lo] == x[above][lo]) lo++;
					if (lo == pivot) { // switch to binary search if end not found
						int hi = n;
						lt = hz[k] < h[above][k] ? MOD - h[above][k] + hz[k] : hz[k] - h[above][k];

						while (lo < hi) {
							int g = (lo + hi + 1) >> 1;
							rt = hz[g] < h[above][g] ? MOD - h[above][g] + hz[g] : hz[g] - h[above][g];

							if (lt * xp[g - k] % MOD == rt) {
								lo = g;
							} else {
								hi = g - 1;
							}
						}
					}
					up_end[lo]++;

					// refine to single-site match [start, end)
					int start = 0, end = N, e_idx = lo;
					if (B == 64) {
						if (up_beg > -1) start = (up_beg + 1) * B - __builtin_ctzll(x[M][up_beg] ^ x[above][up_beg]);
						if (e_idx < n) end = e_idx * B + __builtin_clzll(x[M][e_idx] ^ x[above][e_idx]);
					} else {
						if (up_beg > -1) {
							unsigned long long right = x[M][up_beg] ^ x[above][up_beg];
							if (right == 0) {
								start = (up_beg + 1) * B - 64 - __builtin_ctzll((x[M][up_beg] >> 64) ^ (x[above][up_beg] >> 64));
							} else {
								start = (up_beg + 1) * B - __builtin_ctzll(right);
							}
						}
						if (e_idx < n) {
							unsigned long long left = (x[M][e_idx] >> 64) ^ (x[above][e_idx] >> 64);
							if (left == 0) {
								end = e_idx * B + 64 + __builtin_clzll(x[M][e_idx] ^ x[above][e_idx]);
							} else {
								end = e_idx * B + __builtin_clzll(left);
							}
						}
					}

					// report match if long
					U len = locs[end - 1] - locs[start];
					if (len >= L) {
						out << start << '\t' << end << '\t' << IDs[above] << '\t' << len << '\n';
					}
				}
				up_on = t - p; // add new matches to ongoing match count
			}

			// analogous process for below query

			touch = false;
			if (dn_on == 0 && t < M) {
				int below = a[kp1][t];
				if (up_on > 0) touch = d[kp1][t] <= req_idx;
				else if (x[below][req_idx] == x[M][req_idx]) {
					lt = hz[req_idx] < h[below][req_idx] ? MOD - h[below][req_idx] + hz[req_idx] : hz[req_idx] - h[below][req_idx];
					rt = hz[kp1] < h[below][kp1] ? MOD - h[below][kp1] + hz[kp1] : hz[kp1] - h[below][kp1];
					touch = lt * xp[kp1 - req_idx] % MOD == rt;
				}

				if (touch) {
					int dn_beg_ = req_idx;
					while (--dn_beg_ > dn_beg && x[M][dn_beg_] == x[below][dn_beg_]) {}
					dn_beg = dn_beg_;
				}
			}

			if (dn_on > 0 || touch) {
				int p = t + dn_on;
				while ((p < M && d[kp1][p] <= req_idx) || touch) {
					if (touch) touch = false;
					else dn_beg = max(dn_beg, d[kp1][p] - 1);
					int below = a[kp1][p++];

					int lo = kp1, pivot = min(n, lo + 10);
					while (lo < pivot && x[M][lo] == x[below][lo]) lo++;
					if (lo == pivot) {
						int hi = n;
						lt = hz[k] < h[below][k] ? MOD - h[below][k] + hz[k] : hz[k] - h[below][k];

						while (lo < hi) {
							int g = (lo + hi + 1) >> 1;
							rt = hz[g] < h[below][g] ? MOD - h[below][g] + hz[g] : hz[g] - h[below][g];

							if (lt * xp[g - k] % MOD == rt) {
								lo = g;
							} else {
								hi = g - 1;
							}
						}
					} 
					dn_end[lo]++;

					int start = 0, end = N, e_idx = lo;
					if (B == 64) {
						if (dn_beg > -1) start = (dn_beg + 1) * B - __builtin_ctzll(x[M][dn_beg] ^ x[below][dn_beg]);
						if (e_idx < n) end = e_idx * B + __builtin_clzll(x[M][e_idx] ^ x[below][e_idx]);
					} else {
						if (dn_beg > -1) {
							unsigned long long right = x[M][dn_beg] ^ x[below][dn_beg];
							if (right == 0) {
								start = (dn_beg + 1) * B - 64 - __builtin_ctzll((x[M][dn_beg] >> 64) ^ (x[below][dn_beg] >> 64));
							} else {
								start = (dn_beg + 1) * B - __builtin_ctzll(right);
							}
						}
						if (e_idx < n) {
							unsigned long long left = (x[M][e_idx] >> 64) ^ (x[below][e_idx] >> 64);
							if (left == 0) {
								end = e_idx * B + 64 + __builtin_clzll(x[M][e_idx] ^ x[below][e_idx]);
							} else {
								end = e_idx * B + __builtin_clzll(left);
							}
						}
					}

					U len = locs[end - 1] - locs[start];
					if (len >= L) {
						out << start << '\t' << end << '\t' << IDs[below] << '\t' << len << '\n';
					}
				}
				dn_on = p - t;
			}
		}
	}

	return 0;
}

template<class T>
int SparseQuery<T>::query(const char* query_file, const char* output_dir, const int L) { // 1-length sites
	if (L < minSiteL) return 4;
	if (access(output_dir, W_OK) != 0) return 2; 
	ifstream in(query_file);
        if (in.fail()) return 1;

	const int l = (L - B + 1) / B; // min # sparse sites that must be covered by an L-site match
	string line;
	while (getline(in, line)) {
		if (line.size() < 2u) return 5;
		if (line[0] != '#' || line[1] != '#') break;
	}
	int Q = -9;
	stringstream ss(line);
	while (getline(ss, line, '\t')) Q++;
	if (Q < 1) return 5;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	vector<string> qIDs(Q);
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < Q; i++) getline(ss, qIDs[i], '\t');
	Q <<= 1;
	vector<vector<T>> z(n, vector<T>(Q));

	for (int K = 0; K < N; K++) {
		if (!getline(in, line)) return 5;
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int k = K / B, i = 0;
		while (getline(ss, line, '\t')) {
			if (i == Q || line.size() < 3u) return 5;
			z[k][i] = (z[k][i] << 1) | (line[0] != '0'), i++;
                        z[k][i] = (z[k][i] << 1) | (line[2] != '0'), i++;
		}
	}
	in.close();

	unsigned __int128 lt, rt;

	for (int q = 0; q < Q; q++) {
		ofstream out(string(output_dir) + "/" + qIDs[q >> 1] + "-" + to_string(q & 1) + ".txt");

		z[n - 1][q] <<= n * B - N;
		for (int k = 0, kp1 = 1; k < n; k++, kp1++) {
			x[M][k] = z[k][q];
			hz[kp1] = (unsigned __int128) hz[k] * xp[1] % MOD;
			unsigned __int128 xp3b = 1;
			T x_copy = x[M][k];
			for (int b = 0; b < B; b++) {
				hz[kp1] = (hz[kp1] + ((x_copy & 1) + 1) * xp3b) % MOD;
				xp3b = xp3b * 3 % MOD;
				x_copy >>= 1;
			}
		}
		memset(&up_end[1], 0, n * sizeof(int));
		memset(&dn_end[1], 0, n * sizeof(int));

		for (int k = 0, kp1 = 1, t = 0, up_on = 0, dn_on = 0; k < n; k++, kp1++) {
			inv_a[k][M] = t;
			t = lower_bound(a[kp1].begin(), a[kp1].end(), M, [&, k](int i, int j) {
				if (x[i][k] != x[j][k]) return x[i][k] < x[j][k];
				return inv_a[k][i] < inv_a[k][j];
			}) - a[kp1].begin();
			if (kp1 < l) continue;
			int s_idx = k - l;

			up_on -= up_end[k];
			dn_on -= dn_end[k];
			bool touch = false;
			if (up_on == 0 && t > 0) {
				int above = a[kp1][t - 1];
				if (dn_on > 0) touch = d[kp1][t] <= kp1 - l;
				else if (x[M][kp1 - l] == x[above][kp1 - l]) {
					lt = hz[kp1 - l] < h[above][kp1 - l] ? MOD - h[above][kp1 - l] + hz[kp1 - l] : hz[kp1 - l] - h[above][kp1 - l];
					rt = hz[kp1] < h[above][kp1] ? MOD - h[above][kp1] + hz[kp1] : hz[kp1] - h[above][kp1];
					touch = lt * xp[l] % MOD == rt;
				}
			}

			if (up_on > 0 || touch) {
				int p = t - up_on;
				while (d[kp1][p] <= kp1 - l || touch) {
					touch = false;
					int above = a[kp1][--p];

					int lo = kp1, pivot = min(n, lo + 10);
					while (lo < pivot && x[M][lo] == x[above][lo]) lo++;
					if (lo == pivot) {
						int hi = n;
						lt = hz[k] < h[above][k] ? MOD - h[above][k] + hz[k] : hz[k] - h[above][k];

						while (lo < hi) {
							int g = (lo + hi + 1) >> 1;
							rt = hz[g] < h[above][g] ? MOD - h[above][g] + hz[g] : hz[g] - h[above][g];

							if (lt * xp[g - k] % MOD == rt) {
								lo = g;
							} else {
								hi = g - 1;
							}
						}
					}
					up_end[lo]++;

					int start = 0, end = N, e_idx = lo;
					if (B == 64) {
						if (s_idx > -1) start = (s_idx + 1) * B - __builtin_ctzll(x[M][s_idx] ^ x[above][s_idx]);
						if (e_idx < n) end = e_idx * B + __builtin_clzll(x[M][e_idx] ^ x[above][e_idx]);
					} else {
						if (s_idx > -1) {
							unsigned long long right = x[M][s_idx] ^ x[above][s_idx];
							if (right == 0) {
								start = (s_idx + 1) * B - 64 - __builtin_ctzll((x[M][s_idx] >> 64) ^ (x[above][s_idx] >> 64));
							} else {
								start = (s_idx + 1) * B - __builtin_ctzll(right);
							}
						}
						if (e_idx < n) {
							unsigned long long left = (x[M][e_idx] >> 64) ^ (x[above][e_idx] >> 64);
							if (left == 0) {
								end = e_idx * B + 64 + __builtin_clzll(x[M][e_idx] ^ x[above][e_idx]);
							} else {
								end = e_idx * B + __builtin_clzll(left);
							}
						}
					}

					if (end - start >= L) {
						out << start << '\t' << end << '\t' << IDs[above] << '\n';
					}
				}
				up_on = t - p;
			}

			touch = false;
			if (dn_on == 0 && t < M) {
				int below = a[kp1][t];
				if (up_on > 0) touch = d[kp1][t] <= kp1 - l;
				else if (x[below][kp1 - l] == x[M][kp1 - l]) {
					lt = hz[kp1 - l] < h[below][kp1 - l] ? MOD - h[below][kp1 - l] + hz[kp1 - l] : hz[kp1 - l] - h[below][kp1 - l];
					rt = hz[kp1] < h[below][kp1] ? MOD - h[below][kp1] + hz[kp1] : hz[kp1] - h[below][kp1];
					touch = lt * xp[l] % MOD == rt;
				}
			}

			if (dn_on > 0 || touch) {
				int p = t + dn_on;
				while ((p < M && d[kp1][p] <= kp1 - l) || touch) {
					touch = false;
					int below = a[kp1][p++];

					int lo = kp1, pivot = min(n, lo + 10);
					while (lo < pivot && x[M][lo] == x[below][lo]) lo++;
					if (lo == pivot) {
						int hi = n;
						lt = hz[k] < h[below][k] ? MOD - h[below][k] + hz[k] : hz[k] - h[below][k];

						while (lo < hi) {
							int g = (lo + hi + 1) >> 1;
							rt = hz[g] < h[below][g] ? MOD - h[below][g] + hz[g] : hz[g] - h[below][g];

							if (lt * xp[g - k] % MOD == rt) {
								lo = g;
							} else {
								hi = g - 1;
							}
						}
					} 
					dn_end[lo]++;

					int start = 0, end = N, e_idx = lo;
					if (B == 64) {
						if (s_idx > -1) start = (s_idx + 1) * B - __builtin_ctzll(x[M][s_idx] ^ x[below][s_idx]);
						if (e_idx < n) end = e_idx * B + __builtin_clzll(x[M][e_idx] ^ x[below][e_idx]);
					} else {
						if (s_idx > -1) {
							unsigned long long right = x[M][s_idx] ^ x[below][s_idx];
							if (right == 0) {
								start = (s_idx + 1) * B - 64 - __builtin_ctzll((x[M][s_idx] >> 64) ^ (x[below][s_idx] >> 64));
							} else {
								start = (s_idx + 1) * B - __builtin_ctzll(right);
							}
						}
						if (e_idx < n) {
							unsigned long long left = (x[M][e_idx] >> 64) ^ (x[below][e_idx] >> 64);
							if (left == 0) {
								end = e_idx * B + 64 + __builtin_clzll(x[M][e_idx] ^ x[below][e_idx]);
							} else {
								end = e_idx * B + __builtin_clzll(left);
							}
						}
					}

					if (end - start >= L) {
						out << start << '\t' << end << '\t' << IDs[below] << '\n';
					}
				}
				dn_on = p - t;
			}
		}
	}

	return 0;
}
