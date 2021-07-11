#include "Reporter.h"

template<class T>
int Reporter<T>::precompute(const char* input_file) {
	ifstream in(input_file);
	if (in.fail()) return 1;

	// get M = # haplotypes, N = # sites, n = # merged sites, sampleIDs
	string line;
	while (getline(in, line)) {
		if (line.size() < 2u) return 2;
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	M = -9;
	while (getline(ss, line, '\t')) M++;
	if (M < 1) return 2;
	sampleIDs.resize(M);
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
	for (int i = 0; i < M; i++) getline(ss, sampleIDs[i], '\t');
	M <<= 1;

	// resize data structures to panel with merged sites
	x.resize(M + 1, vector<T>(n));
	h.resize(M, vector<unsigned long long>(n + 1));
	inv_a.resize(n, vector<int>(M + 1));
	a.resize(n + 1, vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0);
	d.resize(n + 1, vector<int>(M));
	xp.resize(n + 1), hz.resize(n + 1), up_end.resize(n), dn_end.resize(n);

	// recurring variables for hash comparison and binary search
	unsigned __int128 lt, rt;
	int lo, hi, g;

	// K = current site, k = current merged site, kp1 = k+1
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
		// update merged panel with current site
		int index = 0;
		while (getline(ss, line, '\t')) {
			if (index == M || line.size() < 3u) return 2;
			x[index][k] = (x[index][k] << 1) | (line[0] != '0'), index++;
			x[index][k] = (x[index][k] << 1) | (line[2] != '0'), index++;
		}
		if (index != M) return 2;

		if (K == N - 1) { // if last site, pad last merged site with 0s
			int pad = n * B - N;
			for (int i = 0; i < M; i++) x[i][k] <<= pad;
			K += pad;
		}
		if ((K + 1) % B > 0) continue; // skip if not end of merged site

		// min query length = max distance across any 2 adjacent merged sites
		// xp[i] = pow(2, i * B) % MOD
		if (k > 0) {
			minPhysL = max(minPhysL, physLocs[min(K, N - 1)] - physLocs[(k - 1) * B]);
			xp[kp1] = (unsigned __int128) xp[k] * xp[1] % MOD;
		} else {
			xp[0] = xp[1] = 1;
			for (int i = 0; i < B; i++) {
				xp[1] <<= 1;
				if (xp[1] >= MOD) xp[1] -= MOD;
			}
		}

		// extend prefix hashes; h[i][k+1] = hash of x[i][0, k]
		for (int i = 0; i < M; i++) {
			h[i][kp1] = ((unsigned __int128) h[i][k] * xp[1] + x[i][k] % MOD) % MOD;
		}

		// inv_a[k][i] = j such that a[k][j] = i
		for (int i = 0; i < M; i++) {
			inv_a[k][a[k][i]] = i;
		}

		// sort haplotypes in a[k+1] by (value at site k, position in a[k])
		iota(a[kp1].begin(), a[kp1].end(), 0);
		sort(a[kp1].begin(), a[kp1].end(), [&, k](int i, int j) {
			if (x[i][k] != x[j][k]) return x[i][k] < x[j][k];
			return inv_a[k][i] < inv_a[k][j];
		});

		d[kp1][0] = kp1;
		for (int i = 1; i < M; i++) {
			int x_cur = a[kp1][i], x_above = a[kp1][i - 1], a_prev = inv_a[k][x_cur];

			if (x[x_cur][k] != x[x_above][k]) { // no match
				d[kp1][i] = kp1;
			} else if (a_prev > 0 && x_above == a[k][a_prev - 1]) { // continued match
				d[kp1][i] = d[k][a_prev];
			} else { // different match, binary search for divergence value using hashes
				lo = d[k][a_prev], hi = k;
				rt = h[x_cur][kp1] < h[x_above][kp1] ? MOD - h[x_above][kp1] + h[x_cur][kp1] : h[x_cur][kp1] - h[x_above][kp1];

				while (lo < hi) {
					g = (lo + hi) >> 1;
					lt = h[x_cur][g] < h[x_above][g] ? MOD - h[x_above][g] + h[x_cur][g] : h[x_cur][g] - h[x_above][g];

					if (lt * xp[kp1 - g] % MOD == rt) {
						hi = g;
					} else {
						lo = g + 1;
					}
				}

				d[kp1][i] = lo;
			}
		}

		k++, kp1++;
	}

	return 0;
}

template<class T>
int Reporter<T>::save(const char* save_file) {
	ofstream out(save_file, ios::binary);
	if (out.fail()) return 1;

	out.write((char*) &M, sizeof(M));
	out.write((char*) &N, sizeof(N));
	out.write((char*) &B, sizeof(B));
	out.write((char*) &minPhysL, sizeof(minPhysL));
	for (int loc : physLocs) out.write((char*) &loc, sizeof(loc));
	for (string ID : sampleIDs) {
		size_t sz = ID.size();
		out.write((char*) &sz, sizeof(sz));
		out.write(ID.c_str(), sz);
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < n; j++) {
			out.write((char*) &x[i][j], sizeof(T));
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
int Reporter<T>::load(const char* load_file) {
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
	sampleIDs.resize(M / 2);
	for (int i = 0; i < M / 2; i++) {
		size_t sz;
		if (!in.read((char*) &sz, sizeof(sz)) || sz < 1) return 2;
		sampleIDs[i].resize(sz);
		if (!in.read(&sampleIDs[i][0], sz)) return 2;
	}
	n = (N + B - 1) / B, minSiteL = B * 2 - 1;
	xp.resize(n + 1), hz.resize(n + 1), up_end.resize(n), dn_end.resize(n);
	xp[0] = xp[1] = 1;
	for (int i = 0; i < B; i++) {
		xp[1] <<= 1;
		if (xp[1] >= MOD) xp[1] -= MOD;
	}
	for (int i = 1; i < n; i++) xp[i + 1] = (unsigned __int128) xp[i] * xp[1] % MOD;
	x.resize(M + 1, vector<T>(n));
	h.resize(M, vector<unsigned long long>(n + 1));
	for (int i = 0; i < M; i++) {
		for (int j = 0, jp1 = 1; j < n; j++, jp1++) {
			if (!in.read((char*) &x[i][j], sizeof(T))) return 2;
			h[i][jp1] = ((unsigned __int128) h[i][j] * xp[1] + x[i][j] % MOD) % MOD;
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
void Reporter<T>::show_attributes() {
	cout << "\nData attributes:\n" <<
		"\tM = " << M << " haplotypes, N = " << N << " sites, B = " << B << '\n';
        if (minGenL == -1) {
		cout << "\tYour minimum allowable query lengths are " << minSiteL << " sites or " << minPhysL << " bps.\n" <<
                        "\tProvide a genetic map to query in units of cM." << endl;
	} else {
                cout << "\tYour minimum allowable query lengths are " << minGenL << " cM, " << minSiteL << " sites, or " << minPhysL << " bps." << endl;
        }
}

template<class T>
int Reporter<T>::set_gen_map(const char* map_file) {
	ifstream in(map_file);
	if (in.fail()) return 1;

	if ((int) genLocs.size() != N) genLocs.resize(N);
	minGenL = -1;
	string line;
	for (int i = 0; i < N; i++) {
		if (!getline(in, line)) return 2;
		stringstream ss(line);
		while (getline(ss, line, '\t')) {}
		try { genLocs[i] = stod(line); }
		catch (out_of_range &e) {}
		catch (exception &e) { return 3; }
	}
	// min query length = max distance across any 2 adjacent merged sites
	for (int i = 1; i < n; i++) {
		minGenL = max(minGenL, genLocs[min(i * B, N) - 1] - genLocs[(i - 1) * B]);
	}

	return 0;
}

template<class T>
template<class U>
int Reporter<T>::query(const char* query_file, const char* output_dir, const U L, const vector<U> &locs) { // variably distributed site locations
	if (sizeof(U) == sizeof(double)) {
		if (minGenL == -1) return 3;
		if (L < minGenL) return 4;
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
	vector<string> queryIDs(Q);
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < Q; i++) getline(ss, queryIDs[i], '\t');
	Q <<= 1;
	vector<vector<T>> z(n, vector<T>(Q));

	// read query panel; "site" will now refer to "merged site"
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
        int lo, hi, g;
        const T ones = -1;

	for (int q = 0; q < Q; q++) {
		ofstream out(string(output_dir) + "/" + queryIDs[q >> 1] + "-" + to_string(q & 1) + ".txt");

		z[n - 1][q] <<= n * B - N; // pad last site with 0s
		for (int i = 0, ip1 = 1; i < n; i++, ip1++) {
			x[M][i] = z[i][q]; // copy query into panel
			hz[ip1] = ((unsigned __int128) hz[i] * xp[1] + x[M][i] % MOD) % MOD; // build prefix hashes of query
		}
		// up_end[i] = # matches (above the query) that end inclusively at site i, dn_end[i] = ... below ... (above/below refer to positions in a[i])
		memset(&up_end[0], 0, n * sizeof(int));
		memset(&dn_end[0], 0, n * sizeof(int));

		/* i = current site, ip1 = i+1, req_idx = rightmost site that must fully match, up_idx = non-inclusive start site of above match,
		   dn_idx = ... below ..., up_ct = # ongoing above matches, dn_ct = ... below ... */
		for (int i = 0, ip1 = 1, t = 0, req_idx = 0, up_idx = -1, dn_idx = -1, up_ct = 0, dn_ct = 0; ip1 < n; i++, ip1++) {
			// binary search for query's position in a[i+1], according to sorting by (value at site i, position in a[i])
			inv_a[i][M] = t;
			t = lower_bound(a[ip1].begin(), a[ip1].end(), M, [&, i](int i1, int i2) { if (x[i1][i] != x[i2][i]) return x[i1][i] < x[i2][i]; return inv_a[i][i1] < inv_a[i][i2]; }) - a[ip1].begin();
			// increment req_idx to the rightmost index such that matching over [req_idx, i+1] is not long enough
			/* the last individual site can be omitted from site i+1 here since a full match of site i+1
			   means the match will have a chance to be detected next time, unless this is the last site */
			U loc = locs[ip1 == n - 1 ? N - 1 : (ip1 + 1) * B - 2];
			while (req_idx < i && loc - locs[req_idx * B] >= L) req_idx++;
			if (req_idx == 0) continue; // no possible matches

			// check for match immediately above query
			bool touch = false;
			if (up_ct == 0 && t > 0) {
				int above = a[ip1][t - 1];
				if (dn_ct > 0) touch = d[ip1][t] <= req_idx; // can use divergence value if there are ongoing matches below query
				else if (x[M][i] == x[above][i]) { // otherwise use hashes to check for match over sites [req_idx, i]
					lt = hz[req_idx] < h[above][req_idx] ? MOD - h[above][req_idx] + hz[req_idx] : hz[req_idx] - h[above][req_idx];
					rt = hz[ip1] < h[above][ip1] ? MOD - h[above][ip1] + hz[ip1] : hz[ip1] - h[above][ip1];
					touch = lt * xp[ip1 - req_idx] % MOD == rt;
				}

				if (touch) {
					// find beginning of match; total # increments bounded by n for a given query
					hi = req_idx - 1;
					int pivot = max(up_idx, hi - 10);
					while (hi > pivot && x[M][hi] == x[above][hi]) hi--;
					if (hi == pivot) {
						lo = up_idx + 1, hi++;
						rt = hz[ip1] < h[above][ip1] ? MOD - h[above][ip1] + hz[ip1] : hz[ip1] - h[above][ip1];

						while (lo < hi) {
							g = (lo + hi) >> 1;
							lt = hz[g] < h[above][g] ? MOD - h[above][g] + hz[g] : hz[g] - h[above][g];

							if (lt * xp[ip1 - g] % MOD == rt) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						up_idx = hi - 1;
					} else up_idx = hi;
				}
			}

			if (up_ct > 0 || touch) {
				int p = t - up_ct; // a[i+1] pointer for finding matches above query
				while (d[ip1][p] <= req_idx || touch) { // while there's another match
					if (touch) touch = false;
					else up_idx = max(up_idx, d[ip1][p] - 1);
					int above = a[ip1][--p];

					// briefly linear search for end of match
					lo = ip1;
					int pivot = min(n, lo + 10);
					while (lo < pivot && x[M][lo] == x[above][lo]) lo++;
					if (lo == pivot) { // switch to binary search if end not found
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

					// binary search for individual site resolution ends of match
					int start = 0, end = N, e_idx = lo;
					if (e_idx < n) {
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][e_idx] >> g) == (x[above][e_idx] >> g)) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						end = (e_idx + 1) * B - lo;
					}
					if (up_idx > -1) {
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][up_idx] & (ones >> g)) == (x[above][up_idx] & (ones >> g))) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						start = up_idx * B + lo;
					}

					// report match if long
					U len = locs[end - 1] - locs[start];
					if (len >= L) {
						out << start << '\t' << end << '\t' << sampleIDs[above >> 1] << '-' << (above & 1) << '\t' << len << '\n';
					}
				}
				up_ct = t - p - up_end[i]; // update # ongoing matches
			}

			// analogous process for below query

			touch = false;
			if (dn_ct == 0 && t < M) {
				int below = a[ip1][t];
				if (up_ct > 0) touch = d[ip1][t] <= req_idx;
				else if (x[below][i] == x[M][i]) {
					lt = hz[req_idx] < h[below][req_idx] ? MOD - h[below][req_idx] + hz[req_idx] : hz[req_idx] - h[below][req_idx];
					rt = hz[ip1] < h[below][ip1] ? MOD - h[below][ip1] + hz[ip1] : hz[ip1] - h[below][ip1];
					touch = lt * xp[ip1 - req_idx] % MOD == rt;
				}

				if (touch) {
					hi = req_idx - 1;
					int pivot = max(dn_idx, hi - 10);
					while (hi > pivot && x[M][hi] == x[below][hi]) hi--;
					if (hi == pivot) {
						lo = dn_idx + 1, hi++;
						rt = hz[ip1] < h[below][ip1] ? MOD - h[below][ip1] + hz[ip1] : hz[ip1] - h[below][ip1];

						while (lo < hi) {
							g = (lo + hi) >> 1;
							lt = hz[g] < h[below][g] ? MOD - h[below][g] + hz[g] : hz[g] - h[below][g];

							if (lt * xp[ip1 - g] % MOD == rt) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						dn_idx = hi - 1;
					} else dn_idx = hi;
				}
			}

			if (dn_ct > 0 || touch) {
				int p = t + dn_ct;
				while ((p < M && d[ip1][p] <= req_idx) || touch) {
					if (touch) touch = false;
					else dn_idx = max(dn_idx, d[ip1][p] - 1);
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
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][e_idx] >> g) == (x[below][e_idx] >> g)) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						end = (e_idx + 1) * B - lo;
					}
					if (dn_idx > -1) {
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][dn_idx] & (ones >> g)) == (x[below][dn_idx] & (ones >> g))) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						start = dn_idx * B + lo;
					}

					U len = locs[end - 1] - locs[start];
					if (len >= L) {
						out << start << '\t' << end << '\t' << sampleIDs[below >> 1] << '-' << (below & 1) << '\t' << len << '\n';
					}
				}
				dn_ct = p - t - dn_end[i];
			}
		}
	}

	return 0;
}

template<class T>
int Reporter<T>::query(const char* query_file, const char* output_dir, const int L) { // 1-length sites
	if (L < minSiteL) return 4;
	if (access(output_dir, W_OK) != 0) return 2; 
	ifstream in(query_file);
        if (in.fail()) return 1;

	const int l = (L - B + 1) / B; // min # merged sites that must be covered by an L-site match
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
	vector<string> queryIDs(Q);
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < Q; i++) getline(ss, queryIDs[i], '\t');
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
        int lo, hi, g;
        const T ones = -1;

	for (int q = 0; q < Q; q++) {
		ofstream out(string(output_dir) + "/" + queryIDs[q >> 1] + "-" + to_string(q & 1) + ".txt");

		z[n - 1][q] <<= n * B - N;
                for (int i = 0, ip1 = 1; i < n; i++, ip1++) {
                        x[M][i] = z[i][q];
                        hz[ip1] = ((unsigned __int128) hz[i] * xp[1] + x[M][i] % MOD) % MOD;
                }
		memset(&up_end[0], 0, n * sizeof(int));
		memset(&dn_end[0], 0, n * sizeof(int));

		for (int i = 0, ip1 = 1, t = 0, up_ct = 0, dn_ct = 0; i < n; i++, ip1++) {
			inv_a[i][M] = t;
			t = lower_bound(a[ip1].begin(), a[ip1].end(), M, [&, i](int i1, int i2) { if (x[i1][i] != x[i2][i]) return x[i1][i] < x[i2][i]; return inv_a[i][i1] < inv_a[i][i2]; }) - a[ip1].begin();
			if (ip1 < l) continue;

			int s_idx = i - l, min_start = max(0, s_idx * B + 1);
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
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][e_idx] >> g) == (x[above][e_idx] >> g)) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						end = (e_idx + 1) * B - lo;
					}
					if (end - min_start < L) continue;
					if (s_idx > -1) {
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][s_idx] & (ones >> g)) == (x[above][s_idx] & (ones >> g))) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						start = s_idx * B + lo;
					}

					if (end - start >= L) {
						out << start << '\t' << end << '\t' << sampleIDs[above >> 1] << '-' << (above & 1) << '\n';

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
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][e_idx] >> g) == (x[below][e_idx] >> g)) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						end = (e_idx + 1) * B - lo;
					}
					if (end - min_start < L) continue;
					if (s_idx > -1) {
						lo = 1, hi = B;
						while (lo < hi) {
							g = (lo + hi) >> 1;
							if ((x[M][s_idx] & (ones >> g)) == (x[below][s_idx] & (ones >> g))) {
								hi = g;
							} else {
								lo = g + 1;
							}
						}
						start = s_idx * B + lo;
					}

					if (end - start >= L) {
						out << start << '\t' << end << '\t' << sampleIDs[below >> 1] << '-' << (below & 1) << '\n';
					}
				}
				dn_ct = p - t - dn_end[i];
			}
		}
	}

	return 0;
}
