#include "SyllableQuery.h"
// This code is licensed under MIT license (see LICENSE for details)

template<class T>
int SyllableQuery<T>::precompute(const char* input_file) {
	ifstream in(input_file);
	if (in.fail()) return 1;

	// get M = # haplotypes, N = # sites, n = # syllables, IDs
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

	// resize data structures for syllabic panel
	x.resize(M, vector<int>(n));
	r.resize(n);
	a.resize(n + 1, vector<int>(M));
	iota(a[0].begin(), a[0].end(), 0);

	{ vector<T> x_(M);

	// K = current site, k = current syllable, kp1 = k+1
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
		// update current syllable with current site
		int index = 0;
		while (getline(ss, line, '\t')) {
			if (index == M || line.size() < 3u) return 2;
			x_[index] = (x_[index] << 1) | (line[0] != '0'), index++;
			x_[index] = (x_[index] << 1) | (line[2] != '0'), index++;
		}
		if (index != M) return 2;

		if (K == N - 1) { // if last site, pad last syllable with 0s
			int pad = n * B - N;
			for (int i = 0; i < M; i++) x_[i] <<= pad;
			K += pad;
		}
		if ((K + 1) % B > 0) continue; // skip if not end of syllable

		// coordinate-compress syllable values
		vector<T> x_copy = x_;
		sort(x_copy.begin(), x_copy.end());
		x_copy.resize(unique(x_copy.begin(), x_copy.end()) - x_copy.begin());
		r[k] = x_copy;
		for (int i = 0; i < M; i++) {
			x[i][k] = lower_bound(r[k].begin(), r[k].end(), x_[i]) - r[k].begin();
		}
		memset(&x_[0], 0, M * sizeof(T));

		// counting sort to build a[k+1] with a[k] and x[][k]
		int r_sz = r[k].size();
		vector<vector<int>> counter(r_sz);
		for (int i = 0; i < M; i++) {
			int a_ = a[k][i];
			counter[x[a_][k]].push_back(a_);
		}
		for (int i = 0, p = 0; i < r_sz; i++) {
			for (int v : counter[i]) {
				a[kp1][p++] = v;
			}
		}

		// min query length = max distance across any 2 adjacent syllables
		if (k > 0) minPhysL = max(minPhysL, (K >= N - 1 ? physLocs[N - 1] + 1 : physLocs[K]) - physLocs[(k - 1) * B]);

		k++, kp1++;
	}}

	xp.resize(n), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1), z_.resize(n), z.resize(n);
	// xp[i] = pow(BASE, i) % MOD
	xp[0] = 1;
	for (int k = 1; k < n; k++) {
		xp[k] = xp[k-1] * BASE % MOD;
	}

	// build prefix hashes; h[i][k+1] = hash of x[i][0, k]
	h.resize(M, vector<unsigned long long>(n + 1));
	for (int i = 0; i < M; i++) {
		for (int k = 0; k < n; k++) {
			h[i][k+1] = (h[i][k] + xp[k] * (x[i][k] + 1)) % MOD;
		}
	}

	return 0;
}

template<class T>
int SyllableQuery<T>::save(const char* save_file) {
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
			out.write((char*) &x[i][j], sizeof(int));
		}
	}
	for (int i = 0; i < n; i++) {
		size_t sz = r[i].size();
		out.write((char*) &sz, sizeof(sz));
		for (T v : r[i]) {
			out.write((char*) &v, sizeof(T));
		}
	}
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			out.write((char*) &a[i][j], sizeof(int));
		}
	}

	if (minGeneL != -1) {
		out.write((char*) &minGeneL, sizeof(double));
		for (int i = 0; i < N; i++) {
			out.write((char*) &geneLocs[i], sizeof(double));
		}
	}

	return 0;
}

template<class T>
int SyllableQuery<T>::load(const char* load_file) {
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
	xp.resize(n), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1), z_.resize(n), z.resize(n);
	xp[0] = 1;
	for (int i = 1; i < n; i++) {
		xp[i] = xp[i - 1] * BASE % MOD;
	}
	x.resize(M + 1, vector<int>(n));
	h.resize(M, vector<unsigned long long>(n + 1));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < n; j++) {
			if (!in.read((char*) &x[i][j], sizeof(int))) return 2;
			h[i][j+1] = (h[i][j] + xp[j] * (x[i][j] + 1)) % MOD;
		}
	}
	r.resize(n);
	for (int i = 0; i < n; i++) {
		size_t sz;
		if (!in.read((char*) &sz, sizeof(sz)) || sz < 1) return 2;
		r[i].resize(sz);
		for (int j = 0; j < (int) sz; j++) {
			if (!in.read((char*) &r[i][j], sizeof(T))) return 2;
		}
	}
	a.resize(n + 1, vector<int>(M));
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < M; j++) {
			if (!in.read((char*) &a[i][j], sizeof(int)) || !(-1 < a[i][j] && a[i][j] < M)) return 2;
		}
	}

	if (!in.read((char*) &minGeneL, sizeof(double))) {
		minGeneL = -1;
		return 0;
	}
	geneLocs.resize(N);
	for (int i = 0; i < N; i++) {
		if (!in.read((char*) &geneLocs[i], sizeof(double))) {
			minGeneL = -1;
			return 3;
		}
	}

	return 0;
}

template<class T>
string SyllableQuery<T>::show_attributes() {
	string attr = "\nData attributes:\n"
		"\tM = " + to_string(M) + " haplotypes, N = " + to_string(N) + " sites, B = " + to_string(B) + " sites\n";
	if (minGeneL == -1) {
		attr += "\tYour minimum allowable query lengths are " + to_string(minSiteL) + " sites or " + to_string(minPhysL) + " bps.\n"
		"\tThe server program must have provided a genetic map to allow queries in cM.";
	} else {
		attr += "\tYour minimum allowable query lengths are " + to_string(minGeneL) + " cM, " + to_string(minSiteL) + " sites, or " + to_string(minPhysL) + " bps.";
	}
	return attr;
}

template<class T>
int SyllableQuery<T>::set_gen_map(const char* map_file) {
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
		catch (exception &e) { i--; }
	}
	if (getline(in, line)) return 2;
	// min query length = max distance across any 2 adjacent syllables
	for (int i = 1; i < n; i++) {
		minGeneL = max(minGeneL, (i == n - 1 ? geneLocs[N - 1] + 1e-6 : geneLocs[(i + 1) * B - 1]) - geneLocs[(i - 1) * B]);
	}

	return 0;
}

template <>
inline void SyllableQuery<unsigned long long>::refine(int s_idx, int e_idx, int xi, int *start, int *end) {
	if (s_idx > -1) {
		*start = (s_idx + 1) * B - __builtin_ctzll(z_[s_idx] ^ r[s_idx][x[xi][s_idx]]);
	} else *start = 0;
	if (e_idx < n) {
		*end = e_idx * B + __builtin_clzll(z_[e_idx] ^ r[e_idx][x[xi][e_idx]]);
	} else *end = N;
}

template <>
inline void SyllableQuery<unsigned __int128>::refine(int s_idx, int e_idx, int xi, int *start, int *end) {
	if (s_idx > -1) {
		unsigned __int128 x_val = r[s_idx][x[xi][s_idx]];
		unsigned long long right = z_[s_idx] ^ x_val;
		if (right == 0) {
			*start = (s_idx + 1) * B - 64 - __builtin_ctzll((z_[s_idx] >> 64) ^ (x_val >> 64));
		} else {
			*start = (s_idx + 1) * B - __builtin_ctzll(right);
		}
	} else *start = 0;
	if (e_idx < n) {
		unsigned __int128 x_val = r[e_idx][x[xi][e_idx]];
		unsigned long long left = (z_[e_idx] >> 64) ^ (x_val >> 64);
		if (left == 0) {
			*end = e_idx * B + 64 + __builtin_clzll(z_[e_idx] ^ x_val);
		} else {
			*end = e_idx * B + __builtin_clzll(left);
		}
	} else *end = N;
}

template<class T>
template<class U>
int SyllableQuery<T>::query(const char* query_file, const char* output_file, const U L, const vector<U> &locs, bool inclusive) { // variably distributed site locations
	if (sizeof(U) == sizeof(double)) {
		if (minGeneL == -1) return 3;
		if (L < minGeneL) return 4;
	} else {
		if (L < minPhysL) return 4;
	}
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
	vector<vector<T>> Z(n, vector<T>(Q));

	// read query panel
	for (int K = 0; K < N; K++) {
		if (!getline(in, line)) return 5;
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int k = K / B, i = 0;
		while (getline(ss, line, '\t')) {
			if (i == Q || line.size() < 3u) return 5;
			Z[k][i] = (Z[k][i] << 1) | (line[0] != '0'), i++;
			Z[k][i] = (Z[k][i] << 1) | (line[2] != '0'), i++;
		}
	}
	in.close();
	for (int q = 0; q < Q; q++) Z[n - 1][q] <<= n * B - N; // pad last syllable with 0s

	ofstream out(output_file);
	if (out.fail()) return 2;
	for (int q = 0; q < Q; q++) {
		string qID = qIDs[q >> 1] + "-" + to_string(q & 1);

		for (int k = 0; k < n; k++) {
			z_[k] = Z[k][q];
			z[k] = lower_bound(r[k].begin(), r[k].end(), z_[k]) - r[k].begin(); // binary search for compressed syllable value
			if (z[k] < (int) r[k].size() && z_[k] != r[k][z[k]]) z[k] = r[k].size(); // assign unique value if unique raw value
			hz[k+1] = (hz[k] + xp[k] * (z[k] + 1)) % MOD;
		}
		// up_end[k] = # matches (above the query) that end at syllable k, dn_end[k] = ... below ... (above/below refer to positions in a[k])
		memset(&up_end[1], 0, n * sizeof(int));
		memset(&dn_end[1], 0, n * sizeof(int));

		/* k = current syllable, kp1 = k+1, req_idx = rightmost syllable that must fully match,
		   up_on = # ongoing above matches, dn_on = ... below ... */
		for (int k = 0, kp1 = 1, req_idx = 0, up_on = 0, dn_on = 0, no_beg = 0; kp1 < n; k++, kp1++) {
			// binary search for virtual position t of z with x[][k] and hashes
			int t;
			if (z[k] == (int) r[k].size()) {
				t = M;
				no_beg = k + 1;
			} else {
				int lo_i = 0, hi_i = M, tu_beg = k, td_beg = k;
				while (lo_i < hi_i) {
					int g_i = (lo_i + hi_i) >> 1, xi = a[kp1][g_i];

					if (z[k] < x[xi][k]) {
						hi_i = g_i;
						td_beg = k;
					} else if (z[k] > x[xi][k]) {
						lo_i = g_i + 1;
						tu_beg = k;
					} else {
						// find match start
						int lo_k;
						if (k - no_beg < 10) {
							lo_k = k;
							while (lo_k-- > no_beg && z[lo_k] == x[xi][lo_k]) {}
						} else {
							lo_k = no_beg;
							int hi_k = k;
							unsigned long long rt = hz[kp1] < h[xi][kp1] ? MOD - h[xi][kp1] + hz[kp1] : hz[kp1] - h[xi][kp1];
							while (lo_k < hi_k) {
								int g_k = (lo_k + hi_k) >> 1;

								if ((hz[g_k] < h[xi][g_k] ? MOD - h[xi][g_k] + hz[g_k] : hz[g_k] - h[xi][g_k]) == rt) {
									hi_k = g_k;
								} else {
									lo_k = g_k + 1;
								}
							}
							lo_k--;
						}

						if (lo_k > -1 && z[lo_k] < x[xi][lo_k]) {
							hi_i = g_i;
							td_beg = lo_k;
						} else {
							lo_i = g_i + 1;
							tu_beg = lo_k;
						}
					}
				}

				t = lo_i;
				no_beg = max(no_beg, min(tu_beg, td_beg) + 1);
			}

			/* increment req_idx until matching over syllables [req_idx, k+1] is not enough.
			   the last site can be omitted from syllable k+1 here since a full match of syllable k+1 means the
			   match will have a chance to be detected at the next syllable, unless this is the last syllable */
			U loc = locs[kp1 == n - 1 ? N - 1 : (kp1 + 1) * B - 2];
			while (req_idx < k && loc - locs[req_idx * B] >= L) req_idx++;

			// subtract terminated matches from ongoing match counts
			up_on -= up_end[k];
			dn_on -= dn_end[k];

			if (req_idx < no_beg) continue; // no possible matches

			unsigned long long zh = hz[kp1] < hz[req_idx] ? MOD - hz[req_idx] + hz[kp1] : hz[kp1] - hz[req_idx];

			// binary search for block of long matches above z
			int p = t - up_on, lo_p = 0, hi_p = p;
			while (lo_p < hi_p) {
				int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
				if ((h[xi][kp1] < h[xi][req_idx] ? MOD - h[xi][req_idx] + h[xi][kp1] : h[xi][kp1] - h[xi][req_idx]) == zh) {
					hi_p = g_p;
				} else {
					lo_p = g_p + 1;
				}
			}
			up_on = t - lo_p; // add new matches in block to ongoing match count

			for (int i = lo_p, s_idx = req_idx; i < p; i++) { // iterate through new matches above
				int xi = a[kp1][i];
				// find match start; total # iterations bounded by n
				while (s_idx > -1 && z[s_idx] == x[xi][s_idx]) s_idx--;

				// briefly linear search for end of match
				int lo = kp1, pivot = min(n, lo + 10);
				while (lo < pivot && z[lo] == x[xi][lo]) lo++;
				if (lo == pivot) { // switch to binary search if end not found
					int hi = n;
					unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

					while (lo < hi) {
						int g = (lo + hi + 1) >> 1;

						if (lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g])) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				}
				up_end[lo]++;

				// refine site boundaries [start, end)
				int start, end;
				refine(s_idx, lo, xi, &start, &end);
				// report match if long
				U len = locs[end - 1] - locs[start] + inclusive;
				if (len >= L) {
					out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\t' << len << '\n';
				}
			}

			// analogous process for below query

			p = t + dn_on, lo_p = p, hi_p = M;
			while (lo_p < hi_p) {
				int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
				if ((h[xi][kp1] < h[xi][req_idx] ? MOD - h[xi][req_idx] + h[xi][kp1] : h[xi][kp1] - h[xi][req_idx]) == zh) {
					lo_p = g_p + 1;
				} else {
					hi_p = g_p;
				}
			}
			dn_on = lo_p - t;

			for (int i = lo_p - 1, s_idx = req_idx; i >= p; i--) {
				int xi = a[kp1][i];
				while (s_idx > -1 && z[s_idx] == x[xi][s_idx]) s_idx--;

				int lo = kp1, pivot = min(n, lo + 10);
				while (lo < pivot && z[lo] == x[xi][lo]) lo++;
				if (lo == pivot) {
					int hi = n;
					unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

					while (lo < hi) {
						int g = (lo + hi + 1) >> 1;

						if (lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g])) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				} 
				dn_end[lo]++;

				int start, end;
				refine(s_idx, lo, xi, &start, &end);
				U len = locs[end - 1] - locs[start] + inclusive;
				if (len >= L) {
					out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\t' << len << '\n';
				}
			}
		}
	}

	return 0;
}

template<class T>
int SyllableQuery<T>::query(const char* query_file, const char* output_file, const int L) { // 1-length sites
	if (L < minSiteL) return 4;
	ifstream in(query_file);
	if (in.fail()) return 1;

	const int l = (L - B + 1) / B; // min # syllables that must be covered by an L-site match
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
	vector<vector<T>> Z(n, vector<T>(Q));

	for (int K = 0; K < N; K++) {
		if (!getline(in, line)) return 5;
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int k = K / B, i = 0;
		while (getline(ss, line, '\t')) {
			if (i == Q || line.size() < 3u) return 5;
			Z[k][i] = (Z[k][i] << 1) | (line[0] != '0'), i++;
			Z[k][i] = (Z[k][i] << 1) | (line[2] != '0'), i++;
		}
	}
	in.close();
	for (int q = 0; q < Q; q++) Z[n - 1][q] <<= n * B - N;

	ofstream out(output_file);
	if (out.fail()) return 2;
	for (int q = 0; q < Q; q++) {
		string qID = qIDs[q >> 1] + "-" + to_string(q & 1);

		for (int k = 0; k < n; k++) {
			z_[k] = Z[k][q];
			z[k] = lower_bound(r[k].begin(), r[k].end(), z_[k]) - r[k].begin();
			if (z[k] < (int) r[k].size() && z_[k] != r[k][z[k]]) z[k] = r[k].size();
			hz[k+1] = (hz[k] + xp[k] * (z[k] + 1)) % MOD;
		}
		memset(&up_end[1], 0, n * sizeof(int));
		memset(&dn_end[1], 0, n * sizeof(int));

		for (int k = l-1, kp1 = l, t = 0, up_on = 0, dn_on = 0, no_beg = 0; k < n; k = kp1-1) {
			if (z[k] == (int) r[k].size()) {
				t = M;
				no_beg = k + 1;
			} else {
				int lo_i = 0, hi_i = M, tu_beg = k, td_beg = k;
				while (lo_i < hi_i) {
					int g_i = (lo_i + hi_i) >> 1, xi = a[kp1][g_i];

					if (z[k] < x[xi][k]) {
						hi_i = g_i;
						td_beg = k;
					} else if (z[k] > x[xi][k]) {
						lo_i = g_i + 1;
						tu_beg = k;
					} else {
						int lo_k;
						if (k - no_beg < 10) {
							lo_k = k;
							while (lo_k-- > no_beg && z[lo_k] == x[xi][lo_k]) {}
						} else {
							lo_k = no_beg;
							int hi_k = k;
							unsigned long long rt = hz[kp1] < h[xi][kp1] ? MOD - h[xi][kp1] + hz[kp1] : hz[kp1] - h[xi][kp1];
							while (lo_k < hi_k) {
								int g_k = (lo_k + hi_k) >> 1;

								if ((hz[g_k] < h[xi][g_k] ? MOD - h[xi][g_k] + hz[g_k] : hz[g_k] - h[xi][g_k]) == rt) {
									hi_k = g_k;
								} else {
									lo_k = g_k + 1;
								}
							}
							lo_k--;
						}

						if (lo_k > -1 && z[lo_k] < x[xi][lo_k]) {
							hi_i = g_i;
							td_beg = lo_k;
						} else {
							lo_i = g_i + 1;
							tu_beg = lo_k;
						}
					}
				}

				t = lo_i;
				no_beg = max(no_beg, min(tu_beg, td_beg) + 1);
			}

			up_on -= up_end[k];
			dn_on -= dn_end[k];
			if (kp1 < no_beg + l) {
				kp1 = no_beg + l;
				continue;
			}
			int s_idx = k - l;

			unsigned long long zh = hz[kp1] < hz[kp1 - l] ? MOD - hz[kp1 - l] + hz[kp1] : hz[kp1] - hz[kp1 - l];

			int p = t - up_on, lo_p = 0, hi_p = p;
			while (lo_p < hi_p) {
				int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
				if ((h[xi][kp1] < h[xi][kp1 - l] ? MOD - h[xi][kp1 - l] + h[xi][kp1] : h[xi][kp1] - h[xi][kp1 - l]) == zh) {
					hi_p = g_p;
				} else {
					lo_p = g_p + 1;
				}
			}
				
			while (p > lo_p) {
				int xi = a[kp1][--p];

				int lo = kp1, pivot = min(n, lo + 10);
				while (lo < pivot && z[lo] == x[xi][lo]) lo++;
				if (lo == pivot) {
					int hi = n;
					unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

					while (lo < hi) {
						int g = (lo + hi + 1) >> 1;

						if (lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g])) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				}
				up_end[lo]++;

				int start, end;
				refine(s_idx, lo, xi, &start, &end);
				if (end - start >= L) {
					out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\n';
				}
			}
			up_on = t - p;

			p = t + dn_on, lo_p = p, hi_p = M;
			while (lo_p < hi_p) {
				int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
				if ((h[xi][kp1] < h[xi][kp1 - l] ? MOD - h[xi][kp1 - l] + h[xi][kp1] : h[xi][kp1] - h[xi][kp1 - l]) == zh) {
					lo_p = g_p + 1;
				} else {
					hi_p = g_p;
				}
			}

			while (p < lo_p) {
				int xi = a[kp1][p++];

				int lo = kp1, pivot = min(n, lo + 10);
				while (lo < pivot && z[lo] == x[xi][lo]) lo++;
				if (lo == pivot) {
					int hi = n;
					unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

					while (lo < hi) {
						int g = (lo + hi + 1) >> 1;

						if (lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g])) {
							lo = g;
						} else {
							hi = g - 1;
						}
					}
				} 
				dn_end[lo]++;

				int start, end;
				refine(s_idx, lo, xi, &start, &end);
				if (end - start >= L) {
					out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\n';
				}
			}
			dn_on = p - t;

			kp1++;
		}
	}

	return 0;
}
