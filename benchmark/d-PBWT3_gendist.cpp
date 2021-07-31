#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
	ifstream in(argv[1]), in_map(argv[2]);

        string line;
        while (getline(in, line)) {
                if (line[0] != '#' || line[1] != '#') break;
        }
        stringstream ss(line);
        int M = -9, N = 0;
        while (getline(ss, line, '\t')) M++;
        M <<= 1;
	while (getline(in, line)) N++;
	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	vector<string> IDs(M);
	for (int i = 0; i < M; i += 2) {
		getline(ss, IDs[i], '\t');
		IDs[i + 1] = IDs[i] + "-1";
		IDs[i] += "-0";
	}

	vector<vector<bool>> x(N, vector<bool>(M));
	vector<vector<int>> a(N, vector<int>(M)), d(N, vector<int>(M)), u(N, vector<int>(M + 1)), v(N, vector<int>(M + 1));
	vector<int> a0(M), a1(M), d0(M), d1(M);
	vector<double> locs(N);

	for (int k = 0; k < N; k++) {
		getline(in, line);
                ss = stringstream(line);
                for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int index = 0;
		while (getline(ss, line, '\t')) {
			x[k][index++] = line[0] != '0';
			x[k][index++] = line[2] != '0';
		}

                int u_ = 0, v_ = 0, p = k + 1, q = k + 1;
                for (int i = 0; i < M; i++) {
			int d_ = k > 0 ? d[k - 1][i] : 0;
			int a_ = k > 0 ? a[k - 1][i] : i;
                        p = max(p, d_);
                        q = max(q, d_);
			u[k][i] = u_;
			v[k][i] = v_;

                        if (x[k][a_]) {
                                a1[v_] = a_;
                                d1[v_] = q;
                                v_++;
                                q = 0;
                        } else {
                                a0[u_] = a_;
                                d0[u_] = p;
                                u_++;
                                p = 0;
                        }
                }
		u[k][M] = u_;
		v[k][M] = M;

                for (int i = 0; i < M; i++) {
			v[k][i] += u_;
                        if (i < u_) {
                                a[k][i] = a0[i];
                                d[k][i] = d0[i];
                        } else {
                                a[k][i] = a1[i - u_];
                                d[k][i] = d1[i - u_];
                        }
                }

		getline(in_map, line);
		ss = stringstream(line);
		while (getline(ss, line, '\t')) {}
		try {
			locs[k] = stod(line);
		} catch (exception& e) {}
        }
	in.close();
	in_map.close();

	vector<int> t(N), zd(N), bd(N), dz(M);

	clock_t START = clock();
	ifstream q_in(argv[3]);
	while (getline(q_in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	int Q = -9;
	ss = stringstream(line);
	while (getline(ss, line, '\t')) Q++;
	q_in.clear(), q_in.seekg(0);
	while (getline(q_in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	vector<string> qIDs(Q);
	ss = stringstream(line);
	for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
	for (int i = 0; i < Q; i++) getline(ss, qIDs[i], '\t');
	Q <<= 1;
	vector<vector<bool>> z(N, vector<bool>(Q));
	for (int k = 0; k < N; k++) {
		getline(q_in, line);
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) getline(ss, line, '\t');
		int i = 0;
		while (getline(ss, line, '\t')) {
			z[k][i++] = line[0] != '0';
			z[k][i++] = line[2] != '0';
		}
	}
	q_in.close();

	double L = stod(argv[4]);
	ofstream out(argv[5]);
	for (int q = 0; q < Q; q++) {
		string qID = qIDs[q >> 1] + "-" +  to_string(q & 1);

		t[0] = z[0][q] ? M : 0;
		for (int i = 1; i < N; i++) {
			if (z[i][q]) {
				t[i] = v[i][t[i - 1]];
			} else {
				t[i] = u[i][t[i - 1]];
			}
		}

		int z_idx = N, b_idx = N;
		for (int i = N - 1; i > -1; i--) {
			z_idx = min(z_idx, i + 1);
			b_idx = min(b_idx, i + 1);
			if (t[i] > 0) {
				while (z_idx > 0 && z[z_idx - 1][q] == x[z_idx - 1][a[i][t[i] - 1]]) z_idx--;
				zd[i] = z_idx;
			} else zd[i] = i + 1;
			if (t[i] < M) {
				while (b_idx > 0 && z[b_idx - 1][q] == x[b_idx - 1][a[i][t[i]]]) b_idx--;
				bd[i] = b_idx;
			} else bd[i] = i + 1;
		}

		int f = t[0], g = t[0], f_end, g_end;
		for (int i = 1, req_idx = 0, up_idx = 0, dn_idx = 0; i < N; i++) {
			if (z[i][q]) {
				f_end = u[i][f];
				g_end = u[i][g];
				f = v[i][f];
				g = v[i][g];
			} else {
				f_end = v[i][f];
				g_end = v[i][g];
				f = u[i][f];
				g = u[i][g];
			}

			while (req_idx < i && locs[i] - locs[req_idx] >= L) req_idx++;
			if (req_idx == 0) continue;

			while (f_end < g_end) {
				out << qID << '\t' << IDs[a[i][f_end]] << '\t' << dz[a[i][f_end]] << '\t' << i << '\t' << locs[i - 1] - locs[dz[a[i][f_end]]] << '\n';
				f_end++;
			}

			if (f == g) {
				if (f > 0 && zd[i] < req_idx) {
					up_idx = zd[i];
					dz[a[i][--f]] = up_idx;
				}
				if (g < M && bd[i] < req_idx) {
					dn_idx = bd[i];
					dz[a[i][g++]] = dn_idx;
				}
			}

			if (f < g) {
				while (f > 0 && d[i][f] < req_idx) {
					up_idx = max(up_idx, d[i][f]);
					dz[a[i][--f]] = up_idx;
				}
				while (g < M && d[i][g] < req_idx) {
					dn_idx = max(dn_idx, d[i][g]);
					dz[a[i][g++]] = dn_idx;
				}
			}
		}

		while (f < g) {
			out << qID << '\t' << IDs[a[N - 1][f]] << '\t' << dz[a[N - 1][f]] << '\t' << N << '\t' << locs[N - 1] - locs[dz[a[N - 1][f]]] << '\n';
			f++;
		}
	}

	cout << "d-PBWT3 " << M << ' ' << N << ' ' << L << '\n';
	cout << clock() - START << " time\n";

	return 0;
}
