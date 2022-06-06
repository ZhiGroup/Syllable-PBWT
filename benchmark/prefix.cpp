#include <bits/stdc++.h>

using namespace std;

int main(int argc, char **argv) {
	ifstream in(argv[1]);

	string line;
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}
	stringstream ss(line);
	int M = -9;
	while (getline(ss, line, '\t')) M++;

	vector<int> Ms(8);
	for (int i = 0; i < 8; i++) {
		Ms[i] = M * (i + 1) / 8;
	}

	vector<set<int>> queries(8);
	for (int i = 0; i < 8; i++) {
		vector<int> rand_ord(Ms[i]);
		iota(rand_ord.begin(), rand_ord.end(), 0);
		mt19937 g(17);
		shuffle(rand_ord.begin(), rand_ord.end(), g);

		for (int j = 0; j < 50; j++) {
			queries[i].insert(rand_ord[j]);
		}
	}

	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}

	vector<ofstream> panel_outs(8), query_outs(8);
	for (int i = 0; i < 8; i++) {
		panel_outs[i] = ofstream(string("panelM-") + to_string(i) + string(".vcf"));
		query_outs[i] = ofstream(string("queryM-") + to_string(i) + string(".vcf"));
	}

	while (true) {
		ss = stringstream(line);

		for (int _ = 0; _ < 9; _++) {
			getline(ss, line, '\t');
			for (int i = 0; i < 8; i++) {
				if (_ > 0) {
					panel_outs[i] << '\t';
					query_outs[i] << '\t';
				}
				panel_outs[i] << line;
				query_outs[i] << line;
			}
		}

		for (int i = 0; i < M; i++) {
			getline(ss, line, '\t');
			for (int j = 0; j < 8; j++) {
				if (i < Ms[j]) {
					if (queries[j].count(i)) {
						query_outs[j] << '\t' << line;
					} else {
						panel_outs[j] << '\t' << line;
					}
				}
			}
		}
		for (int j = 0; j < 8; j++) {
			query_outs[j] << '\n';
			panel_outs[j] << '\n';
		}

		if (!getline(in, line)) break;
	}

	return 0;
}
