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

	vector<int> rand_ord(M);
	iota(rand_ord.begin(), rand_ord.end(), 0);
	mt19937 g(17);
	shuffle(rand_ord.begin(), rand_ord.end(), g);
	
	set<int> queries;
	for (int i = 0; i < 50; i++) {
		queries.insert(rand_ord[i]);
	}

	in.clear(), in.seekg(0);
	while (getline(in, line)) {
		if (line[0] != '#' || line[1] != '#') break;
	}

	vector<ofstream> panel_outs(6), query_outs(6);
	for (int i = 0; i < 6; i++) {
		panel_outs[i] = ofstream(string("panel-") + to_string(i) + string(".vcf"));
		query_outs[i] = ofstream(string("query-") + to_string(i) + string(".vcf"));
	}

	int N = 0;
	while (true) {
		ss = stringstream(line);

		for (int j = N % 6; j < 6; j++) {
			for (int _ = 0; _ < 9; _++) {
				getline(ss, line, '\t');
				if (_ > 0) {
					panel_outs[j] << '\t';
					query_outs[j] << '\t';
				}
				panel_outs[j] << line;
				query_outs[j] << line;
			}

			for (int i = 0; i < M; i++) {
				getline(ss, line, '\t');
				if (queries.count(i)) {
					query_outs[j] << '\t' << line;
				} else {
					panel_outs[j] << '\t' << line;
				}
			}
			query_outs[j] << '\n';
			panel_outs[j] << '\n';
		}

		if (!getline(in, line)) break;
		N++;
	}

	return 0;
}
