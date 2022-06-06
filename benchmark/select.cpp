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
	ofstream out_panel(argv[2] + string("_panel.vcf")), out_query(argv[2] + string("_query.vcf"));
	while (true) {
		ss = stringstream(line);
		for (int _ = 0; _ < 9; _++) {
			getline(ss, line, '\t');
			if (_ > 0) {
				out_panel << '\t';
				out_query << '\t';
			}
			out_panel << line;
			out_query << line;
		}

		for (int i = 0; i < M; i++) {
			getline(ss, line, '\t');
			if (queries.count(i)) {
				out_query << '\t' << line;
			} else {
				out_panel << '\t' << line;
			}
		}
		out_query << '\n';
		out_panel << '\n';

		if (!getline(in, line)) break;
	}

	return 0;
}
