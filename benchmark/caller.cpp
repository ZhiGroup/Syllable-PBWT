#include <bits/stdc++.h>
#include "../SyllableQuery.cpp"

using namespace std;

typedef unsigned __int128 T;
// typedef unsigned long long T;

int main(int argc, char **argv) {
	SyllableQuery<T> sq;
	string in_file = argv[1], si_load = argv[2], query_file = argv[3], out_file = argv[4], unit = argv[5];
	cout << in_file << endl;
	if (si_load == "load") {
		sq.load(in_file.c_str());
	} else {
		sq.precompute(in_file.c_str());
		// sq.save(si_load.c_str());
	}

	if (unit == "cM") {
		double L = stod(argv[6]);
		if (si_load != "load") {
			string gen_file = argv[7];
			sq.set_gen_map(gen_file.c_str());
		}
		clock_t start = clock();
		sq.query(query_file.c_str(), out_file.c_str(), L, sq.geneLocs, false);
		cout << clock() - start << " us\n";
	} else {
		int L = stoi(argv[6]);
		clock_t start = clock();
		sq.query(query_file.c_str(), out_file.c_str(), L);
		cout << clock() - start << " us\n";
	}

	return 0;
}
