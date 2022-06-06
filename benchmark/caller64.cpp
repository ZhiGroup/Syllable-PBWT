#include <time.h>
#include "../SyllableQuery.cpp"

using namespace std;

typedef unsigned long long T;

int main(int argc, char **argv) {
	SyllableQuery<T> sq;
	string panel_file = argv[1], query_file = argv[2], out_file = argv[3];
	assert(sq.precompute(panel_file.c_str()) == 0);
	int L = stoi(argv[4]);
	const clock_t START = clock();
	assert(sq.query(query_file.c_str(), out_file.c_str(), L) == 0);
	cout << "in: " << panel_file << "\nL: " << L << "\ntime: " << double(clock() - START) / CLOCKS_PER_SEC << '\n';

	return 0;
}
	
	
