#pragma once
#pragma GCC diagnostic ignored "-Wshift-count-overflow"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unistd.h>
#include <vector>

using namespace std;

template<class T>
struct SparseQuery {
	static const unsigned long long MOD{-1ull - 58};
	const int B{sizeof(T) * 8};
	int M{0}, N{0}, n{0};
	int minSiteL{B * 2 - 1}, minPhysL{0};
	double minGeneL{-1};
	vector<double> geneLocs;
	vector<int> physLocs;
	vector<string> IDs;
	vector<vector<T>> x;
	vector<vector<unsigned long long>> h;
	vector<vector<int>> a, d, inv_a;
	vector<unsigned long long> xp, hz;
	vector<int> up_end, dn_end;

	int precompute(const char* input_file);
	int save(const char* save_file);
	int load(const char* load_file);
	void show_attributes();
	int set_gen_map(const char* map_file);
	template<class U>
	int query(const char* query_file, const char* output_dir, const U L, const vector<U> &locs);
	int query(const char* query_file, const char* output_dir, const int L);
};
