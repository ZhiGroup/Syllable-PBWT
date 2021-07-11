#pragma once

#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <chrono>
#include <unistd.h>

using namespace std;

template<class T>
class Reporter {
public:
	const unsigned long long MOD = (1ull << 63) - 25;
	int M{0}, N{0}, n{0}, B{sizeof(T) * 8};
	int minSiteL{B * 2 - 1}, minPhysL{0};
	double minGenL{-1};
	vector<double> genLocs;
	vector<int> physLocs;
	vector<string> sampleIDs;
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
