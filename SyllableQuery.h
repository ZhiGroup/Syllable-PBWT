// This code is licensed under MIT license (see LICENSE for details)
#pragma once

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
struct SyllableQuery {
	static const unsigned long long MOD{-1ull - 58}, BASE{(1ull<<48) - 59};
	const int B{sizeof(T) * 8};
	int M{0}, N{0}, n{0};
	int minSiteL{B * 2 - 1}, minPhysL{0};
	double minGeneL{-1};
	vector<double> geneLocs;
	vector<int> physLocs;
	vector<string> IDs;
	vector<vector<int>> x, a;
	vector<vector<T>> r;
	vector<vector<unsigned long long>> h;
	vector<T> z_;
	vector<unsigned __int128> xp;
	vector<unsigned long long> hz;
	vector<int> z, up_end, dn_end;

	int precompute(const char* input_file);
	int save(const char* save_file);
	int load(const char* load_file);
	string show_attributes();
	int set_gen_map(const char* map_file);
	inline void refine(int s_idx, int e_idx, int xi, int *start, int *end);
	template<class U>
	int query(const char* query_file, const char* output_file, const U L, const vector<U> &locs, bool inclusive);
	int query(const char* query_file, const char* output_file, const int L);
};
