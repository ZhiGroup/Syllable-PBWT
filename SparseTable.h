#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

using namespace std;

struct SparseTable {
	static const int B = 30;
	int N, blocks;
	vector<int> v, mask;
	vector<vector<int>> table;

	SparseTable(int _N) {
		N = _N;
		blocks = N / B;
		mask = vector<int>(N);
		table = vector<vector<int>>(blocks, vector<int>(msb(blocks) + 1));
	}

	void build(vector<int>& _v) {
		v = _v;
		assert(N == (int) v.size());

		int cur = 0; // sliding mask
		for (int i = 0; i < N; ++i) {
			cur = (cur << 1) & ((1 << B) - 1);
			while (cur > 0 && max(v[i], v[i - msb(lsb(cur))]) == v[i]) cur ^= lsb(cur);
			cur |= 1;
			mask[i] = cur;
		}

		for (int i = 0; i < blocks; ++i) table[i][0] = max_query(B * i + B - 1);
		for (int j = 1; (1 << j) <= blocks; ++j) {
			for (int i = 0; i + (1 << j) - 1 < blocks; ++i) {
				table[i][j] = max(table[i][j - 1], table[i + (1 << (j - 1))][j - 1]);
			}
		}
	}

	// least significant set bit
	int lsb(int num) {return num & -num;}

	// index of most significant set bit
	int msb(int num) {return __builtin_clz(1) - __builtin_clz(num);}

	int max_query(int r, int len = B) {
		return v[r - msb(mask[r] & ((1 << len) - 1))];
	}

	int query(int l, int r) {
		if (r - l + 1 <= B) return max_query(r, r - l + 1);
		int ret = max(max_query(l + B - 1), max_query(r));
		int blockL = l / B + 1, blockR = r / B - 1;
		if (blockL <= blockR) {
			int j = msb(blockR - blockL + 1);
			ret = max({ret, table[blockL][j], table[blockR - (1 << j) + 1][j]});
		}
		return ret;
	}
};
