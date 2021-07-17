all: main

main: main.cpp SparseQuery.cpp SparseQuery.h SparseTable.h
	g++ -std=c++17 -Wshadow -Wall -o main main.cpp -O2 -Wno-unused-result

clean:
	rm -f main SparseQuery

.PHONY: all clean
