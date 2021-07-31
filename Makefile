all: main

main: server.cpp client.cpp StrideQuery.cpp StrideQuery.h SparseTable.h
	g++ -std=c++17 -Wshadow -Wall -o server server.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o client client.cpp -O2 -Wno-unused-result

clean:
	rm -f server client StrideQuery

.PHONY: all clean
