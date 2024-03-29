all: main

main: server.cpp client.cpp SyllableQuery.cpp SyllableQuery.h
	g++ -std=c++17 -Wshadow -Wall -o server server.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o client client.cpp -O2 -Wno-unused-result

clean:
	rm -f server client

.PHONY: all clean
