all: fullmem merge merge2

fullmem: fullmem.cpp
	g++ -std=c++17 -Wshadow -Wall -o fullmem fullmem.cpp -O2 -Wno-unused-result

merge: merge.cpp
	g++ -std=c++17 -Wshadow -Wall -o merge merge.cpp -O2 -Wno-unused-result

merge2: merge2.cpp
	g++ -std=c++17 -Wshadow -Wall -o merge2 merge2.cpp -O2 -Wno-unused-result

clean:
	rm -f fullmem merge merge2

.PHONY: all clean
