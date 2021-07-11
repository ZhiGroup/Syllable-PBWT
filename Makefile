all: main fullmem

main: main.cpp
	g++ -std=c++17 -Wshadow -Wall -o main main.cpp -O2 -Wno-unused-result

fullmem: fullmem.cpp
	g++ -std=c++17 -Wshadow -Wall -o fullmem fullmem.cpp -O2 -Wno-unused-result

clean:
	rm -f main Reporter fullmem

.PHONY: all clean
