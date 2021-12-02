all: d-PBWT3 caller

d-PBWT3: d-PBWT3.cpp d-PBWT3_gendist.cpp
	g++ -std=c++17 -Wshadow -Wall -o d-PBWT3 d-PBWT3.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o d-PBWT3_gendist d-PBWT3_gendist.cpp -O2 -Wno-unused-result

caller: caller.cpp
	g++ -std=c++17 -Wshadow -Wall -o caller caller.cpp -O2 -Wno-unused-result

clean:
	rm -f d-PBWT3 d-PBWT3_gendist

.PHONY: all clean
