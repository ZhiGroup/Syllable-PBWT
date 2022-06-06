all: d-PBWT3 caller

d-PBWT3: d-PBWT3.cpp d-PBWT3_gendist.cpp
	g++ -std=c++17 -Wshadow -Wall -o d-PBWT3 d-PBWT3.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o d-PBWT3_gendist d-PBWT3_gendist.cpp -O2 -Wno-unused-result

test: select.cpp stride.cpp prefix.cpp
	g++ -std=c++17 -Wshadow -Wall -o select select.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o stride stride.cpp -O2 -Wno-unused-result
	g++ -std=c++17 -Wshadow -Wall -o prefix prefix.cpp -O2 -Wno-unused-result

clean:
	rm -f d-PBWT3 d-PBWT3_gendist select stride prefix

.PHONY: all clean
