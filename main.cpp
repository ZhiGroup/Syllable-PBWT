#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>

#include "Reporter.cpp"

using namespace std;

bool input_file_set, save_file_set, load_file_set;
char *input_file, *save_file, *load_file;
int B{128};

void show_usage(string name) {
	cout << "\nSpace-efficient long-match query algorithm.\n" <<
		"Usage: " << name << " <flags>\n" <<
		"Flags for precomputation stage:\n" <<
		"\t-h,--help                            Show this help message\n" <<
		"\t-i,--input_panel <FILE>              Path to VCF input file\n" <<
		"\t-b,--bits <VALUE>                    Value of B (64 or 128; default: 128), the number of individual sites to be combined into one merged site\n" <<
		"\t                                     A higher B leads to stricter minimum query lengths but improved space efficiency.\n" <<
		"\t-s,--save <FILE>                     Path to output file to save panel data, for fast loading when rerunning this program with --load\n" <<
		"\t-l,--load <FILE>                     Path to input file to load panel data (must have been the save file of a previous successful run)\n" <<
		"Example: " << name << " -i sample.vcf -s save/sample.txt\n" <<
		"Usage note: The --save flag is optional. Exactly one of the flags --input_panel and --load should be used." << endl;
}

void show_query_usage() {
	cout << "\nYou may enter queries using the following flags:\n" <<
		"\t-q,--query_panel <FILE>              Path to VCF query file\n" <<
		"\t-l,--length <VALUE>                  Minimum length of desired matches\n" <<
		"\t-u,--unit <STRING>                   Unit of query length: cM, sites, or bps\n" <<
		"\t-o,--output <DIRECTORY>              Directory to output matches (default: directory of query panel)\n" <<
		"Example: -q query.vcf -l 700 -u sites -o results\n" <<
		"Below are other actions you can perform:\n" <<
		"\t-h,--help                            Show this help message\n" <<
		"\t-a,--attr                            Show data attributes of panel\n" <<
		"\t-g,--gen_map <FILE>                  Path to genetic map file, to be used when query length unit is cM (overrides previous map, if present)\n" <<
		"\t                                     Format: Sites described by one line each with genetic location as the last tab-delimited field\n" <<
		"\t-s,--save <FILE>                     Path to output file to save panel data, for fast loading when rerunning this program with --load\n" <<
		"\t-E,--exit                            Exit the program\n" <<
		"Format of reported matches: Haplotypes will be referred to with its sample ID in the VCF file, followed by -0 or -1, indicating the first or second haplotype, respectively, of a given individual. For each query haplotype, a file titled the haplotype's name will be written to, with one match reported per line, consisting of at least 3 tab-delimited fields $1 $2 $3, indicating a match over sites [$1, $2) with panel haplotype $3. If the query unit is cM or bps, there will be a fourth field indicating the length of the match in terms of the unit specified." << endl;
}

string getDir(string &s) {
	size_t i = s.rfind('/', s.length());
	if (i != string::npos) return s.substr(0, i);
	return "";
}

template<class T>
void listen() {
	Reporter<T> rep;

	if (input_file_set) {
		cout << "Precomputing on " << input_file << endl;
		int code = rep.precompute(input_file);
		if (code == 1) {
			cerr << "Input file does not exist or cannot be read." << endl;
			return;
		}
		if (code == 2) {
			cerr << "Input file not in VCF format." << endl;
			return;
		}
		cout << "Precomputation complete" << endl;
	} else {
		cout << "Loading data from " << load_file << endl;
		int code = rep.load(load_file);
		if (code == 1) {
			cerr << "Load file does not exist or cannot be read." << endl;
			return;
		}
		if (code == 2) {
			cerr << "Load file is incomplete." << endl;
			return;
		}
		if (code == 3) {
			cerr << "The specified --bits value, " << B << ", does not coincide with the one in the save file, " << rep.B << "." << endl;
			return;
		}
		cout << "Loading data complete" << endl;
	}

	rep.show_attributes();

	if (save_file_set) {
		cout << "Saving data to " << save_file << endl;
		int code = rep.save(save_file);
		if (code == 1) {
			cerr << "Save file cannot be written to." << endl;
		} else {
			cout << "Saving data complete" << endl;
		}
	}

	show_query_usage();

	while (true) {
		cout << "\nEnter a query or action: " << flush;
		string line, flag, val, query_file, unit, output_dir;
		double L = 0;
		bool skip = false;
		getline(cin, line);
		stringstream ss(line);

		while (getline(ss, flag, ' ')) {
			if (flag == "-h" || flag == "--help") {
				show_query_usage();
				skip = true; break;
			}
			if (flag == "-a" || flag == "--attr") {
				rep.show_attributes();
				skip = true; break;
			}
			if (flag == "-E" || flag == "--exit") {
				return;
			}

			if (!getline(ss, val, ' ')) {
				cerr << "No value specified for your last flag." << endl;
				skip = true; break;
			}

			if (flag == "-g" || flag == "--gen_map") {
				cout << "Loading genetic map from " << val << endl;
				int code = rep.set_gen_map(val.c_str());
				if (code == 1) {
					cerr << "Genetic map file does not exist or cannot be read." << endl;
				} else if (code == 2) {
					cerr << "Genetic map file does not have one line per site." << endl;
				} else if (code == 3) {
					cerr << "Genetic map file does not contain a number in the last tab-delimited field of every line." << endl;
				} else {
					cout << "Genetic map loaded." << endl;
					rep.show_attributes();
				}
				skip = true; break;
			}
			if (flag == "-s" || flag == "--save") {
				cout << "Saving data to " << val << endl;
				int code = rep.save(val.c_str());
				if (code == 1) {
					cerr << "Save file cannot be written to." << endl;
				} else {
					cout << "Saving data complete" << endl;
				}
				skip = true; break;
			}

			if (flag == "-q" || flag == "--query_panel") {
				query_file = val;
			} else if (flag == "-l" || flag == "--length") {
				try { L = stod(val); }
				catch (exception &e) {
					cerr << "Length must be a number." << endl;
					skip = true; break;
				}
			} else if (flag == "-u" || flag == "--unit") {
				unit = val;
			} else if (flag == "-o" || flag == "--output") {
				output_dir = val;
			} else {
				cerr << "Invalid flag: " << flag << endl;
				show_query_usage();
				skip = true; break;
			}
		}

		if (skip) continue;
		if (query_file.size() == 0u) {
			cerr << "Query file must be specified." << endl;
			continue;
		}
		if (L == 0) {
			cerr << "Match length must be specified." << endl;
			continue;
		}
		if (unit != "cM" && unit != "sites" && unit != "bps") {
			cerr << "Unit must be cM, sites, or bps." << endl;
			continue;
		}
		if (output_dir.size() == 0u) output_dir = getDir(query_file);

		cout << "Performing query on " << query_file << endl;
		const clock_t START = clock();
		int code;
		if (unit == "sites") {
			code = rep.query(query_file.c_str(), output_dir.c_str(), (int) L);
		} else if (unit == "cM") {
			code = rep.query(query_file.c_str(), output_dir.c_str(), L, rep.genLocs);
		} else { // "bps"
			code = rep.query(query_file.c_str(), output_dir.c_str(), (int) L, rep.physLocs);
		}

		if (code == 1) {
			cerr << "Query file does not exist or cannot be read." << endl;
		} else if (code == 2) {
			cerr << "Output directory " << output_dir << " does not exist or cannot be written to." << endl;
		} else if (code == 3) {
			cerr << "A genetic map must have been specified to query with genetic distance." << endl;
		} else if (code == 4) {
			cerr << "The query length, " << L << " " << unit << ", is less than the minimum allowable query length," << (unit == "sites" ? rep.minSiteL : unit == "cM" ? rep.minGenL : rep.minPhysL) << " " << unit << "." << endl;
		} else if (code == 5) {
			cerr << "Query file not in VCF format." << endl;
		} else {
			cout << "Query results outputted to " << output_dir << "\nElapsed CPU time (s): " << double(clock() - START) / CLOCKS_PER_SEC << endl;
		}
	}
}

int main(int argc, char** argv) {
	for (int i = 1; i < argc; i++) {
		string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			show_usage(argv[0]);
			return 0;
		}
		if (++i == argc) {
			cerr << "No value specified for your last flag." << endl;
			return 0;
		}

		if (arg == "-i" || arg == "-input_panel") {
			input_file = argv[i];
			input_file_set = true;
		} else if (arg == "-b" || arg == "--bits") {
                        try { B = stoi(argv[i]); }
			catch (exception &e) {
				cerr << "--bits must be set to 64 or 128 (128 by default)." << endl;
				return 0;
			} 
		} else if (arg == "-s" || arg == "--save") {
			save_file = argv[i];
			save_file_set = true;
		} else if (arg == "-l" || arg == "--load") {
			load_file = argv[i];
			load_file_set = true;
		} else {
			cerr << "Invalid flag: " << arg << endl;
			show_usage(argv[0]);
			return 0;
		}
	}

	if (!(input_file_set ^ load_file_set)) {
		cerr << "Out of the flags --input_panel and --load, exactly one must be specified." << endl;
		return 0;
	}
	if (B != 64 && B != 128) {
		cerr << "--bits must be set to 64 or 128 (128 by default)." << endl;
		return 0;
	}

	if (B == 64) listen<unsigned long long>();
	else listen<unsigned __int128>();

	return 0;
}
