#include <iostream>
#include <fstream>
#include <math.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/stat.h>
#include <time.h>

#include "SyllableQuery.cpp"

using namespace std;

long long save_ct = 0, out_ct = 0;
bool fifo_query_set, input_file_set, load_file_set, gen_map_file_set, si_save;
char *input_file, *load_file, *gen_map_file;
int B{128};

char *fifo_query;

void signal_handler(int signum) {
	unlink(fifo_query);
	exit(signum);
}

void show_server_usage(string name) {
	cout << "Syllable-Query: space-efficient long-match query algorithm.\n" <<
		"Server usage: " << name << " <flags>\n" <<
		"Flags:\n" <<
		"\t-h,--help                            Show this help message\n" <<
		"\t-f,--fifo <FILE>                     File path to create named pipe for communication with clients\n" <<
		"\t-i,--input_panel <FILE>              Path to VCF input file\n" <<
		"\t-b,--bits <VALUE>                    Value of B (64 or 128; default: 128), the number of sites to be grouped into one syllable\n" <<
		"\t                                     A higher B leads to stricter minimum query lengths but improved space efficiency.\n" <<
		"\t-s,--save                            Save panel data for fast loading when rerunning this program with --load\n" <<
		"\t                                     You will be told the save file upon creation.\n" <<
		"\t-l,--load <FILE>                     Path to load file (must have been the save file of a previous successful run)\n" <<
		"\t-g,--gen_map <FILE>                  Path to genetic map file, to be used when query length unit is cM\n" <<
		"\t                                     Format: Sites described by one line each with genetic location as the last tab-delimited field.\n" <<
		"\t                                             To accommodate header lines, lines with a non-numeric last field will be ignored.\n" <<
		"Example: " << name << " -f myfifo -i sample.vcf -s\n" <<
		"Usage note: The --save and --gen_map flags are optional. Exactly one of the flags --input_panel and --load should be used." << endl;
}

string show_client_usage() {
	return "\nYou may enter queries using the following flags:\n"
		"\t-q,--query_panel <FILE>              Path to VCF query file\n"
		"\t-l,--length <VALUE>                  Minimum length of desired matches\n"
		"\t-u,--unit <STRING>                   Unit of query length: sites (default), cM, or bps\n"
		"Example: -q query.vcf -l 700 -u sites\n"
		"Format of reported matches: Haplotypes will be referred to with the sample ID in the VCF file, followed by -0 or -1, denoting the first or second haplotype, respectively, of the individual. You will be told the location of the txt file containing one match reported per line, consisting of at least 4 tab-delimited fields $1 $2 $3 $4, indicating a match between query haplotype $1 and panel haplotype $2 over sites [$3, $4). If the query unit is cM or bps, there will be a fifth field indicating the length of the match in terms of the unit specified.";
}

void output(int id, string s) {
	size_t sz = s.length();
	if (write(id, (char*) &sz, sizeof(sz)) != (int) sizeof(sz)) cout << "Unsuccessful write to user" << endl;
	if (write(id, s.c_str(), sz) != (int) sz) cout << "Unsuccessful write to user" << endl;
}

template<class T>
void listen() {
	SyllableQuery<T> sq;

	char *wdc = getcwd(NULL, 0);
	string swd = string(wdc) + "/";
	free(wdc);

	if (input_file_set) {
		cout << "Precomputing on " << input_file << endl;
		int code = sq.precompute(input_file);
		if (code == 1) {
			cerr << "Input file does not exist or cannot be read." << endl;
			return;
		}
		if (code == 2) {
			cerr << "Input file not in VCF format." << endl;
			return;
		}
		cout << "Precomputation complete" << endl;

		if (gen_map_file_set) {
			cout << "Loading genetic map from " << gen_map_file << endl;
			code = sq.set_gen_map(gen_map_file);
			if (code == 1) {
				cerr << "Genetic map file does not exist or cannot be read." << endl;
			} else if (code == 2) {
				cerr << "Genetic map not in specified format." << endl;
			} else {
				cout << "Genetic map loaded." << endl;
			}
		}
	} else {
		cout << "Loading data from " << load_file << endl;
		int code = sq.load(load_file);
		if (code == 1) {
			cerr << "Load file does not exist or cannot be read." << endl;
			return;
		}
		if (code == 2) {
			cerr << "Load file is incomplete." << endl;
			return;
		}
		if (code == 3) {
			cerr << "The specified --bits value, " << B << ", does not coincide with the one in the load file, " << sq.B << "." << endl;
			return;
		}
		cout << "Loading data complete" << endl;
	}

	if (si_save) {
		string save_file;
		while (true) {
			save_file = string(fifo_query) + "-save" + to_string(save_ct++) + ".txt";
			ifstream tmp(save_file);
			if (tmp.fail()) break;
		}
		save_file = swd + save_file;
		cout << "Saving data to " << save_file << endl;
		int code = sq.save(save_file.c_str());
		if (code == 1) {
			cerr << "Save file cannot be written to." << endl;
		} else {
			cout << "Saving data complete" << endl;
		}
	}

	string query;
	size_t sz;
	unlink(fifo_query);
	umask(0);
	if (mkfifo(fifo_query, 0666) == -1) {
		cout << "Unable to create named pipe " << fifo_query << endl;
		return;
	}
	cout << "Ready for queries. Named pipe for server: " << swd << fifo_query << endl;
	int fdq = open(fifo_query, O_RDONLY);
	if (fdq < 0) {
		cerr << "Cannot read from named pipe " << fifo_query << endl;
		return;
	}
	while (true) {
		if (!read(fdq, (char*) &sz, sizeof(sz))) {
			fdq = open(fifo_query, O_RDONLY);
			if (fdq < 0) {
				cerr << "Cannot read from named pipe " << fifo_query << endl;
				return;
			}
			continue;
		}
		query.resize(sz);
		if (!read(fdq, &query[0], sz)) {
			cerr << "Failed to read client message" << endl;
			continue;
		}
		cout << "Query: " << query << endl;

		stringstream ss(query);
		string fifo_feedback;
		getline(ss, fifo_feedback, ' ');
		int fdf = open(fifo_feedback.c_str(), O_WRONLY);
		if (fdf < 0) {
			cerr << "Cannot write to named pipe " << fifo_feedback << endl;
			continue;
		} 

		string qwd;
		if (!getline(ss, qwd, ' ')) {
			output(fdf, sq.show_attributes() + "\n" + show_client_usage() + "\n\nEnter query:");
			continue;
		}

		string flag, val, query_file, unit = "sites";
		double L = 0;
		bool skip = false;

		while (getline(ss, flag, ' ')) {
			if (!getline(ss, val, ' ')) {
				output(fdf, "No value specified for your last flag.");
				skip = true; break;
			}

			if (flag == "-q" || flag == "--query_panel") {
				query_file = val;
			} else if (flag == "-l" || flag == "--length") {
				try { L = stod(val); }
				catch (exception &e) {
					output(fdf, "Length must be a number.");
					skip = true; break;
				}
			} else if (flag == "-u" || flag == "--unit") {
				unit = val;
			} else {
				output(fdf, "Invalid flag: " + flag);
				output(fdf, show_client_usage());
				skip = true; break;
			}
		}

		if (skip) continue;
		if (query_file.size() == 0u) {
			output(fdf, "Query file must be specified.");
			continue;
		}
		if (L == 0) {
			output(fdf, "Match length must be specified.");
			continue;
		}
		if (unit != "sites" && unit != "cM" && unit != "bps") {
			output(fdf, "Unit must be sites, cM, or bps.");
			continue;
		}
		if (query_file[0] != '/') query_file = qwd + query_file;
		string output_file;
		while (true) {
			output_file = string(fifo_query) + "-results" + to_string(out_ct++) + ".txt";
			ifstream tmp(output_file);
			if (tmp.fail()) break;
		}
		output_file = swd + output_file;

		output(fdf, "\tPerforming query on " + query_file);
		const clock_t START = clock();
		int code;
		if (unit == "sites") {
			code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L));
		} else if (unit == "cM") {
			code = sq.query(query_file.c_str(), output_file.c_str(), L, sq.geneLocs, false);
		} else { // "bps"
			code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L), sq.physLocs, true);
		}

		if (code == 1) {
			output(fdf, "Query file does not exist or cannot be read.");
		} else if (code == 2) {
			output(fdf, "Output file " + output_file + " cannot be written to.");
		} else if (code == 3) {
			output(fdf, "A genetic map must have been specified to query with genetic distance.");
		} else if (code == 4) {
			output(fdf, "The query length, " + to_string(L) + " " + unit + ", is less than the minimum allowable query length, " + to_string(unit == "sites" ? sq.minSiteL : unit == "cM" ? sq.minGeneL : sq.minPhysL) + " " + unit + ".");
		} else if (code == 5) {
			output(fdf, "Query file not in VCF format or has a different number of sites from the initial panel.");
		} else {
			output(fdf, "Query results outputted to " + output_file + "\nElapsed CPU time (s): " + to_string(double(clock() - START) / CLOCKS_PER_SEC));
		}
	}
}

int main(int argc, char** argv) {
	signal(SIGINT, signal_handler);

	if (argc == 1) {
		show_server_usage(argv[0]);
		return 0;
	}

	for (int i = 1; i < argc; i++) {
		string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			show_server_usage(argv[0]);
			return 0;
		}
		if (arg == "-s" || arg == "--save") {
			si_save = true;
			continue;
		}
		if (++i == argc) {
			cerr << "No value specified for your last flag." << endl;
			return 0;
		}

		if (arg == "-f" || arg == "--fifo") {
			fifo_query = argv[i];
			fifo_query_set = true;
		} else if (arg == "-i" || arg == "-input_panel") {
			input_file = argv[i];
			input_file_set = true;
		} else if (arg == "-b" || arg == "--bits") {
                        try { B = stoi(argv[i]); }
			catch (exception &e) {
				cerr << "--bits must be set to 64 or 128 (128 by default)." << endl;
				return 0;
			} 
		} else if (arg == "-l" || arg == "--load") {
			load_file = argv[i];
			load_file_set = true;
		} else if (arg == "-g" || arg == "--gen_map") {
			gen_map_file = argv[i];
			gen_map_file_set = true;
		} else {
			cerr << "Invalid flag: " << arg << endl;
			show_server_usage(argv[0]);
			return 0;
		}
	}

	if (!fifo_query_set) {
		cerr << "A named pipe must be specified with the --pipe flag." << endl;
		return 0;
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
