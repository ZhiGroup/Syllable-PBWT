// This code is licensed under MIT license (see LICENSE for details)

#include <iostream>
#include <fcntl.h>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

char *feed_ch;
char *fifo_query;
string fifo_feedback;
string ffs, wd;

void send(string msg) {
	cout << "Waiting on server" << endl;
	int fdq = open(fifo_query, O_WRONLY);
	if (fdq < 0) {
		cout << "Invalid named pipe for writing to server" << endl;
		return;
	}
	size_t sz = msg.length();
	if (write(fdq, (char*) &sz, sizeof(sz)) != (int) sizeof(sz)) {
		cout << "Unsuccessful write to server" << endl;
		return;
	}
	if (write(fdq, msg.c_str(), sz) != (int) sz) {
		cout << "Unsuccessful write to server" << endl;
		return;
	}
	int fdf = open(fifo_feedback.c_str(), O_RDONLY);
	if (fdf < 0) {
		cout << "Invalid named pipe for reading from server" << endl;
		return;
	}

	string feed;
	while (true) {
		if (!read(fdf, (char*) &sz, sizeof(sz))) continue;
		feed.resize(sz);
		if (!read(fdf, &feed[0], sz)) {
			cout << "Failed to read server message" << endl;
			continue;
		}

		if (feed[0] == '\n') {
			cout << feed << endl;
			string line; getline(cin, line);
			string query = ffs + " " + wd + " " + line;
			sz = query.length();
			if (write(fdq, (char*) &sz, sizeof(sz)) != (int) sizeof(sz)) {
				cout << "Unsuccessful write to server" << endl;
				break;
			}
			if (write(fdq, query.c_str(), sz) != (int) sz) {
				cout << "Unsuccessful write to server" << endl;
				break;
			}
		} else if (feed[0] == '\t') {
			cout << feed.substr(1, feed.length() - 1) << endl;
		} else {
			cout << feed << endl;
			break;
		}
	}
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "The first parameter should be the named pipe for writing to the server." << endl;
		return 1;
	}
	fifo_query = argv[1];
	int ct = 0;
	while (true) {
		fifo_feedback = "fifo-" + to_string(ct++);
		if (access(fifo_feedback.c_str(), F_OK) == -1) break;
	}
	unlink(fifo_feedback.c_str());
	umask(0);
	if (mkfifo(fifo_feedback.c_str(), 0666) == -1) {
		cout << "Unable to create named pipe " << fifo_feedback << endl;
		return 0;
	}

	char *wdc = getcwd(NULL, 0);
	wd = string(wdc) + "/";
	free(wdc);
	string fqs;
	if (fifo_query[0] != '/') {
		fqs = wd + string(fifo_query);
		fifo_query = &fqs[0];
	}
	ffs = wd + string(fifo_feedback);
	fifo_feedback = &ffs[0];

	if (argc == 2) {
		send(ffs);
	} else {
		ifstream query_list(argv[2]);
		if (query_list.fail()) {
			cout << "File of queries does not exist or cannot be read." << endl;
			return 0;
		}

		string line;
		while (getline(query_list, line)) {
			cout << '\n' << line << endl;
			send(ffs + " " + wd + " " + line);
		}
	}

	unlink(fifo_feedback.c_str());

	return 0;
}
