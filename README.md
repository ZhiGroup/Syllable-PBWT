# Syllable-PBWT
## Introduction
Syllable-PBWT is a space-efficient variation of the PBWT data structure introduced by Richard Durbin. In this repository, you'll find an implementation of Syllable-Query, an algorithm for querying long haplotype matches using the lightweight Syllable-PBWT data structures.

## Dependencies
- Linux
- C++ (at least GCC 5)

## Installation
Clone the repository:

`git clone https://github.com/ZhiGroup/Syllable-PBWT`

Enter the repository and compile the code:

`cd Syllable-PBWT`

`make`

## Usage Instructions
Our code allows multiple client users to query from one server program as long as they are on the same host. The server help page can be viewed by running
`./server`:
|              Flag             |      Description         |        Details         |
|:-----------------------------:|:------------------------:|:---------------------:|
| `-h` or `--help`  | Show this help message   |  |
| `-f <FILE>` or `--fifo <FILE>`  | File path to create named pipe | For communication with clients |
| `-i <FILE>` or `--input_panel <FILE>` | Path to VCF input file | |
| `-b <VALUE>` or `--bits <VALUE>` | Value of B | B must be 64 or 128, by default 128. B is the number of sites to be grouped into one syllable. |
| `-s` or `--save` | Will save panel data | For fast loading when rerunning this program with `--load`. You will be told the save file upon creation. |
| `-l <FILE>` or `--load <FILE>` | Path to load file | Must be the save file of a previous successful run |
| `-g <FILE>` or `--gen_map <FILE>` | Path to genetic map file | The genetic map will be used when the query length unit is cM. Format: Sites described by one line each with genetic location as the last tab-delimited field. To accommodate header lines, lines with a non-numeric last field will be ignored. |

The client program should be run with the server's named pipe (specified with `--fifo`) as the only parameter (you will be told so if you run `./client` with no parameters). When the server is done preprocessing, the client will be prompted with directions on how to query as well as the output format. Below are the flags used to query:
|             Flag              |     Description         |
|:-----------------------------:|:-----------------------:|
| `-q <FILE>` or `--query_panel <FILE>` | Path to VCF query file |
| `-l <VALUE>` or `--length <VALUE>` | Minimum length of desired matches |
| `-u <STRING>` or `--units <STRING>` | Unit of query length: sites (default), cM, or bps |

## Sample Use
The server and client programs are to be run from different terminals on the same host. In this example, we will use the data found in the `sample_data` folder. On terminal 1, run

`./server -f fifo -i sample_data/sample_panel.vcf`

On terminal 2, run

`./client fifo`

and when prompted to enter a query, type

`-q sample_data/sample_query.vcf -l 255`

and then you will be told the file location of your match results.

## Benchmarking
For reference, the `benchmark` folder includes our implementation for Algorithm 3 of [Sanaullah et al.](https://academic.oup.com/bioinformatics/article/37/16/2390/6149123) against which we benchmarked Syllable-Query.

## Citation
If you found our work useful in your research, please consider citing the following paper:
```
@article{Syllable-PBWT,
  author = {Wang, Victor and Naseri, Ardalan and Zhang, Shaojie and Zhi, Degui},
  title = {{Syllable-PBWT for space-efficient haplotype long-match query}},
  year = {2022},
  doi = {10.1101/2022.01.31.478234}
}
```
