#ifndef COMMAND_LINE_ARGS_HPP
#define COMMAND_LINE_ARGS_HPP

#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


struct parsed_options {
    string mcs_fname;
    string compr_rxn_fname;
    unsigned int num_uncompressed_rxns = 0;
    unsigned int max_d = 0;
    unsigned int threads = 1;
    double p = 1e-3;
    bool use_cache=true;
};


/**
 * split string at whitespaces
 */
vector<string> string2words(const string& str) {
    string word;
    stringstream sstream(str);
    vector<string> words;

    while (sstream >> word) {
        words.push_back(word);
    }
    return words;
}


/**
 * wrap potentially long string into a field of a given length (i.e. into lines
 * with a maximum length)
 */
void wrap_in_field(string str, unsigned int field_width,
                   unsigned int offset = 0, unsigned int initial_offset = 0) {
    printf("%*s", initial_offset, "");
    if (str.size() <= field_width) {
        cout << setw(field_width) << left << str;
    } else {
        unsigned int curr_width = 0;
        vector<string> words = string2words(str);
        for (const string& word : words) {
            if (word.size() > field_width) {
                cout << endl << word << endl; // single word larger than field
                continue;
            }
            if (curr_width + word.size() <= field_width) {
                cout << word << " ";
                curr_width += word.size() + 1;
            } else {
                printf("\n%*s", offset + initial_offset, "");
                cout << word << " ";
                curr_width = word.size();
            }
        }
    }
}


/**
 * print formatted help string
 */
void print_help() {
    string header{
        "Calculate the probability of failure (PoF) for a metabolic network "
        "given low-cardinality minimal cut sets (MCS). Capable of using MCS "
        "obtained from a compressed network to yield the PoF of the original "
        "uncompressed one."};
    vector<vector<string>> options{
        {{"-m, --mcs"}, {"file with binary-encoded MCS (e.g. '100010100...')"}},
         {{"-c, --compr"}, {"file with space-separated numbers of linearly "
          "compressed "
          "reactions per column in MCS file (e.g. '1 1 3 1 2...'). "
          "If not provided, assumes uncompressed network."}},
        {{"-r, --rxns"},
         {"number of reactions in the uncompressed network; not "
          "required for uncompressed case. "
          "[default=sum of numbers in compr. rxn file"}},
        {{"-d, --d_max"},
         {"maximum cardinality up to which PoF should be "
          "calculated. [default=number of uncrompr. rxns]"}},
        {{"-t, --threads"}, {"number of threads. [default=1]"}},
        {{"-p, --prob"},
         {"estimated probability of a loss-of-function mutation. "
          "[default=1e-3]"}},
        {{"-n, --no_cache"},
          "provide this flag to disable caching results when resolving "
          "compressed cutsets"},
        {{"-h, --help"}, {"print this message"}}};
    wrap_in_field(header, 75);
    cout << endl << endl;
    cout << "Usage:" << endl;
    wrap_in_field({"PoF_calculator -m MCS_file [OPTIONS...]"}, 50, 0, 3);
    cout << endl << endl;
    for (const auto& opt : options) {
        wrap_in_field(opt[0], 15, 0, 3);
        wrap_in_field(opt[1], 50, 18, 0);
        cout << endl << endl;
    }
}


/**
 * function to parse command line options
 */
parsed_options parse_cmd_line(int argc, char* argv[]) {
    parsed_options parsed_options;
    if (argc == 1) {
        print_help();
        exit(0);
    }

    for (size_t i = 1; i < argc; i++) {
        string argument(argv[i]);
        if ((argument == "-h") || (argument == "--help")) {
            print_help();
            exit(0);
        } else if ((argument == "-m") || (argument == "--mcs")) {
            parsed_options.mcs_fname = argv[i + 1];
            i++;
        } else if ((argument == "-c") || (argument == "--compr")) {
            parsed_options.compr_rxn_fname = argv[i + 1];
            i++;
        } else if ((argument == "-r") || (argument == "--rxns")) {
            parsed_options.num_uncompressed_rxns = atoi(argv[i + 1]);
            i++;
        } else if ((argument == "-d") || (argument == "--d_max")) {
            parsed_options.max_d = atoi(argv[i + 1]);
            i++;
        } else if ((argument == "-p") || (argument == "--prob")) {
            double p = atof(argv[i + 1]);
            if ((p >= 1) || (p < 0)){
                cout << "ERROR: p should be between 0 and 1\n" << endl;
                print_help();
                exit(1);
            }
            parsed_options.p = p;
            i++;
        } else if ((argument == "-t") || (argument == "--threads")) {
            parsed_options.threads = atoi(argv[i + 1]);
            i++;
        } else if ((argument == "-n") || (argument == "--no_cache")) {
            parsed_options.use_cache = false;
        } else {
            cout << argv[i] << endl;
            print_help();
            exit(2);
        }
    }
    return parsed_options;
}

#endif /* COMMAND_LINE_ARGS_HPP */
