#include "cxxopts.hpp"
#include <iostream>
using namespace std;

struct parsed_options{
    string mcs_fname;
    string compr_rxn_fname;
    unsigned int num_uncompressed_rxns;
    unsigned int max_d;
    unsigned int threads;
    double lambd;
};


parsed_options parse_cmd_line(int argc, char* argv[]){
    parsed_options parsed_options;

    cxxopts::Options options("PoF_calculator",
        "Calculate failure probability (PoF) of metabolic networks given low- "
        "cardinality minimal cut sets (MCS). Capable of using MCS obtained from "
        "a compressed network to yield the PoF of the original uncompressed one.\n");

    options.custom_help("-m MCS_file [OPTIONS...]");
	options.add_options()
        ("m,mcs", "file with binary-encoded MCS (e.g. '100010100...')",
            cxxopts::value<string>(), " ")
        ("c,compr", "file with space-separated numbers of linearly compressed "
            "reactions per column in MCS file (e.g. '1 1 3 1 2...'). "
            "If not provided, assumes uncompressed network.",
            cxxopts::value<string>(), " ")
        ("r,num_rxns", "number of reactions in the uncompressed network; not "
            "required for uncompressed case. "
            "[default=sum of numbers in compr. rxn file",
            cxxopts::value<unsigned int>(), " ")
        ("d,d_max", "maximum cardinality up to which PoF should be calculated. "
            "[default=number of uncrompr. rxns]",
            cxxopts::value<unsigned int>(), " ")
        ("t,threads", "number of threads.",
            cxxopts::value<unsigned int>(), " ")
        ("l,lambda", "lambda value used in poisson distribution modeling "
            "mutation frequency. [default=0.5]",
            cxxopts::value<double>(), " ")
        ("h,help", "Print this message");

	auto result = options.parse(argc, argv);

    if (result.count("help")) {
        cout << options.help() << endl;
        exit(0);
    }

    if (result.count("mcs")) {
        parsed_options.mcs_fname = result["mcs"].as<string>();
    } else {
        cout << "Please provide at least an MCS file\n" << endl;
        cout << options.help() << endl;
        exit(0);
    }

    parsed_options.compr_rxn_fname = (result.count("compr")) ?
        result["compr"].as<string>() : "";

    parsed_options.num_uncompressed_rxns = (result.count("num_rxns")) ?
        result["num_rxns"].as<unsigned int>() : 0;

    parsed_options.max_d = (result.count("d_max")) ?
        result["d_max"].as<unsigned int>() : 0;

    parsed_options.lambd = (result.count("lambd")) ?
        result["lambd"].as<double>() : 0.5;

    parsed_options.threads = (result.count("threads")) ?
        result["threads"].as<unsigned int>() : 1;

    return parsed_options;
}
