#include "command_line_args.hpp"
#include <iomanip>
#include <ios>
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
    parsed_options cmd_opts = parse_cmd_line(argc, argv);
    cout << "mcs: " << cmd_opts.mcs_fname << endl;
    cout << "compr: " << cmd_opts.compr_rxn_fname << endl;
    cout << "rxns: " << cmd_opts.num_uncompressed_rxns << endl;
    cout << "d: " << cmd_opts.max_d << endl;
    cout << "l: " << cmd_opts.lambd << endl;
    cout << "t: " << cmd_opts.threads << endl;

}
