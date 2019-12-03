#include <iostream>
#include "command_line_args.cpp"
using namespace std;



int main(int argc, char* argv[]) {
    parsed_options opts = parse_cmd_line(argc, argv);
    cout << "mcs file: " << opts.mcs_fname << endl;
    cout << "comp rxn file: " << opts.compr_rxn_fname << endl;
    printf("r: %d\n", opts.num_uncompressed_rxns);
    printf("max_d: %d\n", opts.max_d);
    printf("lambda: %f\n", opts.lambd);
    printf("threads: %d\n", opts.threads);
}
