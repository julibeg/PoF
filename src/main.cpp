#include "PoF_calculator.hpp"
#include "command_line_args.hpp"

using namespace std;


int main(int argc, char* argv[]){
	// parse command line arguments
	parsed_options cmd_opts = parse_cmd_line(argc, argv);

	// print command line args
	cout << "Running PoFcalc:\n";
	cmd_opts.print();
	cout << endl;

	// instantiate calculator class for compressed or uncompressed case
	PoF_calculator calc;
	if (cmd_opts.compr_rxn_fname.size() == 0) {                 // uncompressed
		calc = PoF_calculator(cmd_opts.mcs_fname);
	} else {
		calc = PoF_calculator(cmd_opts.mcs_fname,               // compressed
		                      cmd_opts.compr_rxn_fname,
		                      cmd_opts.num_uncompressed_rxns);
	}

	// perform recursive cutset search
	calc.get_cardinalities(cmd_opts.max_d, cmd_opts.threads, cmd_opts.use_cache);

	// print result
	calc.print_results(cmd_opts.p, cmd_opts.dm);

	// print result table for debugging
	// cout << endl;
	// calc.print_cd_table();

	return 0;
}
