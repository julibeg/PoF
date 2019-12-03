#include "PoF_calculator.cpp"
#include "command_line_args.cpp"

using namespace std;

int main(int argc, char* argv[]){
	// parse command line arguments
	parsed_options cmd_opts = parse_cmd_line(argc, argv);
	// logic for command line arguments (validated input and compr. vs uncompr.)
	PoF_calculator calc;
	if (cmd_opts.compr_rxn_fname.size() == 0) {
		calc = PoF_calculator(cmd_opts.mcs_fname);
	} else {
		calc = PoF_calculator(cmd_opts.mcs_fname,
		                      cmd_opts.compr_rxn_fname,
		                      cmd_opts.num_uncompressed_rxns);
	}

    calc.get_cardinalities(cmd_opts.max_d, cmd_opts.threads);
    for (size_t i = 1; i <= calc.m_max_d; i++) {
         printf("d=%li\t%.18f\n", i, calc.score_cd_table(i));
    }
    calc.print_cd_table_2d_vec();
    cout << calc.m_r << endl;
	return 0;
}
