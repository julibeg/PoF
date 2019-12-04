#include "PoF_calculator.cpp"
#include "command_line_args.cpp"
#include "table.cpp"
#include <boost/math/distributions/poisson.hpp>

using namespace std;


int main(int argc, char* argv[]){
	// parse command line arguments
	parsed_options cmd_opts = parse_cmd_line(argc, argv);

	// instantiate calculator class for compressed or uncompressed case
	PoF_calculator calc;
	if (cmd_opts.compr_rxn_fname.size() == 0) {                     // uncompressed
		calc = PoF_calculator(cmd_opts.mcs_fname);
	} else {
		calc = PoF_calculator(cmd_opts.mcs_fname,               // compressed
		                      cmd_opts.compr_rxn_fname,
		                      cmd_opts.num_uncompressed_rxns);
	}

	// perform recursive cutset search
	calc.get_cardinalities(cmd_opts.max_d, cmd_opts.threads);

	// get poisson distribution for weights
	boost::math::poisson_distribution<double> weight_dist(cmd_opts.lambd);
	double score, weight, weighted_score, acc_weighted_score = 0, found_cutsets,
	       all_cutsets;

	// initialize table t0 print results
	Table table { {"d", "weight", "F(d)", "weighted F(d)", "acc. weighted F(d)",
		       "lethal CS", "possible CS"},
		      {5, 20, 20, 20, 20, 20, 20},
		      {"%.5g", "%.10g", "%.10g", "%.10g", "%.10g",
		       "%.10g", "%.10g"} };
	table.print_header();

	for (size_t d = 1; d <= calc.m_max_d; d++) {
		score = calc.score_cd_table(d);
		weight = boost::math::pdf(weight_dist, d);
		weighted_score = score * weight;
		acc_weighted_score += weighted_score;
		all_cutsets = binom((double) calc.m_r, (double) d);
		found_cutsets = score * all_cutsets;

		// print all the numbers
		vector<double> numbers {d, weight, score, weighted_score,
			                acc_weighted_score, found_cutsets, all_cutsets};
		table.print_row(numbers);
	}

	// calc.print_cd_table_2d_vec();
	// cout << calc.m_r << endl;

	return 0;
}
