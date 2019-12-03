#include "VariadicTable.h"
#include "cutset_compressed.cpp"
#include "types.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
using namespace std;

class PoF_calculator {
public:
	bool m_compressed = false;
	vector<unsigned int> m_compr_rxn_counts;
	vector<size_t> m_mcs1_rxns;
	size_t m_r, m_nbytes, m_nMCS, m_r_reduced, m_nbytes_reduced, m_nMCS_reduced,
	       m_last_d1, m_num_mcs1 = 0, m_num_mcs1_uncompressed = 0;
	bool m_MCS_d1_present = false;
	vector<Cutset> m_MCSs;
	Matrix<long> m_cd_table;
	unsigned int m_max_d;

	// default constructor
	PoF_calculator(){
	}

	// constructor for compressed networks --> delegating to constructor for
	// uncomprressed (constructor delegation is new in C++11)
	PoF_calculator(const string &mcs_input_fname, const string &comp_rxn_fname,
		size_t r=0)
		: PoF_calculator(mcs_input_fname){    // delegate to other constructor
		read_comp_rxn_file(comp_rxn_fname);
		m_compressed = true;
		// if # of uncompr. rxns has not been provided, take sum of compr. rxns
		if (r == 0) {
			m_r = 0;
			for (auto rxn_count : m_compr_rxn_counts) {
				m_r += rxn_count;
			}
		} else {
			m_r = r;
		}
		// find the number of uncompressed mcs1 rxns
		for (size_t rxn_id : m_mcs1_rxns) {
			m_num_mcs1_uncompressed += m_compr_rxn_counts[rxn_id];
		}
		// reduce the vector of compressed rxns
		vector<unsigned int> temp;
		temp.reserve(m_compr_rxn_counts.size() - m_mcs1_rxns.size());
		for (size_t i = 0; i < m_compr_rxn_counts.size(); i++) {
			if (!in_vec(m_mcs1_rxns, i)) {
				temp.push_back(m_compr_rxn_counts[i]);
			}
		}
		m_compr_rxn_counts = temp;
	}

	// constructor for uncompressed case
	PoF_calculator(const string &mcs_input_fname){
		vector<Cutset> MCSs = read_MCS_file(mcs_input_fname);
		m_r = MCSs[0].m_len;
		m_nbytes = MCSs[0].m_nbytes;
		m_nMCS = MCSs.size();
		// only reduce matrix if MCS with d=1 and MCS with d>1 are present
		if ((MCSs[0].CARDINALITY() == 1) && MCSs.back().CARDINALITY() > 1) {
			m_MCS_d1_present = true;
			cout << "Reducing MCS matrix..." << endl;
			MCSs = reduce_MCS_arr(MCSs);
		} else { // no MCS with d=1 or d>1 present --> no reduction
			m_r_reduced = m_r;
			m_nbytes_reduced = m_nbytes;
			m_nMCS_reduced = m_nMCS;
		}
		m_MCSs = MCSs;
	}

	void read_comp_rxn_file(const string &fname){
		string line;
		ifstream file(fname);
		if (file.is_open()) {
			cout << "Reading comp. rxns file..." << endl;
			getline(file, line);
			file.close();
		} else {
			cout << "Error opening compr. rxns file" << endl;
			exit(EXIT_FAILURE);
		}
		// split line at whitespace
		string comp_rxn;
		istringstream iss(line);
		while (getline(iss, comp_rxn, ' ')) {
			m_compr_rxn_counts.push_back(stoi(comp_rxn));
		}
	}

	vector<Cutset> read_MCS_file(const string &fname){
		string line;
		vector<Cutset> MCSs;
		ifstream file(fname);
		if (file.is_open()) {
			cout << "Reading MCS file..." << endl;
			while (getline(file, line)) {
				Cutset cs(line);
				MCSs.push_back(cs);
			}
			file.close();
		} else {
			cout << "Error opening MCS file" << endl;
			exit(EXIT_FAILURE);
		}
		return MCSs;
	}

	vector<Cutset> reduce_MCS_arr(const vector<Cutset> &MCSs){
		// find number of MCS(d=1) reactions and the corresponding reactions
		for (size_t i = 0; i < MCSs.size(); i++) {
			if (MCSs[i].CARDINALITY() > 1) {
				m_last_d1 = i - 1;
				m_num_mcs1 = i;
				break;
			}
			m_mcs1_rxns.push_back(MCSs[i].get_first_active_rxn());
		}
		m_r_reduced = m_r - (m_last_d1 + 1);
		m_nbytes_reduced = int_div_ceil(m_r_reduced, (size_t)CHAR_BIT);
		m_nMCS_reduced = m_nMCS - (m_last_d1 + 1);
		// initialize return vector
		vector<Cutset> reduced_arr;
		reduced_arr.reserve(m_nMCS_reduced);
		for (size_t i = m_last_d1 + 1; i < m_nMCS; i++) {
			Cutset reduced_cs = MCSs[i].remove_rxns(m_mcs1_rxns);
			reduced_arr.push_back(reduced_cs);
		}
		return reduced_arr;
	}

	void add_MCS1_to_table(){
		cout << "adding MCS(d=1) to table..." << endl;
		unsigned int num_mcs1 =
			(m_compressed) ? m_num_mcs1_uncompressed : m_num_mcs1;
		size_t Mj = 1; // cardinality of MCS1 = 1
		for (size_t a = 0; a < num_mcs1; a++) {
			m_cd_table[Mj - 1][a]++;
		}
	}

	void print_cd_table_2d_vec() const {
		using cell_type = decltype(m_cd_table)::value_type::value_type;
		VariadicTable<size_t, size_t, cell_type> vt({"|Mj|", "a", "#"});
		for (size_t Mj = 0; Mj < m_cd_table.size(); Mj++) {
			for (size_t a = 0; a < m_cd_table[Mj].size(); a++) {
				auto count = m_cd_table[Mj][a];
				if (count != 0) {
					// actual Mj is larger (0-based indexing)
					tuple<size_t, size_t, cell_type> vt_row {Mj + 1, a, count};
					vt.addRow(vt_row);
				}
			}
		}
		vt.print(cout);
	}

	template <typename T>

	void fill_2d_vec(vector<vector<T> > &vec, T new_value=0) const {
		for (auto &row : vec) {
			fill(row.begin(), row.end(), new_value);
		}
	}

	void add_to_cd_table(const Counter &temp_table){
		for (auto elem : temp_table) {
			auto Mj = get<0>(elem.first);
			auto a = get<1>(elem.first);
			auto count = elem.second;
			m_cd_table[Mj - 1][a] += count;
		}
	}

	void get_cardinalities(unsigned int max_d, unsigned int num_threads=1){
		// check if max_d is greater than the number of reactions
		if ((max_d > m_r) || (max_d == 0)) {
			max_d = m_r;
		}
		m_max_d = max_d;
		m_cd_table = Matrix<long>(max_d, vector<long>(m_r, 0));
		size_t last_MCS_to_consider = m_MCSs.size();
		if (m_MCS_d1_present) {
			add_MCS1_to_table();
		}

		cout << "Starting recursion..." << endl;
		// check if there are MCS with d > max_d
		if (m_MCSs.back().CARDINALITY() > max_d) {
			// get the last element with d <= max_d
			for (size_t i = 0; i < m_MCSs.size(); i++) {
				if (m_MCSs[i].CARDINALITY() > max_d) {
					last_MCS_to_consider = i;
					break;
				}
			}
		}
		// initialize openMP for loop
	#pragma omp parallel for num_threads(num_threads)
		for (size_t i = 0; i < last_MCS_to_consider; i++) {
			unsigned int mcs_card = m_MCSs[i].CARDINALITY();
			GET_CARDINALITIES(i, m_MCSs[i], mcs_card, max_d, 1,
			                  Cutset(m_r_reduced));
		}
	}

	void GET_CARDINALITIES(size_t index, const Cutset &Cs, unsigned int Cd,
	                       unsigned int max_d, unsigned int depth,
	                       Cutset stored){
		tuple<bool, bool, size_t> plus1_rxn_result;
		vector<size_t> still_to_check;
		size_t plus1_rxns = 0;
		Cutset testCs(Cs.m_len);
		unsigned int testCd;
		cout << "depth: " << depth << endl;
		// check for plus 1 rxns first
		for (size_t i = 0; i < index; i++) {
			if (!(m_MCSs[i] && stored)) {
				plus1_rxn_result = Cs.find_plus1_rxn(m_MCSs[i]);
				if (get<0>(plus1_rxn_result)) {
					stored.set_bit(get<2>(plus1_rxn_result));
				} else if (get<1>(plus1_rxn_result)) {
					still_to_check.push_back(i);
				} else { // m_MCSs[i] is a subset
					return;
				}
			}
		}
		// get number of plus 1 rxns
		if (m_compressed) {
			for (size_t rxn_id : stored.get_active_rxns()) {
				plus1_rxns += m_compr_rxn_counts[rxn_id];
			}
			plus1_rxns += m_num_mcs1_uncompressed;          // add MCS1 rxns
		} else {
			plus1_rxns += m_num_mcs1;
		}
		// perform additional recursions if required
		if (Cd < max_d) {
			for (size_t j : still_to_check) {
				if (!(m_MCSs[j] && stored)) {
					testCs = Cs | m_MCSs[j];
					testCd = testCs.CARDINALITY();
					if (testCd <= max_d) {
						// no need to check for testCd > Cd, since testCs must
						// have at least 2 extra rxns
						GET_CARDINALITIES(j, testCs, testCd, max_d, depth + 1,
						                  stored);
					}
				}
			}
		}
		// get active rxns of current cutset
		if (m_compressed) {
			vector<size_t> Cs_rxns = Cs.get_active_rxns();
			vector<unsigned int> NCRs;
			NCRs.reserve(Cs_rxns.size());
			for (size_t rxn_id : Cs.get_active_rxns()) {
				NCRs.push_back(m_compr_rxn_counts[rxn_id]);
			}
			map<size_t, int> table = resolve_compressed_cutset(NCRs, max_d, depth);
			// single threaded now
	#pragma omp critical
			{
				for (auto elem : table) {
					size_t Mj = elem.first;
					int count = elem.second;
					m_cd_table[Mj - 1][plus1_rxns] += count;
				}
			}
		} else {
			int sign = pow(-1, depth - 1);
	#pragma omp critical
			{
				m_cd_table[Cd - 1][plus1_rxns] += sign;
			}
		}
	}

	double score_cd_table(unsigned int d) const {
		long count;
		double score = 0;
		for (size_t Mj = 0; Mj < m_cd_table.size(); Mj++) {
			for (size_t a = 0; a < m_cd_table[Mj].size(); a++) {
				count = m_cd_table[Mj][a];
				if (count != 0) {
					// the actual Mj is larger (0 indexing in 2d-vector)
					score += count * SCORE(m_r, a, Mj + 1, d);
				}
			}
		}
		return score;
	}

	static double SCORE(unsigned int r, unsigned int a, unsigned int Mj,
	                    unsigned int d){
		double prod1 = 1, prod2 = 1;
		for (size_t mj = 0; mj < Mj; mj++) {
			prod1 *= (d - mj) / (double)(r - a - mj);
		}
		for (size_t i = 0; i < a; i++) {
			prod2 *= (r - d - i) / (double)(r - i);
		}
		return prod1 * prod2;
	}

};
