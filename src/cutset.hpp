#ifndef CUTSET_HPP
#define CUTSET_HPP

#include "combinatorics.hpp"
#include <algorithm>
#include <iostream>
#include <limits.h>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

typedef map<pair<size_t, size_t>, int> Counter;

// make '10000000' mask
unsigned char MASK = pow(2, CHAR_BIT - 1);

// integer division with ceiling
template <typename T> T int_div_ceil(T a, T b){
	return (a + b - 1) / b;
}


template <typename T> bool in_vec(const vector<T> &vec, const T &x){
	if (find(vec.begin(), vec.end(), x) != vec.end()) {
		return true;
	} else {
		return false;
	}
}


vector<unsigned int> generate_lookup_table(size_t nbits){
	vector<unsigned int> lookup_table;
	unsigned int max = pow(2, nbits);
	for (size_t i = 0; i < max; i++) {
		unsigned int count = 0;
		for (size_t j = 0; j < nbits; j++) {
			if (i & (1 << j)) {
				count++;
			}
		}
		lookup_table.push_back(count);
	}
	return lookup_table;
}


// initialize lookup table
vector<unsigned int> lookup_table = generate_lookup_table(CHAR_BIT);

// use lookup table
inline unsigned int count_byte(const unsigned char &byte){
	return lookup_table[byte];
}


class Cutset {
public:
	size_t m_len;
	size_t m_nbytes;
	vector<unsigned char> m_bitarr;

	Cutset(size_t);                 // initialize empty bitarr with number of bytes
	Cutset(const string &); // initialize from string

	void print(bool) const;
	unsigned int CARDINALITY() const;
	Cutset operator | (const Cutset &) const;
	bool operator && (const Cutset &) const;
	vector<size_t> get_active_rxns() const;
	size_t get_first_active_rxn() const;
	Cutset remove_rxns(const vector<size_t> &) const;
	tuple<bool, bool, size_t> find_plus1_rxn(const Cutset &) const;
	void set_bit(size_t);
	static Counter resolve_compressed_cutset(const vector<unsigned int> &,
	                                         unsigned int, unsigned int);

private:
	bool read_bit(size_t) const;
	static void print_byte(const unsigned char);
	void string_to_bitarr(const string &);
};

Cutset::Cutset(size_t num_rxns) : m_len(num_rxns){
	m_nbytes = int_div_ceil(m_len, (size_t)CHAR_BIT); // int. div. with ceiling
	m_bitarr.resize(m_nbytes);
}


Cutset::Cutset(const string &cs){
	m_len = cs.length();
	m_nbytes = int_div_ceil(m_len, (size_t)CHAR_BIT); // int. div. with ceiling
	m_bitarr.resize(m_nbytes);
	string_to_bitarr(cs);
}


inline void Cutset::set_bit(size_t pos){
	size_t byte_pos = pos / CHAR_BIT;
	size_t bit_pos = pos % CHAR_BIT;
	m_bitarr[byte_pos] |= (MASK >> bit_pos);
}


void Cutset::string_to_bitarr(const string &cs){
	for (size_t i = 0; i < cs.length(); i++) {
		if (cs[i] == '1') {
			set_bit(i);
		}
	}
}


void Cutset::print_byte(unsigned char x){
	for (size_t i = 0; i < CHAR_BIT; i++) {
		putchar((x & (MASK >> i)) ? '1' : '0');
	}
}


void Cutset::print(bool new_line=true) const {
	for (size_t i = 0; i < m_nbytes; i++) {
		print_byte(m_bitarr[i]);
	}
	if (new_line) {
		cout << endl;
	}
}


inline unsigned int Cutset::CARDINALITY() const {
	unsigned int card = 0;
	for (size_t i = 0; i < m_nbytes; i++) {
		card += count_byte(m_bitarr[i]);
	}
	return card;
}


inline Cutset Cutset::operator | (const Cutset &other_CS) const {
	Cutset new_cs(m_len);
	for (size_t i = 0; i < new_cs.m_nbytes; i++) {
		new_cs.m_bitarr[i] = m_bitarr[i] | other_CS.m_bitarr[i];
	}
	return new_cs;
}


inline bool Cutset::operator && (const Cutset &other_CS) const {
	for (size_t i = 0; i < m_nbytes; i++) {
		if (m_bitarr[i] & other_CS.m_bitarr[i]) {
			return true;
		}
	}
	return false;
}


vector<size_t> Cutset::get_active_rxns() const {
	vector<size_t> active_rxns;
	// reserve some space to minimize frequent reallocations
	active_rxns.reserve(m_len / 100 > 10 ? m_len / 100 : 10);
	size_t pos;
	for (size_t byte = 0; byte < m_nbytes; byte++) {
		if (m_bitarr[byte]) {
			for (size_t bit = 0; bit < CHAR_BIT; bit++) {
				if (m_bitarr[byte] & MASK >> bit) {
					pos = byte * CHAR_BIT + bit;
					active_rxns.push_back(pos);
				}
			}
		}
	}
	return active_rxns;
}


size_t Cutset::get_first_active_rxn() const {
	size_t pos = 0;
	for (size_t byte = 0; byte < m_nbytes; byte++) {
		if (m_bitarr[byte]) {
			for (size_t bit = 0; bit < CHAR_BIT; bit++) {
				if (m_bitarr[byte] & MASK >> bit) {
					pos = byte * CHAR_BIT + bit;
					return pos;
				}
			}
		}
	}
	throw invalid_argument("Cutset does not have an active reaction");
}


bool Cutset::read_bit(size_t pos) const {
	size_t byte_pos = pos / CHAR_BIT;
	size_t bit_pos = pos % CHAR_BIT;
	if (MASK >> bit_pos & m_bitarr[byte_pos]) {
		return true;
	} else {
		return false;
	}
}


Cutset Cutset::remove_rxns(const vector<size_t> &del_rxns) const {
	Cutset new_cs(m_len - del_rxns.size());
	size_t old_pos = 0, new_pos = 0;
	for (size_t byte = 0; byte < m_nbytes; byte++) {
		for (size_t bit = 0; bit < CHAR_BIT; bit++) {
			old_pos = byte * CHAR_BIT + bit;
			if (!in_vec(del_rxns, old_pos)) {
				if (read_bit(old_pos)) {
					new_cs.set_bit(new_pos);
				}
				new_pos++;
			}
		}
	}
	return new_cs;
}


/**
 * loops over another cutset and returns whether
 *      -) there is only a single additional reaction active (i.e. not present
 *      in this CS)
 *           --> first bool == true
 *      -) there are more than 1 extra reactions in the other cutset
 *           --> second bool == true
 *      in the first case, the third element in the tuple is the index of the
 *      single extra reaction.
 */
tuple<bool, bool, size_t> Cutset::find_plus1_rxn(const Cutset &other_CS) const {
	size_t plus1_rxn_idx;
	unsigned char b1, b2, plus1_rxns;
	unsigned int plus1_rxn_count = 0;
	for (size_t byte = 0; byte < m_nbytes; byte++) {
		b1 = m_bitarr[byte];
		b2 = other_CS.m_bitarr[byte];
		plus1_rxns = b2 & ~b1;
		plus1_rxn_count += count_byte(plus1_rxns);
		// end if there's already more than 1 extra rxn
		if (plus1_rxn_count > 1) {
			return tuple<bool, bool, size_t>{false, true, 0};
		}
		// if there's only 1 extra rxn, save its index
		if (count_byte(plus1_rxns) == 1) {
			for (size_t bit = 0; bit < CHAR_BIT; bit++) {
				if (!(b1 & MASK >> bit) && (b2 & MASK >> bit)) {
					plus1_rxn_idx = byte * CHAR_BIT + bit;
					break;
				}
			}
		}
	}
	// loop completed; check if there was an extra rxn
	if (plus1_rxn_count == 1) {
		return tuple<bool, bool, size_t>{true, false, plus1_rxn_idx};
	} else {
		return tuple<bool, bool, size_t>{false, false, 0};
	}
}


// template <typename T>
// map<size_t, int> resolve_compressed_cutset(const vector<T> &NCRs,
//                                            unsigned int max_d,
//                                            unsigned int depth=1){
// 	map<size_t, int> table;
// 	unsigned int Mj, J, m = NCRs.size();
// 	char sign;
// 	Matrix<T> NSRs = get_NSRs(NCRs);
// 	vector<T> counts = get_combs(NCRs, NSRs);
// 	for (size_t i = 0; i < NSRs.size(); i++) {
// 		Mj = sum_vec(NSRs[i]);
// 		J = depth + Mj - m;
// 		sign = (J % 2) ? 1 : -1;
// 		if (Mj > max_d) {
// 			continue;
// 		}
// 		table[Mj] += round((int) counts[i] * sign);
// 	}
// 	return table;
// }


template <typename T>
map<size_t, int> resolve_compressed_cutset(const vector<T> &NCRs,
                                           unsigned int max_d,
                                           unsigned int depth=1){
	map<size_t, int> table;
	unsigned int Mj, J, m = NCRs.size();
	char sign;
	auto Mjs_counts = get_Mjs_and_counts(NCRs);
	auto Mjs = Mjs_counts.first;
	auto counts = Mjs_counts.second;
	for (size_t i = 0; i < Mjs.size(); i++) {
		Mj = Mjs[i];
		J = depth + Mj - m;
		sign = (J % 2) ? 1 : -1;
		if (Mj > max_d) {
			continue;
		}
		table[Mj] += round((int) counts[i] * sign);
	}
	return table;
}

#endif
