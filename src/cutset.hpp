#ifndef CUTSET_HPP
#define CUTSET_HPP

#include "combinatorics.hpp"
#include "types.hpp"
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


// define datatype to be used to store the indices of deletions in a cut set.
typedef unsigned int rxn_idx;


// make '10000000' bitmask
unsigned char MASK = pow(2, CHAR_BIT - 1);


/*
 * integer division with ceiling
 */
template <typename T> T int_div_ceil(T a, T b){
	return (a + b - 1) / b;
}


/**
 * check whether x is in vec
 */
template <typename T> bool in_vec(const vector<T> &vec, const T &x){
	if (find(vec.begin(), vec.end(), x) != vec.end()) {
		return true;
	} else {
		return false;
	}
}


/*
 * class for representing a single cut set. it stores the indices of deleted
 * reactions in a vector. given the sparsity of cut sets in most situations
 * (e.g. max. 10 deletions in a network of 2000 reactions) this is more
 * efficient than most generic binary representations. in the context of cut
 * sets "active rxns" means deleted reactions.
 */
class Cutset {
public:
	size_t m_len;						// total number of reactions
	vector<rxn_idx> m_active_rxns;		// indices of deleted reactions

	Cutset(size_t);             // init. empty instance with number of rxns
	Cutset(const string &);     // init. from string of '0's and '1's

	void print(bool) const;
	unsigned int CARDINALITY() const;
	Cutset operator | (const Cutset &) const;
	bool operator && (const Cutset &) const;
	vector<rxn_idx> get_active_rxns() const;
	size_t get_first_active_rxn() const;
	void add_reaction(rxn_idx);
	Cutset remove_rxns(const vector<rxn_idx> &) const;
	tuple<bool, bool, rxn_idx> find_plus1_rxn(const Cutset &) const;

private:
	void extract_active_rxns_from_string(const string &);
};


/*
 * construct empty instance with a given number of reactions.
 */
Cutset::Cutset(size_t num_rxns) : m_len(num_rxns){
}


/*
 * construct from string of the form "0001010001..."
 */
Cutset::Cutset(const string &cs){
	m_len = cs.length();
	extract_active_rxns_from_string(cs);
}


/*
 * get positions of deletions in a binary string
 */
void Cutset::extract_active_rxns_from_string(const string &cs){
	for (size_t i = 0; i < cs.length(); i++) {
		if (cs[i] == '1') {
			m_active_rxns.push_back(i);
		}
	}
}


/*
 * print cut set for debugging
 */
void Cutset::print(bool new_line=true) const {
	for (rxn_idx i = 0; i < m_len; i++) {
		cout << (in_vec(m_active_rxns, i) ? "1" : "0");
	}
	if (new_line) {
		cout << endl;
	}
}


/*
 * return cardinality of cut set
 */
inline unsigned int Cutset::CARDINALITY() const {
	return m_active_rxns.size();
}


/*
 * overload | operator to give the union of two cut sets
 */
inline Cutset Cutset::operator | (const Cutset &other_CS) const {
	Cutset new_cs(m_len);
	new_cs.m_active_rxns.reserve(m_active_rxns.size() +
	                             other_CS.m_active_rxns.size());
	set_union(begin(m_active_rxns), end(m_active_rxns),
	          begin(other_CS.m_active_rxns), end(other_CS.m_active_rxns),
	          back_inserter(new_cs.m_active_rxns));
	return new_cs;
}


/*
 * overload && operator to check wheter two cut sets have at least one deletion
 * in common
 */
inline bool Cutset::operator && (const Cutset &other_CS) const {
	auto first1 = m_active_rxns.begin();
	auto last1 = m_active_rxns.end();
	auto first2 = other_CS.m_active_rxns.begin();
	auto last2 = other_CS.m_active_rxns.end();
	while (first1 != last1) {
		if (first2 == last2)
			return false;
		if (*first2 < *first1) {
			first2++;
		} else {
			if (*first1 == *first2) {
				return true;
			}
			++first1;
		}
	}
	return false;
}


/*
 * return indices of deletions
 */
vector<rxn_idx> Cutset::get_active_rxns() const {
	return m_active_rxns;
}


/*
 * get index of first deletion
 */
size_t Cutset::get_first_active_rxn() const {
	return m_active_rxns[0];
}


/*
 * return new Cutset with the reactions at the positions in del_rxns removed.
 * requires del_rxns to be sorted
 */
Cutset Cutset::remove_rxns(const vector<rxn_idx> &del_rxns) const {
	Cutset new_cs(m_len - del_rxns.size());
	new_cs.m_active_rxns.reserve(m_active_rxns.size());

	auto first1 = m_active_rxns.cbegin();
	auto last1 = m_active_rxns.cend();
	auto first2 = del_rxns.cbegin();
	auto last2 = del_rxns.cend();

	size_t del_rxn_count = 0;
	while (first1 != last1) {
		if (first2 == last2) {
			for (; first1 != last1; first1++) {
				new_cs.m_active_rxns.push_back(*first1 - del_rxn_count);
			}
			break;
		}
		if (*first2 <= *first1) {
			del_rxn_count++;
			if (*first1 == *first2) {
				first1++;
			}
			first2++;
		} else {
			new_cs.m_active_rxns.push_back(*first1++ - del_rxn_count);
		}
	}
	return new_cs;
}


/*
 * add deletion to cut set while keeping the indices in order
 */
void Cutset::add_reaction(rxn_idx new_rxn){
	auto itr = m_active_rxns.begin();
	for (; itr != m_active_rxns.end(); ++itr) {
		if (*itr >= new_rxn) {
			break;
		}
	}
	m_active_rxns.insert(itr, new_rxn);
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
 *      in the second case, the third element in the tuple is the number of
 *      extra reactions
 */
tuple<bool, bool, rxn_idx> Cutset::find_plus1_rxn(const Cutset &other_CS) const {
	typedef tuple<bool, bool, rxn_idx> result;
	auto first1 = m_active_rxns.begin();
	auto last1 = m_active_rxns.end();
	auto first2 = other_CS.m_active_rxns.begin();
	auto last2 = other_CS.m_active_rxns.end();
	unsigned int plus1_rxn_count = 0;
	rxn_idx plus1_rxn_idx;
	while (first1 != last1) {
		if (first2 == last2)
			break;
		if (*first2 < *first1) {
			plus1_rxn_idx = *first2++;
			plus1_rxn_count++;
			// if (plus1_rxn_count > 1) {
			// 	return result {false, true, 0};
			// }
		} else {
			if (*first1 == *first2) {
				*first2++;
			}
			first1++;
		}
	}
	if (last2 - first2) {
		plus1_rxn_count += last2 - first2;
		plus1_rxn_idx = *(last2 - 1);
		// if (plus1_rxn_count > 1) {
		// 	return result {false, true, 0};
		// }
	}
	if (plus1_rxn_count == 1) {
		return result {true, false, plus1_rxn_idx};
	} else if (plus1_rxn_count == 0) {
		return result {false, false, 0};
	} else {
		return result {false, true, plus1_rxn_count};
	}
}

#endif
