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
	vector<size_t> m_active_rxns;

	Cutset(size_t);                 // initialize empty bitarr with number of bytes
	Cutset(const string &); 		// initialize from string

	void print(bool) const;
	unsigned int CARDINALITY() const;
	Cutset operator | (const Cutset &) const;
	bool operator && (const Cutset &) const;
	vector<size_t> get_active_rxns() const;
	size_t get_first_active_rxn() const;
	void add_reaction(size_t);
	Cutset remove_rxns(const vector<size_t> &) const;
	tuple<bool, bool, size_t> find_plus1_rxn(const Cutset &) const;
	static Counter resolve_compressed_cutset(const vector<unsigned int> &,
	                                         unsigned int, unsigned int);

private:
	void extract_active_rxns_from_string(const string &);
};

Cutset::Cutset(size_t num_rxns) : m_len(num_rxns){}


Cutset::Cutset(const string &cs){
	m_len = cs.length();
	extract_active_rxns_from_string(cs);
}


void Cutset::extract_active_rxns_from_string(const string &cs){
	for (size_t i = 0; i < cs.length(); i++) {
		if (cs[i] == '1') {
			m_active_rxns.push_back(i);
		}
	}
}


void Cutset::print(bool new_line=true) const {
	for (size_t i = 0; i < m_len; i++) {
		cout << (in_vec(m_active_rxns, i) ? "1" : "0");
	}
	if (new_line){
		cout << endl;
	}
}


inline unsigned int Cutset::CARDINALITY() const {
	return m_active_rxns.size();
}


inline Cutset Cutset::operator | (const Cutset &other_CS) const {
	Cutset new_cs(m_len);
	new_cs.m_active_rxns.reserve(m_active_rxns.size() +
		other_CS.m_active_rxns.size());
	set_union(begin(m_active_rxns), end(m_active_rxns),
		begin(other_CS.m_active_rxns), end(other_CS.m_active_rxns),
        back_inserter(new_cs.m_active_rxns));
	return new_cs;
}


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


vector<size_t> Cutset::get_active_rxns() const {
	return m_active_rxns;
}


size_t Cutset::get_first_active_rxn() const {
	return m_active_rxns[0];
}

/*
 * requires del_rxns to be sorted!
 */
Cutset Cutset::remove_rxns(const vector<size_t>& del_rxns) const{
	Cutset new_cs(m_len - del_rxns.size());
	new_cs.m_active_rxns.reserve(m_active_rxns.size());

	auto first1 = m_active_rxns.cbegin();
    auto last1 = m_active_rxns.cend();
    auto first2 = del_rxns.cbegin();
    auto last2 = del_rxns.cend();

    size_t del_rxn_count = 0;
    while (first1 != last1) {
        if (first2 == last2){
            for (; first1 != last1; first1++){
                new_cs.m_active_rxns.push_back(*first1 - del_rxn_count);
            }
            break;
        }
        if (*first2 <= *first1) {
            del_rxn_count++;
            if (*first1 == *first2){
                first1++;
            }
            first2++;
        } else {
            new_cs.m_active_rxns.push_back(*first1++ - del_rxn_count);
        }
    }
    return new_cs;
}


void Cutset::add_reaction(size_t new_rxn){
	auto itr = m_active_rxns.begin();
	for (; itr!=m_active_rxns.end(); ++itr){
		if (*itr >= new_rxn){
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
 */
tuple<bool, bool, size_t> Cutset::find_plus1_rxn(const Cutset &other_CS) const {
	auto first1 = m_active_rxns.begin();
	auto last1 = m_active_rxns.end();
	auto first2 = other_CS.m_active_rxns.begin();
	auto last2 = other_CS.m_active_rxns.end();
	unsigned int plus1_rxn_count = 0;
	size_t plus1_rxn_idx;
	while (first1 != last1) {
		if (first2 == last2)
			break;
		if (*first2 < *first1) {
			plus1_rxn_idx = *first2++;
			plus1_rxn_count++;
			if (plus1_rxn_count > 1) {
				return tuple<bool, bool, size_t>(false, true, 0);
			}
		} else {
			if (*first1 == *first2){
				*first2++;
			}
			first1++;
		}
	}
	if (last2 - first2){
		plus1_rxn_count += last2 - first2;
		plus1_rxn_idx = *(last2 - 1);
		if (plus1_rxn_count > 1) {
			return tuple<bool, bool, size_t>(false, true, 0);
		}
	}
	if (plus1_rxn_count == 1) {
		return tuple<bool, bool, size_t>(true, false, plus1_rxn_idx);
	} else {
		return tuple<bool, bool, size_t>(false, false, 0);
	}
}


template <typename T>
map<size_t, int> resolve_compressed_cutset(const vector<T> &NCRs,
                                           unsigned int max_d,
                                           unsigned int depth=1){
	map<size_t, int> table;
	unsigned int Mj, J, m = NCRs.size();
	char sign;
	Matrix<T> NSRs = get_NSRs(NCRs);
	vector<T> counts = get_combs(NCRs, NSRs);
	for (size_t i = 0; i < NSRs.size(); i++) {
		Mj = sum_vec(NSRs[i]);
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
