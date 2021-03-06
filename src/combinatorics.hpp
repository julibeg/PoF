#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include "types.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include <vector>

using namespace std;

/**
 * type-casting wrapper for boost library function to calculate
 * binomial coefficients. might overflow when used with integer types
 */
template <typename T> T binom(T n, T k) {
    T result = boost::math::binomial_coefficient<double>(n, k);
    return result;
}

/**
 * print Matrix for debugging
 */
template <typename T> void print_matrix(const Matrix<T> m) {
    for (auto const& row : m) {
        for (auto const& elem : row) {
            printf("%d ", elem);
        }
        printf("\n");
    }
}

/**
 * recursively get all possible combinations of numbers {1 ... n} for each n in
 * Ns.
 */
template <typename T>

void rec_combine_ranges(const vector<T>& Ns, size_t col, vector<T> row_so_far,
                        Matrix<T>& final_result) {
    T n = Ns[col];
    for (size_t i = 1; i <= n; i++) {
        row_so_far[col] = i;
        // if we are in the last column add to result
        if (col == Ns.size() - 1) {
            final_result.push_back(row_so_far);
        } else {
            rec_combine_ranges(Ns, col + 1, row_so_far, final_result);
        }
    }
}

/**
 * get combinations of numbers of selected rxns (NSRs) from a vector with
 * number of compressed rxns (NCRs) in the respective cutset.
 */
template <typename T> Matrix<T> get_NSRs(vector<T> NCRs) {
    // determine the dimensions of final Matrix and reserve space
    size_t num_rows = 1;
    for (auto const& n : NCRs) {
        num_rows *= n;
    }
    Matrix<T> NSRs;
    NSRs.reserve(num_rows);
    rec_combine_ranges(NCRs, 0, vector<T>(NCRs.size()), NSRs);
    return NSRs;
}

/**
 * get the number of uncompressed cutsets that can be built from a given
 * set of NSRs and the corresponding NCRs.
 * example:
 * NCRs={2, 2, 3}, select 1 rxn each (--> NSRs={1, 1, 1}) yields 12
 * possible combinations.
 */
template <typename T>
vector<T> get_combs(const vector<T>& NCRs, const Matrix<T>& NSRs) {
    vector<T> result;
    result.reserve(NSRs.size());
    unsigned int prod;
    for (auto const& row : NSRs) {
        prod = 1;
        for (size_t i = 0; i < NCRs.size(); i++) {
            T n = NCRs[i];
            prod *= binom(n, row[i]);
        }
        result.push_back(prod);
    }
    return result;
}

/**
 * get sum of a vector
 */
template <typename T> inline T sum_vec(vector<T> vec) {
    T sum = 0;
    for (auto const& elem : vec) {
        sum += elem;
    }
    return sum;
}

/**
 * binomial distribution as weight function
 */
template <typename T1, typename T2>
double binom_dist_weight(T1 n, T2 k, double p = 0.5) {
    return binom(static_cast<double>(n), static_cast<double>(k)) * pow(p, k) *
           pow(1 - p, n - k);
}

/**
 * resolve compressed cutsets from vector with number of compressed rxns (NCRs).
 * the function goes through all possible combinations of deletions in the
 * uncompressed model that would lead to the cut set in the compressed case
 * and returns the respective Mjs and counts.
 * can cache results.
 */
template <typename T>
pair<vector<T>, vector<T>> get_Mjs_and_counts(const vector<T>& NCRs,
                                              bool use_cache = true) {
    static map<vector<T>, pair<vector<T>, vector<T>>> cache;
    if (use_cache) {
        // search the cache table
        auto search = cache.find(NCRs);
        if (search != cache.end()) {
            return search->second;
        }
    }
    // result not found --> calculate it and add to cache
    Matrix<T> NSRs = get_NSRs(NCRs);
    vector<T> counts = get_combs(NCRs, NSRs);
    vector<T> Mjs;
    Mjs.reserve(NSRs.size());
    for (size_t i = 0; i < NSRs.size(); i++) {
        Mjs.push_back(sum_vec(NSRs[i]));
    }
    if (use_cache) {
#pragma omp critical
        { cache[NCRs] = pair<vector<T>, vector<T>>{Mjs, counts}; }
    }
    return pair<vector<T>, vector<T>>{Mjs, counts};
}

/*
 * resolve compressed cut set returning a table that can be added to the main
 * result in PoF_calculator.hpp
 */
template <typename T>
map<size_t, int>
resolve_compressed_cutset(const vector<T>& NCRs, unsigned int max_d,
                          unsigned int depth = 1, bool use_cache = true) {
    map<size_t, int> table;
    unsigned int Mj, J, m = NCRs.size();
    char sign;
    auto Mjs_counts = get_Mjs_and_counts(NCRs, use_cache);
    auto Mjs = Mjs_counts.first;
    auto counts = Mjs_counts.second;
    for (size_t i = 0; i < Mjs.size(); i++) {
        Mj = Mjs[i];
        J = depth + Mj - m;
        sign = (J % 2) ? 1 : -1;
        if (Mj > max_d) {
            continue;
        }
        table[Mj] += round((int)counts[i] * sign);
    }
    return table;
}

#endif /*  COMBINATORICS_HPP */
