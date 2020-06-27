#ifndef POF_CALCULATOR_HPP
#define POF_CALCULATOR_HPP

#include "cutset.hpp"
#include "table.hpp"
#include "types.hpp"
// Luigi Pertoldi's progress bar from https://github.com/gipert/progressbar
#include "../include/progressbar/progressbar.hpp"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
using namespace std;

const double WEIGHT_LIMIT = 1e-20;
// const double WEIGHT_LIMIT = numeric_limits<double>::min();

/*
 * class for MCS import, pre-processing (i.e. reduction) and recursion.
 */
class PoF_calculator {
  public:
    // whether in compressed or uncompressed case
    bool m_compressed = false;
    // numb. of compr. rxns on each position if in compressed case
    vector<unsigned int> m_compr_rxn_counts;
    // positions of essential reactions
    vector<rxn_idx> m_mcs1_rxns;
    // numb. of rxns and MCSs for original MCS matrix and reduced form.
    // numb. of essential reactions in compressed and uncompressed network
    size_t m_r, m_nMCS, m_r_reduced, m_nMCS_reduced,
        m_num_mcs1 = 0, m_num_mcs1_uncompressed = 0;
    // whether essential rxns are present --> no recursion necessery otherwise
    bool m_MCS_d1_present = false;
    // cut sets
    vector<Cutset> m_MCSs;
    // result table
    Matrix<long> m_cd_table;
    // d0
    unsigned int m_max_d;

    // default constructor
    PoF_calculator() {
    }

    /*
     * constructor for compressed networks --> delegating to constructor for
     * uncomprressed case (constructor delegation is new in C++11)
     */
    PoF_calculator(const string& mcs_input_fname, const string& comp_rxn_fname,
                   size_t r = 0)
        : PoF_calculator(mcs_input_fname) { // delegate to other constructor
        read_comp_rxn_file(comp_rxn_fname);
        m_compressed = true;
        // if the numb. of uncompr. rxns has not been provided, r is the sum of
        // compr. rxns
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
        for (rxn_idx i = 0; i < m_compr_rxn_counts.size(); i++) {
            if (!in_vec(m_mcs1_rxns, i)) {
                temp.push_back(m_compr_rxn_counts[i]);
            }
        }
        m_compr_rxn_counts = temp;
    }

    /*
     * constructor for uncompressed case
     */
    PoF_calculator(const string& mcs_input_fname) {
        vector<Cutset> MCSs = read_MCS_file(mcs_input_fname);
        m_r = MCSs[0].m_len;
        m_nMCS = MCSs.size();
        // set 'reduced' variables (are overwritten in reduce_MCS_arr() later)
        m_r_reduced = m_r;
        m_nMCS_reduced = m_nMCS;
        // only reduce matrix if MCS with d=1 are present
        if ((MCSs[0].CARDINALITY() == 1)) {
            m_MCS_d1_present = true;
            cout << "Reducing MCS matrix...\n" << endl;
            MCSs = reduce_MCS_arr(MCSs);
        }
        m_MCSs = MCSs;
    }

    /*
     * read file with number of compressed rxns for each rxn in the compressed
     * network (file in the format of "3 4 2 1 4 3 6 ...")
     */
    void read_comp_rxn_file(const string& fname) {
        string line;
        ifstream file(fname);
        if (file.is_open()) {
            cout << "Reading comp. rxns file...\n" << endl;
            getline(file, line);
            file.close();
        } else {
            cout << "Error opening compr. rxns file" << endl;
            exit(EXIT_FAILURE);
        }
        // split line at whitespaces
        string comp_rxn;
        istringstream iss(line);
        while (getline(iss, comp_rxn, ' ')) {
            m_compr_rxn_counts.push_back(stoi(comp_rxn));
        }
    }

    /*
     * read file with MCSs in binary format
     * file should be sorted by MCS cardinality and look like:
     * '001000000...
     *  000000100...'
     */
    vector<Cutset> read_MCS_file(const string& fname) {
        string line;
        vector<Cutset> MCSs;
        ifstream file(fname);
        if (file.is_open()) {
            cout << "Reading MCS file...\n" << endl;
            while (getline(file, line)) {
                Cutset cs(line); // instantiate Cutset object for every line
                MCSs.push_back(cs);
            }
            file.close();
        } else {
            cout << "Error opening MCS file" << endl;
            exit(EXIT_FAILURE);
        }
        return MCSs;
    }

    /*
     * reduce MCS matrix by removing all essential reactions and the MCSs
     * containing those.
     */
    vector<Cutset> reduce_MCS_arr(const vector<Cutset>& MCSs) {
        // find number of essential rxns and the corresponding indices
        if (MCSs.back().CARDINALITY() == 1) {
            // matrix only has MCS with d=1 and is reduced to nothing
            m_num_mcs1 = MCSs.size();
            for (const Cutset& cs : MCSs) {
                m_mcs1_rxns.push_back(cs.get_first_active_rxn());
            }
            m_r_reduced = 0;
            m_nMCS_reduced = 0;
            return vector<Cutset>();

        } else {
            // matrix has MCS with d=1 and d>1 --> remove all MCS with d=1 and
            // the corresponding reactions
            for (size_t i = 0; i < MCSs.size(); i++) {
                if (MCSs[i].CARDINALITY() > 1) {
                    m_num_mcs1 = i;
                    break;
                }
                // add index of essential rxn to m_mcs1_rxns
                m_mcs1_rxns.push_back(MCSs[i].get_first_active_rxn());
            }
            // sort mcs1 rxns --> required for Cutset::remove_rxns
            sort(m_mcs1_rxns.begin(), m_mcs1_rxns.end());
            m_r_reduced = m_r - m_num_mcs1;
            m_nMCS_reduced = m_nMCS - m_num_mcs1;
            // initialize reduced MCS matrix
            vector<Cutset> reduced_arr;
            reduced_arr.reserve(m_nMCS_reduced);
            // add MCSs with non-essential rxns to reduced MCS matrix
            for (size_t i = m_num_mcs1; i < m_nMCS; i++) {
                Cutset reduced_cs = MCSs[i].remove_rxns(m_mcs1_rxns);
                reduced_arr.push_back(reduced_cs);
            }
            return reduced_arr;
        }
    }

    /*
     * pre-populate result table based on number of essential reactions
     */
    void add_MCS1_to_table() {
        unsigned int num_mcs1 =
            (m_compressed) ? m_num_mcs1_uncompressed : m_num_mcs1;
        size_t Mj = 1; // cardinality of MCS1 = 1
        for (size_t a = 0; a < num_mcs1; a++) {
            m_cd_table[Mj - 1][a]++;
        }
    }

    /*
     * pretty-print result table (e.g. for debugging)
     */
    void print_cd_table() const {
        Table table({"|Mj|", "a", "#"}, {15, 15, 15}, {"%.5g", "%g", "%g"}
                    // uses scientific fmt in case the numbers
                    ); // get too large
        table.print_header();
        for (size_t Mj = 0; Mj < m_cd_table.size(); Mj++) {
            for (size_t a = 0; a < m_cd_table[Mj].size(); a++) {
                auto count = m_cd_table[Mj][a];
                if (count != 0) {
                    // actual Mj is larger by 1 (0-based indexing)
                    table.print_row(vector<double>{Mj + 1, a, count});
                }
            }
        }
    }

    /*
     * main function preparing for and then invoking recursion
     */
    void get_cardinalities(unsigned int max_d, unsigned int num_threads = 1,
                           bool use_cache = true) {
        // check if d0 supplied at cmd line is greater than the number of rxns
        if ((max_d > m_r) || (max_d == 0)) {
            max_d = m_r;
        }
        m_max_d = max_d;
        // initialize result table with rows for every 0 < d <= d0 and columns
        // for every reaction (representing plus 1 rxns)
        m_cd_table = Matrix<long>(max_d, vector<long>(m_r, 0));
        size_t last_MCS_to_consider = m_MCSs.size();
        if (m_MCS_d1_present) {
            // add MCS1 to table
            cout << "adding MCS(d=1) to table...\n" << endl;
            add_MCS1_to_table();
            if (m_MCSs.size() == 0) {
                // MCS matrix was reduced to nothing --> only MCS with d=1 in
                // original matrix
                cout << "no MCS with d>1 present --> no recursion required\n"
                     << endl;
                return;
            }
        }
        cout << "Starting recursion...\n" << endl;
        // check if there are MCS with d > d0 (max_d)
        if (m_MCSs.back().CARDINALITY() > max_d) {
            // get the last element with d <= max_d
            for (size_t i = 0; i < m_MCSs.size(); i++) {
                if (m_MCSs[i].CARDINALITY() > max_d) {
                    last_MCS_to_consider = i;
                    break;
                }
            }
        }
        // setup progress bar
        progressbar prog_bar(last_MCS_to_consider);
        size_t j;
// initialize openMP for loop
#pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < last_MCS_to_consider; i++) {
            // start with high cardinality MCSs first --> loop speeds up towards
            // the end instead of slowing down --> nicer.
            size_t j = last_MCS_to_consider - i - 1;
#pragma omp task
            {
                unsigned int mcs_card = m_MCSs[j].CARDINALITY();
                // start recursion
                GET_CARDINALITIES(j, m_MCSs[j], mcs_card, max_d, 1,
                                  Cutset(m_r_reduced), use_cache);
#pragma omp critical
                { prog_bar.update(); }
            }
        }
        // add new lines after progress bar
        cout << "\n\n" << endl;
    }

    /*
     * implement the recursive algorithm
     */
    void GET_CARDINALITIES(size_t index, const Cutset& Cs, unsigned int Cd,
                           unsigned int max_d, unsigned int depth,
                           Cutset stored, bool use_cache) {
        tuple<bool, bool, size_t> plus1_rxn_result;
        vector<size_t> still_to_check;
        still_to_check.reserve(index);
        size_t plus1_rxns = 0;
        unsigned int testCd;
        // check for plus 1 rxns first
        for (size_t i = 0; i < index; i++) {
            if (!(m_MCSs[i] && stored)) {
                plus1_rxn_result = Cs.find_plus1_rxn(m_MCSs[i]);
                if (get<0>(plus1_rxn_result)) {
                    stored.add_reaction(get<2>(plus1_rxn_result));
                } else if (get<1>(plus1_rxn_result)) {
                    // check if the union of Cs and m_MCSs[i] would have too
                    // high cardinality --> don't check later
                    if (Cd + get<2>(plus1_rxn_result) > max_d) {
                        continue;
                    }
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
            plus1_rxns += m_num_mcs1_uncompressed; // add MCS1 rxns
        } else {
            plus1_rxns += stored.CARDINALITY() + m_num_mcs1;
        }
        // perform additional/deeper recursions if required
        if (Cd < max_d) {
            for (size_t j : still_to_check) {
                if (!(m_MCSs[j] && stored)) {
                    Cutset testCs = Cs | m_MCSs[j]; // looks inefficient to
                    // create a new Cutset here in every iteration but is
                    // actually not done due to compiler optimization
                    testCd = testCs.CARDINALITY();
                    if (testCd <= max_d) {
                        // no need to check for testCd > Cd, since testCs must
                        // have at least 2 extra rxns
                        GET_CARDINALITIES(j, testCs, testCd, max_d, depth + 1,
                                          stored, use_cache);
                    }
                }
            }
        }
        // get active rxns of current cutset
        if (m_compressed) {
            vector<rxn_idx> Cs_rxns = Cs.get_active_rxns();
            vector<unsigned int> NCRs;
            // get NCRs to resolved to uncompressed case later
            NCRs.reserve(Cs_rxns.size());
            for (auto rxn_id : Cs_rxns) {
                NCRs.push_back(m_compr_rxn_counts[rxn_id]);
            }
            sort(NCRs.begin(), NCRs.end());
            // resolve compressed cut set
            map<size_t, int> table =
                resolve_compressed_cutset(NCRs, max_d, depth, use_cache);
// single threaded now while adding to result table
#pragma omp critical
            {
                for (const auto& elem : table) {
                    size_t Mj = elem.first;
                    int count = elem.second;
                    m_cd_table[Mj - 1][plus1_rxns] += count;
                }
            }
        } else {
            int sign = (depth % 2) ? 1 : -1;
#pragma omp critical
            { m_cd_table[Cd - 1][plus1_rxns] += sign; }
        }
    }

    /*
     * evaluate result table to give f(d, m0)
     */
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

    /*
     * evaluate result table to give f(d, m1) --> i.e. use only the essential
     * reactions
     */
    double get_f1(unsigned int d) const {
        long count;
        double f1 = 0;
        // only iterate over the first row in table --> all MCSs with d=1
        for (size_t a = 0; a < m_cd_table[0].size(); a++) {
            count = m_cd_table[0][a];
            if (count != 0) {
                // the actual Mj 1
                f1 += count * SCORE(m_r, a, 1, d);
            }
        }
        return f1;
    }

    /*
     * evaluate result table to give f(d, m0)
     */
    tuple<double, double> score_cd_table2(unsigned int d) const {
        long count;
        double f1 = get_f1(d);
        double score = f1;
        if (m_cd_table.size() > 1) {
            for (size_t Mj = 1; Mj < m_cd_table.size(); Mj++) {
                for (size_t a = 0; a < m_cd_table[Mj].size(); a++) {
                    count = m_cd_table[Mj][a];
                    if (count != 0) {
                        // the actual Mj is larger (0 indexing in 2d-vector)
                        score += count * SCORE(m_r, a, Mj + 1, d);
                    }
                }
            }
        }
        return make_tuple(score, f1);
    }

    /*
     * scoring function (eq. 3 in paper) in a form that avoids binomial
     * coefficients
     */
    static double SCORE(unsigned int r, unsigned int a, unsigned int Mj,
                        unsigned int d) {
        double prod1 = 1, prod2 = 1;
        for (size_t mj = 0; mj < Mj; mj++) {
            prod1 *= (d - mj) / (double)(r - a - mj);
        }
        for (size_t i = 0; i < a; i++) {
            prod2 *= (r - d - i) / (double)(r - i);
        }
        return prod1 * prod2;
    }

    /*
     * converts result table of Mjs and plus1_rxns into the vector that would
     * have been generated if the canonical result table (i.e. of all Mjs and
     * Js - nothing left out due to the plus1_rxns-trick) had been condensed.
     * the vector can then be used to calculate the final PoF(d0=r) via p**Mj
     * for the Mj of every superset found in recursion.
     * to avoid overflows in binom or numbers that are so large that they are
     * coerced to `inf` when printed later, stop at Mj=100 (assuming that p^100
     * is essentially 0 anyway).
     */
    template <typename T> vector<double> convert_table(const Matrix<T>& table) {
        vector<double> result(table.size() * table[0].size());
        size_t Mj_limit = 100;
        for (size_t row = 0; row < table.size(); row++) {
            size_t Mj = row + 1;
            if (Mj > Mj_limit) {
                break;
            }
            // add original result
            for (size_t plus1_rxns = 0; plus1_rxns < table[row].size();
                 plus1_rxns++) {
                double count = table[row][plus1_rxns];
                if (count == 0) {
                    continue;
                }
                result[Mj - 1] += count;
                for (size_t p1rx = 1; p1rx < plus1_rxns + 1; p1rx++) {
                    int sign = (p1rx % 2) ? -1 : 1;
                    size_t Mj_new = Mj + p1rx;
                    if (Mj_new > Mj_limit) {
                        break;
                    }
                    try {
                        double count_new =
                            binom(static_cast<double>(plus1_rxns),
                                  static_cast<double>(p1rx));
                        result[Mj_new - 1] += sign * count_new * count;
                    } catch (std::overflow_error) {
                        break;
                    }
                }
            }
        }
        return result;
    }

    /*
     * get final PoF according to eq. 5 of the paper (i.e. d --> r). The
     * result is not 100% accurate, however, and can overestimate the PoF
     * slightly. The result of F(d=d0), on the other hand, represents a
     * definite lower bound. Can also return a string representation of
     * the polynomial.
     */
    template <typename T>
    tuple<double, string> get_final_PoF(const vector<T>& Mjs, double p,
                                        bool get_poly = false) {
        double final_PoF = 0;
        stringstream poly_buf;
        string polynomial;
        for (size_t i = 0; i < Mjs.size(); i++) {
            T count = Mjs[i];
            if (count == 0) {
                continue;
            }
            T Mj = i + 1;
            final_PoF += count * pow(p, Mj);
            if (get_poly) {
                if (i == 0) {
                    poly_buf << count << "p";
                } else {
                    poly_buf << ((count < 0) ? " - " : " + ") << abs(count)
                             << "p^" << Mj;
                }
            }
        }
        if (get_poly) {
            poly_buf << endl;
            polynomial = poly_buf.str();
            if (polynomial.substr(0, 3) == " + ") {
                polynomial.erase(0, 3);
            }
        }
        return make_tuple(final_PoF, polynomial);
    }

    /*
     * print results for d --> d0 and then the result of the polynomial
     * (according to the final equation in the paper). for d <= d0 the
     * actual F(d) is used. as F(d) cannot decrease with increasing d,
     * it can be carried over to d > d0 for d --> r. additionally, F1(d)
     * (i.e. the contribution of the essential reactions) can be calculated
     * without recursion up to d = r. to obtain the highest-possible lower
     * bound for every d > d0 the larger one of F(d0) or F1(d) is selected.
     */
    void print_results(double p, unsigned int dm = 0, bool print_poly = false) {
        double score, f1, weight, weighted_score,
            acc_weighted_score = 0, found_CS, possible_CS, error = 0;
        // initialize table to print results
        Table table{
            {"d", "weight", "F(d)", "weighted F(d)", "acc. weighted F(d)",
             "lethal CS", "possible CS"},
            {5, 20, 20, 20, 20, 20, 20},
            {"%.5g", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g"}};
        table.print_header();

        // get score for each d and print the corresponding numbers
        for (size_t d = 1; d <= m_r; d++) {
            if (d <= m_max_d) {
                // update score while d < d0
                tie(score, f1) = score_cd_table2(d);
            } else {
                // only update f1, but keep old score (i.e. f(d0))
                f1 = get_f1(d);
                // use f1 instead of f(d0) as soon as it's larger
                if (f1 > score) {
                    score = f1;
                }
            }
            weight = binom_dist_weight(m_r, d, p);
            if (weight < WEIGHT_LIMIT) {
                printf("\nFor d>%d the weight is smaller than the threshold "
                       "(%.2e).\n",
                       d, WEIGHT_LIMIT);
                break;
            }
            weighted_score = score * weight;
            acc_weighted_score += weighted_score;
            possible_CS =
                binom(static_cast<double>(m_r), static_cast<double>(d));
            found_CS = score * possible_CS;
            if (d > dm) {
                error += weight - weighted_score;
            }
            vector<double> numbers{
                d,        weight,     score, weighted_score, acc_weighted_score,
                found_CS, possible_CS};
            table.print_row(numbers);
        }
        cout << "\n\nFinal results:" << endl;
        cout << string(22, '-') << endl;
        printf("Iterative PoF(d0=r)\t= %.15e\t\t--> lower bound\n",
               acc_weighted_score);
        double final_PoF;
        string polynomial;
        tie(final_PoF, polynomial) =
            get_final_PoF(convert_table(m_cd_table), p, print_poly);
        printf("Polynomial PoF(d0=r)\t= %.15e\t\t--> best guess, but might "
               "overshoot\n",
               final_PoF);
        printf("It. PoF(d0=r) + error\t= %.15e\t\t--> upper bound\n",
               acc_weighted_score + error);
        cout << string(22, '-') << endl;

        if (print_poly) {
            cout << "Polynomial: " << polynomial << endl;
        }
    }
};

#endif /* POF_CALCULATOR_HPP */
