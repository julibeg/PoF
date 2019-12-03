#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

vector<unsigned int> read_comp_rxn_file(const string &fname){
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
	vector<unsigned int> comp_rxn_counts;
	string comp_rxn;
	istringstream iss(line);
	while (getline(iss, comp_rxn, ' ')) {
		comp_rxn_counts.push_back(stoi(comp_rxn));
	}
	return comp_rxn_counts;
}


int main(){
	auto vec = read_comp_rxn_file("../example_files/"
	                              "d4.mcs.binary.compressed.num_rxns");
	for (auto v : vec) {
		cout << v << " ";
	}
	cout << endl;
}
