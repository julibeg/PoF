#include "cxxopts.hpp"
#include <iostream>
using namespace std;

cxxopts::Options options("Test", "bla bli blu die kuh macht muh");

int main(int argc, char* argv[]) {
    options.add_options()("t,threads", "number of threads",
                          cxxopts::value<int>())("f,file", "input file name",
                                                 cxxopts::value<string>());

    auto result = options.parse(argc, argv);

    int threads;
    string file;
    if (result.count("threads")) {
        threads = result["threads"].as<int>();
    } else {
        threads = 1;
        cout << "1 thread selected" << endl;
    }
    if (result.count("file")){
        string file = result["file"].as<string>();
    } else {
        cout << "please provide file name" << endl;
        cout << options.help() << endl;
        exit(0);
    }

    cout << threads << " threads" << endl;
    cout << "filename: " << file << endl;
}
