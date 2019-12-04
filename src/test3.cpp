#include "PoF_calculator.cpp"
#include <iostream>
#include <vector>
#include <boost/math/distributions/poisson.hpp>
using namespace std;

int main() {
    VariadicTable<size_t, double, double, double, double, double, double>
        result_table({"d", "weight", "F(d)", "weighted F(d)",
        "acc. weighted F(d)", "lethal CS", "possible CS"});

}
