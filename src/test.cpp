#include "cutset_compressed.cpp"
#include <iostream>
#include <vector>
#include <boost/math/distributions/poisson.hpp>
using namespace std;

int main() {
    boost::math::poisson_distribution<double> dist(0.5);
double res = boost::math::pdf(dist, 5);
cout << res << endl;
}
