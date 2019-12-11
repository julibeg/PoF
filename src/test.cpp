#include <iostream>
#include "combinatorics.hpp"
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>

using namespace std;



int main(){
	unsigned int n = 82;
	unsigned int d_max = 50;
	double p = 6e-3;
	// double p = 0.5/n;

	boost::math::poisson pois(0.5);
	boost::math::binomial dist(n, p);
	// boost::math::poisson pois(p * n);

	for (size_t d = 1; d < d_max + 1; d++) {
		printf("%.18g\n", binom_dist_weight(n, d, p));
		printf("%.18g\n", boost::math::pdf(dist, d));
		printf("%.18g\n", boost::math::pdf(pois, d));
		// double diff = binom_prob(n, d, p) - boost::math::pdf(dist, d);
		cout << "===========" << endl;
	}
	return 0;
}
