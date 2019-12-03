#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <map>
using namespace std;

template <typename T> using Matrix = vector<vector<T>>;
typedef map<pair<size_t, size_t>, int> Counter;

#endif /* TYPES_H */
