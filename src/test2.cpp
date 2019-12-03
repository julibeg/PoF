#include <boost/multi_array.hpp>
#include <iostream>
using namespace std;

int main() {
    int bla = 6;
    boost::multi_array<char, 1> a{boost::extents[bla]};

    a[0] = 'B';
    a[1] = 'o';
    a[2] = 'o';
    a[3] = 's';
    a[4] = 't';
    a[5] = '\0';
    for (auto it : a){
        cout << it << " ";
    }
    cout << endl;
}
