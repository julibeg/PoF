#include <iostream>
#include <iomanip>
#include <ios>

using namespace std;

int main() {

    // Initializing the integer
    int num = 500000000;

    cout << "Before setting the width: \n" << num << endl;

    // Using setw()
    cout << "Setting the width"
         << " using setw to 5: \n"
         << setw(5);

    cout << num << endl;

    return 0;
}
