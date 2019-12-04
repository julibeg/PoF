#include <iostream>
#include <vector>
#include <string>
#include <tuple>
using namespace std;

void centerText(string text, unsigned int fieldWidth){
	unsigned int padLen_left = (fieldWidth - text.size()) / 2;
	unsigned int padLen_right = padLen_left;
	if (2 * padLen_left + text.size() < fieldWidth) {padLen_left++;}
	printf("%*s%s%*s\n", padLen_left, "", text.c_str(), padLen_right, "");
}


template<typename T> void centerNumber(T numb, unsigned int fieldWidth,
                                       string fmt="%f"){
	char char_arr[fieldWidth + 1];
	sprintf(char_arr, fmt.c_str(), numb);
	string text(char_arr);
	unsigned int padLen_left = (fieldWidth - text.size()) / 2;
	unsigned int padLen_right = padLen_left;
	if (2 * padLen_left + text.size() < fieldWidth) {
		padLen_left++;
	}
	printf("%*s%s%*s", padLen_left, "", text.c_str(), padLen_right, "");
}


int main(int argc, char* argv[]){
	// vector<string> bla {"seas", "a;sfdlj", "f", "asdfa", "aaaaaaaaaa"};
	// vector<double> bla {0.1234, 124.12445, 1.523e84, 0.000000000000000545};
	// for (auto b : bla){
	//     centerNumber(b, 20, "%20.10g");
	// }
	// double a = 123e50;
	// char t[20];
	// sprintf(t, "%5.5g", a);
	// printf("%s\n", t);
	// printf("%5.5g\n", a);
    int a = 1;
    int b = 10;
    tuple <int, int> bla {a, b};
    for (auto h : bla) {
        cout << h << " ";
    }
    coud << endl;
}
