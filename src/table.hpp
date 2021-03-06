#ifndef TABLE_HPP
#define TABLE_HPP

#include <iostream>
#include <string>
#include <vector>
using namespace std;

/**
 * class used to print formatted tables for debugging
 */
class Table {
  public:
    vector<string> header;
    vector<unsigned int> field_widths;
    vector<string> number_formats;

    Table(vector<string> h, vector<unsigned int> fw, vector<string> nf)
        : header(h), field_widths(fw), number_formats(nf) {
    }

    void print_border_row() {
        cout << string(m_line_length, '-') << endl;
    }

    void print_header() {
        m_line_length = 0;
        for (auto const& fw : field_widths) {
            m_line_length += fw + 1;
        }
        for (size_t i = 0; i < header.size(); i++) {
            print_centered_text(header[i], field_widths[i]);
            if (i < header.size() - 1) { // no border at the end
                cout << "|";
            }
        }
        cout << endl;
        print_border_row();
    }

    template <typename T> void print_row(vector<T> numbers) {
        for (size_t i = 0; i < numbers.size(); i++) {
            print_centered_number(numbers[i], field_widths[i],
                                  number_formats[i]);
            if (i < header.size() - 1) { // no border at the end
                cout << "|";
            }
        }
        cout << endl;
    }

  private:
    unsigned int m_line_length;

    void print_centered_text(string text, unsigned int fieldWidth) {
        unsigned int padLen_left = (fieldWidth - text.size()) / 2;
        unsigned int padLen_right = padLen_left;
        if (2 * padLen_left + text.size() < fieldWidth) {
            padLen_left++;
        }
        printf("%*s%s%*s", padLen_left, "", text.c_str(), padLen_right, "");
    }

    template <typename T>
    void print_centered_number(T numb, unsigned int fieldWidth,
                               string fmt = "%f") {
        char char_arr[fieldWidth + 1];
        sprintf(char_arr, fmt.c_str(), numb);
        string text(char_arr);
        print_centered_text(text, fieldWidth);
    }
};

#endif /* TABLE_HPP */
