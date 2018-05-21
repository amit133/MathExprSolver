#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

const int NUMBER = 0, OPERATION = 1, BRACKET = 2, WHITESPACE =3, OTHER = 4;

map<string, string> equationParams;
map<string, string>::iterator paramsItr;

int get_type(char ch) {
    if((ch >= '0' && ch <= '9') || ch == '.') { return NUMBER; }

    // OPs + - * / ^
    if(ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^' ) { return OPERATION; }

    // bracket
    if(ch == '(') { return BRACKET; }

    // whitespace
    if(ch == ' ') { return WHITESPACE; }
    return OTHER;
}

void clearConsecutiveWhitespaces() {

}

int main() {
    equationParams["govUtilW1"] = "0.0";
    equationParams["targRNW"] = "9.5";
    equationParams["integcost"] = "0.5";
    equationParams["gamma"] = "1.0";
    equationParams["alpha"] = "0.32";

    string equ = "PRNW-   23/gamma/ 44+alpha*0.45*  POSUB";
    //cout << "Enter an expression" << endl;
    //cin >> equ;

    cout << "Expression to be parsed: " << equ << endl;
    size_t pos = 0;
    size_t len = equ.length();
    cout << "Size of string: " << len << endl;
    char ch;// = equ[pos];
    string const_buffer;
    string var_buffer;
    string other;
    size_t whiteSpaceCount = 0;
    while(pos < len) {
        ch = equ[pos];
        int type = get_type(ch);
        //cout << ch << endl;
        switch(type) {
            case WHITESPACE:
                whiteSpaceCount = 1;
                ++pos;
                // Find and remove all continuous whitespaces
                while(get_type(equ[pos]) == WHITESPACE) {
                    ++whiteSpaceCount;
                    ++pos;
                }

                pos = pos-whiteSpaceCount;

                // Erase all whitespaces found in this iteration
                equ = equ.erase(pos, whiteSpaceCount);

                // Reset the white space counter
                whiteSpaceCount = 0;

                //cout << "Equation after erasing whitespace: " << equ << endl;
                //cout << "Pos = " << pos << " Char at pos: " << equ[pos] << endl;

            break;

            case NUMBER:
            break;

            case OPERATION:
                cout << "An operator found: " << ch << endl;

                // Omit any whitespace if present.
                ++pos;

                // Start with zero white space count
                whiteSpaceCount = 0;
                while(get_type(equ[pos]) == WHITESPACE) {
                    ++whiteSpaceCount;
                    ++pos;
                }

                pos = pos-whiteSpaceCount-1;

                // Erase all whitespaces found in this iteration
                equ = equ.erase(pos+1, whiteSpaceCount);

                cout << "Equation after erasing whitespace: " << equ << endl;
                cout << "Pos = " << pos << " Char at pos: " << equ[pos] << endl;
                // Reset the white space counter
                whiteSpaceCount = 0;

                // Reset the length of the string
                len = equ.length();

                // Check if consecutive operators are present in the expression.
                if(get_type(equ[pos+1]) == OPERATION) {
                    cout << "Invalid expression as back to back to operator is not allowed" << endl;
                    exit(-1);
                }

                // Write logic if a param or variable name is found
                if (other.size() > 0) {
                    paramsItr = equationParams.find(other);
                    if(paramsItr != equationParams.end()) {
                        // Replace the param name with the actual value in the expression string
                        //cout << equ.substr(0, pos - other.size()) << paramsItr->second << equ.substr(pos) << endl;
                        
                        equ = equ.substr(0, pos - other.size()) + paramsItr->second + equ.substr(pos);
                        // Reset the length of the string
                        len = equ.length();
                        cout << "New equation: " << equ << endl;
                        pos = pos - other.size() + (paramsItr->second).size();

                    }
                    other.clear();
                }
            break;

            case BRACKET:
            break;

            case OTHER:
                other += ch;
                //cout << "other: " << other << endl;
            break;
            
        }
        //cout << ch << endl;
        ++pos;
    }

    cout << "Equation after parse: " << equ << endl;
    return 0;
}