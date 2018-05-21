#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

const int NUMBER = 0, OPERATION = 1, BRACKET = 2, WHITESPACE =3, OTHER = 4;

map<string, string> equationParams;
map<string, string>::iterator isParam;

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

void clearConsecutiveWhitespaces(string& equ, size_t &pos) {
    size_t whiteSpaceCount = 0;

    // Find and remove all continuous whitespaces
    while(get_type(equ[pos + whiteSpaceCount]) == WHITESPACE) {
        ++whiteSpaceCount;
    }

    // Erase all whitespaces found in this iteration
    equ = equ.erase(pos, whiteSpaceCount);
}

int main() {
    equationParams["govUtilW1"] = "0.0";
    equationParams["targRNW"] = "9.5";
    equationParams["integcost"] = "0.5";
    equationParams["gamma"] = "1.0";
    equationParams["alpha"] = "0.32";

    string equ = "PRNW-   23/gamma/ 44+alpha*0.45*  POSUB  ";
    //cout << "Enter an expression" << endl;
    //cin >> equ;

    cout << "Expression to be parsed: " << equ << endl;
    size_t pos = 0;
    size_t len = equ.length();
    cout << "Size of string: " << len << endl;
    char ch; // = equ[pos];
    string const_buffer;
    string var_buffer;
    string other;

    while(pos < len) {
        ch = equ[pos];
        int type = get_type(ch);
        //cout << ch << endl;
        switch(type) {
            case WHITESPACE:
                // Omit any whitespace till the next char or end of the string.
                clearConsecutiveWhitespaces(equ, pos);

                // Reset the length of the equation string
                len = equ.length();

            break;

            case NUMBER:
            break;

            case OPERATION:
                cout << "An operator found: " << ch << endl;

                // Omit any whitespace till the next char or end of the string.
                ++pos;
                clearConsecutiveWhitespaces(equ, pos);

                // Reposition the index
                --pos;

                // Reset the length of the equation string
                len = equ.length();

                // Check if consecutive operators are present in the expression.
                if(get_type(equ[pos+1]) == OPERATION) {
                    cout << "Invalid expression as back to back to operator is not allowed" << endl;
                    exit(-1);
                }

                // Write logic if a param or a variable name is found
                if (other.size() > 0) {
                    isParam = equationParams.find(other);
                    if(isParam != equationParams.end()) {
                        //cout << equ.substr(0, pos - other.size()) << isParam->second << equ.substr(pos) << endl;
                        
                        // Replace the param name in the expression string with its numerical value
                        equ = equ.substr(0, pos - other.size()) + isParam->second + equ.substr(pos);

                        // Reset the length of the string
                        len = equ.length();
                        cout << "New equation: " << equ << endl;
                        pos = pos - other.size() + (isParam->second).size();
                    }

                    // Reset the buffer
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