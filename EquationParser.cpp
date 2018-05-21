#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

const int NUMBER = 0, OPERATION = 1, BRACKET = 2, OTHER = 3;

map<string, string> equationParams;
map<string, string>::iterator paramsItr;

int get_type(char ch) {
    if(ch >= '0' && ch <= '9') { return NUMBER;}

    // OPs + - * / ^
    if(ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^' ) { return OPERATION; }

    // bracket
    if(ch == '(') { return BRACKET;}
    return OTHER;
}

int main() {
    equationParams["govUtilW1"] = "0.0";
    equationParams["targRNW"] = "9.5";
    equationParams["integcost"] = "0.5";
    equationParams["gamma"] = "1.0";
    equationParams["alpha"] = "0.32";

    string equ = "PRNW-gamma/alpha*POSUB";
    //cout << "Enter an expression" << endl;
    //cin >> equ;

    cout << "Expression to be parsed: " << equ << endl;
    size_t pos = 0;
    size_t len = equ.length();
    cout << "Size of string: " << len << endl;
    char ch;// = equ[pos];
    string const_buffer;
    string var_buffer;
    size_t found_op = 0;
    bool found_other = false;
    string other;
    while(pos < len) {
        ch = equ[pos];
        int type = get_type(ch);
        switch(type) {
            case NUMBER:
                // reset the operation coutner
                found_op = 0;
                
            break;
            case OPERATION:
                cout << "Operation found" << endl;
                ++found_op;
                if(found_op > 1) {
                    cout << "Invalid expression as back to back to operator is not allowed" << endl;
                    exit(-1);
                }


                // Write logic if a param or variable name is found
                if (other.size() > 0) {
                    paramsItr = equationParams.find(other);
                    if(paramsItr != equationParams.end()) {
                        // Replace the param name with the actual value in the expression string
                        cout << equ.substr(0, pos - other.size()) << paramsItr->second << equ.substr(pos) << endl;
                        
                        equ = equ.substr(0, pos - other.size()) + paramsItr->second + equ.substr(pos);
                        // Reset the length of the string
                        len = equ.length();
                        cout << "New equation: " << equ << endl;
                        pos = pos - other.size() + (paramsItr->second).size();

                    }
                    found_other = false;
                    other.clear();
                }
            break;

            case BRACKET:
            break;

            case OTHER:
                other += ch;
                // reset the operation coutner
                found_op = 0;
                found_other = true;
                cout << "other: " << other << endl;
            break;
            
        }
        //cout << ch << endl;
        ++pos;
    }
    return 0;
}