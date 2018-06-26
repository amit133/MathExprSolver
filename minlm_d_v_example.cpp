#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "exprtk.hpp"
#include <vector>
#include<cassert>
#include "xml/xmlParser.h"

using namespace std;

using namespace alglib;

typedef exprtk::expression<double>     expression_t;
std::vector<expression_t> expressions;

typedef exprtk::symbol_table<double> symbol_table_t;
symbol_table_t symbol_table;


typedef std::string policyName;
typedef double policyValue;
typedef std::map<policyName, policyValue> policies;
policies policyWeights;
std::vector<policyName> listOfPolicyNames;

typedef std::string paramName;
typedef double paramValue;
typedef std::map<paramName, paramValue> parameters;
parameters params;

//typedef std::vector<std::string> variables;
//variables vars;// = {"OE","E","G","RNW","PE","S","OS","PS","N","KP","SF","W","R","SH","C","RSTAR","TR","PO","B","O","KG","Y","POSUB","PG","PRNW"};

// An intermediate vector z is required because input parameter x in the callback function of optimization algorithm is a const which causes error in compilation
// Size of this vector needs to be fixed. Keeping it dynamic causes runtime trouble (neither error nor crash) with symbol table of exprtk library.
// I think it is due to change in location of vector in memory if the vector undergoes a change in its size due to push back operation.
vector<double> z; // Rename z to some meaningful name

typedef std::string variableName;
typedef double variableValue;
typedef std::map<variableName, variableValue> valuesOfVars;
valuesOfVars initValuesOfVars;

typedef std::string varName;
std::vector<varName> listOfVars;

typedef exprtk::parser<double> parser_t;

string xmlFile("sampleinput.xml"); 
KXml xmlParser(xmlFile);

void setInitValuesOfVars() {
	listOfVars = xmlParser.getVariables();
	auto initValExpressions = xmlParser.getInitComputeExpressions();

    // The map would be populated in a loop when xml reading comes into play
    for(auto var : listOfVars) {
        initValuesOfVars[var] = 0.0;
    }

    symbol_table_t symbol_table_initvals;
    symbol_table_initvals.add_constants();

    //typedef exprtk::details::variable_node<double> exprtk_var_ptr;
 
    //auto params = xmlParser.getParameters();
    for(auto param : params) {
        symbol_table_initvals.add_constant(param.first, param.second);
        //exprtk_var_ptr *v = symbol_table_initvals.get_variable(param.first);
        //cout << "Symbol name: " << param.first << ", symbol value: " << v->value() << endl;
    }

    for(auto policy : policyWeights) {
        symbol_table_initvals.add_constant(policy.first, policy.second);
        //exprtk_var_ptr *v = symbol_table_initvals.get_variable(policy.first);
        //cout << "Symbol name: " << policy.first << ", symbol value: " << v->value() << endl;
    }
    
    parser_t parser;

    expression_t expr;
    expr.register_symbol_table(symbol_table_initvals);

    for(auto var : listOfVars) {
        //cout << "Variable: " << var << ", Expression: " << initValExpressions[var] << endl;
        parser.compile(initValExpressions[var] , expr);

        // evaluate the expression and store in the map for the values of each variable
        initValuesOfVars[var] = expr.value();

        // Add each constant name with its value in the symbol table
        symbol_table_initvals.add_constant(var, initValuesOfVars.at(var));
        //cout << "Variable: " << var << " , Init Value: " << initValuesOfVars[var] << endl;
    }
}

void setPolicyWeights() {
	string POSUB0     = "POSUB0";
    string PG0        = "PG0";
    string PRNW0      = "PRNW0";
    string igexg      = "igexg";
    string laborShare = "laborShare";

	string POSUB0_value     = "0.5";
    string PG0_value        = "0.6563";
    string PRNW0_value      = "1.5625";
    string igexg_value      = "0";
    string laborShare_value = "0.5";

    // Push the variables in listOfVars. Make sure to call setInitValuesOfVars() before this function
    // Else the order of the variables would be different and the output of the program would be wrong
    listOfPolicyNames.push_back( POSUB0 );
    listOfPolicyNames.push_back( PG0    );
    listOfPolicyNames.push_back( PRNW0  );

    policyWeights.insert(policies::value_type( POSUB0     , std::stod( POSUB0_value     )));
    policyWeights.insert(policies::value_type( PG0        , std::stod( PG0_value        )));
    policyWeights.insert(policies::value_type( PRNW0      , std::stod( PRNW0_value      )));
    policyWeights.insert(policies::value_type( igexg      , std::stod( igexg_value      )));
    policyWeights.insert(policies::value_type( laborShare , std::stod( laborShare_value )));

    for(auto policy : policyWeights) {
        symbol_table.add_constant(policy.first, policy.second);
    }
}

// Make sure to call setInitValuesOfVars() and setPolicyWeights() before calling this function
std::string getInitValues() {
    assert(!listOfVars.empty());

    string initVals ("[");
    for(auto var : listOfVars) { // This option of for loop has to rechecked
        try {
            initVals += std::to_string(initValuesOfVars.at(var)) + ",";
        } catch(...) {
            // Note: The key is not present in the map. Do nothing and let the loop continue
        }
    }

    for(auto policyName : listOfPolicyNames) {
        initVals += std::to_string(policyWeights[policyName]) + ",";
    }

    initVals.back() = ']';

    //cout << "Init Vals: " << initVals << endl;
    return initVals;
}

void setParameters() {
    params = xmlParser.getParameters();
    for(auto param : params) {
        //cout << param.first << ": " << param.second << endl;
        symbol_table.add_constant(param.first, param.second);
    }
}

void setMathExpressions() {
    typedef std::string equation;
    typedef std::vector<equation> equations;
    equations systemOfEquations;

    // The push_back() calls would run in a loop when xml reading comes into play
    systemOfEquations.push_back("E-alpha*OE-beta*G-gamma*RNW");
    systemOfEquations.push_back("PE-1/alpha*POSUB");
    systemOfEquations.push_back("PG-beta/alpha*POSUB");
    systemOfEquations.push_back("PRNW-gamma/alpha*POSUB");
    systemOfEquations.push_back("S-(a*E^lambda+(1-a)*OS^lambda)^(1/lambda)");
    systemOfEquations.push_back("PS-POSUB*(OS^(1-lambda))/(1-a)*(a*E^lambda+(1-a)*OS^lambda)^((lambda-1)/lambda)");
    systemOfEquations.push_back("PS-PE*(E^(1-lambda))/a*(a*E^lambda+(1-a)*OS^lambda)^((lambda-1)/lambda)");
    systemOfEquations.push_back("W-theta*N^(theta-1)*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu)");
    systemOfEquations.push_back("R-(1-theta)*(1-b)*KP^(niu-1)*N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu-1)");
    systemOfEquations.push_back("PS-(1-theta)*b*SF^(niu-1)*N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu-1)");
    systemOfEquations.push_back("N-1");
    systemOfEquations.push_back("PS-d*(SH/C)^(sigmac-1)");
    systemOfEquations.push_back("(B-(1+RSTAR)*B)-bss*(Y+PO*(O-OE-OS))");
    systemOfEquations.push_back("1-betadisc*(1-delta+R)");
    systemOfEquations.push_back("B+KP-(1-delta)*KP+C+PS*SH-W*N-R*KP-(1+RSTAR)*B-TR");
    systemOfEquations.push_back("PRNW*RNW+POSUB*(OS+OE)+PO*(O-OS-OE)+PG*G-TR-igexg");
    systemOfEquations.push_back("RNW-A/(1+integcost*RNW/E)*KG");
    systemOfEquations.push_back("S-SF-SH");
    systemOfEquations.push_back("PO-poaverage");
    systemOfEquations.push_back("RSTAR-RBAR");
    systemOfEquations.push_back("nu*KG-igexg");
    systemOfEquations.push_back("POSUB-posubexg");
    systemOfEquations.push_back("O-oexg");
    systemOfEquations.push_back("G-gexg");
    systemOfEquations.push_back("Y-N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu)");

    z.resize(listOfVars.size() + 3);
    int i=0;
    for( ; i < z.size()-3 ; ++i) {
        cout << "z is set to: " << listOfVars[i] << " " << endl;
        symbol_table.add_variable(listOfVars[i].substr(0, listOfVars[i].length() - 1),     z[i]);
    }
    symbol_table.add_variable("POSUB",     z[i]); ++i;
    symbol_table.add_variable("PG",     z[i]); ++i;
    symbol_table.add_variable("PRNW",     z[i]);

    symbol_table.add_constants();

    // typedef exprtk::parser<double> parser_t;
    parser_t parser;

    size_t equationsCount =  systemOfEquations.size();

    for(int i = 0; i < equationsCount /* Note: count of equations under <optimize> tag of xml*/; ++i) {
        expressions.push_back(expression_t());
        expressions[i].register_symbol_table(symbol_table);
        parser.compile(systemOfEquations[i] , expressions[i]);
    }
}

void  function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{

    //
    // this callback calculates the values of expressions as per the inputs passed by optimization algorithm
    //

    // Confirm that the length of the number of expressions is equal to the length of input real_1d_array x
    assert(expressions.size() == x.length());

    // A non-const vector z is required because input parameter x is a const which causes c++ error in code compilation
    // Size of this vector needs to be fixed and equal to that of the number of variables ( which is also equal to number of equations
    // and the size of the input real_1d_array x in this callback method.
    // Keeping it dynamic causes runtime trouble (neither error nor crash) with symbol table of exprtk library.
    // I think it is due to change in location of vector in memory if the vector undergoes a change in its size due to push back operation.
    // We can't use push_back() method as the referenece to individual members of z are given to the symbol table.
    for(int i=0 ; i < x.length() ; ++i) {
        z[i] = x[i]; // Values of z are used to evaluate the expression values
    }

    size_t equationsCount =  expressions.size();
    for(size_t i = 0; i < equationsCount /* Note: count of equations under <optimize> tag of xml*/; ++i) {
        fi[i] = expressions[i].value(); // this evaluation uses the vector z which was mapped in the symbol table of expressions
        cout << "fi[" << i <<  "] = " << fi[i] << endl;
    }
    cout << endl;
}

int main(int argc, char **argv)
{
    //
    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2+...+fn^2
    //
    // using "V" mode of the Levenberg-Marquardt optimizer.
    //
    // Optimization algorithm uses:
    // * function vector f[] = {f1,f2,...,fn}
    //
    // No other information (Jacobian, gradient, etc.) is needed.
    //

    // Read all parameter values
    setParameters();

    setPolicyWeights();

    setInitValuesOfVars();

    string initValuesOfVars = getInitValues();
    cout << "Init Vals: " << initValuesOfVars << endl;

    setMathExpressions();

    //OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0
    //real_1d_array x = "[0.0305, 0.0156, 0.0140, 0, 1.5625, 0.5000, 0.6563, 1.5625, 0.0389, 0.0673, 1.4932, 1.0000, \
    //5.6246, 0.0262, 1.1806, 0.1400, 0.0127, 0.9363, 0.0400, 0.6868, 1.6371, -8.1109, 0.4819, 0, 2.0277]";

    //OE0,E0,G0,RNW0,PE0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0,POSUB0,PG0,PRNW0
    //real_1d_array x = "[0.0305, 0.0156, 0.0140, 0, 1.5625, 0.0389, 0.0673, 1.4932, 1.0000, \
    //5.6246, 0.0262, 1.1806, 0.1400, 0.0127, 0.9363, 0.0400, 0.6868, 1.6371, -8.1109, 0.4819, 0, 2.0277, 0.5000, 0.6563, 1.5625]";
                      //  "[0.030500,1.000000,5.624600,0.936300,0.015632,0.013980,0.000000,1.562500,
                      //  0.067350,0.038908,1.493231,0.026229,1.180557,0.140000,0.012680,0.040000,0.686829,
                      //  1.637100,-8.110867,0.481900,0.000000,2.027717,0.500000,0.656300,1.562500]"
    //OE0,N0,KP0,C0,E0,G0,RNW0,PE0,OS0,S0,PS0,SF0,W0,R0,SH0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0,POSUB0,PG0,PRNW0

    real_1d_array x = initValuesOfVars.c_str();
    double epsx = 0.0000000001;
    ae_int_t maxits = 0;
    minlmstate state;
    minlmreport rep;

    minlmcreatev(25, 25, x, 0.000001, state);
    minlmsetcond(state, epsx, maxits);

    alglib::minlmoptimize(state, function1_fvec);

    minlmresults(state, x, rep);

    //cout << "x.length " << x.length() << endl;
    cout << "Solution set: ";

    for(auto var: listOfVars) {
        cout << var << "  ";
    }
    cout << endl;
    for(int index=0; index < x.length(); ++index) {
        printf("%.4f\n", x[index]);
    }

    return 0;
}