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

typedef std::vector<std::string> variables;
variables vars;// = {"OE","E","G","RNW","PE","S","OS","PS","N","KP","SF","W","R","SH","C","RSTAR","TR","PO","B","O","KG","Y","POSUB","PG","PRNW"};

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

void setInitValuesOfVars() {
	KXml xmlParser(xmlFile);
	listOfVars = xmlParser.getVariables();
	auto initValExpressions = xmlParser.getInitComputeExpressions();

    // The map would be populated in a loop when xml reading comes into play
    for(auto var : listOfVars) {
        initValuesOfVars[var] = 0.0;
    }

    symbol_table_t symbol_table_initvals;
    symbol_table_initvals.add_constants();

    typedef exprtk::details::variable_node<double> exprtk_var_ptr;
 
    for(auto param : params) {
        symbol_table_initvals.add_constant(param.first, param.second);
        exprtk_var_ptr *v = symbol_table_initvals.get_variable(param.first);
        //cout << "Symbol name: " << param.first << ", symbol value: " << v->value() << endl;
    }

    for(auto policy : policyWeights) {
        symbol_table_initvals.add_constant(policy.first, policy.second);
        exprtk_var_ptr *v = symbol_table_initvals.get_variable(policy.first);
        //cout << "Symbol name: " << policy.first << ", symbol value: " << v->value() << endl;
    }
    
    // typedef exprtk::parser<double> parser_t;
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
    assert(!vars.empty());

    string initVals ("[");
    for(auto var : listOfVars) { // This option of for loop has to rechecked
    //for(auto var : vars) {
        //initVals += std::to_string(initValuesOfVars[varName + string("0")]) + ","; // This code looks more apt instead of that in try-catch block
        try {
            // initVals += std::to_string(initValuesOfVars.at(var + string("0"))) + ",";
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
    // Each parameter's name and value would come from the xml
    string alpha        =    "alpha";
    string beta         =    "beta";
    string gamma        =    "gamma";
    string a            =    "a";
    string lambda       =    "lambda";
    string niu          =    "niu";
    string theta        =    "theta";
    string b            =    "b";
    string d            =    "d";
    string sigmac       =    "sigmac";
    string omega        =    "omega";
    string sigma        =    "sigma";
    string tau          =    "tau";
    string phi          =    "phi";
    string betadisc     =    "betadisc";
    string delta        =    "delta";
    string nu           =    "nu";
    string poaverage    =    "poaverage";
    string A            =    "A";
    string rho          =    "rho";
    string RBAR         =    "RBAR";
    string bss          =    "bss";
    string psi          =    "psi";
    string gexg         =    "gexg";
    string oexg         =    "oexg";
    string posubexg     =    "posubexg";
    string integcost    =    "integcost";
    string igexg        =    "igexg";

    string alpha_Value        =    "0.32";
    string beta_Value         =    "0.42";
    string gamma_Value        =    "1";
    string a_Value            =    "0.33230";
    string lambda_Value       =    "-0.25786";
    string niu_Value          =    "-1.381";
    string theta_Value        =    "0.58221";
    string b_Value            =    "0.00003";
    string d_Value            =    "0.0028";
    string sigmac_Value       =    "-0.33333"; //"-1.0/3.0"; // *****Note 1*****: 1/3 results into 0 because of integer division rule. Also
                                               // *****Note 2*****: 1.0/3.0 has to be evaluated. Is it possible to put a constraint to provide a final value instead of expression
                                               // It would make programming less complex. Else every value has to be checked and evaluated before the actual use.
    string omega_Value        =    "0.0735";
    string sigma_Value        =    "2";
    string tau_Value          =    "1";
    string phi_Value          =    "0";
    string betadisc_Value     =    "0.96";
    string delta_Value        =    "0.1";
    string nu_Value           =    "0.05";
    string poaverage_Value    =    "1.6371";
    string A_Value            =    "0.0328";
    string rho_Value          =    "0.9";
    string RBAR_Value         =    "0.04";
    string bss_Value          =    "0.16";
    string psi_Value          =    "0";
    string gexg_Value         =    "0.01398";
    string oexg_Value         =    "0.4819";
    string posubexg_Value     =    "0.5";
    string integcost_Value    =    "0.5";

    params.insert(parameters::value_type( alpha    , std::stod( alpha_Value     )));
    params.insert(parameters::value_type( beta     , std::stod( beta_Value      )));
    params.insert(parameters::value_type( gamma    , std::stod( gamma_Value     )));
    params.insert(parameters::value_type( a        , std::stod( a_Value         )));
    params.insert(parameters::value_type( lambda   , std::stod( lambda_Value    )));
    params.insert(parameters::value_type( niu      , std::stod( niu_Value       )));
    params.insert(parameters::value_type( theta    , std::stod( theta_Value     )));
    params.insert(parameters::value_type( b        , std::stod( b_Value         )));
    params.insert(parameters::value_type( d        , std::stod( d_Value         )));
    params.insert(parameters::value_type( sigmac   , std::stod( sigmac_Value    )));
    params.insert(parameters::value_type( omega    , std::stod( omega_Value     )));
    params.insert(parameters::value_type( sigma    , std::stod( sigma_Value     )));
    params.insert(parameters::value_type( tau      , std::stod( tau_Value       )));
    params.insert(parameters::value_type( phi      , std::stod( phi_Value       )));
    params.insert(parameters::value_type( betadisc , std::stod( betadisc_Value  )));
    params.insert(parameters::value_type( delta    , std::stod( delta_Value     )));
    params.insert(parameters::value_type( nu       , std::stod( nu_Value        )));
    params.insert(parameters::value_type( poaverage, std::stod( poaverage_Value )));
    params.insert(parameters::value_type( A        , std::stod( A_Value         )));
    params.insert(parameters::value_type( rho      , std::stod( rho_Value       )));
    params.insert(parameters::value_type( RBAR     , std::stod( RBAR_Value      )));
    params.insert(parameters::value_type( bss      , std::stod( bss_Value       )));
    params.insert(parameters::value_type( psi      , std::stod( psi_Value       )));
    params.insert(parameters::value_type( gexg     , std::stod( gexg_Value      )));
    params.insert(parameters::value_type( oexg     , std::stod( oexg_Value      )));
    params.insert(parameters::value_type( posubexg , std::stod( posubexg_Value  )));
    params.insert(parameters::value_type( integcost, std::stod( integcost_Value )));

    for(auto param : params) {
        symbol_table.add_constant(param.first, param.second);
    }
}

void setVarNames() {
    // Variables' names would be read from the xml file. The count of these variables would be unknown.
    std::string OE = "OE";
    std::string E = "E";
    std::string G = "G";
    std::string RNW = "RNW";
    std::string PE = "PE";
    std::string POSUB = "POSUB";
    std::string PG = "PG";
    std::string PRNW = "PRNW";
    std::string S = "S";
    std::string OS = "OS";
    std::string PS = "PS";
    std::string N = "N";
    std::string KP = "KP";
    std::string SF = "SF";
    std::string W = "W";
    std::string R = "R";
    std::string SH = "SH";
    std::string C = "C";
    std::string RSTAR = "RSTAR";
    std::string TR = "TR";
    std::string PO = "PO";
    std::string B = "B";
    std::string O = "O";
    std::string KG = "KG";
    std::string Y = "Y";

    // The push_back() calls would run in a loop when xml reading comes into play
    //OE0,N0,KP0,C0,E0,G0,RNW0,PE0,OS0,S0,PS0,SF0,W0,R0,SH0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0,POSUB0,PG0,PRNW0
    vars.push_back(OE);
    vars.push_back(E);
    vars.push_back(G);
    vars.push_back(RNW);
    vars.push_back(PE);
    vars.push_back(S);
    vars.push_back(OS);
    vars.push_back(PS);
    vars.push_back(N);
    vars.push_back(KP);
    vars.push_back(SF);
    vars.push_back(W);
    vars.push_back(R);
    vars.push_back(SH);
    vars.push_back(C);
    vars.push_back(RSTAR);
    vars.push_back(TR);
    vars.push_back(PO);
    vars.push_back(B);
    vars.push_back(O);
    vars.push_back(KG);
    vars.push_back(Y);

    // Policy variables. Note that the sequence is different from that in the matlab code
    vars.push_back(POSUB);
    vars.push_back(PG);
    vars.push_back(PRNW);
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

    // z.resize(vars.size());

    // for(int i=0 ; i < vars.size() ; ++i) {
    //     symbol_table.add_variable(vars[i],     z[i]);
    // }
 
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
    
    setVarNames();
   
    // Read all parameter values
    setParameters();

    setPolicyWeights();

    setInitValuesOfVars();

    string initValuesOfVars = getInitValues();
    //cout << "Init Vals: " << initValuesOfVars << endl;

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