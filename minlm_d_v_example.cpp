#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "exprtk.hpp"
#include <vector>
#include<cassert>
using namespace std;

using namespace alglib;

    typedef std::string paramName;
    typedef double paramValue;
    typedef std::map<paramName, paramValue> parameters;
    parameters params;

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
    string igexg_Value        =    "0";

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
    params.insert(parameters::value_type( igexg    , std::stod( igexg_Value     )));
}

void  function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{

    //
    // this callback calculates the values of expressions as per the inputs passed by optimization algorithm
    //

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

    typedef std::vector<std::string> variables;
    variables vars;// = {"OE","E","G","RNW","PE","POSUB","PG","PRNW","S","OS","PS","N","KP","SF","W","R","SH","C","RSTAR","TR","PO","B","O","KG","Y"};

    // The push_back() calls would run in a loop when xml reading comes into play
    vars.push_back(OE);
    vars.push_back(E);
    vars.push_back(G);
    vars.push_back(RNW);
    vars.push_back(PE);
    vars.push_back(POSUB);
    vars.push_back(PG);
    vars.push_back(PRNW);
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

    // Confirm that the length of the vector vars is equal to the length of input real_1d_array x
    assert(vars.size() == x.length());

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>     expression_t;
    typedef exprtk::parser<double>             parser_t;

    //OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0
    symbol_table_t symbol_table;

    // A temporary vector z is required because input parameter x is a const which causes error in compilation
    // Size of this vector needs to be fixed. Keeping it dynamic causes runtime trouble (neither error nor crash) with symbol table of exprtk library.
    // I think it is due to change in location of vector in memory if the vector undergoes a change in its size due to push back operation.
    vector<double> z(x.length(), 0.0);

    for(int i=0 ; i < x.length() /* or vars.size() */; ++i) {
        // z.push_back(x[i]); // push_back() can't be used because symbol_table.add_variable() takes a reference to double type of variable
        z[i] = x[i];
        symbol_table.add_variable(vars[i],     z[i]);
    }
 
    parameters::iterator paramsItr = params.begin();
    parameters::iterator paramsEndItr = params.end();
    for(; paramsItr != paramsEndItr; ++paramsItr) {
        symbol_table.add_constant(paramsItr->first, paramsItr->second);
    }

    symbol_table.add_constants();

    parser_t parser;
    std::vector<expression_t> expressions;
    size_t equationsCount =  systemOfEquations.size();

    for(int i = 0; i < equationsCount /* Note: count of equations under <optimize> tag of xml*/; ++i) {
        expressions.push_back(expression_t());
        expressions[i].register_symbol_table(symbol_table);
        parser.compile(systemOfEquations[i] , expressions[i]);
    }

    // expression_t expression;
    // expression.register_symbol_table(symbol_table);

    for(int i = 0; i < equationsCount /* Note: count of equations under <optimize> tag of xml*/; ++i) {
        fi[i] = expressions[i].value();
    }

    for(int k = 0; k < 25; ++k) {
        cout << "fi[" << k <<  "] = " << fi[k] << endl;
    }
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
    
    //OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0
    real_1d_array x = "[0.0305, 0.0156, 0.0140, 0, 1.5625, 0.5000, 0.6563, 1.5625, 0.0389, 0.0673, 1.4932, 1.0000, \
    5.6246, 0.0262, 1.1806, 0.1400, 0.0127, 0.9363, 0.0400, 0.6868, 1.6371, -8.1109, 0.4819, 0, 2.0277]";
    double epsx = 0.0000000001;
    ae_int_t maxits = 0;
    minlmstate state;
    minlmreport rep;

    minlmcreatev(25, 25, x, 0.000001, state);
    minlmsetcond(state, epsx, maxits);

    alglib::minlmoptimize(state, function1_fvec);

    minlmresults(state, x, rep);

    cout << "x.length " << x.length() << endl;
    for(int index=0; index < x.length(); ++index) {
        printf("%.4f\t", x[index]);
    }

    return 0;
}