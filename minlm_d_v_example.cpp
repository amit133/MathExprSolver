#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "exprtk.hpp"
#include <vector>
using namespace std;

using namespace alglib;
void  function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
    //
    // this callback calculates
    // f0(x0,x1) = 100*(x0+3)^4,
    // f1(x0,x1) = (x1-3)^4
    //

    //fi[0] = 10*pow(x[0]+3,2);
    //fi[1] = pow(x[1]-3,2);
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>     expression_t;
    typedef exprtk::parser<double>             parser_t;

    //std::string expression_string1 = "10*pow(x0+3,2)";
    //std::string expression_string2 = "pow(x1-3,2)";


    std::string expression_string_0="Y-N^0.58221*((1-0.00003)*KP^-1.381+0.00003*SF^-1.381)^((1-0.58221)/-1.381)";                                   
    std::string expression_string_1="E-0.32*OE-0.42*G-1*RNW";
    std::string expression_string_2="PE-1/0.32*POSUB";
    std::string expression_string_3="PG-0.42/0.32*POSUB";
    std::string expression_string_4="PRNW-1/0.32*POSUB";
    std::string expression_string_5="S-(0.33230*E^-0.25786+(1-0.33230)*OS^-0.25786)^(1/-0.25786)";
    std::string expression_string_6="PS-POSUB*(OS^(1--0.25786))/(1-0.33230)*(0.33230*E^-0.25786+(1-0.33230)*OS^-0.25786)^((-0.25786-1)/-0.25786)";
    std::string expression_string_7="PS-PE*(E^(1--0.25786))/0.33230*(0.33230*E^-0.25786+(1-0.33230)*OS^-0.25786)^((-0.25786-1)/-0.25786)";
    std::string expression_string_8="W-0.58221*N^(0.58221-1)*((1-0.00003)*KP^-1.381+0.00003*SF^-1.381)^((1-0.58221)/-1.381)";
    std::string expression_string_9="R-(1-0.58221)*(1-0.00003)*KP^(-1.381-1)*N^0.58221*((1-0.00003)*KP^-1.381+0.00003*SF^-1.381)^((1-0.58221)/-1.381-1)";
    std::string expression_string_10="PS-(1-0.58221)*0.00003*SF^(-1.381-1)*N^0.58221*((1-0.00003)*KP^-1.381+0.00003*SF^-1.381)^((1-0.58221)/-1.381-1)";
    std::string expression_string_11="N-1";
    std::string expression_string_12="PS-0.0028*(SH/C)^((-1/3)-1)";
    std::string expression_string_13="(B-(1+RSTAR)*B)-0.16*(Y+PO*(O-OE-OS))";
    std::string expression_string_14="1-0.96*(1-0.1+R)";
    std::string expression_string_15="B+KP-(1-0.1)*KP+C+PS*SH-W*N-R*KP-(1+RSTAR)*B-TR";
    std::string expression_string_16="PRNW*RNW+POSUB*(OS+OE)+PO*(O-OS-OE)+PG*G-TR-0";
    std::string expression_string_17="RNW-.0328/(1+0.5*RNW/E)*KG";
    std::string expression_string_18="S-SF-SH";
    std::string expression_string_19="PO-1.6371";
    std::string expression_string_20="RSTAR-0.04";
    std::string expression_string_21="0.05*KG-0";
    std::string expression_string_22="POSUB-0.5";
    std::string expression_string_23="O-0.4819";
    std::string expression_string_24="G-0.01398";


    vector<double> z;
    for(int i=0 ; i < x.length(); ++i) {
        z.push_back(x[i]);
    }

    //double x0 = x[0];
    //double x1 = x[1];
    
    //OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0
    symbol_table_t symbol_table;
    symbol_table.add_variable("OE",    z[0]);
    symbol_table.add_variable("E",     z[1]);
    symbol_table.add_variable("G",     z[2]);
    symbol_table.add_variable("RNW",   z[3]);
    symbol_table.add_variable("PE",    z[4]);
    symbol_table.add_variable("POSUB", z[5]);
    symbol_table.add_variable("PG",    z[6]);
    symbol_table.add_variable("PRNW",  z[7]);
    symbol_table.add_variable("S",     z[8]);
    symbol_table.add_variable("OS",    z[9]);
    symbol_table.add_variable("PS",    z[10]);
    symbol_table.add_variable("N",     z[11]);
    symbol_table.add_variable("KP",    z[12]);
    symbol_table.add_variable("SF",    z[13]);
    symbol_table.add_variable("W",     z[14]);
    symbol_table.add_variable("R",     z[15]);
    symbol_table.add_variable("SH",    z[16]);
    symbol_table.add_variable("C",     z[17]);
    symbol_table.add_variable("RSTAR", z[18]);
    symbol_table.add_variable("TR",    z[19]);
    symbol_table.add_variable("PO",    z[20]);
    symbol_table.add_variable("B",     z[21]);
    symbol_table.add_variable("O",     z[22]);
    symbol_table.add_variable("KG",    z[23]);
    symbol_table.add_variable("Y",     z[24]);

    symbol_table.add_constants();

    expression_t expression0;
    expression_t expression1;
    expression_t expression2;
    expression_t expression3;
    expression_t expression4;
    expression_t expression5;
    expression_t expression6;
    expression_t expression7;
    expression_t expression8;
    expression_t expression9;
    expression_t expression10;
    expression_t expression11;
    expression_t expression12;
    expression_t expression13;
    expression_t expression14;
    expression_t expression15;
    expression_t expression16;
    expression_t expression17;
    expression_t expression18;
    expression_t expression19;
    expression_t expression20;
    expression_t expression21;
    expression_t expression22;
    expression_t expression23;
    expression_t expression24;	

	expression0.register_symbol_table(symbol_table);
    expression1.register_symbol_table(symbol_table);
    expression2.register_symbol_table(symbol_table);
    expression3.register_symbol_table(symbol_table);
    expression4.register_symbol_table(symbol_table);
    expression5.register_symbol_table(symbol_table);
    expression6.register_symbol_table(symbol_table);
    expression7.register_symbol_table(symbol_table);
    expression8.register_symbol_table(symbol_table);
    expression9.register_symbol_table(symbol_table);
    expression10.register_symbol_table(symbol_table);
    expression11.register_symbol_table(symbol_table);
    expression12.register_symbol_table(symbol_table);
    expression13.register_symbol_table(symbol_table);
    expression14.register_symbol_table(symbol_table);
    expression15.register_symbol_table(symbol_table);
    expression16.register_symbol_table(symbol_table);
    expression17.register_symbol_table(symbol_table);
    expression18.register_symbol_table(symbol_table);
    expression19.register_symbol_table(symbol_table);
    expression20.register_symbol_table(symbol_table);
    expression21.register_symbol_table(symbol_table);
    expression22.register_symbol_table(symbol_table);
    expression23.register_symbol_table(symbol_table);
    expression24.register_symbol_table(symbol_table);


    parser_t parser;
    parser.compile(expression_string_0 ,expression0);
    parser.compile(expression_string_1 ,expression1);
    parser.compile(expression_string_2 ,expression2);
    parser.compile(expression_string_3 ,expression3);
    parser.compile(expression_string_4 ,expression4);
    parser.compile(expression_string_5 ,expression5);
    parser.compile(expression_string_6 ,expression6);
    parser.compile(expression_string_7 ,expression7);
    parser.compile(expression_string_8 ,expression8);
    parser.compile(expression_string_9  ,expression9);
    parser.compile(expression_string_10 ,expression10);
    parser.compile(expression_string_11 ,expression11);
    parser.compile(expression_string_12 ,expression12);
    parser.compile(expression_string_13 ,expression13);
    parser.compile(expression_string_14 ,expression14);
    parser.compile(expression_string_15 ,expression15);
    parser.compile(expression_string_16 ,expression16);
    parser.compile(expression_string_17 ,expression17);
    parser.compile(expression_string_18 ,expression18);
    parser.compile(expression_string_19 ,expression19);
    parser.compile(expression_string_20 ,expression20);
    parser.compile(expression_string_21 ,expression21);
    parser.compile(expression_string_22 ,expression22);
    parser.compile(expression_string_23 ,expression23);
    parser.compile(expression_string_24 ,expression24);

    fi[0 ] = expression0.value();
    fi[1 ] = expression1.value();
    fi[2 ] = expression2.value();
    fi[3 ] = expression3.value();
    fi[4 ] = expression4.value();
    fi[5 ] = expression5.value();
    fi[6 ] = expression6.value();
    fi[7 ] = expression7.value();
    fi[8 ] = expression8.value();
    fi[9 ] = expression9.value();
    fi[10] = expression10.value();
    fi[11] = expression11.value();
    fi[12] = expression12.value();
    fi[13] = expression13.value();
    fi[14] = expression14.value();
    fi[15] = expression15.value();
    fi[16] = expression16.value();
    fi[17] = expression17.value();
    fi[18] = expression18.value();
    fi[19] = expression19.value();
    fi[20] = expression20.value();
    fi[21] = expression21.value();
    fi[22] = expression22.value();
    fi[23] = expression23.value();
    fi[24] = expression24.value();
    for(int k = 0; k < 25; ++k) {
        cout << "fi[" << k <<  "] = " << fi[k] << endl;
    }

    cout << endl << endl;
}

int main(int argc, char **argv)
{
    //
    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
    //
    //     f0(x0,x1) = 10*(x0+3)^2
    //     f1(x0,x1) = (x1-3)^2
    //
    // using "V" mode of the Levenberg-Marquardt optimizer.
    //
    // Optimization algorithm uses:
    // * function vector f[] = {f1,f2}
    //
    // No other information (Jacobian, gradient, etc.) is needed.
    //
    
    //OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0
    //real_1d_array x = "[0.030500,0.015632,0.013980,0.000000,1.562500,0.500000,0.656300,1.562500,0.038908,0.067350,1.493231,1.000000, 
    // 5.624600,0.026229,1.180557,0.140000,0.012680,0.936300,0.040000,0.686829,1.637100,-8.110867,0.481900,0.000000,2.027717]";
    //real_1d_array x = "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]";
    real_1d_array x = "[0.0305, 0.0156, 0.0140, 0, 1.5625, 0.5000, 0.6563, 1.5625, 0.0389, 0.0673, 1.4932, 1.0000, \
    5.6246, 0.0262, 1.1806, 0.1400, 0.0127, 0.9363, 0.0400, 0.6868, 1.6371, -8.1109, 0.4819, 0, 2.0277]";
    double epsx = 0.0000000001;
    ae_int_t maxits = 0;
    minlmstate state;
    minlmreport rep;

    //minlmcreatev(25, x, 0.000001, state);
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