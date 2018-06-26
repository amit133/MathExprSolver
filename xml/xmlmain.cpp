#include "xmlParser.h"
#include <iostream>
using namespace std;

int main() {
	string xmlFile("sampleinput.xml"); 
	KXml xmlParser(xmlFile);
	auto listOfVars = xmlParser.getVariables();
	auto initValExpressions = xmlParser.getInitComputeExpressions();

	cout << "Expressions for initial values: " << endl;
	for(auto var : listOfVars) {
		cout << var << ": " << initValExpressions.at(var) << endl;
	}

	cout << "Parameters: " << endl;
	auto params = xmlParser.getParameters();
	for(auto param : params) {
		cout << param.first << ": " << param.second << endl;
	}
	cout  << endl << endl;

	return 0;
}
