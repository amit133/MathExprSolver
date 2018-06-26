#include <iostream>
#include <cassert>
#include "xmlParser.h"

using namespace std;
using namespace tinyxml2;

KXml::KXml(string & xmlFileName) {
	xmlFileName_m = xmlFileName;
	bool fileOpened = LoadXmlFile(xmlFileName);
	if(!fileOpened) {
		cerr << "XML File load failure" << endl;
		return;
	}

	setActors();
	setVariables();
	setPolicies();
	setParameters();
	setOptimizeFunctions();
	setEquations();
}

bool KXml::LoadXmlFile(string & xmlFileName) {
	XMLError eResult = doc.LoadFile(xmlFileName_m.c_str());
	if(eResult != XML_SUCCESS) {
		cerr << "Could not open xml file: " << xmlFileName_m << endl; 
		return false;
	}
	rootElement = doc.RootElement();

	return true;
}

void KXml::setActors() {
	auto getChildText = [](XMLElement *parentElement, const char * childElementName) {
		auto childElement = parentElement->FirstChildElement(childElementName);
		assert(childElement != 0);
		auto childText = childElement->GetText();
		return childText;
	};

	auto getSiblingText = [](XMLElement *currentElement, const char * siblingElementName){
		auto *siblingElement = currentElement->NextSiblingElement(siblingElementName);
		assert(siblingElement != 0);
		auto text = siblingElement->GetText();
		assert(text != 0);
		return text;
	};

	XMLElement *actorsList = rootElement->FirstChildElement("actors");
	XMLElement *actor = actorsList->FirstChildElement("actor");
	while(actor) {
		auto *nameElement = actor->FirstChildElement("name");
		auto *actorName = nameElement->GetText();

		// auto actorName = getChildText(actor, "name");

		auto *descriptionElement = nameElement->NextSiblingElement("description");
		auto description = descriptionElement->GetText();

		auto *influenceElement = descriptionElement->NextSiblingElement("influence");
		auto influence = influenceElement->GetText();

		auto *utilityElement = influenceElement->NextSiblingElement("utility");
		auto utility = utilityElement->GetText();

		//auto *influence = actor->FirstChildElement("influence");
		//cout << actorName << ", " << description << ", Influence: " << influence << ", utility: " << utility << endl;

		actor = actor->NextSiblingElement("actor");
	}
}

void KXml::setVariables() {
	XMLElement *variables = rootElement->FirstChildElement("variables");
	auto *varElem = variables->FirstChildElement("variable");
	while(varElem) {
		auto *varNameElem = varElem->FirstChildElement("name");
		assert (varNameElem != nullptr);

		auto varName = varNameElem->GetText();
		assert(varName != 0);
		listOfVars.push_back(string(varName));

		auto varInitValElem = varElem->FirstChildElement("initial");
		auto varInitVal = varInitValElem->GetText();

		// If an initial value is not there, the formula to compute the initial value should be present
		if(varInitVal == 0) {
			auto varInitValComputeElem = varElem->FirstChildElement("compute");
			assert(varInitValComputeElem != 0);

			auto varInitValComputeExpr = varInitValComputeElem->GetText();
			//cout << varName << ": " << varInitValComputeExpr << endl;

			initValExpressions[listOfVars.back()] = string(varInitValComputeExpr);
		} else {
			//cout << varName << ": " << varInitVal << endl;
			initValExpressions[listOfVars.back()] = string(varInitVal);
		}

		// Go to the next variable
		varElem = varElem->NextSiblingElement();
	}

	// cout << "Expressions for init values: ";
	// for(auto var : listOfVars) {
	// 	cout << var << ": " << initValExpressions.at(var) << endl;
	// }
	// cout  << endl << endl;
}

void KXml::setPolicies() {
	XMLElement *policies = rootElement->FirstChildElement("policy");
	auto *policyElem = policies->FirstChildElement();
	while(policyElem) {
		auto *policyNameElem = policyElem->FirstChildElement("name");
		assert (policyNameElem != nullptr);

		auto policyName = policyNameElem->GetText();

		auto policyBaseElem = policyElem->FirstChildElement("base");
		auto policyBase = policyBaseElem->GetText();

		string policyType = policyElem->Name();
		//cout << "policy Name " << policyElem->Name() << endl;
		if(policyType == "variable") {
			//cout << "policy variable, " << policyName << ": " << policyBase << endl;
		} else if(policyType == "constant") {
			//cout << "policy constant, " << policyName << ": " << policyBase << endl;
		}

		// Go to the next policy
		policyElem = policyElem->NextSiblingElement();
	}
}

void KXml::setParameters() {
	XMLElement *parameters = rootElement->FirstChildElement("parameters");
	XMLElement *paramElement = parameters->FirstChildElement("parameter");
	while(paramElement) {
		auto paramNameElem = paramElement->FirstChildElement("name");
		auto paramName = paramNameElem->GetText();

		auto paramValueElem = paramElement->FirstChildElement("value");
		auto paramValue = paramValueElem->DoubleText();

		params[paramName] = paramValue;

		//cout << paramName << ": " << paramValue << endl;

		paramElement = paramElement->NextSiblingElement("parameter");
	}
}

KXml::parameters KXml::getParameters() {
	return params;
}

void KXml::setOptimizeFunctions() {
	XMLElement *optimizeFunctions = rootElement->FirstChildElement("optimize");
	XMLElement *optimizeEquation = optimizeFunctions->FirstChildElement("equation");
	while(optimizeEquation) {
		auto optimizeEquationElem = optimizeEquation->FirstChildElement("function");
		auto optimizeFunction = optimizeEquationElem->GetText();
		//cout << "optimize function: " << optimizeFunction << endl;
		optimizeEquation = optimizeEquation->NextSiblingElement();
	}
}

void KXml::setEquations() {
	XMLElement *equations = rootElement->FirstChildElement("equations");
	XMLElement *equationElement = equations->FirstChildElement("equation");
	while(equationElement) {
		auto equationNameElem = equationElement->FirstChildElement("name");
		auto equationName = equationNameElem->GetText();

		auto equationFunctionElem = equationElement->FirstChildElement("function");
		auto equationFunction = equationFunctionElem->GetText();

		//cout << equationName << ": " << equationFunction << endl;

		equationElement = equationElement->NextSiblingElement();
	}
}

KXml::listOfVariables KXml::getVariables() {
	return listOfVars;
}

KXml::initValExprMap KXml::getInitComputeExpressions() {
	return initValExpressions;
}


/* int main() {
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
 */