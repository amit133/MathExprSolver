#include <iostream>
#include <cassert>
#include "xmlParser.h"

using namespace std;
using namespace tinyxml2;

bool KXml::LoadXmlFile(string & xmlFileName) {
	xmlFileName_m = xmlFileName;
	XMLError eResult = doc.LoadFile(xmlFileName_m.c_str());
	if(eResult != XML_SUCCESS) {
		cerr << "Could not open xml file: " << xmlFileName_m << endl; 
		return false;
	}
	rootElement = doc.RootElement();

	return true;
}

std::vector<string> KXml::getActors() {
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
		cout << actorName << ", " << description << ", Influence: " << influence << ", utility: " << utility << endl;

		actor = actor->NextSiblingElement("actor");
	}

	return vector<string>();
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

		// for(auto initvar : initValExpressions) {
		// 	cout << initvar.first << initvar.second << endl;
		// }

		// Go to the next variable
		varElem = varElem->NextSiblingElement();
	}

	cout << "Expressions for init values: ";
	for(auto var : listOfVars) {
		cout << var << ": " << initValExpressions.at(var) << endl;
	}
	cout  << endl << endl;
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
			cout << "policy variable, " << policyName << ": " << policyBase << endl;
		} else if(policyType == "constant") {
			cout << "policy constant, " << policyName << ": " << policyBase << endl;
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
		auto paramValue = paramValueElem->GetText();

		cout << paramName << ": " << paramValue << endl;

		paramElement = paramElement->NextSiblingElement("parameter");
	}
}

void KXml::setOptimizeFunctions() {
	XMLElement *optimizeFunctions = rootElement->FirstChildElement("optimize");
	XMLElement *optimizeEquation = optimizeFunctions->FirstChildElement("equation");
	while(optimizeEquation) {
		auto optimizeEquationElem = optimizeEquation->FirstChildElement("function");
		auto optimizeFunction = optimizeEquationElem->GetText();
		cout << "optimize function: " << optimizeFunction << endl;
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

		cout << equationName << ": " << equationFunction << endl;

		equationElement = equationElement->NextSiblingElement();
	}
}

int main() {
	KXml xmlParser;
	string xmlFile("sampleinput.xml"); 
	bool fileOpened = xmlParser.LoadXmlFile(xmlFile);
	if(!fileOpened) {
		cerr << "XML File load failure" << endl;
		return -1;
	}

	xmlParser.getActors();
	xmlParser.setVariables();
	xmlParser.setPolicies();
	xmlParser.setParameters();
	xmlParser.setOptimizeFunctions();
	xmlParser.setEquations();

/* 	const char *pFilename = "sampleinput.xml";
	XMLDocument doc;
	XMLError eResult = doc.LoadFile(pFilename);
	if(eResult != XML_SUCCESS) {
		cout << "Could not open xml file: " << pFilename << endl; 
		return -1;
	}

	cout << "Loaded xml file: " << pFilename << endl;

	// Get the root element
	//XMLElement *modelElement = doc.FirstChildElement( "model" );
	XMLElement *modelElement = doc.RootElement();
 */
/* 
	XMLElement *actorsList = modelElement->FirstChildElement("actors");
	XMLElement *actor = actorsList->FirstChildElement("actor");
	while(actor) {
		auto *nameElement = actor->FirstChildElement("name");
		auto *actorName = nameElement->GetText();

		auto *descriptionElement = nameElement->NextSiblingElement("description");
		auto description = descriptionElement->GetText();

		auto *influenceElement = descriptionElement->NextSiblingElement("influence");
		auto influence = influenceElement->GetText();

		auto *utilityElement = influenceElement->NextSiblingElement("utility");
		auto utility = utilityElement->GetText();

		//auto *influence = actor->FirstChildElement("influence");
		cout << actorName << ", " << description << ", Influence: " << influence << ", utility: " << utility << endl;

		actor = actor->NextSiblingElement();
	}
 */	

/* 	XMLElement *variables = modelElement->FirstChildElement("variables");
	auto *varElem = variables->FirstChildElement("variable");
	while(varElem) {
		auto *varNameElem = varElem->FirstChildElement("name");
		assert (varNameElem != nullptr);

		auto varName = varNameElem->GetText();

		auto varInitValElem = varElem->FirstChildElement("initial");
		auto varInitVal = varInitValElem->GetText();

		// If an initial value is not there, the formula to compute the initial value should be present
		if(varInitVal == 0) {
			auto varInitValComputeElem = varElem->FirstChildElement("compute");
			assert(varInitValComputeElem != 0);

			auto varInitValComputeExpr = varInitValComputeElem->GetText();
			cout << varName << ": " << varInitValComputeExpr << endl;
		} else {
			cout << varName << ": " << varInitVal << endl;
		}

		// Go to the next variable
		varElem = varElem->NextSiblingElement();
	}
 */

/* 	XMLElement *policies = modelElement->FirstChildElement("policy");
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
			cout << "policy variable, " << policyName << ": " << policyBase << endl;
		} else if(policyType == "constant") {
			cout << "policy constant, " << policyName << ": " << policyBase << endl;
		}

		// Go to the next policy
		policyElem = policyElem->NextSiblingElement();
	}
 */	

/* 	XMLElement *parameters = modelElement->FirstChildElement("parameters");
	XMLElement *paramElement = parameters->FirstChildElement("parameter");
	while(paramElement) {
		auto paramNameElem = paramElement->FirstChildElement("name");
		auto paramName = paramNameElem->GetText();

		auto paramValueElem = paramElement->FirstChildElement("value");
		auto paramValue = paramValueElem->GetText();

		cout << paramName << ": " << paramValue << endl;

		paramElement = paramElement->NextSiblingElement("parameter");
	}
 */

/* 
	XMLElement *optimizeFunctions = modelElement->FirstChildElement("optimize");
	XMLElement *optimizeEquation = optimizeFunctions->FirstChildElement("equation");
	while(optimizeEquation) {
		auto optimizeEquationElem = optimizeEquation->FirstChildElement("function");
		auto optimizeFunction = optimizeEquationElem->GetText();
		cout << "optimize function: " << optimizeFunction << endl;
		optimizeEquation = optimizeEquation->NextSiblingElement();
	}
 */	

/* 	XMLElement *equations = modelElement->FirstChildElement("equations");
	XMLElement *equationElement = equations->FirstChildElement("equation");
	while(equationElement) {
		auto equationNameElem = equationElement->FirstChildElement("name");
		auto equationName = equationNameElem->GetText();

		auto equationFunctionElem = equationElement->FirstChildElement("function");
		auto equationFunction = equationFunctionElem->GetText();

		cout << equationName << ": " << equationFunction << endl;

		equationElement = equationElement->NextSiblingElement();
	}

 */	//doc.SaveFile("sample.xml");
	return 0;
}
