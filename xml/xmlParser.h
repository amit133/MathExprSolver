#include <string>
#include <vector>
#include <map>
#include "tinyxml2.h"

class KXml {
private:
	std:: string xmlFileName_m;
	tinyxml2::XMLDocument doc;
	tinyxml2::XMLElement *rootElement = nullptr;

	typedef std::string variableName;
	typedef std::vector<variableName> listOfVariables;
	listOfVariables listOfVars;

    typedef std::string initComputeExpr;
    typedef std::map<variableName, initComputeExpr> initValExprMap;
    initValExprMap initValExpressions;

	// Policies
	typedef std::string policyName;
	typedef double policyValue;
	typedef std::map<policyName, policyValue> policies;
	policies policyVariables;
	policies policyConstants;
	std::vector<policyName> listOfPolicyVariableNames;

	// parameters
	typedef std::string paramName;
	typedef double paramValue;
	typedef std::map<paramName, paramValue> parameters;
	parameters params;

	// System of Optimization Equations
    typedef std::string equation;
    typedef std::vector<equation> optimizingFunctions;
    optimizingFunctions systemOfFunctions;

private:
	bool LoadXmlFile(std::string & xmlFileName);
	void setActors();
	void setVariables();
	void setPolicies();
	void setParameters();
	void setOptimizeFunctions();
	void setEquations();

public:
	KXml() = delete;
	KXml(const KXml&) = delete;
	KXml& operator=(const KXml&) = delete;
	explicit KXml(std::string & xmlFileName);

	listOfVariables getVariables();
	initValExprMap getInitComputeExpressions();
	parameters getParameters();
	optimizingFunctions getOptimizeFunctions();
	policies getPolicyVariables();
	policies getPolicyConstants();
	std::vector<policyName> getPolicyVariableNames();
};

