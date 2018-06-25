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
};

