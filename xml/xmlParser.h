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
	std::vector<variableName> listOfVars;
    typedef std::string initComputeExpr;
    typedef std::map<variableName, initComputeExpr> initValExprMap;
    initValExprMap initValExpressions;

public:
	bool LoadXmlFile(std::string & xmlFileName);
	std::vector<std::string> getActors();
	void setVariables();
	void setPolicies();
	void setParameters();
	void setOptimizeFunctions();
	void setEquations();
};

