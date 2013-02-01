#ifndef __QUAD_H__
#define __QUAD_H__

#include "Model.h"

namespace madai {

class Model;

class QuadHandler{
public:
	QuadHandler(parameterMap *parmap, Model * m_Model);
	~QuadHandler();

	void QueryQuad(std::vector<double> Theta,vector<double> &Means, vector<double> &Errors);

private:

  Model*      m_Model;
  std::string m_EmulatedParams;
	std::string m_QuadScriptHome;
	std::string m_EmInputFile;
	std::string m_EmOutputFile;
	std::string m_EmErrorFile;
	char*       m_pPath;
};

} // end namespace madai

#endif // end __QUAD_H__