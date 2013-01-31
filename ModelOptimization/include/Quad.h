#ifndef __Quad_h__
#define __Quad_h__

#include "Model.h"

namespace madai {

class Model;

/** \class QuadHandler
 *
 * \todo Add documentation for this class. */
class QuadHandler {
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
  char*       m_Path;
};

} // end namespace madai

#endif // end __Quad_h__
