#include "Model.h"
#include "Distribution.h"

madai::PriorDistribution_Cosmo::PriorDistribution_Cosmo(madai::Model * in_Model){
  m_Model = in_Model;
	m_SepMap = parameter::getB(m_Model->m_ParameterMap, "PRIOR_PARAMETER_MAP", false);
	
	if(m_SepMap){
		std::string parmapfile = m_Model->m_DirectoryName + "/parameters/prior.param";
		m_ParameterMap = new parameterMap;
		parameter::ReadParsFromFile(*m_ParameterMap, parmapfile);
		//parameter::ReadParsFromFile(parmap, parameter_file_name);
	}else{
		m_ParameterMap = &(m_Model->m_ParameterMap);
	}
}

double madai::PriorDistribution_Cosmo::Evaluate(std::vector<double> Theta){
	return 1.0;
}