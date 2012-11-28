#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "parametermap.h"
#include <set>
#include "Quad.h"

class emulator;

namespace madai {
    
class Model;
//class QuadHandler;

class Distribution{
public:
	Distribution();
	virtual ~Distribution();
	
	double Normal(double x, double mu, double sigma);
	double IntegratedNormal(double x, double mu, double sigma, double data_sigma);
	double Gaussian(double x, double mu, double sigma);
	double Gaussian(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double Gaussian(gsl_vector x, gsl_vector mu, gsl_matrix sigma, gsl_matrix data_sigma);
	double Log_MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double LogNormal(double x, double mu, double sigma);

protected:

  Model*        m_Model;
  gsl_rng*      m_RandNumGen;
  bool          m_SepMap;
  bool          m_Timing;
  bool          m_Verbose;
  bool          m_Debug;
  parameterMap* m_ParameterMap;
};

/** ---------------------------------------- */

class ProposalDistribution:public Distribution {
public:
	ProposalDistribution(Model *m_Model);
  std::vector<double> Iterate(std::vector<double>& current,
                              float& scale, 
                              std::set<std::string>& activeParameters);
	virtual double Evaluate(std::vector<double> Theta1, 
                          std::vector<double> Theta2, 
                          float scale);

	protected:

	bool                m_SymmetricProposal;
	bool                m_RescaledMethod;
	double              m_Scale;
	double              m_Offset;
  std::vector<double> m_MixingStdDev;
	double*             m_Min_Ranges;
	double*             m_Max_Ranges;

	int FindParam(std::string param_name);
};

class LikelihoodDistribution:public Distribution{
public:
	LikelihoodDistribution();
	virtual ~LikelihoodDistribution();
	virtual double Evaluate(std::vector<double> Theta);
	
  protected:

	virtual std::vector<double> GetData();

  std::vector<double> m_Data;
	bool                m_UseEmulator;
  std::ofstream       m_EmulatorTest;
};

class PriorDistribution:public Distribution {
public:
	PriorDistribution();
	virtual ~PriorDistribution();
	virtual double Evaluate(std::vector<double> Theta);
};

/** ---------------------------------------- */

class PriorDistribution_RHIC:public PriorDistribution{
public:
  PriorDistribution_RHIC(Model *in_Model);
  double Evaluate(std::vector<double> Theta); //ParameterSet Theta

  private:

  std::string              m_Prior;
  bool                     m_Scaled;
  std::vector<double>      m_GaussianMeans;
  std::vector<double>      m_GaussianSTDVS;
  std::vector<double>      m_StepMeans;
  std::vector<std::string> m_StepSide;
};

class PriorDistribution_RHIC_PCA:public PriorDistribution{
public:
  PriorDistribution_RHIC_PCA(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class PriorDistribution_Cosmo:public PriorDistribution {
public:
  PriorDistribution_Cosmo(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class PriorDistribution_Test:public PriorDistribution {
public:
  PriorDistribution_Test(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class LikelihoodDistribution_RHIC:public LikelihoodDistribution{
public:
  LikelihoodDistribution_RHIC(Model *in_Model);
  ~LikelihoodDistribution_RHIC();
  double Evaluate(std::vector<double> Theta);//ParameterSet Theta

  private:

  std::vector<double>  GetFakeData();
  std::vector<double>  GetRealData();
  std::vector<double>  GetRealError();

  double*              m_DataMean;
  double*              m_DataError;
  std::vector<double>  m_Data;
  std::vector<double>  m_Error;
  bool                 m_FakeData;
  bool                 m_SuppressErrors;
  emulator*            m_Emulator;
  parameterMap         m_ObservablesParamMap;

  int FindParam(std::string param_name, std::vector<std::string> PNames);
};

class LikelihoodDistribution_RHIC_PCA:public LikelihoodDistribution{
public:
  LikelihoodDistribution_RHIC_PCA(Model *in_Model);
  ~LikelihoodDistribution_RHIC_PCA();
  double Evaluate(std::vector<double> Theta);

private:

  vector<double> GetRealData(); 


  double*             m_DataMean;
  double*             m_DataError;
  std::vector<double> m_Data;
  std::vector<double> m_Error;
  QuadHandler*        m_Quad;
  parameterMap        m_ObservablesParamMap;
};

class LikelihoodDistribution_Cosmo:public LikelihoodDistribution {
public:
  LikelihoodDistribution_Cosmo(Model *in_Model);
  ~LikelihoodDistribution_Cosmo();
  double Evaluate(std::vector<double> Theta);

private:

  std::vector<double> GetData();

  std::vector<double> m_Data;
  std::vector<int>    m_intData;
};

class LikelihoodDistribution_Test:public LikelihoodDistribution {
public:
  LikelihoodDistribution_Test(Model *in_Model);
  ~LikelihoodDistribution_Test();
  double Evaluate(std::vector<double> Theta);

private:

  std::vector<double> GetData();

  std::vector<double> m_Data;
};
    
}

#endif
