#ifndef __Distribution_h__
#define __Distribution_h__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "parametermap.h"
#include "Model.h"

#include <set>

namespace madai {
    
class Model;

/** \class Distribution
 *
 * Base class for distributions. */
class Distribution {
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

// Proposal isn't currently being used (transfered functions to MCMCRun)
class ProposalDistribution : public Distribution {
public:
  ProposalDistribution(Model *m_Model);
  std::vector<double> Iterate(std::vector<double>& current,
                              double& scale);
  virtual double Evaluate(std::vector<double> Theta1, 
                          std::vector<double> Theta2, 
                          double scale);
    
protected:
    
  bool                  m_SymmetricProposal;
  bool                  m_RescaledMethod;
  double                m_MinScale;
  double                m_MaxScale;
  double                m_Prefactor;
  double                m_Scale;
  std::vector<double>   m_MixingStdDev;
  std::set<std::string> m_ActiveParameters;
    
  int FindParam(std::string param_name);
};

/** ---------------------------------------- */

class LikelihoodDistribution : public Distribution {
public:
  LikelihoodDistribution();
  virtual ~LikelihoodDistribution();
  virtual double Evaluate(std::vector<double> ModelMeans,
                          std::vector<double> ModelErrors);
  
protected:
  virtual std::vector<double> GetData();

  std::vector<double> m_Data;
  bool                m_UseEmulator;
  bool                m_ProcessPipe;
  std::ofstream       m_EmulatorTest;
};

/** ---------------------------------------- */

class PriorDistribution : public Distribution {
public:
  PriorDistribution();
  virtual ~PriorDistribution();
  virtual double Evaluate(std::vector<double> Theta);
};

/** ---------------------------------------- */

class RHIC_PriorDistribution : public PriorDistribution {
public:
  RHIC_PriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta); //ParameterSet Theta

private:
  std::string              m_Prior;
  bool                     m_Scaled;
  std::vector<double>      m_GaussianMeans;
  std::vector<double>      m_GaussianSTDVS;
  std::vector<double>      m_StepMeans;
  std::vector<std::string> m_StepSide;
};

/** ---------------------------------------- */

class InterpolatorPriorDistribution : public PriorDistribution{
public:
  InterpolatorPriorDistribution(Model * in_Model);
  double Evaluate(std::vector<double> Theta);
  
private:

  std::string              m_Prior;
  bool                     m_Scaled;
  std::vector<double>      m_GaussianMeans;
  std::vector<double>      m_GaussianSTDVS;
  std::vector<double>      m_StepMeans;
  std::vector<std::string> m_StepSide;
};
  
class RHIC_PCA_PriorDistribution : public PriorDistribution {
public:
  RHIC_PCA_PriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class CosmoPriorDistribution : public PriorDistribution {
public:
  CosmoPriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class TestPriorDistribution : public PriorDistribution {
public:
  TestPriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

class RHIC_LikelihoodDistribution : public LikelihoodDistribution{
public:
  RHIC_LikelihoodDistribution(Model *in_Model);
  ~RHIC_LikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);

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
  parameterMap         m_ObservablesParamMap;

  int FindParam(std::string param_name, std::vector<std::string> PNames);
};

class RHIC_PCA_LikelihoodDistribution : public LikelihoodDistribution{
public:
  RHIC_PCA_LikelihoodDistribution(Model *in_Model);
  ~RHIC_PCA_LikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);
    
private:
  std::vector<double> GetRealData(); 
  
  double*             m_DataMean;
  double*             m_DataError;
  std::vector<double> m_Data;
  std::vector<double> m_Error;
  parameterMap        m_ObservablesParamMap;
};

class CosmoLikelihoodDistribution : public LikelihoodDistribution {
public:
  CosmoLikelihoodDistribution(Model *in_Model);
  ~CosmoLikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);

private:
  std::vector<double> GetData();
  
  std::vector<double> m_Data;
  std::vector<int>    m_intData;
};

class TestLikelihoodDistribution : public LikelihoodDistribution {
public:
  TestLikelihoodDistribution(Model *in_Model);
  ~TestLikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);
  
private:
  std::vector<double> GetData();
  
  std::vector<double> m_Data;
};

} // namespace madai

#endif
