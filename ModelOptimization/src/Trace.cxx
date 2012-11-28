/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
	Cory Quammen <cquammen AT cs.unc.edu>
	Russell Taylor <taylorr AT cs.unc.edu>
	Scott Pratt <pratt AT nscl.msu.edu>
	Kevin Novak <novakkev AT msu.edu>
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/

#include "Trace.h"
#include "MCMCRun.h"

madai::TraceElement::TraceElement(const std::vector< double > & parameter_values,
                                  const std::vector< double > & output_values ) :
  m_ParameterValues(parameter_values),
  m_OutputValues(output_values) { }

madai::TraceElement::TraceElement(const std::vector< double > & parameter_values) :
	m_ParameterValues(parameter_values),
  m_Used(true) { }

void madai::TraceElement::Reset()
{
  m_ParameterValues.clear();
  m_Used=false;
  m_InTrace = false;
}

madai::TraceElement::TraceElement(madai::Trace *list)
{
  m_Used=false;
  m_LocalTrace = list;
}

void madai::TraceElement::VizTrace()
{
  if(!m_InTrace){
    m_InTrace=true;
  }
}

void madai::TraceElement::Initialize(madai::TraceElement TElement)
{
  m_ParameterValues = TElement.m_ParameterValues;
  m_Used=true;
}

std::vector<std::string> madai::Trace::GetParNames(){
  return (this->m_ParameterNames);
}

unsigned int madai::Trace::length() const
{
	return this->m_TraceElements.size();
}

void madai::Trace::add(const std::vector< double > & parameterValues,
                       const std::vector< double > & OutputValues )
{
	this->m_TraceElements.push_back(
            TraceElement(parameterValues,OutputValues));
}

void madai::Trace::add(const std::vector< double > & parameterValues)
{
  if(m_CurrentIteration >= m_MCMC->m_Writeout){
    std::cerr << "Error: Trace class out of bounds (Greater than WRITEOUT).\n\n";
    exit(1);
  }else{
    for( int i = 0; i<m_MCMC->GetModel()->GetNumberOfParameters();i++)
      m_TraceElements[m_CurrentIteration].m_ParameterValues.push_back(parameterValues[i]);
    m_TraceElements[m_CurrentIteration].m_Used = true;
    m_CurrentIteration++;
  }
}

//Not Working?
void madai::Trace::add(madai::TraceElement TElement)
{
  if(m_CurrentIteration >= m_MCMC->m_Writeout){
    std::cerr << "Error: Trace class out of bounds (Greater than WRITEOUT).\n\n";
    exit(1);
  }else{
    m_TraceElements[m_CurrentIteration]=TElement;
    m_CurrentIteration++;
  }
}

madai::TraceElement & madai::Trace::operator[](unsigned int idx)
{
	return this->m_TraceElements[idx];
}

const madai::TraceElement & madai::Trace::operator[](unsigned int idx) const
{
	return this->m_TraceElements[idx];
}

// Not Working
void madai::TraceElement::Print()
{
  if(m_Used){
    std::cout << "This parameter set is alive." << std::endl;
  }else{
    std::cout << "This parameter set isn't alive." << std::endl;
  }
  std::cerr << "This parameter set contains " << m_ParameterValues.size() << " parameters." << std::endl;
  std::vector<std::string> temp = this->m_LocalTrace->GetParNames();
    
  for(int i=0;i<m_ParameterValues.size();i++)
    std::cout << temp[i] << ":\t" << m_ParameterValues[i] << std::endl;
}

madai::Trace::Trace(madai::MCMCRun *mcmc_in)
{
  m_MCMC = mcmc_in;
  std::string temp;
  int N = m_MCMC->GetModel()->GetNumberOfParameters();
  for(int i=0; i<N;i++){
    m_ParameterNames.push_back((m_MCMC->GetModel()->GetParameters())[i].m_Name);
  }
  m_WriteOutCounter=0;
  m_CurrentIteration=0;
  m_TraceElements.reserve(m_MCMC->m_Writeout+1);
  for(int i=0; i<m_MCMC->m_Writeout;i++){//127
    m_TraceElements.push_back(TraceElement(this));
  }
}

template <class T>
void write_vector(std::ostream& o, std::vector< T > const & v, char delim) {
	if (! v.empty()) {
		typename std::vector< T >::const_iterator itr = v.begin();
		o << *(itr++);
		while (itr < v.end())
			o << delim << *(itr++);
	}
}

void madai::Trace::write(std::ostream & out) const {
	unsigned int N = this->length();
	for (unsigned int i = 0; i < N; i ++) {
		write_vector(out, (*this)[i].m_ParameterValues, ',');
		out << ',';
    write_vector(out, (*this)[i].m_OutputValues, ',');
		if ((*this)[i].m_Comments.size() > 0) {
			out << ",\"";
			write_vector(out, (*this)[i].m_Comments, ';');
			out << '"';
		}
		out << '\n';
	}
	out.flush();
}

	/*
	  Assert:
	    FOR ALL i < this->m_TraceElements.size():
	      this->m_TraceElements[i].m_ParameterValues.size() == params.size()
	      this->m_TraceElements[i].m_OutputValues.size() == outputs.size()
	*/

void madai::Trace::writeHead(std::ostream & o, 
                             const std::vector<madai::Parameter> & params,
                             const std::vector<std::string> & outputs) const
{
	if (! params.empty()) {
		std::vector< madai::Parameter >::const_iterator itr = params.begin();
		o << '"' << itr->m_Name << '"';
		for (itr++; itr < params.end(); itr++)
			o << ',' << '"' << itr->m_Name << '"';
		if (! outputs.empty())
			o << ',';
	}
	if (! outputs.empty()) {
		std::vector<std::string>::const_iterator itr = outputs.begin();
		o << '"' << *itr << '"';
		for (itr++; itr < outputs.end(); itr++)
			o << ',' << '"' << *itr << '"';
		o << '\n';
	}
}

void madai::Trace::writeHead(std::ostream & o, 
                             const std::vector<madai::Parameter> & params) const
{
	if (! params.empty()) {
		std::vector< madai::Parameter >::const_iterator itr = params.begin();
		o << '"' << itr->m_Name << '"';
		for (itr++; itr < params.end(); itr++)
			o << ',' << '"' << itr->m_Name << '"';
	}
}

void madai::Trace::PrintDataToFile()
{
  std::cout << "Printing data to file." << std::endl;
  std::ofstream outputfile;
  std::stringstream ss;
  ss << m_WriteOutCounter+1;
  std::string out_file = m_MCMC->m_TraceDirectory+"/output"+ss.str()+".dat";
    
  if(m_TraceElements[0].m_Used){
    outputfile.open(out_file.c_str());
    std::cout << "Writing out to: " << out_file << std::endl;
    if(outputfile){
      outputfile << "#ITERATION,";
      if(!m_ParameterNames.empty()){
        std::vector<std::string> temp_names = m_ParameterNames;
        for(unsigned int i=0;i<temp_names.size();i++){
          outputfile << temp_names[i] << ",";
        }
      }
      outputfile << std::endl;
            
      for(int i=0;i<m_MCMC->m_Writeout;i++){
        if(m_TraceElements[i].m_Used){
          outputfile << i+m_WriteOutCounter*m_MCMC->m_Writeout << ",";
          if(!m_TraceElements[i].m_ParameterValues.empty()){
            if(m_MCMC->m_RescaledTrace){
              double * range = new double[2]();
              for(int j=0;j<m_TraceElements[i].m_ParameterValues.size();j++){
                m_MCMC->GetModel()->GetRange(j,range);
                outputfile << (m_TraceElements[i].m_ParameterValues[j]-range[0])/(range[1]-range[0]);
                //outputfile << mcmc->RescaleTheta(i,j);
                if(j!=m_TraceElements[i].m_ParameterValues.size()-1){
                  outputfile << ",";
                }
              }
            }else{
              for(int j=0;j<m_TraceElements[i].m_ParameterValues.size();j++){
                outputfile << m_TraceElements[i].m_ParameterValues[j];
                if(j!=m_TraceElements[i].m_ParameterValues.size()-1){
                  outputfile<< ",";
                }
              }
            }
          }else{
            std::cout << "Error: Accessing empty element." << std::endl;
            exit(1);
          }
        }
        outputfile << std::endl;
      }
      outputfile.close();
    }else{
      std::cout << "Error: Couldn't open the output file" << std::endl;
      exit(1);
    }
  }else{
    std::cerr << "The first element of the list is not used. ERROR\n\n";
    exit(1);
  }
}

void madai::Trace::WriteOut()
{
  this->PrintDataToFile();
  for(int i=0;i<m_MCMC->m_Writeout;i++){
    m_TraceElements[i].Reset();
  }
  m_WriteOutCounter++;
  m_CurrentIteration=0;
}

void madai::Trace::MakeTrace()
{
  std::stringstream ss;
  ss << "cat ";
    
  for(int i=1;i<=ceil((double)(m_MCMC->m_MaxIterations)/(double)(m_MCMC->m_Writeout));i++){
    std::cout << "Parsing " << m_MCMC->m_TraceDirectory << "/output" << i << ".dat" << std::endl;
    ss << m_MCMC->m_TraceDirectory << "/output" << i << ".dat ";
  }
  ss << "> " << m_MCMC->m_TraceDirectory << "/trace.dat" << std::endl;

  std::string command = ss.str();
  std::system(command.c_str());
    
  ss.str(string());
}