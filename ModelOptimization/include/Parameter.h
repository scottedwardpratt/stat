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
#ifndef __Parameter_h_
#define __Parameter_h_

#include <string>
#include <vector>

namespace madai {
	using namespace std;
	
	class CParameterInfo {
	public:
			//ParameterInfo() {};
			//virtual ~ParameterInfo() {};
		string Name;
		unsigned int index;
		double MinimumPossibleValue;
		double MaximumPossibleValue;
		double BestValue;
		
	}; // end class ParameterInfo

	class Parameter {
	public:
		Parameter(std::string nm, 
							double mn=0.0,
							double mx=1.0 ) :
			m_Name(nm),
			m_MinimumPossibleValue(mn), 
			m_MaximumPossibleValue(mx) { }
		virtual ~Parameter() { }
		std::string m_Name;
		double      m_MinimumPossibleValue;
		double      m_MaximumPossibleValue;
	}; // end class Parameter
	
	class CParameterList{
	public:
		vector<CParameterInfo> *parinfo;
		CParameterInfo *GetParInfo(unsigned int parindex);
		CParameterInfo *GetParInfo(string parname);
		void AddParameter(string parname,double maxvalue,double minvalue);
		void DeleteParameter(string parname);
		void DeleteParameter(int parindex);
		int Length;
	};

} // end namespace madai

#endif // __Parameter_h_
