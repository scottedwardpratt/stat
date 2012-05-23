#ifndef __Parameter_h_
#define __Parameter_h_

#include <string>

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
	
	class CParameterList{
	public:
		vector<CParameterInfo> *parinfo;
		CParameterInfo *GetParInfo(unsigned int parindex);
		CParameterInfo *GetParInfo(string parname);
		void AddParameter(string parname,double maxvalue,double minvalue);
		void DeleteParameter(string parname);
		void DeleteParameter(int parindex);
		int Length;
	}
	
} // end namespace madai

#endif // __Parameter_h_
