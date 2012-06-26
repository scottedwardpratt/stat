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

madai::TraceElement::TraceElement(
		const std::vector< double > & parameter_values,
		const std::vector< double > & output_values ) :
	m_ParameterValues(parameter_values),
	m_OutputValues(output_values) { }

unsigned int madai::Trace::length() const
{
	return this->m_TraceElements.size();
}
void madai::Trace::add(
		const std::vector< double > & parameterValues,
		const std::vector< double > & OutputValues )
{
	this->m_TraceElements.push_back(
		TraceElement(parameterValues,OutputValues));
}
madai::TraceElement & madai::Trace::operator[](unsigned int idx)
{
	return this->m_TraceElements[idx];
}

const madai::TraceElement & madai::Trace::operator[](unsigned int idx) const
{
	return this->m_TraceElements[idx];
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

void madai::Trace::writeHead(
	std::ostream & o, 
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
