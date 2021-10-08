//!  SimulationParameters.cpp
/*!
Container class for all simulation related parameters that need to get passed between classes and functions

Author: Kai Herz <kai.herz@tuebingen.mpg.de>

Copyright 2020 Kai Herz

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "BMSimFit.h"

//! Get fit data
/*!	\return pointer to fit data */
std::vector<double>* SimFitParameters::GetFitData()
{
	return &fitData;
}

//! Get fit data weights
/*!	\return pointer to fit data */
std::vector<double>* SimFitParameters::GetFitDataWeights()
{
	return &weights;
}

//! Get ADC positions
/*!	\return pointer to vector with indices of adc events */
std::vector<long>* SimFitParameters::GetADCPositions()
{
	return &adcPos;
}


//! Get the scaling of the initial M
/*!	\return scale of magnetization vector */
double SimFitParameters::GetMagnetizationVectorScale()
{
	return scale;
}

//! Set the scaling of the initial M
void SimFitParameters::SetMagnetizationVectorScale(double s)
{
	scale = s;
}


//! Register a fit parameter
/*!
    We can fit water and cest and b0 shift so far
	\param name name of the parameter
	\param start start value
	\param upper upper boundary
	\param lower lower boundary
	\return true if registration worked
*/
bool SimFitParameters::RegisterFitParameter(std::string name, double start, double upper, double lower)
{
	using std::placeholders::_1; // need this for std::bind
	
	// check if values are valid
	if (lower > start || upper < start) {
		std::cout << "Start value of " << name << " is not between boundaries " << std::endl;
		return false;
	}

	// init the fit parameter	
	FitParameter fp;
	fp.name = name;
	fp.lower = lower;
	fp.upper = upper;
	
	// divide the name string by underline
	std::stringstream namestr(name);
	std::vector<std::string> seglist; // seg list contains the string parts divided by _
	std::string segment;
	while (std::getline(namestr, segment, '_')) 
		seglist.push_back(segment);
	
	// we can choose the fit parameters dynamically by binding the corresponding get and set funtions 
	// in the SimulationParameter object
	if (seglist.size() == 1) { // single word ->  b0 shift
		if (seglist.at(0).compare("b0shift") == 0) {
			fp.set = std::bind(&SimFitParameters::SetScannerB0Inhom, this, _1);
			fp.get = std::bind(&SimFitParameters::GetScannerB0Inhom, this);
		}
		else if(seglist.at(0).compare("scale") == 0) {
			fp.set = std::bind(&SimFitParameters::SetMagnetizationVectorScale, this, _1);
			fp.get = std::bind(&SimFitParameters::GetMagnetizationVectorScale, this);
		}
		else {
			std::cout << "ERROR: " << name << " is not a valid name for a fit parameter! " << std::endl;
		    return false;
		}
	}
	else if (seglist.size() == 2) { // two words -> has to be water
		if (seglist.at(0).compare("water") == 0)	{
			if (seglist.at(1).compare("t1") == 0) {
				fp.set = std::bind(&WaterPool::SetT1, this->GetWaterPool(), _1);
				fp.get = std::bind(&WaterPool::GetT1, this->GetWaterPool());
			}
			else if (seglist.at(1).compare("t2") == 0) {
				fp.set = std::bind(&WaterPool::SetT2, this->GetWaterPool(), _1);
				fp.get = std::bind(&WaterPool::GetT2, this->GetWaterPool());
			}
			else {
				std::cout << "ERROR: Can only fit T1 and T2 of water, but not " << seglist.at(1) << std::endl;
				return false;
			}
		}
		else {
			std::cout << "ERROR: " << name << " is not a valid name for a fit parameter! " << std::endl;
			return false;
		}
	}
	else if (seglist.size() == 3) { //three words -> has to be cest
		if (seglist.at(0).compare("cest") == 0) {
			int npool = atoi(seglist.at(1).c_str()) - 1;
			if (npool < 0 || npool + 1 > this->GetNumberOfCESTPools()) {
				std::cout << "ERROR: Can not fit params of cest pool " << npool+1 << "! There are only " << this->GetNumberOfCESTPools() << " pools in the param file!" << std::endl;
				return false;
			}
			if (seglist.at(2).compare("t1") == 0) {
				fp.set = std::bind(&CESTPool::SetT1, this->GetCESTPool(npool), _1);
				fp.get = std::bind(&CESTPool::GetT1, this->GetCESTPool(npool));
			}
			else if (seglist.at(2).compare("t2") == 0) {
				fp.set = std::bind(&CESTPool::SetT2, this->GetCESTPool(npool), _1);
				fp.get = std::bind(&CESTPool::GetT2, this->GetCESTPool(npool));
			}
			else if (seglist.at(2).compare("k") == 0) {
				fp.set = std::bind(&CESTPool::SetExchangeRateInHz, this->GetCESTPool(npool), _1);
				fp.get = std::bind(&CESTPool::GetExchangeRateInHz, this->GetCESTPool(npool));
			}
			else if (seglist.at(2).compare("dw") == 0) {
				fp.set = std::bind(&CESTPool::SetShiftinPPM, this->GetCESTPool(npool), _1);
				fp.get = std::bind(&CESTPool::GetShiftinPPM, this->GetCESTPool(npool));
			}
			else if (seglist.at(2).compare("f") == 0) {
				fp.set = std::bind(&CESTPool::SetFraction, this->GetCESTPool(npool), _1);
				fp.get = std::bind(&CESTPool::GetFraction, this->GetCESTPool(npool));
			}
			else {
				std::cout << "ERROR: Can only fit T1, T2, k, f and dw of cest pool, but not " << seglist.at(2) << std::endl;
				return false;
			}
		}
		else {
			std::cout << "ERROR: " << name << " is not a valid name for a fit parameter! " << std::endl;
			return false;
		}
	}
	else{
		std::cout << "ERROR: " << name << " is not a valid name for a fit parameter! " << std::endl;
		return false;
	}

	// try to access the parameter via the set function
	try
	{
		fp.set(start);
		fitParams.push_back(fp); // register param
	}
	catch (...)
	{
		std::cout << "Unspecified ERROR: Could not register " << name << std::endl;
		return false;
	}
	return true;
}

//! Get fit parameters
/*!
	\return a pointer to avector for the fit parameters, or NULL if none were registered
*/
std::vector<FitParameter>* SimFitParameters::GetFitParams()
{
	return fitParams.size() == 0 ? NULL : &fitParams;
}


//! Constructor
/*!	\param simFitParams initial SimFitParameters object */
FitFramework::FitFramework(SimFitParameters& simFitParams) : BMCSim(simFitParams) {
	sfp = &simFitParams;
}

//! Copy Constructor
/*!	
    The copy constructor is needed to spread the solver across threads
    \param ff FitFramework object that should be copied
*/
FitFramework::FitFramework(FitFramework& ff) : BMCSim(*(ff.GetSimFitParameters()))
{
	// copy all parameters
	sfp = ff.sfp;
	seq = ff.seq;
	sequenceLoaded = ff.sequenceLoaded;
	uniquePulses = ff.uniquePulses;
	numberOfADCBlocks = ff.numberOfADCBlocks;
	Mvec = ff.Mvec;
}

//! Get the Simulation and fit parameters
/*!	\return pointer to SimFitParameters object */
SimFitParameters* FitFramework::GetSimFitParameters()
{
	return sfp;
}

//! Get the Exxternal Sequence
/*!	\return pointer to ExternalSequence object */
ExternalSequence* FitFramework::GetExternalSequence()
{
	return &seq;
}

//! Run a specific event block
/*!
	Wrapper for protected BMCSim function
    \param M reference to current magnetization vector
    \param accumPhase reference to current accumulated phase
    \param seqBlock pointer to specific ebent block
*/
void FitFramework::RunNonADCEventBlock(Eigen::VectorXd &M, float &accumPhase, SeqBlock* seqBlock)
{
	RunEventBlock(M, accumPhase, seqBlock);
}

//! Update the Bloch matrix matrix
/*! Wrapper for protected solver method*/
void FitFramework::UpdateSolver()
{
	solver->UpdateSimulationParameters(*sfp);
}

//! Load external Pulseq sequence
/*!
    Wraps BMCSim load function and decodes adc positions for parallel evaluation
	\param path full filepath of the .seq-file
	\return true if sequence could be loaded and contains adc events
*/
bool FitFramework::SetExternalSequence(std::string path)
{
	ExternalSequence extSeq;
	bool seqLoaded = extSeq.load(path);
	if (seqLoaded) {
		for (unsigned int nSample = 0; nSample < extSeq.GetNumberOfBlocks(); nSample++)
		{
			SeqBlock* seqBlock = extSeq.GetBlock(nSample);
			if (seqBlock->isADC()) {
				sfp->GetADCPositions()->push_back(nSample);
			}
			delete seqBlock;
		}
		sequenceLoaded = LoadExternalSequence(path);
	}
	return sequenceLoaded;
}