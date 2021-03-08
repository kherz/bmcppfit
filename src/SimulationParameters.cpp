//!  SimulationParameters.cpp
/*!
Container class for all simulation related parameters that need to get passed between classes and functions

kai.herz@tuebingen.mpg.de

Copyright 2020 Kai Herz

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "SimulationParameters.h"


// Water Pool Function Definitions ////

//! Default Constructor
WaterPool::WaterPool() : R1(0), R2(0), f(0) {}

//! Constructor
/*!
  \param nR1 R1 of water pool [1/s]
  \param nR2 R2 of water pool [1/s]
  \param nf fraction of water pool
*/
WaterPool::WaterPool(double nR1, double nR2, double nf) : R1(nR1), R2(nR2), f(nf) {}

//! Default destructor
WaterPool::~WaterPool() {}

//! Get R1
/*! \return 1/T1 of pool */
const double WaterPool::GetR1() { return R1; }

//! Get R2
/*! \return 1/T2 of pool */
const double WaterPool::GetR2() { return R2; }

//! Get T1
/*! \return T1 of pool */
const double WaterPool::GetT1() { return 1.0 / R1; }

//! Get T2
/*! \return T2 of pool */
const double WaterPool::GetT2() { return 1.0 / R2; }

//! Get f
/*! \return fraction of pool */
const double WaterPool::GetFraction() { return f; }

//! Set R1
/*! \param new 1/T1 of pool */
void WaterPool::SetR1(double nR1) { R1 = nR1; }

//! Set R2
/*! \param new 1/T2 of pool */
void WaterPool::SetR2(double nR2) { R2 = nR2; }

//! Set T1
/*! \param new T1 of pool */
void WaterPool::SetT1(double nT1) { R1 = 1.0 / nT1; }

//! Set T2
/*! \param new T2 of pool */
void WaterPool::SetT2(double nT2) { R2 = 1.0 / nT2; }

//! Set f
/*! \param new fraction of pool */
void WaterPool::SetFraction(double nf) { f = nf; }


// CEST Pool Function Definitions ////

//! Default Constructor
CESTPool::CESTPool() : dw(0), k(0) {}

//! Constructor
/*!
  \param nR1 R1 of CEST pool [1/s]
  \param nR2 R2 of CEST pool [1/s]
  \param nf fraction of CEST pool
  \param ndw chemical shift of CEST pool [ppm]
  \param nk exchange rate of CEST pool [Hz]
*/
CESTPool::CESTPool(double nR1, double nR2, double nf, double ndw, double nk) : WaterPool(nR1, nR2, nf), dw(ndw), k(nk) {}

//! Copy Constructor
/*!
  \param c pointer to existing CESTPool class object
*/
CESTPool::CESTPool(CESTPool* c) : WaterPool(c->R1, c->R2, c->f), dw(c->dw), k(c->k) {}

//! Default destructor
CESTPool::~CESTPool() {}

//! Get chemical shift
/*! \return chemical shift of pool in ppm*/
const double CESTPool::GetShiftinPPM() { return dw; }

//! Get exchange rate
/*! \return exchage rate of pool in Hz*/
const double CESTPool::GetExchangeRateInHz() { return k; }

//! Set shift
/*! \param ndw new chemical shift of pool in ppm*/
void CESTPool::SetShiftinPPM(double ndw) { dw = ndw; }

//! Set exchange rate
/*! \param nk new exchage rate of pool in Hz*/
void CESTPool::SetExchangeRateInHz(double nk) { k = nk; }


// MT Pool Function Definitions ////

//! Default Constructor
MTPool::MTPool() : ls(None) {}

//! Constructor
/*!
  \param nR1 R1 of MT pool [1/s]
  \param nR2 R2 of MT pool [1/s]
  \param nf fraction of MT pool
  \param ndw chemical shift of MT pool [ppm]
  \param nk exchange rate of MT pool [Hz]
  \param nls lineshape of MT pool
*/
MTPool::MTPool(double nR1, double nR2, double nf, double ndw, double nk, MTLineshape nls) : CESTPool(nR1, nR2, nf, ndw, nk), ls(nls) {}

//! Default destructor
MTPool::~MTPool() {}

//! Get the lineshape
/*! \return lineshape of MT pool*/
MTLineshape MTPool::GetMTLineShape() { return ls; }

//! Set the lineshape
/*! \param nls new lineshape of MT pool*/
void MTPool::SetMTLineShape(MTLineshape nls) { ls = nls; }

//! Get the MT parameter at the current offset
/*!
	The return value considers the lineshape of the MT pool
	It's defined in doi:10.1088/0031-9155/58/22/R221 as Rrfb
	The return value is Rrfb/(w1^2)
	\param offset frequency offset of rf pulse
	\param omega0 larmor frequency
	\return Rrfb/(w1^2) of MT pool
*/
double MTPool::GetMTLineAtCurrentOffset(double offset, double omega0)
{
	double mtLine = 0.0;
	switch (ls)
	{
	case None: // 0.0 stays
	{
		break;
	}
	case Lorentzian:
	{
		double T2 = 1 / R2;
		mtLine = T2 / (1 + pow((offset - dw * omega0)*T2, 2.0));
		break;
	}
	case SuperLorentzian: // integrated SL lineshape
	{
		double dw0 = dw * omega0;
		double dwPool = (offset - dw0);
		if (abs(dwPool) >= omega0) // empirical cutoff is 1 ppm
		{
			mtLine = InterpolateSuperLorentzianShape(dwPool);
		}
		else // return lorentian lineshape if we are in a pol region
		{
			std::vector<double> px{ -300 - omega0, -100 - omega0, 100 + omega0, 300 + omega0 };
			std::vector<double> py(px.size(),0);
			for (int i = 0; i < px.size(); i++) {
				py[i] = InterpolateSuperLorentzianShape(px[i]);
			}
			mtLine = CubicHermiteSplineInterpolation(dwPool, px, py);
		}
		break;
	}
	default:
		break;
	}
	return mtLine;
}

//! Calculate the SuperLorentzian Lineshape
/*!
	\param dw frequency offset between rf pulse and offset of MT pool [rad]
	\return lineshape of MT pool
*/
double MTPool::InterpolateSuperLorentzianShape(double dw)
{
	double mtLine = 0.0;
	double T2 = 1 / R2;
	int numberOfIntegrationSamples = 101; // number of
	double integrationSampleStepSize = 0.01;
	double sqrt_2_pi = sqrt(2.0 / M_PI);
	for (int i = 0; i < numberOfIntegrationSamples; i++)
	{
		double powcu2 = abs(3.0 * pow(integrationSampleStepSize*double(i), 2.0) - 1.0); // helper variable
		mtLine += sqrt_2_pi * T2 / powcu2 * exp(-2.0 * pow(dw * T2 / powcu2, 2.0)); // add to integrate
	}
	return mtLine * (M_PI * integrationSampleStepSize); // final line
}

//! Spline interpolation to avoid pol in superlorentzian lineshape function
/*!
	\param px_int position at which the splin eposition should be calculated
	\param px vector with the x postion of the 4 grid points of the spline
	\param py vector with the y postion of the 4 grid points of the spline
	\return lineshape of MT pool, 0.0 if non-valid input
*/
double MTPool::CubicHermiteSplineInterpolation(double px_int, std::vector<double> px, std::vector<double> py)
{
	if (px.size() != 4 || py.size() != 4)
		return 0.0;

	//y values
	double p0y = py[1]; //points
	double p1y = py[2];

	double tangentWeight = 30; // empirically chosen 
	double d0y = tangentWeight * (p0y - py[0]); //tangents
	double d1y = tangentWeight * (py[3] - p1y);

	//calculate the interpolation points
	double cStep = abs((px_int - px[1] + 1.0) / (px[2] - px[1] + 1.0));//
	double c3 = cStep * cStep*cStep;
	double c2 = cStep * cStep;

	// hermite spline
	double h0 = 2 * c3 - 3 * c2 + 1;
	double h1 = -2 * c3 + 3 * (c2);
	double h2 = c3 - 2 * c2 + cStep;
	double h3 = c3 - c2;
	//calculate y value
	return h0 * p0y + h1 * p1y + h2 * d0y + h3 * d1y;
}



// Simulation Parameters Function Definitions ////

//! Constructor
SimulationParameters::SimulationParameters()
{
	numberOfCESTPools = 0;
	simulateMTPool = false;
	cestMemAllocated = false;
	verboseMode = false;
	useInitMagnetization = true;
	maxNumberOfPulseSamples = 100;
	scale = 1.0;
	InitScanner(0.0);
}

//! Destructor
SimulationParameters::~SimulationParameters()
{
	if (numberOfCESTPools > 0 && cestMemAllocated)
		delete[] cestPools;
}

//! Set external sequence object
/*!	\param seq ExternalSequence object that should be simulated */
void SimulationParameters::SetExternalSequence(std::string filename)
{
	sequence.load(filename);
	this->DecodeSeqInfo();
}

//! Get external sequence object
/*!	\return ExternalSequence object that should be simulated */
ExternalSequence* SimulationParameters::GetExternalSequence()
{
	return &sequence;
}

//! Init Magnetitazion Vector Array
/*!
	Replicate the initial Magnetization vector for output
	\param M initial magnetization vector after ADC
	\param numOutput number of ADC events in external sequence
*/
void SimulationParameters::InitMagnetizationVectors(VectorXd &M, unsigned int numOutput)
{
	Mvec = M.rowwise().replicate(numOutput);
}

//! Get Magnetization vectors
/*!	\return Magnetization vectors at each ADC event */
MatrixXd* SimulationParameters::GetMagnetizationVectors()
{
	return &Mvec;
}

//! Set Water Pool
/*!	\param wp new water pool */
void SimulationParameters::SetWaterPool(WaterPool wp)
{
	waterPool = wp;
}

//! Get Water Pool
/*!	\return pointer to water pool */
WaterPool* SimulationParameters::GetWaterPool()
{
	return &waterPool;
}

//! Init CEST pool memory
/*!
	Allocate heap memory for CEST pools
	\param numPools number of CEST pools that should be simulated
*/
void  SimulationParameters::InitCESTPoolMemory(unsigned int numPools)
{
	numberOfCESTPools = numPools;
	cestPools = new CESTPool[numberOfCESTPools];
	cestMemAllocated = true;
}

//! Set CEST Pool
/*!
	\param cp new CEST pool
	\param poolIdx id of cest pool [0 ... numberOfCESTPools-1]
*/
void  SimulationParameters::SetCESTPool(CESTPool cp, unsigned int poolIdx)
{
	if (poolIdx < numberOfCESTPools)
		cestPools[poolIdx] = cp;
}

//! Get CEST Pool
/*!
	\param poolIdx id of cest pool [0 ... numberOfCESTPools-1]
	\return pointer to cest pool at poolIdx
*/
CESTPool* SimulationParameters::GetCESTPool(unsigned int poolIdx)
{
	return poolIdx < numberOfCESTPools ? &cestPools[poolIdx] : NULL;
}

//! Set MT Pool
/*!	\param mp new mt pool */
void SimulationParameters::SetMTPool(MTPool mp)
{
	mtPool = mp;
	simulateMTPool = true;
}

//! Get MT Pool
/*!	\return pointer to mt pool */
MTPool* SimulationParameters::GetMTPool()
{
	return simulateMTPool ? &mtPool : NULL;
}

//! Set Scanner related info
/*!
    \param b0 static field [T]
    \param relB1 relative B1
    \param B0Inhomogeneity field inhomogeneity [ppm]
    \param Gamma gyromagnetic ratio [rad/uT]
*/
void SimulationParameters::InitScanner(double b0, double b1, double b0Inh, double gamma)
{
	scanner.B0 = b0;
	scanner.relB1 = b1;
	scanner.B0Inhomogeneity = b0Inh;
	scanner.Gamma = gamma;
}

//! Get Scanner B0
/*!	\return static field of scanner [T] */
double SimulationParameters::GetScannerB0()
{
	return scanner.B0;
}

//! Get Scanner relative B1
/*!	\return relative b1 of scanner */
double SimulationParameters::GetScannerRelB1()
{
	return scanner.relB1;
}

//! Set Scanner relative B1
/*!	\param relative b1 of scanner */
void SimulationParameters::SetScannerRelB1(double b1)
{
	scanner.relB1 = b1;
}

//! Get Scanner B0 inhomogeneity
/*!	\return field inhomogeneity [ppm] of scanner */
double SimulationParameters::GetScannerB0Inhom()
{
	return scanner.B0Inhomogeneity;
}

//! Get Scanner B0 inhomogeneity
/*!	\return field inhomogeneity [ppm] of scanner */
void SimulationParameters::SetScannerB0Inhom(double inho)
{
	scanner.B0Inhomogeneity = inho;
}

//! Get Scanner Gamma
/*!	\return gyromagnetic ratio of simulated nucleus [rad/uT] */
double SimulationParameters::GetScannerGamma()
{
	return scanner.Gamma;
}

//! Get bool if MT should be simulated
/*!	\return true, if an MT pool should be simulated */
bool SimulationParameters::IsMTActive()
{
	return simulateMTPool;
}

//! Get Number of CEST Pools
/*!	\return total number of CEST pools that should be simulates */
unsigned int SimulationParameters::GetNumberOfCESTPools()
{
	return numberOfCESTPools;
}

//! Set Verbose mode
/*!	\param v info about verbosity (true/false) */
void SimulationParameters::SetVerbose(bool v)
{
	verboseMode = v;
}

//! Get verbose mode
/*!	\return info about verbosity (true/false) */
bool SimulationParameters::IsVerbose()
{
	return verboseMode;
}

//! Set Use init magnetization
/*!
	True, if the magnetization should be reset to the initial M vector after each ADC
	This can be set to false if the readout is simulated as well
	e.g. for MRF sequences where a transient magnetization is important
	\param initMag true if magnetizytion should be initialized after ADC
*/
void SimulationParameters::SetUseInitMagnetization(bool initMag)
{
	useInitMagnetization = initMag;
}

//! Get Use init magnetization
/*!	\return info about freshly initialized magnetization after ADC */
bool SimulationParameters::GetUseInitMagnetization()
{
	return useInitMagnetization;
}

//! Set number of max pulse samples
/*!
    pulseq samples pulses at each us. For simulation these pulses are 
	resampled to <= maxNumberOfPulseSamples (default = 100)
	\param numSamples max samples for shaped pulses
*/
void SimulationParameters::SetMaxNumberOfPulseSamples(unsigned int numSamples)
{
	maxNumberOfPulseSamples = numSamples;
}

//! Get number of max pulse samples
/*!	\return max samples of shaped pulses */
unsigned int SimulationParameters::GetMaxNumberOfPulseSamples()
{
	return maxNumberOfPulseSamples;
}

//! Get fit data
/*!	\return pointer to fit data */
std::vector<double>* SimulationParameters::GetFitData()
{
	return &fitData;
}

//! Get fit data weights
/*!	\return pointer to fit data */
std::vector<double>* SimulationParameters::GetFitDataWeights()
{
	return &weights;
}

//! Get ADC positions
/*!	\return pointer to vector with indices of adc events */
std::vector<long>* SimulationParameters::GetADCPositions()
{
	return &adcPos;
}


//! Get the scaling of the initial M
/*!	\return scale of magnetization vector */
double SimulationParameters::GetMagnetizationVectroScale()
{
	return scale;
}

//! Set the scaling of the initial M
void SimulationParameters::SetMagnetizationVectroScale(double s)
{
	scale = s;
}

//! Decode the unique pulses from the seq file
void SimulationParameters::DecodeSeqInfo()
{
	std::vector<std::pair<int, int>> uniquePairs;
	for (unsigned int nSample = 0; nSample < sequence.GetNumberOfBlocks(); nSample++)
	{
		SeqBlock* seqBlock = sequence.GetBlock(nSample);
		if (seqBlock->isRF())
		{
			RFEvent rf = seqBlock->GetRFEvent();
			// ake unique mag phase pair 
			auto p = std::make_pair(rf.magShape, rf.phaseShape);
			if (!(std::find(uniquePairs.begin(), uniquePairs.end(), p) != uniquePairs.end())) {
				// register pulse
				std::vector<PulseSample> uniqueSamples;
				// get rf and check its length
				sequence.decodeBlock(seqBlock);
				unsigned int rfLength = seqBlock->GetRFLength();
				// check arrays of uncompresed shape
				std::vector<float> amplitudeArray(seqBlock->GetRFAmplitudePtr(), seqBlock->GetRFAmplitudePtr() + rfLength);
				std::vector<float> phaseArray(seqBlock->GetRFPhasePtr(), seqBlock->GetRFPhasePtr() + rfLength);
				// rfDeadTime is usually zeros at the end of the pulse, we search for them here
				int nEnd;
				int delayAfterPulse = 0;
				for (nEnd = rfLength; nEnd > 0; --nEnd) {
					if (fabs(amplitudeArray[nEnd - 1]) > 1e-6)// because of the round-up errors in the ascii and derivative/integral reconstructuion
						break;
				}
				delayAfterPulse = rfLength - nEnd;
				rfLength = nEnd;

				amplitudeArray.erase(amplitudeArray.end() - delayAfterPulse, amplitudeArray.end());
				phaseArray.erase(phaseArray.end() - delayAfterPulse, phaseArray.end());
				// search for unique samples in amplitude and phase
				std::vector<float> amplitudeArrayUnique(rfLength);
				std::vector<float>::iterator it_amplitude = std::unique_copy(amplitudeArray.begin(), amplitudeArray.end(), amplitudeArrayUnique.begin());
				amplitudeArrayUnique.resize(std::distance(amplitudeArrayUnique.begin(), it_amplitude));
				std::vector<float> phaseArrayUnique(rfLength);
				std::vector<float>::iterator it_phase = std::unique_copy(phaseArray.begin(), phaseArray.end(), phaseArrayUnique.begin());
				phaseArrayUnique.resize(std::distance(phaseArrayUnique.begin(), it_phase));
				//
				float rfAmplitude = 0.0;
				float rfPhase = 0.0;
				float rfFrequency = seqBlock->GetRFEvent().freqOffset;
				float timestep;
				// need to resample pulse
				unsigned int max_samples = std::max(amplitudeArrayUnique.size(), phaseArrayUnique.size());
				if (max_samples > maxNumberOfPulseSamples) {
					int sampleFactor = ceil(float(rfLength) / maxNumberOfPulseSamples);
					float pulseSamples = rfLength / sampleFactor;
					timestep = float(sampleFactor) * 1e-6;
					// resmaple the original pulse with max ssamples and run the simulation
					uniqueSamples.resize(pulseSamples);
					for (int i = 0; i < pulseSamples; i++) {
						uniqueSamples[i].magnitude = seqBlock->GetRFAmplitudePtr()[i*sampleFactor] * seqBlock->GetRFEvent().amplitude;
						uniqueSamples[i].phase = seqBlock->GetRFPhasePtr()[i*sampleFactor] + seqBlock->GetRFEvent().phaseOffset;
						uniqueSamples[i].timestep = timestep;
					}
				}
				else {
					std::vector<unsigned int>samplePositions(max_samples + 1);
					unsigned int sample_idx = 0;
					if (amplitudeArrayUnique.size() >= phaseArrayUnique.size()) {
						std::vector<float>::iterator it = amplitudeArray.begin();
						for (it_amplitude = amplitudeArrayUnique.begin(); it_amplitude != amplitudeArrayUnique.end(); ++it_amplitude) {
							it = std::find(it, amplitudeArray.end(), *it_amplitude);
							samplePositions[sample_idx++] = std::distance(amplitudeArray.begin(), it);
						}
					}
					else {
						std::vector<float>::iterator it = phaseArray.begin();
						for (it_phase = phaseArrayUnique.begin(); it_phase != phaseArrayUnique.end(); ++it_phase) {
							it = std::find(it, phaseArray.end(), *it_phase);
							samplePositions[sample_idx++] = std::distance(phaseArray.begin(), it);
						}
					}
					uniqueSamples.resize(max_samples);
					samplePositions[max_samples] = rfLength;
					// now we have the duration of the single samples -> simulate it
					for (int i = 0; i < max_samples; i++) {
						uniqueSamples[i].magnitude = seqBlock->GetRFAmplitudePtr()[samplePositions[i]] * seqBlock->GetRFEvent().amplitude;
						uniqueSamples[i].phase = seqBlock->GetRFPhasePtr()[samplePositions[i]] + seqBlock->GetRFEvent().phaseOffset;
						uniqueSamples[i].timestep = (samplePositions[i + 1] - samplePositions[i])*1e-6;
					}
				}
				uniquePairs.push_back(p);
				uniquePulses.insert(std::make_pair(p, uniqueSamples));
			}
		}
		else if (seqBlock->isADC()) {
			adcPos.push_back(nSample);
		}
		delete seqBlock;
	}
}

//! Get a unique pulse
/*!	
    \param pair a pair containing the magnitude and phase id of the seq file
    \return pointer to vector containing the pulse samples of a unique pulse
*/
std::vector<PulseSample>* SimulationParameters::GetUniquePulse(std::pair<int, int> pair)
{
	std::map<std::pair<int, int>, std::vector<PulseSample>>::iterator it;
	it = uniquePulses.find(pair);
	return &(it->second);
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
bool SimulationParameters::RegisterFitParameter(std::string name, double start, double upper, double lower)
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
			fp.set = std::bind(&SimulationParameters::SetScannerB0Inhom, this, _1);
			fp.get = std::bind(&SimulationParameters::GetScannerB0Inhom, this);
		}
		else if(seglist.at(0).compare("scale") == 0) {
			fp.set = std::bind(&SimulationParameters::SetMagnetizationVectroScale, this, _1);
			fp.get = std::bind(&SimulationParameters::GetMagnetizationVectroScale, this);
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
std::vector<FitParameter>* SimulationParameters::GetFitParams()
{
	return fitParams.size() == 0 ? NULL : &fitParams;
}