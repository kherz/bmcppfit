//!  SimulationParameters.h
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

#pragma once

#include "SimulationParameters.h"
#include "BMCSim.h"
#include <functional>

//!  SimulationParameters class. 
/*!
	Container for all the relevant simulation parameters

*/

struct FitParameter
{
	std::string name;
	std::function<void(double)> set;
	std::function<double()> get;
	double lower;
	double upper;
};

class SimFitParameters : public SimulationParameters
{
public: // TODO: write get and set methods for member variables and make them private

	//! Get fit data
	std::vector<double>* GetFitData();

	//! Get weights
	std::vector<double>* GetFitDataWeights();

	//! Get ADC Positions
	std::vector<long>* GetADCPositions();

	// OpenMP Thrads
	int numberOfThreads;

	//! Register the fit parameters
	bool RegisterFitParameter(std::string name, double start, double lower, double upper);

	//! Get the registered fit parameters
	std::vector<FitParameter>* GetFitParams();

	//! Get the scaling of the initial M
	double GetMagnetizationVectorScale();

	//! Set the scaling of the initial M
	void SetMagnetizationVectorScale(double s);


protected:
	std::vector<long> adcPos;            /*!< vector with index of ADC events in seq file */
	double scale;                        /*!< scaling of initial magnetization */
	std::vector<double> fitData;         /*!< the Z-spec that shoud be fitted */
	std::vector<double> weights;         /*!< weights for the fit data */
	std::vector<FitParameter> fitParams; /*!< vector containg parameters to fit */


};


class FitFramework : public BMCSim {
public:
	FitFramework(SimFitParameters& sfp);
	FitFramework(FitFramework& ff);
	SimFitParameters* GetSimFitParameters();
	ExternalSequence* GetExternalSequence();
	void RunNonADCEventBlock(Eigen::VectorXd &M, float &accumPhase, SeqBlock* seqBlock);
	bool SetExternalSequence(std::string path);
	void UpdateSolver();

private:
	SimFitParameters* sfp;

};

struct CostFunctor {
	FitFramework* fitFramework;
	bool operator()(double const* const* x, double* residual) const
	{
		for (int p = 0; p < fitFramework->GetSimFitParameters()->GetFitParams()->size(); p++) {
			fitFramework->GetSimFitParameters()->GetFitParams()->at(p).set((*x)[p]);
		}
		fitFramework->UpdateSolver();
#pragma omp parallel for
		for (int i = 0; i < fitFramework->GetSimFitParameters()->GetADCPositions()->size(); i++)
		{
			FitFramework localSimFramework(*fitFramework);
			if (fitFramework->GetSimFitParameters()->GetFitDataWeights()->at(i) < 1e-12) { // only run the stuff if we wnat to consider this offset
				residual[i] = 0.0;
			}
			else {
				// spread solver across thread
				//FitFramework localFitFramework = *fitFramework;
				float accumPhase = 0;
				Eigen::VectorXd M = fitFramework->GetMagnetizationVectors()->col(i); // magnetization vector
				M *= fitFramework->GetSimFitParameters()->GetMagnetizationVectorScale(); // use init scale 
				int startIdx = i == 0 ? 0 : fitFramework->GetSimFitParameters()->GetADCPositions()->at(i - 1) + 1; // find the idx in the pulseq file from which the current offset starts
				for (int j = startIdx; j <= fitFramework->GetSimFitParameters()->GetADCPositions()->at(i); j++)
				{
					SeqBlock* seqBlock = fitFramework->GetExternalSequence()->GetBlock(j);
					if (seqBlock->isADC()) { // we calculate the residuals at the ADC event
						double Zdiff = fitFramework->GetSimFitParameters()->GetFitData()->at(i) - M(2 * (fitFramework->GetSimFitParameters()->GetNumberOfCESTPools() + 1));
						residual[i] = Zdiff * fitFramework->GetSimFitParameters()->GetFitDataWeights()->at(i);
					}
					else { // delay or single gradient -> simulated as delay
						localSimFramework.RunNonADCEventBlock(M, accumPhase, seqBlock);
					}
					delete seqBlock; // gets allocated with new during getBlock() call
				}
			}
		}
		return true; // TODO: make error checks that retun false
	}
};

