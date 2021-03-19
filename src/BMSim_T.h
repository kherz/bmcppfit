//!  BMSim_T.h 
/*!
Implementation of the Z-Spectrum Bloch-McConnell simulation for N pools as a CostFunctor for the seres solver

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

#include "BlochMcConnellSolver.h"

//! Templated cost functor for the ceres fit
template <int size> struct CostFunctor {
	SimulationParameters* sp;
	//! Runs the Z-spectrum simulation and calculates the residuals
    /*!
       \param x ceres parameters
	   \param residual ceres residuals
	   \return operator expects a return value
    */
	bool operator()(double const* const* x, double* residual) const
	{
		// update the parameters
		for (int p = 0; p < sp->GetFitParams()->size(); p++) {
			sp->GetFitParams()->at(p).set((*x)[p]);
		}

		// IMPORTATNT!!!
		// we expect a magnetization vector at the beginning of each offset that is independen from the previous readout. 
		// By doing so, we cann parallelize the the different offsets
		// TODO: tweak openmp call
#pragma omp parallel for
		for (int i = 0; i < sp->GetADCPositions()->size(); i++)
		{
			if (sp->GetFitDataWeights()->at(i) < 1e-12) { // only run the stuff if we wnat to consider this offset
				residual[i] = 0.0; 
			}
			else {
				BlochMcConnellSolver<size> bm_solver = BlochMcConnellSolver<size>(*sp);
				double accummPhase = 0;
				Matrix<double, size, 1> M = sp->GetMagnetizationVectors()->col(i); // magnetization vector
				M *= sp->GetMagnetizationVectorScale(); // use init scale 
				int startIdx = i == 0 ? 0 : sp->GetADCPositions()->at(i - 1) + 1; // find the idx in the pulseq file from which the current offset starts
				for (int j = startIdx; j <= sp->GetADCPositions()->at(i); j++)
				{
					SeqBlock* seqBlock = sp->GetExternalSequence()->GetBlock(j);
					if (seqBlock->isADC()) { // we calculate the residuals at the ADC event
						double Zdiff = sp->GetFitData()->at(i) - M(2 * (sp->GetNumberOfCESTPools() + 1));
						residual[i] = Zdiff * sp->GetFitDataWeights()->at(i);
					}
					else if (seqBlock->isTrapGradient(0) && seqBlock->isTrapGradient(1) && seqBlock->isTrapGradient(2)) { // spoil for all 3 gradients
						for (int i = 0; i < (sp->GetNumberOfCESTPools() + 1) * 2; i++)
							M[i] = 0.0;
					}
					else if (seqBlock->isRF()) { // saturation pulse
						auto p = std::make_pair(seqBlock->GetRFEvent().magShape, seqBlock->GetRFEvent().phaseShape); // get the magnitude and phase pair
						std::vector<PulseSample>* pulseSamples = sp->GetUniquePulse(p); // find the unque rf id in the previously decoded seq file library
						double rfFrequency = seqBlock->GetRFEvent().freqOffset;
						for (int p = 0; p < pulseSamples->size(); p++) { // loop through pulse samples
							bm_solver.UpdateBlochMatrix(*sp, pulseSamples->at(p).magnitude, rfFrequency, pulseSamples->at(p).phase + seqBlock->GetRFEvent().phaseOffset - accummPhase);
							bm_solver.SolveBlochEquation(M, pulseSamples->at(p).timestep);
						}
						int phaseDegree = seqBlock->GetDuration() * 1e-6 * 360 * rfFrequency;
						phaseDegree %= 360;
						accummPhase += double(phaseDegree) / 180 * PI;
					}
					else { // delay or single gradient -> simulated as delay
						double timestep = seqBlock->GetDuration()*1e-6;
						bm_solver.UpdateBlochMatrix(*sp, 0, 0, 0);
						bm_solver.SolveBlochEquation(M, timestep);
					}
					delete seqBlock; // gets allocated with new during getBlock() call
				}
			}
		}
		return true; // TODO: make error checks that retun false
	}
};
