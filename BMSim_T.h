//!  BMSim_T.h 
/*!
Implementation of the Z-Spectrum Bloch-McConnell simulation for N pools as a CostFunctor for the seres solver

Kai Herz, 2019
kai.herz@tuebingen.mpg.de

**********************************
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************

*/

#pragma once

#include "BlochMcConnellSolver.h"

//! Templated cost functor for the ceres fit
template <int size> struct CostFunctor {
	SimulationParameters* sp; /*!< Water Pool */
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
			BlochMcConnellSolver<size> bm_solver = BlochMcConnellSolver<size>(*sp);
			double accummPhase = 0;
			Matrix<double, size, 1> M = sp->GetMagnetizationVectors()->col(i); // magnetization vector
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
		return true; // TODO: make error checks that retun false
	}
};
