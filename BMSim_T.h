//!  BMSim_T.h 
/*!
Implementation of the Z-Spectrum Bloch-McConnell simulation for N pools

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

//#include "BlochMatrix.h"
#include "BlochMcConnellSolver.h"

//! Runs the Z-spectrum simulation
/*!
   \param sp SimulationParameters object containing pool and pulse info
*/
template <int size> struct CostFunctor {
	SimulationParameters* sp;
	bool operator()(double const* const* x, double* residual) const
	{

		// update params
		//for (int np = 0; np < sp->GetNumberOfCESTPools(); np++)
		//{
		//	//gsl_vector_set(x, np, max(min(gsl_vector_get(x, np), fp->at(np).upper), fp->at(np).lower));
		//	sp->GetCESTPool(np)->SetExchangeRateInHz(gsl_vector_get(x, np));
		//	//std::cout << sp->GetCESTPool(np)->GetExchangeRateInHz() << std::endl;  
		//}
		sp->GetCESTPool(0)->SetExchangeRateInHz((*x)[0]);
		sp->GetCESTPool(1)->SetExchangeRateInHz((*x)[1]);
		//std::cout << "k1 " << (*x)[0] << ", k2 " << (*x)[1] << std::endl;
		//sp->GetCESTPool(0) =;
		//residual[0] = 0;

		BlochMcConnellSolver<size> bm_solver = BlochMcConnellSolver<size>(*sp);
//#pragma omp parallel for
		for (int i = 0; i < sp->GetADCPositions()->size(); i++)
		{
			double accummPhase = 0;
			Matrix<double, size, 1> M = sp->GetMagnetizationVectors()->col(i);
			// parfor is poosible here
			int startIdx = i == 0 ? 0 : sp->GetADCPositions()->at(i - 1) + 1;
			for (int j = startIdx; j <= sp->GetADCPositions()->at(i); j++)
			{
				SeqBlock* seqBlock = sp->GetExternalSequence()->GetBlock(j);
				if (seqBlock->isADC()) {
					double Zdiff = M(2 * (sp->GetNumberOfCESTPools() + 1)) - sp->GetFitData()->at(i);
					residual[i] = Zdiff;
				}
				else if (seqBlock->isTrapGradient(0) && seqBlock->isTrapGradient(1) && seqBlock->isTrapGradient(2)) {
					for (int i = 0; i < (sp->GetNumberOfCESTPools() + 1) * 2; i++)
						M[i] = 0.0;
				}
				else if (seqBlock->isRF())
				{
					auto p = std::make_pair(seqBlock->GetRFEvent().magShape, seqBlock->GetRFEvent().phaseShape);
					std::vector<PulseSample>* pulseSamples = sp->GetUniquePulse(p);
					double rfFrequency = seqBlock->GetRFEvent().freqOffset;
					for (int p = 0; p < pulseSamples->size(); p++)
					{
						bm_solver.UpdateBlochMatrix(*sp, pulseSamples->at(p).magnitude, rfFrequency, pulseSamples->at(p).phase - accummPhase);
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

				delete seqBlock;
			}
		}

		return true;
	}
};
