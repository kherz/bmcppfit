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

#include "BlochMatrix.h"
#include "BlochMcConnellSolver.h"

//! Runs the Z-spectrum simulation
/*!
   \param sp SimulationParameters object containing pool and pulse info
*/
template <int size> void BMSim_T(SimulationParameters& sp)
{
	typedef Matrix<double, size, 1> VectorNd; 
	typedef Matrix<double, size, size> MatrixNd;

	BlochMatrix<size> ABase(sp); // init base matrix with fixed pool parameters
	VectorNd C; // vector containing relaxation paremeters

	if (size == Dynamic){ // alocate space for dynamic matrices
		C.resize(sp.Mvec.rows());
	}

	// set entries
	C.setConstant(0.0);
	int nP = sp.numberOfCESTPools;
	C((nP + 1) * 2) = sp.waterPool.f * sp.waterPool.R1; // water
	for (int p = 0; p < nP; p++) { // vest
		C((nP + 1) * 2 + (p+1)) = sp.cestPool[p].f*sp.cestPool[p].R1;
	}

	if (sp.simulateMTPool) { // set MT related info
		C(3*(nP + 1)) = sp.mtPool.f*sp.mtPool.R1;
		sp.CalculateMTLineshape();
	}

	double w0 = sp.scanner.B0*sp.scanner.Gamma; // omega0 [rad/s]
	double b0inhomogeneity = sp.scanner.B0Inhomogeneity*w0; //b0 inhomogeneity

#if defined(_OPENMP)
	//set number of threads for openmp
	omp_set_num_threads(sp.numberOfThreads);
#endif

//calculate the different samples of the Z-spectrum parallel
#pragma omp parallel for
	for (int k = 0; k < sp.numPPMSamples; k++) {
		double dw0 = sp.PPMVector[k] * w0; //current frequency offset [rad/s]
		BlochMatrix<size> A(ABase); // get copy of bloch matrix
		
		int l = 0; // loop variable
		while (l < sp.pulseTrain.NumberOfSaturationPulses) { // loop through pulses
			for (int m = 0; m < sp.pulseShape.NumberOfSamples; m++)	{ // loop through pulse samples
				double rfAmplitude  = sp.pulseShape.Amplitude[m]; // current w1 [rad/s]
				double rfTimestep   = sp.pulseShape.Timesteps[m];  // duration of current pulse sample [s]
				double rfFrequency  = sp.pulseShape.hasFrequencyOffset ? sp.pulseShape.Frequency[m] : 0.0; // current pulse frequency offset [rad/s]
				double rfPhase      = sp.pulseShape.hasPhaseOffset     ? sp.pulseShape.Phase[m]     : 0.0; // current pulse phase offset [rad]
				A.UpdateA(sp, rfAmplitude, rfFrequency + dw0, rfPhase, b0inhomogeneity, sp.simulateMTPool ? pow(rfAmplitude, 2)*sp.RrfcContainer[k] : 0.0); // set pulse
				sp.Mvec.col(k) = SolveBlochEquation<size>(sp.Mvec.col(k), A.GetBlochMatrix(), C, rfTimestep); // and solve the equation
			}

			if (++l < sp.pulseTrain.NumberOfSaturationPulses) { // pause between pulses
				if (sp.pulseTrain.Spoiling)	{ // "spoiling"
					for (int p = 0; p < (nP+1)*2; p++) { //spoiling is simulated by setting all the transverse magnetization to 0
						sp.Mvec(p, k) = 0;
					}
				}
				A.UpdateA(sp, 0.0, 0.0, 0.0, b0inhomogeneity); // set pulse to 0 during pause
				sp.Mvec.col(k) = SolveBlochEquation<size>(sp.Mvec.col(k), A.GetBlochMatrix(), C, sp.pulseTrain.PauseBetweenPulses); // and solve the equation
			}
		}
	}
}

