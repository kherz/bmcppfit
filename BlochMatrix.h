//!  BlochMatrix.h 
/*!
Container for pool and pulse related parameters

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
#include "SimulationParameters.h"

//!  BlochMatrix class. 
/*!
    Container for pool and pulse related parameters
*/
template <int size> class BlochMatrix
{
public:
	//! Constructor
	/*!
	  Sets fixed parameters such as fraction, exchange rate etc.
	  \param sp SimulationParamter object containing pool informations
	*/
	BlochMatrix(SimulationParameters &sp)
	{
		int N = sp.numberOfCESTPools;
		if (size == Dynamic)
		{
			A.resize(sp.Mvec.rows(), sp.Mvec.rows()); // allocate space for dynamic matrices
		}
		A.setConstant(0.0); // init A

		// MT
		double k_ac = 0.0; // init with 0 for late
		if (sp.simulateMTPool)
		{
			double k_ca = sp.mtPool.k;
			k_ac = k_ca * sp.mtPool.f;
			A(2 * (N + 1), 3 * (N + 1)) = k_ca;
			A(3 * (N + 1), 2 * (N + 1)) = k_ac;
		}

		//WATER
		double k1a = sp.waterPool.R1 + k_ac;
		double k2a = sp.waterPool.R2;
		for (int i = 0; i < N; i++)
		{
			double k_ai = sp.cestPool[i].f * sp.cestPool[i].k;
			k1a += k_ai;
			k2a += k_ai;
		}
		A(0, 0) = -k2a;
		A(1 + N, 1 + N) = -k2a;
		A(2 + 2 * N, 2 + 2 * N) = -k1a;
		
		//CEST POOLS
		for (int i = 0; i < N; i++)
		{
			double k_ia = sp.cestPool[i].k;
			double k_ai = sp.cestPool[i].f * sp.cestPool[i].k;
			double k1i = sp.cestPool[i].R1 + k_ia;
			double k2i = sp.cestPool[i].R2 + k_ia;
			;
			// 1st submatrix
			A(0, i + 1) = k_ia;
			A(i + 1, 0) = k_ai;
			A(i + 1, i + 1) = -k2i;

			// 2nd Submatrix
			A(1 + N, i + 2 + N) = k_ia;
			A(i + 2 + N, 1 + N) = k_ai;
			A(i + 2 + N, i + 2 + N) = -k2i;

			//3rd submatrix
			A(2 * (N + 1), i + 1 + 2 * (N + 1)) = k_ia;
			A(i + 1 + 2 * (N + 1), 2 * (N + 1)) = k_ai;
			A(i + 1 + 2 * (N + 1), i + 1 + 2 * (N + 1)) = -k1i;
		}
	}

	//! Copy Constructor
	/*!
	  \param nA BlochMatrix object 
	*/
	BlochMatrix(BlochMatrix& nA)
	{
		A = nA.GetBlochMatrix();
	}


	//! Update Matrix with pulse parameters
	/*!
	  \param sp SimulationParameters object
	  \param rfAmplitude Amplitude of pulse sample omega1 [rad/s]
	  \param rfFrequency Offset frequency of current pulse sample [rad/s]
	  \param dw0 B0 inhomogeneity [rad/s]
	  \param Rrfc Rrfb/(w1^2) -> see doi:10.1088/0031-9155/58/22/R221
	*/
	void UpdateA(SimulationParameters &sp, double rfAmplitude, double rfFrequency, double rfPhase, double dw0, double Rrfc = 0.0)
	{
		int N = sp.numberOfCESTPools;

		double rfAmplitudeCosPhi = rfAmplitude * cos(rfPhase);
		double rfAmplitudeSinPhi = rfAmplitude * sin(rfPhase);


		// Set pulse
		//water
		A(0, 2 * (N + 1)) = -rfAmplitudeSinPhi;
		A(2 * (N + 1), 0) = rfAmplitudeSinPhi;
		A(N + 1, 2 * (N + 1)) = rfAmplitudeCosPhi;
		A(2 * (N + 1), N + 1) = -rfAmplitudeCosPhi;

		//CEST 
		for (int i = 1; i <= N; i++)
		{
			A(i, i + 2 * (N + 1)) = -rfAmplitudeSinPhi;
			A(i + 2 * (N + 1), i) = rfAmplitudeSinPhi;
			A(N + 1 + i, i + 2 * (N + 1)) = rfAmplitudeCosPhi;
			A(i + 2 * (N + 1), N + 1 + i) = -rfAmplitudeCosPhi;
		}

		// Set off-resonance
		//water
		A(0, 1 + N) = -dw0;
		A(1 + N, 0) = dw0;
		if (abs(rfAmplitude) > 0) {
			A(0, 1 + N) -= rfFrequency;
			A(1 + N, 0) += rfFrequency;
		}
		//cest
		for (int i = 1; i <= N; i++)
		{
			double dwi = sp.cestPool[i - 1].dw*sp.scanner.Gamma*sp.scanner.B0 - rfFrequency - dw0;
			A(i, i + N + 1) = dwi;
			A(i + N + 1, i) = -dwi;
		}

		// Set MT related paramter
		if (sp.simulateMTPool) {
			A(3 * (N + 1), 3 * (N + 1)) = -sp.mtPool.R1 -sp.mtPool.k - Rrfc;
		}
	}

	//! Get the Current matrix
	/*! \return A*/
	Matrix <double, size, size> GetBlochMatrix() { return A; }

private:
	Matrix <double, size, size> A; /*!< Matrix containing pool and pulse paramters */
};
