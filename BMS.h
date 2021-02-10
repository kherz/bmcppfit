//!  BMS.h
/*!
Definitions of underlying structs and error codes.

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

// MATLAB includes
#include <math.h>

//Eigen include
#include <Eigen/Eigen>

//open mp for parallelization
#if defined(_OPENMP)
	#include <omp.h>
#else
	#pragma message("OPENMP not defined! Parallelization is not possible")

	#if defined(_MSC_VER)
		#pragma message("Add COMPFLAGS=\"$COMPFLAGS /openmp\" to the mex command")
	#elif defined(__GNUG__)
		#pragma message("Add CXXFLAGS=\"$CXXPFLAGS -fopenmp\" to the mex command")
	#endif

#endif

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI // should be in cmath
#define M_PI 3.14159265358979323846
#endif // !M_PI

//! Namespace of the Bloch McConnell Simulation Parameters
namespace BMS
{
	//! Shape of the magnetization transfer pool
	enum MTLineshape
	{
		SuperLorentzian,
		Lorentzian,
		None
	};

	//! B0 field related info
	struct Scanner
	{
		double B0;                /*!< static field [T]*/
		double B0Inhomogeneity;   /*!< field inhomogeneity [ppm] */
		double Gamma;             /*!< gyromagnetic ratio [rad/uT] */
	};

	//! Parameters of the water pool
	struct WaterPool
	{
		double R1; /*!< 1/T1 [Hz] */
		double R2; /*!< 1/T2 [Hz] */
		double dw; /*!< Offset from Water resonance [ppm] */
		double f;  /*!< Pool fraction */
	};

	//! Parameters of the magnetization transfer pool
	struct MTPool
	{
		double R1; /*!< 1/T1 [Hz] */
		double R2; /*!< 1/T2 [Hz] */
		double dw; /*!< Offset from Water resonance [ppm] */
		double f;  /*!< Pool fraction */
		double k;  /*!< Exchange rate [Hz] */
		MTLineshape lineshape; /*!< MT linehape [Hz] */
	};

	//! Parameters of th cest pool
	struct CESTPool
	{
		double R1; /*!< 1/T1 [Hz] */
		double R2; /*!< 1/T2 [Hz] */
		double dw; /*!< Offset from Water resonance [ppm] */
		double f;  /*!< Pool fraction */
		double k;  /*!< Exchange rate [Hz] */
	};
}