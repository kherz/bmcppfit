//!  SimulationParameters.h
/*!
Main interface for the Bloch-McConnell simulation

Kai Herz, 2018
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

#include "BMS.h"
// include pulseq
#include "ExternalSequence.h"

using namespace BMS;
using namespace Eigen;

//!  SimulationParameters class. 
/*!
	Container for all the relevant simulation parameters
*/
class SimulationParameters
{
public: // TODO: write get and set methods for member variables and make them private

	MatrixXd Mvec;  /*!< Matrix containing all magnetization vectors */

	std::vector<double> PPMVector;			 /*!< vector containing ppm values of the z-spectrum  [ppm] */
	unsigned int numPPMSamples;	 /*!<length of the PPMVector */
	std::vector<double> FitData; /*!< Spectra that should be fitted */
	
	WaterPool waterPool; /*!< Water Pool */
	MTPool mtPool;       /*!< MT Pool */
	CESTPool* cestPool;  /*!< CEST Pool(s) */
	
	bool simulateMTPool;    /*!< true if MT should be simulated */
	double *RrfcContainer;  /*!< container for Rrfb / (w1 ^ 2) at each ppm sample ->see doi : 10.1088 / 0031 - 9155 / 58 / 22 / R221 */
	bool mtMemAllocated;    /*!< true if memory forRrfcContainer was allocated*/

	unsigned int numberOfCESTPools; /*!< number of CEST Pools */
	bool cestMemAllocated;          /*!< true if memory for cest pools was allocated*/

	ExternalSequence seq;
	Scanner scanner;            /*!< B0 field related info*/

	bool verboseMode;              /*!< true, if you want to have some output information */

	unsigned int numberOfThreads;  /*!< threads for openmp */

	//! Constructor
	SimulationParameters()
	{
		numberOfCESTPools = 0;
		simulateMTPool = false;
		mtMemAllocated = false;
		cestMemAllocated = false;
	
	}

	//! Destructor
	~SimulationParameters()
	{
		if (numberOfCESTPools>0 && cestMemAllocated)
			delete[] cestPool;
		
		if (simulateMTPool && mtMemAllocated) {
			delete[] RrfcContainer;
		}
	}

	//! Fills the Rrfc container
	void CalculateMTLineshape()
	{
		if (!mtMemAllocated) { // ceck if mem was already allocated
			RrfcContainer = new double[numPPMSamples];
			mtMemAllocated = true;
		}
		double T2c = (1.0 / mtPool.R2);
		double Gamma = scanner.Gamma;
		switch (mtPool.lineshape)
		{
			case None:
			{
				for (int i = 0; i < numPPMSamples; i++) {
					RrfcContainer[i] = 0.0;
				}
				break;

			}
			case Lorentzian: //the easy case
			{
				for (int i = 0; i < numPPMSamples; i++) {
					double dw0c = (PPMVector[i] - mtPool.dw)* Gamma* scanner.B0;
					RrfcContainer[i] = T2c / (1 + pow(dw0c*T2c, 2.0));
				}
				break;
			}
			case SuperLorentzian: //the more complicated case
			{
				//init temporary variables for superlorentzian interpolation
				const int numU = 1001;
				double uStep = 0.001;
				double sqrt2pi = sqrt(2.0 / M_PI);
				double cutoff = 3 * Gamma;

				//calculate the values that are outside the cutoff frequency
				std::vector<int> jContainer;
				for (int j = 0; j < numPPMSamples; j++) {
					double dw0c = (PPMVector[j] - mtPool.dw) * Gamma* scanner.B0 + scanner.B0Inhomogeneity;
					double sum = 0.0;
					if (abs(dw0c) >= cutoff) {
						for (int i = 0; i < numU; i++) {
							double cu = uStep*double(i);
							double powcu2 = pow(cu, 2.0);
							double absPowCu2 = (3.0 * powcu2 - 1.0) >= 0 ? (3.0 * powcu2 - 1.0) : -(3.0 * powcu2 - 1.0);
							sum += sqrt2pi*T2c / absPowCu2 *exp(-2.0 * pow(dw0c * T2c / absPowCu2, 2.0));
						}
						RrfcContainer[j] = M_PI*uStep*sum;
					}
					else { //interpolate these values later
						jContainer.push_back(j);
						RrfcContainer[j] = 0.0;
					}
				}
				//fill the empty part with spline interpolation
				//cubic hermite spline might work ok
				int interpPoints = jContainer.size();
				if (interpPoints > 0) {
					int p0id = jContainer[0] - 1;
					int d0id = jContainer[0] - 2;
					int p1id = jContainer[interpPoints - 1] + 1;
					int d1id = jContainer[interpPoints - 1] + 2;
					//need at least two points on each side of the cutoff for correct interpolation
					if (d0id < 0 || d1id > numPPMSamples - 1) {
					//	std::cerr("Not enough control points for spline interpolation. Choose Lorentzian Lineshape,\
					//			 						simulate PPM values that are further from dw0c or increase the Samples");
					}
					else {
						//y values
						double p0y = RrfcContainer[p0id]; //points
						double p1y = RrfcContainer[p1id];

						double tangentWeight = (cutoff*0.04); // empirically chosen 
						double d0y = tangentWeight*(p0y - RrfcContainer[d0id]); //tangents
						double d1y = tangentWeight*(RrfcContainer[d1id] - p1y);

						//calculate the interpolation points
						for (int i = 0; i < interpPoints; i++) {
							double cStep = double(i + 1.0) / (interpPoints + 1);
							double c3 = cStep*cStep*cStep;
							double c2 = cStep*cStep;

							// hermite spline
							double h0 = 2 * c3 - 3 * c2 + 1;
							double h1 = -2 * c3 + 3 * (c2);
							double h2 = c3 - 2 * c2 + cStep;
							double h3 = c3 - c2;
							//fill the point container
							RrfcContainer[jContainer[i]] = h0 * p0y + h1 * p1y + h2 * d0y + h3 * d1y;
						}
					}
				}
				break;
			}
			default:
			{
				//mexErrMsgTxt("Please choose a valid Lineshape for the MT Pool! (Lorentzian or SuperLorentzian)");
			}
		}
	}
};


