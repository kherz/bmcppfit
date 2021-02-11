//!  BMSim.cpp
/*!
Bloch-McConnell Z-Spectrum simulation

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


#include "YamlIO.h"
#include "NLSFit.h"
#include "BMSim_T.h"


//! 
/*!
    \param nlhs yamlIn Input parameter file 
	\param nlhs yamlOut Output parameter file
*/
bool RunBMSimFit(std::string yamlIn, std::string yamlOut)
{
	//init the simulation interface and read the input
	SimulationParameters sp;
	if (!ParseYamlInputStruct(yamlIn, sp))
		return EXIT_FAILURE;

	// check if fit data has the same length as seq file has adc events
	if (sp.GetADCPositions()->size() != sp.GetFitData()->size())
	{
		std::cout << "Error: pulseq file and fit data have different number of data points!" << std::endl;
		return EXIT_FAILURE;
	}
		
	NLSFit fit(sp.GetADCPositions()->size());
	for (int np = 0; np < sp.GetNumberOfCESTPools(); np++)
	{
		double k = sp.GetCESTPool(np)->GetExchangeRateInHz();
		fit.AddFreeParameter(k, k / 2, 2 * k);
	}
	
	/* For a small number of pools the matrix size can be set at compile time. This ensures allocation on the stack and therefore a faster simulation. 
	   This speed advantade vanishes for more pools and can even result in a stack overflow for very large matrices
	   In this case more than 3 pools are simulated with dynamic matrices, but this could be expanded eventually
	*/
	switch (sp.GetNumberOfCESTPools())
	{
	case 0:
		sp.IsMTActive() ? fit.Init(BMSim_T<4>, sp) : fit.Init(BMSim_T<3>, sp); // only water
		break;
	case 1:
		sp.IsMTActive() ? fit.Init(BMSim_T<7>, sp) : fit.Init(BMSim_T<6>, sp); // one cest pool
		break;
	//case 2:
	//	sp.IsMTActive() ? Sim_pulseqSBB_T<10>(sp) : Sim_pulseqSBB_T<9>(sp); // two cest pools
	//	break;
	//case 3:
	//	sp.IsMTActive() ? Sim_pulseqSBB_T<13>(sp) : Sim_pulseqSBB_T<12>(sp); // three cest pools
	//	break;
	default:
		fit.Init(BMSim_T<Dynamic>, sp); // > three pools
		break;
	}
	return EXIT_SUCCESS;
}
