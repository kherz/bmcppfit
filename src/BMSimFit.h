//!  main.cpp
/*!
Sequential procedure of doing the fit

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
		
	switch (sp.GetNumberOfCESTPools())
	{
	case 0:
	    sp.IsMTActive() ? RunFit<4>(sp) : RunFit<3>(sp); // no cest pool
		break;
	case 1:
		sp.IsMTActive() ? RunFit<7>(sp) : RunFit<6>(sp); // one cest pool
		break;
	case 2:
		sp.IsMTActive() ? RunFit<10>(sp) : RunFit<9>(sp); // two cest pools
		break;
	default:
		RunFit<Dynamic>(sp);
		break;
	}

	// write the result
	if (!WriteFitResult(yamlOut, sp.GetFitParams()))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
