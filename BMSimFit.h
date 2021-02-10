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

#include "YamlIO.h"
#include "BMSim_T.h"

//! Enry point for MATLAB mex function
/*!
    \param nlhs number of output arguments
	\param plhs Array of pointers to the mxArray output arguments
	\param nrhs number of input arguments
	\param prhs Array of pointers to the mxArray input arguments
*/
bool RunBMSimFit(std::string yamlIn, std::string yamlOut)
{
	//init the simulation interface and read the input
	SimulationParameters sp;
	if (!ParseYamlInputStruct(yamlIn, sp))
		return EXIT_FAILURE;

	//// disp info about input
	//if (sp.verboseMode) {
	//	mexPrintf("Read parameters succesfully! \n");
	//	mexPrintf("Found %i CEST pool(s) and %i MT Pool(s) \n", sp.numberOfCESTPools, sp.simulateMTPool ? 1 : 0);
	//}


	///* For a small number of pools the matrix size can be set at compile time. This ensures allocation on the stack and therefore a faster simulation. 
	//   This speed advantade vanishes for more pools and can even result in a stack overflow for very large matrices
	//   In this case more than 3 pools are simulated with dynamic matrices, but this could be expanded eventually
	//*/
	//switch (sp.numberOfCESTPools)
	//{
	//case 0:
	//	sp.simulateMTPool ? BMSim_T<4>(sp) : BMSim_T<3>(sp); // only water
	//	break;
	//case 1:
	//	sp.simulateMTPool ? BMSim_T<7>(sp) : BMSim_T<6>(sp); // one cest pool
	//	break;
	//case 2:
	//	sp.simulateMTPool ? BMSim_T<10>(sp) : BMSim_T<9>(sp); // two cest pools
	//	break;
	//case 3:
	//	sp.simulateMTPool ? BMSim_T<13>(sp) : BMSim_T<12>(sp); // three cest pools
	//	break;
	//default:
	//	BMSim_T<Dynamic>(sp); // > three pools
	//	break;
	//}

	//ReturnResultToMATLAB(plhs, sp.Mvec); // return results after simulation

}
