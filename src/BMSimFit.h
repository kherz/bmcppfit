//!  BMSimFit.h
/*!
Sequential procedure of doing the fit

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


#include "YamlIO.h"
#include "NLSFit.h"
#include "BMSim_T.h"


//! 
/*!
    \param nlhs yamlIn Input parameter file 
	\param nlhs yamlOut Output parameter file
*/
bool RunBMSimFit(std::string yamlIn, std::string yamlOut, std::string yamlFitOptions)
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

	// init solver options
	ceres::Solver::Options opts;
	SetStandardOprionsParameters(opts);
	if (!yamlFitOptions.empty())
		if(!ParseYamlCeresOptions(yamlFitOptions, opts))
			return EXIT_FAILURE;

	switch (sp.GetNumberOfCESTPools())
	{
	case 0:
	    sp.IsMTActive() ? RunFit<4>(sp, opts) : RunFit<3>(sp, opts); // no cest pool
		break;
	case 1:
		sp.IsMTActive() ? RunFit<7>(sp, opts) : RunFit<6>(sp, opts); // one cest pool
		break;
	case 2:
		sp.IsMTActive() ? RunFit<10>(sp, opts) : RunFit<9>(sp, opts); // two cest pools
		break;
	default:
		RunFit<Dynamic>(sp, opts);
		break;
	}

	// write the result
	if (!WriteFitResult(yamlOut, sp.GetFitParams()))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
