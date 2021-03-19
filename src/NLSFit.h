//!  NLSFit.h
/*!
Contains the non-linear least squares ceres fitting stuff

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

#include "BMSim_T.h"
#include "ceres/ceres.h"
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::DynamicNumericDiffCostFunction;
using ceres::Problem;
//using ceres::Solver;


// TODO: Make this a class and give the user the possibility to change fit parameters with some config file!!!

//! Sets some standard options for the fir
/*!
   \param options ceres::Solver::Options object that is used in the fit later
*/
void SetStandardOptionsParameters(ceres::Solver::Options &options)
{
	// todo: add more standards
	options.minimizer_progress_to_stdout = true;
}


//! Runs the fit
/*!
   \param sp SimulationParameters object that contains all the infos
   \return true if success
*/
template <int size> bool RunFit(SimulationParameters &sp, ceres::Solver::Options &options)
{
	Problem problem; // houston? 

	// prepare fit data
	std::vector<FitParameter>* fp = sp.GetFitParams(); // thats the fit parameters
	int num_params = fp->size();
	std::vector<double> x; // lets make them ceres compatible
	for (int i = 0; i < num_params; i++)  
		x.push_back(fp->at(i).get());

	// prepare cost function
	CostFunctor<size>* cf = new CostFunctor<size>(); // the cost functor where the simulation takes place
	cf->sp = &sp; // set the simulation parameters
	DynamicNumericDiffCostFunction<CostFunctor<size>>* cost_function =
		new DynamicNumericDiffCostFunction<CostFunctor<size>>(cf);
	cost_function->AddParameterBlock(x.size()); // number of parametersparameters
	cost_function->SetNumResiduals(sp.GetFitData()->size()); // a single residual for each z-spectrum entry
	problem.AddResidualBlock(cost_function, NULL, x.data());

	// set upper and lower boundaries
	for (int i = 0; i < num_params; i++) {
		problem.SetParameterLowerBound(x.data(), i, fp->at(i).lower);
		problem.SetParameterUpperBound(x.data(), i, fp->at(i).upper);
	}

	// init options
	//Solver::Options options;
	//options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary; 

	// run the fit
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";

	return true; // TODO: make some error checks that reurn false
}
