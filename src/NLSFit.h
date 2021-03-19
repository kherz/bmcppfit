//!  NLSFit.h
/*!
Contains the non-linear least squares ceres fitting stuff

Kai Herz, 2021
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
