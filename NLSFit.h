//!  NLSFit.h
/*!
Contains the gsl fitting stuff

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

#include "BMSim_T.h"
#include "ceres/ceres.h"
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::DynamicNumericDiffCostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


//template <int size> class NLSFit
//{
//public:
//
//	NLSFit(){}
//	~NLSFit(){}
//	bool RunFit(SimulationParameters &sp);
//	std::vector<FitParameter> fitResult;
//};

template <int size> bool RunFit(SimulationParameters &sp)
{
	// Build the problem.
	Problem problem;

	std::vector<FitParameter>* fp = sp.GetFitParams();
	int num_params = fp->size();
	std::vector<double> x;
	for (int i = 0; i < num_params; i++)
		x.push_back(fp->at(i).get());
	// Set up the only cost function (also known as residual). This uses
	// numeric differentiation to obtain the derivative (jacobian).
	CostFunctor<size>* cf = new CostFunctor<size>();
	cf->sp = &sp;
	DynamicNumericDiffCostFunction<CostFunctor<size>>* cost_function =
		new DynamicNumericDiffCostFunction<CostFunctor<size>>(cf);

	cost_function->AddParameterBlock(x.size());
	cost_function->SetNumResiduals(sp.GetFitData()->size());

	problem.AddResidualBlock(cost_function, NULL, x.data());
	for (int i = 0; i < num_params; i++) {
		problem.SetParameterLowerBound(x.data(), i, fp->at(i).lower);
		problem.SetParameterUpperBound(x.data(), i, fp->at(i).upper);
	}

	// Run the solver!
	Solver::Options options;
	options.minimizer_progress_to_stdout = true;
	Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";
	return 0;
}
