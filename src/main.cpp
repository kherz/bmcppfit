/*!
Entry funtion

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


#include <iostream>
#include "argh.h"
#include "BMSimFit.h"
#include "YamlParser.h"
#include "ceres/ceres.h"
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::DynamicNumericDiffCostFunction;
using ceres::Problem;

using namespace std;

int main(int argc, char** argv)
{
	// lets parse the input arguments
	auto cmdl = argh::parser(argc, argv);

	// help wished?
	if(cmdl[{ "-h", "--help" }]) {
		cout << "A command-line based BlochMcConnell fit routine for fast multi-thread fitting." << endl;
		cout << "Please provide the following parameters with -<flag>=<param>" << endl;
		cout << "-p: .yaml-file with all specifications. See example file for help" << endl;
		cout << "-o: output .yaml file with fitted parameters (optional)" << endl;
		cout << "-f: .yaml-file with parameters for the fit algorithm (optional) " << endl;
		return EXIT_SUCCESS;
	}

	// check for yaml params
	string paramFile;
	if (!cmdl( "-p")) {
		cout << "Parameter file not provided" << endl;
		return EXIT_FAILURE;
	}
	else {
		paramFile = cmdl("-p").str();
		cout << "Using parameter file " << paramFile << '\n';
	}

	// output file
	string outFile;
	if (!cmdl("-o")) {
		outFile = "out.yaml";
		cout << "Parameter file not provided. Using standard file " << outFile << endl;
	}
	else {
		outFile = cmdl("-o").str();
		cout << "Writing output to " << outFile << '\n';
	}

	// fit options
	YamlParser ymlParser;
	ceres::Solver::Options options;
	string fitOptionsFile = "";
	if (!cmdl("-f")) {
		cout << "Fit options file not provided, using standard parameters" << endl;
		options.minimizer_progress_to_stdout = true;
	}
	else {
		fitOptionsFile = cmdl("-f").str();
		cout << "Using fit options file " << paramFile << '\n';
		if(!ymlParser.ParseYamlCeresOptions(fitOptionsFile, options))
			return EXIT_FAILURE;
	}
	
	// get parameters from yaml file
	SimFitParameters sfp;
	if (!ymlParser.ParseYamlInputStruct(paramFile, sfp))
		return EXIT_FAILURE;

	// init cost functor
	FitFramework fitFramework(sfp);
	if(!ymlParser.ParseSequenceFileName(paramFile, fitFramework))
		return EXIT_FAILURE;

	// run the fit
	Problem problem; // houston? 

    //prepare fit data
	std::vector<FitParameter>* fp = fitFramework.GetSimFitParameters()->GetFitParams(); // thats the fit parameters
	int num_params = fp->size();
	std::vector<double> x; // lets make them ceres compatible
	for (int i = 0; i < num_params; i++)
		x.push_back(fp->at(i).get());

	// prepare cost function
	auto costFunc = new CostFunctor(); // the cost functor where the simulation takes place
	costFunc->fitFramework = &fitFramework; // set the simulation parameters
	auto costFunctionDiff = new DynamicNumericDiffCostFunction < CostFunctor >(costFunc);
	costFunctionDiff->AddParameterBlock(x.size()); // number of parametersparameters
	costFunctionDiff->SetNumResiduals(fitFramework.GetSimFitParameters()->GetFitData()->size()); // a single residual for each z-spectrum entry
	problem.AddResidualBlock(costFunctionDiff, NULL, x.data());

	// set upper and lower boundaries
	for (int i = 0; i < num_params; i++) {
		problem.SetParameterLowerBound(x.data(), i, fp->at(i).lower);
		problem.SetParameterUpperBound(x.data(), i, fp->at(i).upper);
	}

	// init options
	ceres::Solver::Summary summary;

	// run the fit
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";

	if (!ymlParser.WriteFitResult(outFile, costFunc->fitFramework->GetSimFitParameters()->GetFitParams()))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}