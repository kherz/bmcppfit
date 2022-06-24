/*!
Entry funtion

Author: Kai Herz <kai.herz@tuebingen.mpg.de>

Copyright 2021 Kai Herz

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

	// check if seq filename was read
	if (ymlParser.GetSeqFilename().empty()) {
		cout << "ERROR: Could not find .seq filename in parameter file." << endl;
		return EXIT_FAILURE;
	}

	// load sequence
	FitFramework fitFramework(sfp);
	if(!fitFramework.SetExternalSequence(ymlParser.GetSeqFilename())) {
		cout << "ERROR: Could not load .seq file" << endl;
		return EXIT_FAILURE;
	}

	// check if seq file is valid for fit data
	if (sfp.GetADCPositions()->size() != sfp.GetFitData()->size()) {
		cout << "ERROR: Number of ADC events in .seq file (" << sfp.GetADCPositions()->size()  << ") and fit data points (" << sfp.GetFitData()->size() << ") do not match!" << endl;
		return EXIT_FAILURE;
	}

	// ---------- everything is read, now run the fit ----------//
	//init all variables
	Problem problem;                                                   // houston? 
	ceres::Solver::Summary summary;                                    // summary for output
	CostFunctor* costFunc = new CostFunctor();                         // the cost functor where the simulation takes place
	costFunc->fitFramework = &fitFramework;                            // set the simulation parameters
	DynamicNumericDiffCostFunction < CostFunctor >* costFunctionDiff = 
		new DynamicNumericDiffCostFunction < CostFunctor >(costFunc);  // the numeric diff function

    // save the initial start values in a ceres compatible vector
	std::vector<FitParameter>* fp = fitFramework.GetSimFitParameters()->GetFitParams(); 
	std::vector<double> x;
	for (int i = 0; i < fp->size(); i++) {
		x.push_back(fp->at(i).get());
	}
		
	// set the parameters and residuals for the differentiator
	costFunctionDiff->AddParameterBlock(x.size()); 
	costFunctionDiff->SetNumResiduals(fitFramework.GetSimFitParameters()->GetFitData()->size()); // a single residual for each z-spectrum entry
	
	// add the information to the problem
	problem.AddResidualBlock(costFunctionDiff, NULL, x.data());
	for (int i = 0; i < fp->size(); i++) {
		problem.SetParameterLowerBound(x.data(), i, fp->at(i).lower);
		problem.SetParameterUpperBound(x.data(), i, fp->at(i).upper);
	}

	// run the fit
	ceres::Solve(options, &problem, &summary);

	// output
	std::cout << summary.BriefReport() << "\n";

	// write results to yml file
	if (!ymlParser.WriteFitResult(outFile, costFunc->fitFramework->GetSimFitParameters()->GetFitParams()))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}