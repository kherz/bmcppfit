#include <iostream>
#include "argh.h"
#include "BMSimFit.h"

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
	string param_file;
	if (!cmdl( "-p")) {
		cout << "Parameter file not provided" << endl;
		return EXIT_FAILURE;
	}
	else {
		param_file = cmdl("-p").str();
		cout << "Using parameter file " << param_file << '\n';
	}

	// output file
	string out_file;
	if (!cmdl("-o")) {
		out_file = "out.yaml";
		cout << "Parameter file not provided. Using standard file " << out_file << endl;
	}
	else {
		out_file = cmdl("-o").str();
		cout << "Writing output to " << out_file << '\n';
	}

	// fit options
	string fit_options_file = "";
	if (!cmdl("-f")) {
		cout << "Fit options file not provided, using standard parameters" << endl;
	}
	else {
		fit_options_file = cmdl("-f").str();
		cout << "Using fit options file " << param_file << '\n';
	}
	
	return RunBMSimFit(param_file, out_file, fit_options_file);
}