//!  YamlIO.h
/*!
Read/write to yaml files

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
#include "YamlParser.h"

//! Read the parameter yaml file
/*!
   \param yamlIn file name of the yaml file
   \param sp SimulationParameters object by reference that gets filled
   \return true if success
*/
bool YamlParser::ParseYamlInputStruct(std::string yamlIn, SimFitParameters &sp)
{
	YAML::Node config;
	// try to read the file
	try
	{
		config = YAML::LoadFile(yamlIn);
	}
	catch (...)
	{
		std::cout << "Could not load yaml file " << yamlIn << std::endl;
		return false;
	}

	// get data to fit
	try
	{
		auto fit_data = config["fit_data"];
		auto n = fit_data.size();
		sp.GetFitData()->resize(n);
		for (int sample = 0; sample < n; sample++)
			sp.GetFitData()->at(sample) = (fit_data[sample].as<double>());
	}
	catch (...)
	{
		std::cout << "Could not read fit data" << std::endl;
		return false;
	}
	// init weights
	sp.GetFitDataWeights()->resize(sp.GetFitData()->size());
	std::fill(sp.GetFitDataWeights()->begin(), sp.GetFitDataWeights()->end(), 1.0);
	try
	{
		auto weights = config["weights"];
		if (weights.IsDefined()) {
			auto n = weights.size();
			if (n == sp.GetFitDataWeights()->size()) { // fill the vector with the weights
				for (int sample = 0; sample < n; sample++)
					sp.GetFitDataWeights()->at(sample) = (weights[sample].as<double>());
			}
			else {
				std::cout << "Size of weights (" << n << ") does not match with size of fit data ( " << sp.GetFitDataWeights()->size() << ")" << std::endl;
				std::cout << "Ignoring weights!" << std::endl;
			}

		}
	}
	catch (...)
	{
		std::cout << "Could not read weights" << std::endl;
		return false;
	}

	// water pool
	try
	{
		auto water_pool = config["water_pool"];
		auto f = water_pool["f"].as<double>();
		auto r1 = 1.0 / water_pool["t1"].as<double>();
		auto r2 = 1.0 / water_pool["t2"].as<double>();
		sp.SetWaterPool(WaterPool(r1, r2, f));
	}
	catch (...)
	{
		std::cout << "Could not read water_pool" << std::endl;
		return false;
	}

	// mt pool
	try
	{
		auto mt_pool = config["mt_pool"];
		if (mt_pool.IsDefined()) {
			auto f = mt_pool["f"].as<double>();
			auto r1 = 1.0 / mt_pool["t1"].as<double>();
			auto r2 = 1.0 / mt_pool["t2"].as<double>();
			auto k = mt_pool["k"].as<double>();
			auto dw = mt_pool["dw"].as<double>();
			std::string sring_ls = mt_pool["lineshape"].as<std::string>();
			if (sring_ls.compare("Lorentzian") == 0) {
				sp.SetMTPool(MTPool(r1, r2, f, dw, k, Lorentzian));
			}
			else if (sring_ls.compare("SuperLorentzian") == 0) {
				sp.SetMTPool(MTPool(r1, r2, f, dw, k, SuperLorentzian));
			}
			else if (sring_ls.compare("None") == 0) {
				sp.SetMTPool(MTPool(r1, r2, f, dw, k, None));
			}
			else {
				std::cout << "No valid MT Lineshape! Use None, Lorentzian or SuperLorentzian" << std::endl;
				return false;
			}
		}
	}
	catch (...)
	{
		std::cout << "Could not read mt_pool" << std::endl;
		return false;
	}

	// cest pools
	try
	{
		auto cest_pools = config["cest_pool"];
		if (cest_pools.IsDefined()) {
			sp.SetNumberOfCESTPools(cest_pools.size());
			int i = 0;
			for (YAML::const_iterator it = cest_pools.begin(); it != cest_pools.end(); ++it) {
				auto c_pool = it->second.as<YAML::Node>();
				auto r1 = 1.0 / c_pool["t1"].as<double>();
				auto r2 = 1.0 / c_pool["t2"].as<double>();
				auto f = c_pool["f"].as<double>();
				auto k = c_pool["k"].as<double>();
				auto dw = c_pool["dw"].as<double>();
				sp.SetCESTPool(CESTPool(r1, r2, f, dw, k), i++);
			}
		}
	}
	catch (...)
	{
		std::cout << "Could not read cest_pool" << std::endl;
		return false;
	}

	// as soon as we know how many pools, we can read the fit params
	try
	{
		auto fit_params = config["fit_params"];
		for (YAML::const_iterator it = fit_params.begin(); it != fit_params.end(); ++it) {
			auto name = it->first.as<YAML::Node>();
			auto params = it->second.as<YAML::Node>();
			if (!sp.RegisterFitParameter(name.as<std::string>(), params["start"].as<double>(), params["upper"].as<double>(), params["lower"].as<double>()))
				return false;
		}
	}
	catch (...)
	{
		std::cout << "Could not read fit parameters" << std::endl;
		return false;
	}

	// get scaling of vector
	try
	{
		auto scaleNode = config["scale"];
		if (scaleNode.IsDefined()) {
			sp.SetMagnetizationVectorScale(scaleNode.as<double>());
		}
	}
	catch (...)
	{
		std::cout << "Could not read scale" << std::endl;
		return false;
	}

	// put together magnetization vector
	int ncp = sp.GetNumberOfCESTPools();
	int vecSize = (3 * (ncp + 1)) + (sp.IsMTActive() ? 1 : 0);
	Eigen::VectorXd M(vecSize);
	M.fill(0.0);
	M(2 * (ncp + 1)) = sp.GetWaterPool()->GetFraction();
	for (int i = 0; i < ncp; i++) {
		M(2 * (ncp + 1) + i + 1) = sp.GetCESTPool(i)->GetFraction();
	}
	if (sp.IsMTActive())
		M(vecSize - 1) = sp.GetMTPool()->GetFraction();
	// init with same number as fit data
	sp.SetInitialMagnetizationVector(M);

	// read ther stuff sequence
	try
	{
		auto b0 = config["b0"];
		sp.InitScanner(b0.as<double>());
		auto b1 = config["rel_b1"];
		if (b1.IsDefined()) {
			sp.SetScannerRelB1(b1.as<double>());
		}
	}
	catch (...)
	{
		std::cout << "Could not read scanner stuff" << std::endl;
		return false;
	}

	// get scaling of vector
	int maxSamples = 200;
	try
	{
		auto samplesNode = config["max_pulse_samples"];
		if (samplesNode.IsDefined()) {
			maxSamples = samplesNode.as<double>();
		}
	}
	catch (...)
	{
		std::cout << "Could not read max pulse samples file" << std::endl;
		return false;
	}
	sp.SetMaxNumberOfPulseSamples(maxSamples);


#if defined (_OPENMP)
	//set threads for openmp
	// if there is no specification, the simulation only runs on one thread
	int maxThreads = omp_get_max_threads();
	try {
		auto threads = config["threads"];
		int inputThreads = threads.as<int>();
		if (maxThreads < inputThreads) {
			std::cout << "WARNING: Specified openMP threads " << inputThreads << " are higher than maximum allowed " << maxThreads << "! Using maximum threads instead" << std::endl;
		}
		sp.numberOfThreads = inputThreads <= maxThreads ? inputThreads : maxThreads;
	}
	catch (...) {
		sp.numberOfThreads = 1;
	}
	omp_set_num_threads(sp.numberOfThreads);
#else
	sp.numberOfThreads = 1;
#endif

	// enough for now we come to the rest later
	sp.SetVerbose(false);

	return true;

}

bool YamlParser::ParseSequenceFileName(std::string yamlIn, FitFramework &fitFramework)
{
	YAML::Node config;
	// try to read the file
	try
	{
		config = YAML::LoadFile(yamlIn);
	}
	catch (...)
	{
		std::cout << "Could not load yaml file " << yamlIn << std::endl;
		return false;
	}

	// read pulseq sequence
    // if no path is specified for the pulseq file it is assumed to be in the same folder as the yaml file
	try
	{
		auto seq_fn = config["pulseq_file"];
		std::string seqfile = seq_fn.as<std::string>();
#if __cplusplus >= 201703L
		std::filesystem::path seqpath(seqfile);
		std::filesystem::path yamlpath(yamlIn);
#else
		std::experimental::filesystem::path seqpath(seqfile);
		std::experimental::filesystem::path yamlpath(yamlIn);
#endif
		std::string fullseq;
		if (seqpath.parent_path().empty()) { // get path from yaml file
			fullseq = std::string(yamlpath.parent_path().u8string());
			if (!fullseq.empty()) {
#ifdef WIN32
				fullseq += "\\";
#else
				fullseq += "/";
#endif
			}
			fullseq += seqfile;
		}
		else {
			fullseq = seqfile;
		}
		if (!fitFramework.SetExternalSequence(fullseq)) {
			std::cout << "Could not read pulseq file" << std::endl;
			return false;
		}
	}
	catch (...)
	{
		std::cout << "Could not read pulseq file" << std::endl;
		return false;
	}
	return true;
}


//! Write the results yaml file
/*!
   \param yamlOut file name of the yaml file
   \param fitResult results of the fit
   \return true if success
*/
bool YamlParser::WriteFitResult(std::string yamlOut, std::vector<FitParameter>* fitResult)
{
	try
	{
		YAML::Node node;
		std::ofstream fout(yamlOut);
		for (int i = 0; i < fitResult->size(); i++)
		{
			std::string nodestr = fitResult->at(i).name + " : " + std::to_string(fitResult->at(i).get());
			node = YAML::Load(nodestr);
			fout << node << std::endl;
		}
		fout.close();
	}
	catch (...)
	{
		std::cout << "Could not write results file" << std::endl;
		return false;
	}
	return true;
}

//! Read options yaml file
/*!
   \param yamlOptionsFile file name of the yaml fit options file
   \param opts ceres::Solver::Options by reference
   \return true if success
*/
bool YamlParser::ParseYamlCeresOptions(std::string yamlOptionsFile, ceres::Solver::Options & opts)
{
	YAML::Node config;
	// try to read the file
	try
	{
		config = YAML::LoadFile(yamlOptionsFile);
	}
	catch (...)
	{
		std::cout << "Could not load yaml file " << yamlOptionsFile << std::endl;
		return false;
	}

	try
	{
		auto fit_options = config["fit_options"];
		for (YAML::const_iterator it = fit_options.begin(); it != fit_options.end(); ++it) {
			std::string name = it->first.as<YAML::Node>().as<std::string>();
			// trust region type
			if (name.compare("trust_region_strategy_type") == 0) {
				std::string type = it->second.as<YAML::Node>().as<std::string>();
				if (type.compare("LEVENBERG_MARQUARDT") == 0) {
					opts.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
				}
				else if (type.compare("DOGLEG") == 0) {
					opts.trust_region_strategy_type = ceres::DOGLEG;
				}
				else {
					std::cout << type << " is not a valid trust_region_strategy_type, Using LEVENBERG_MARQUARDT" << std::endl;
				}
			}
			else if (name.compare("max_num_iterations") == 0) {
				opts.max_num_iterations = it->second.as<YAML::Node>().as<int>();
			}
			else if (name.compare("max_solver_time_in_seconds") == 0) {
				opts.max_solver_time_in_seconds = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("initial_trust_region_radius") == 0) {
				opts.initial_trust_region_radius = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("max_trust_region_radius") == 0) {
				opts.max_trust_region_radius = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("min_trust_region_radius") == 0) {
				opts.min_trust_region_radius = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("function_tolerance") == 0) {
				opts.function_tolerance = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("gradient_tolerance") == 0) {
				opts.gradient_tolerance = it->second.as<YAML::Node>().as<double>();
			}
			else if (name.compare("parameter_tolerance") == 0) {
				opts.parameter_tolerance = it->second.as<YAML::Node>().as<double>();
			}
			else {
				std::cout << name << " is not a valid ceres fit option or not implemented as variable yet. Ignring this option" << std::endl;
			}
		}
	}
	catch (...)
	{
		std::cout << "Error during read of ceres fit parameters!" << std::endl;
		return false;
	}
	return true;
}