//!  MatlabIO.h
/*!
Read/return variables to MATLAB

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

#include "yaml-cpp/yaml.h"
#include "SimulationParameters.h"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

//! Read the parameter yaml file
/*!
   \param yamlIn file name of the yaml file
   \param sp SimulationParameters object by reference that gets filled
   \return true if success
*/
bool ParseYamlInputStruct(std::string yamlIn, SimulationParameters &sp)
{
	YAML::Node config;
	// try to read the file
	try
	{
		config = YAML::LoadFile(yamlIn);
	}
	catch (...)
	{
		std::cout << "Could not read load yaml file" << std::endl;
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
			if(sring_ls.compare("Lorentzian") == 0) {
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
			sp.InitCESTPoolMemory(cest_pools.size());
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
	double scale = 1.0;
	try
	{
		auto scaleNode = config["scale"];
		if (scaleNode.IsDefined()) {
			scale = scaleNode.as<double>();
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
	M(2 * (ncp + 1)) = sp.GetWaterPool()->GetFraction();
	for (int i = 0; i < ncp; i++) {
		M(2 * (ncp + 1) + i + 1) = sp.GetCESTPool(i)->GetFraction();
	}
	if (sp.IsMTActive())
		M(vecSize - 1) = sp.GetMTPool()->GetFraction();
	// scale
	M *= scale;
	// init with same number as fit data
	sp.InitMagnetizationVectors(M, sp.GetFitData()->size());

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
#ifdef WIN32
			fullseq += "\\";
#else
			fullseq += "/";
#endif
			fullseq += seqfile;
		}
		else {
			fullseq = seqfile;
		}
		sp.SetExternalSequence(fullseq);
	}
	catch (...)
	{
		std::cout << "Coul not read pulseq file" << std::endl;
		return false;
	}

	// read ther stuff sequence
	try
	{
		auto b0 = config["b0"];
		sp.InitScanner(b0.as<double>());
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
		std::cout << "Could not read pulseq file" << std::endl;
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


//! Write the results yaml file
/*!
   \param yamlOut file name of the yaml file
   \param fitResult results of the fit
   \return true if success
*/
bool WriteFitResult(std::string yamlOut, std::vector<FitParameter>* fitResult)
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