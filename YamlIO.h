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


bool ParseYamlInputStruct(std::string yamlIn, SimulationParameters &sp)
{
	YAML::Node config = YAML::LoadFile(yamlIn);
	// check for offsets
	try
	{
		auto offsets_ppm = config["offsets_ppm"];
		sp.numPPMSamples = offsets_ppm.size();
		sp.PPMVector.resize(sp.numPPMSamples);
		for (int sample = 0; sample < sp.numPPMSamples; sample++)
			sp.PPMVector[sample] = (offsets_ppm[sample].as<double>());
	}
	catch (...)
	{
		std::cout << "Coul not read offset_ppm" << std::endl;
		return false;
	}

	// check for data to fit
	try
	{
		auto fit_data = config["fit_data"];

		if (sp.numPPMSamples != fit_data.size()) {
			std::cout << "offsets and fit data have different no of entries" << std::endl;
			return false;
		}
		sp.FitData.resize(sp.numPPMSamples);
		for (int sample = 0; sample < sp.numPPMSamples; sample++)
			sp.FitData[sample] = (fit_data[sample].as<double>());
	}
	catch (...)
	{
		std::cout << "Could not readfit data" << std::endl;
		return false;
	}

	// water pool
	try
	{
		auto water_pool = config["water_pool"];
		sp.waterPool.f = water_pool["f"].as<double>();
		sp.waterPool.R1 = 1.0 / water_pool["t1"].as<double>();
		sp.waterPool.R2 = 1.0 / water_pool["t2"].as<double>();
		sp.waterPool.dw = 0; // remove later
	}
	catch (...)
	{
		std::cout << "Coul not read water_pool" << std::endl;
		return false;
	}

	// mt pool
	try
	{
		auto mt_pool = config["mt_pool"];
		if (mt_pool.IsDefined()) {
			sp.mtPool.f = mt_pool["f"].as<double>();
			sp.mtPool.R1 = 1.0 / mt_pool["t1"].as<double>();
			sp.mtPool.R2 = 1.0 / mt_pool["t2"].as<double>();
			sp.mtPool.k = mt_pool["k"].as<double>();
			sp.mtPool.f = mt_pool["f"].as<double>();
			std::string sring_ls = mt_pool["lineshape"].as<std::string>();
			if(sring_ls.compare("Lorentzian") == 0) {
				sp.mtPool.lineshape = Lorentzian;
			}
		    else if (sring_ls.compare("SuperLorentzian") == 0) {
			   sp.mtPool.lineshape = SuperLorentzian;
			}
			else if (sring_ls.compare("None") == 0) {
				sp.mtPool.lineshape = None;
			}
			else {
				std::cout << "No valid MT Lineshape! Use None, Lorentzian or SuperLorentzian" << std::endl;
				return false;
			}
		}
		else {
			sp.simulateMTPool = false;
		}
		
	}
	catch (...)
	{
		std::cout << "Coul not read mt_pool" << std::endl;
		return false;
	}

	// cest pools
	try
	{
		auto cest_pools = config["cest_pool"];
		if (cest_pools.IsDefined()) {
			sp.numberOfCESTPools = cest_pools.size();
			sp.cestPool = new CESTPool[sp.numberOfCESTPools];
			sp.cestMemAllocated = true;
			int i = 0;
			for (YAML::const_iterator it = cest_pools.begin(); it != cest_pools.end(); ++it) {
				auto c_pool = it->second.as<YAML::Node>();
				sp.cestPool[i].R1 = 1.0 / c_pool["t1"].as<double>();
				sp.cestPool[i].R2 = 1.0 / c_pool["t2"].as<double>();
				sp.cestPool[i].f = c_pool["f"].as<double>();
				sp.cestPool[i].k = c_pool["k"].as<double>();
				sp.cestPool[i++].dw = c_pool["dw"].as<double>();
			}
		}
		else {
			sp.numberOfCESTPools = 0;
		}

	}
	catch (...)
	{
		std::cout << "Coul not read cest_pool" << std::endl;
		return false;
	}
	
	// read pulseq sequence
	try
	{
		auto seq_fn = config["pulseq_file"];
		sp.seq.load(seq_fn.as<std::string>());
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
		sp.scanner.B0 = b0.as<double>();
	}
	catch (...)
	{
		std::cout << "Could not read pulseq file" << std::endl;
		return false;
	}


#if defined (_OPENMP)
	//set threads for openmp
	// if there is no specification, the simulation only runs on one thread
	int maxThreads = omp_get_max_threads();
	try {
		auto threads = config["threads"];
		int inputThreads = b0.as<int>();
		sp.numberOfThreads = inputThreads <= maxThreads ? inputThreads : maxThreads;
	}
	catch (...) {
			sp.numberOfThreads = 1;
		}
	}
#else
	sp.numberOfThreads = 1;
#endif

	// enough for now we come to the rest later
	sp.verboseMode = false;
	
	return true;

}