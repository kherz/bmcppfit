//!  YamlParser.h
/*!
Read/write to yaml files

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
#pragma once

#include "yaml-cpp/yaml.h"
#include "BMSimFit.h"
#include "BMCSim.h"
#include "ceres/ceres.h"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif


class YamlParser
{
public:
	YamlParser();
	bool ParseYamlInputStruct(std::string yamlIn, SimFitParameters &sp);
	bool WriteFitResult(std::string yamlOut, std::vector<FitParameter>* fitResult);
	bool ParseYamlCeresOptions(std::string yamlOptionsFile, ceres::Solver::Options & opts);
	std::string GetSeqFilename();
private:
	std::string seqFilename;
};
