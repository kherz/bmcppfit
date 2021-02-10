//!  BMSim.cpp
/*!
Bloch-McConnell Z-Spectrum simulation

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

#define MAX_SAMPLES 200

#include "YamlIO.h"
#include "BMSim_T.h"

// decode the unique pulses from the seq file
bool DecodeRFShapes(ExternalSequence &seq, std::map<std::pair<int, int>, std::vector<PulseSample>> &uniquePulses)
{
	std::vector<std::pair<int, int>> uniquePairs;
	for (unsigned int nSample = 0; nSample < seq.GetNumberOfBlocks(); nSample++)
	{
		SeqBlock* seqBlock = seq.GetBlock(nSample);
		if (seqBlock->isRF())
		{
			RFEvent rf = seqBlock->GetRFEvent();
			// ake unique mag phase pair 
			auto p = std::make_pair(rf.magShape, rf.phaseShape);
			if (!(std::find(uniquePairs.begin(), uniquePairs.end(), p) != uniquePairs.end())){
				// register pulse
				std::vector<PulseSample> uniqueSamples;
				// get rf and check its length
				seq.decodeBlock(seqBlock);
				unsigned int rfLength = seqBlock->GetRFLength();
				// check arrays of uncompresed shape
				std::vector<float> amplitudeArray(seqBlock->GetRFAmplitudePtr(), seqBlock->GetRFAmplitudePtr() + rfLength);
				std::vector<float> phaseArray(seqBlock->GetRFPhasePtr(), seqBlock->GetRFPhasePtr() + rfLength);
				// rfDeadTime is usually zeros at the end of the pulse, we search for them here
				int nEnd;
				int delayAfterPulse = 0;
				for (nEnd = rfLength; nEnd > 0; --nEnd) {
					if (fabs(amplitudeArray[nEnd - 1]) > 1e-6)// because of the round-up errors in the ascii and derivative/integral reconstructuion
						break;
				}
				delayAfterPulse = rfLength - nEnd;
				rfLength = nEnd;

				amplitudeArray.erase(amplitudeArray.end() - delayAfterPulse, amplitudeArray.end());
				phaseArray.erase(phaseArray.end() - delayAfterPulse, phaseArray.end());
				// search for unique samples in amplitude and phase
				std::vector<float> amplitudeArrayUnique(rfLength);
				std::vector<float>::iterator it_amplitude = std::unique_copy(amplitudeArray.begin(), amplitudeArray.end(), amplitudeArrayUnique.begin());
				amplitudeArrayUnique.resize(std::distance(amplitudeArrayUnique.begin(), it_amplitude));
				std::vector<float> phaseArrayUnique(rfLength);
				std::vector<float>::iterator it_phase = std::unique_copy(phaseArray.begin(), phaseArray.end(), phaseArrayUnique.begin());
				phaseArrayUnique.resize(std::distance(phaseArrayUnique.begin(), it_phase));
				//
				float rfAmplitude = 0.0;
				float rfPhase = 0.0;
				float rfFrequency = seqBlock->GetRFEvent().freqOffset;
				float timestep;
				// need to resample pulse
				unsigned int max_samples = std::max(amplitudeArrayUnique.size(), phaseArrayUnique.size());
				if (max_samples > MAX_SAMPLES) {
					int sampleFactor = ceil(float(rfLength) / MAX_SAMPLES);
					float pulseSamples = rfLength / sampleFactor;
					timestep = float(sampleFactor) * 1e-6;
					// resmaple the original pulse with max ssamples and run the simulation
					uniqueSamples.resize(pulseSamples);
					for (int i = 0; i < pulseSamples; i++) {
						uniqueSamples[i].magnitude = seqBlock->GetRFAmplitudePtr()[i*sampleFactor] * seqBlock->GetRFEvent().amplitude;
						uniqueSamples[i].phase = seqBlock->GetRFPhasePtr()[i*sampleFactor] + seqBlock->GetRFEvent().phaseOffset;
						uniqueSamples[i].timestep = timestep;
					}
				}
				else {
					std::vector<unsigned int>samplePositions(max_samples + 1);
					unsigned int sample_idx = 0;
					if (amplitudeArrayUnique.size() >= phaseArrayUnique.size()) {
						std::vector<float>::iterator it = amplitudeArray.begin();
						for (it_amplitude = amplitudeArrayUnique.begin(); it_amplitude != amplitudeArrayUnique.end(); ++it_amplitude) {
							it = std::find(it, amplitudeArray.end(), *it_amplitude);
							samplePositions[sample_idx++] = std::distance(amplitudeArray.begin(), it);
						}
					}
					else {
						std::vector<float>::iterator it = phaseArray.begin();
						for (it_phase = phaseArrayUnique.begin(); it_phase != phaseArrayUnique.end(); ++it_phase) {
							it = std::find(it, phaseArray.end(), *it_phase);
							samplePositions[sample_idx++] = std::distance(phaseArray.begin(), it);
						}
					}
					uniqueSamples.resize(max_samples);
					samplePositions[max_samples] = rfLength;
					// now we have the duration of the single samples -> simulate it
					for (int i = 0; i < max_samples; i++) {
						uniqueSamples[i].magnitude = seqBlock->GetRFAmplitudePtr()[samplePositions[i]] * seqBlock->GetRFEvent().amplitude;
						uniqueSamples[i].phase = seqBlock->GetRFPhasePtr()[samplePositions[i]] + seqBlock->GetRFEvent().phaseOffset;
						uniqueSamples[i].timestep = (samplePositions[i + 1] - samplePositions[i])*1e-6;
					}
				}
				uniquePairs.push_back(p);
				uniquePulses.insert(std::make_pair(p, uniqueSamples));
			}
		}
		delete seqBlock;
	}
	return true;
}

//! Enry point for MATLAB mex function
/*!
    \param nlhs number of output arguments
	\param plhs Array of pointers to the mxArray output arguments
	\param nrhs number of input arguments
	\param prhs Array of pointers to the mxArray input arguments
*/
bool RunBMSimFit(std::string yamlIn, std::string yamlOut)
{
	//init the simulation interface and read the input
	SimulationParameters sp;
	if (!ParseYamlInputStruct(yamlIn, sp))
		return EXIT_FAILURE;

	// at first check pulseq for unique pulses, so that we have to decode them only once
	std::map<std::pair<int, int>, std::vector<PulseSample>>  uniquePulses;
	DecodeRFShapes(sp.seq, uniquePulses);

	//// disp info about input
	//if (sp.verboseMode) {
	//	mexPrintf("Read parameters succesfully! \n");
	//	mexPrintf("Found %i CEST pool(s) and %i MT Pool(s) \n", sp.numberOfCESTPools, sp.simulateMTPool ? 1 : 0);
	//}


	///* For a small number of pools the matrix size can be set at compile time. This ensures allocation on the stack and therefore a faster simulation. 
	//   This speed advantade vanishes for more pools and can even result in a stack overflow for very large matrices
	//   In this case more than 3 pools are simulated with dynamic matrices, but this could be expanded eventually
	//*/
	//switch (sp.numberOfCESTPools)
	//{
	//case 0:
	//	sp.simulateMTPool ? BMSim_T<4>(sp) : BMSim_T<3>(sp); // only water
	//	break;
	//case 1:
	//	sp.simulateMTPool ? BMSim_T<7>(sp) : BMSim_T<6>(sp); // one cest pool
	//	break;
	//case 2:
	//	sp.simulateMTPool ? BMSim_T<10>(sp) : BMSim_T<9>(sp); // two cest pools
	//	break;
	//case 3:
	//	sp.simulateMTPool ? BMSim_T<13>(sp) : BMSim_T<12>(sp); // three cest pools
	//	break;
	//default:
	//	BMSim_T<Dynamic>(sp); // > three pools
	//	break;
	//}

	//ReturnResultToMATLAB(plhs, sp.Mvec); // return results after simulation

}
