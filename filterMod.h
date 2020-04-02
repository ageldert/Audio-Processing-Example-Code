#pragma once

#ifndef __FilterMod__
#define __FilterMod__

#include "fxobjects.h"
#include "pluginstructures.h"
#include "envelopeFollowerNew.h"

/**
\struct FilterModParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the FilterMod object.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/

enum class FilterType { LPF, BPF, HPF };
enum class DetectMode { Envelope, LFO };

struct FilterModParameters
{
	FilterModParameters() {}

	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	FilterModParameters& operator=(const FilterModParameters& params)	// need this override for collections to work
	{
		// --- it is possible to try to make the object equal to itself
		//     e.g. thisObject = thisObject; so this code catches that
		//     trivial case and just returns this object
		if (this == &params)
			return *this;

		// --- copy from params (argument) INTO our variables
		filterFc = params.filterFc;
		filterQ = params.filterQ;
		threshold_DB = params.threshold_DB;
		sensitivity = params.sensitivity;
		attackTime_ms = params.attackTime_ms;
		releaseTime_ms = params.releaseTime_ms;
		gainBPF_dB = params.gainBPF_dB;
		dryWet = params.dryWet;
		filterType = params.filterType;
		detectMode = params.detectMode;
		lfoRate = params.lfoRate;
		modDown = params.modDown;
		fcMeter = params.fcMeter;
		coupledQ = params.coupledQ;
		engagePF = params.engagePF;
		parFactor = params.parFactor;
		return *this;
	}

	// --- individual parameters
	double filterFc = 0.0;
	double filterQ = 0.0;
	double threshold_DB = 0.0;
	double sensitivity = 0.0;
	double attackTime_ms = 0.0;
	double releaseTime_ms = 0.0;
	double gainBPF_dB = 0.0;
	double dryWet = 0.0;
	double lfoRate = 0.0;
	bool modDown = false;
	float fcMeter = 0.f;
	bool coupledQ = false;
	bool engagePF = false;
	double parFactor = 1.0;
	FilterType filterType = FilterType::LPF;
	DetectMode detectMode = DetectMode::Envelope;
};


/**
\class FilterMod
\ingroup FX-Objects
\brief
The FilterMod object implements ....

Audio I/O:
- Processes mono input to mono output.
- *** Optionally, process frame *** Modify this according to your object functionality

Control I/F:
- Use FilterModParameters structure to get/set object params.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
class FilterMod : public IAudioSignalProcessor
{
public:
	FilterMod(void) {}	/* C-TOR */
	~FilterMod(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- store the sample rate
		sampleRate = _sampleRate;

		for (int i = 0; i < 4; i++)
			modFilters[i].reset(_sampleRate);

		return true;
	}

	/** process MONO input */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		return xn;
	}

	/** query to see if this object can process frames */
	virtual bool canProcessAudioFrame() { return true; } // <-- change this!

	/** process audio frame: implement this function if you answer "true" to above query */
	virtual bool processAudioFrame(const float* inputFrame,	/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		double xnL = inputFrame[0];
		double xnR = inputFrame[1];

		double ynL;

		if (parameters.engagePF)
			ynL = modFilters[0].processAudioSample(xnL) + modFilters[2].processAudioSample(xnL);
		else
			ynL = modFilters[0].processAudioSample(xnL);

		double gain = 1.0;
		if (parameters.gainBPF_dB != 0.0 && parameters.filterType == FilterType::BPF)
			gain = pow(10.0, parameters.gainBPF_dB / 20.0);
		double dryMix = parameters.dryWet / 100.0;
		double wetMix = 1.0 - (parameters.dryWet / 100.0);

		// mono-mono
		if (inputChannels == 1 && outputChannels == 1)
		{
			outputFrame[0] = ynL*gain*wetMix + xnL*dryMix;
			parameters.fcMeter = modFilters[0].getParameters().fcMeter;
			return true;
		}

		// mono-stereo
		if (inputChannels == 1 && outputChannels == 2)
		{
			outputFrame[0] = ynL*gain*wetMix + xnL*dryMix;
			outputFrame[1] = ynL*gain*wetMix + xnL*dryMix;
			parameters.fcMeter = modFilters[0].getParameters().fcMeter;
			return true;
		}

		// stereo-stereo
		if (inputChannels == 2 && outputChannels == 2)
		{
			double ynR;

			if (parameters.engagePF)
				ynR = modFilters[1].processAudioSample(xnR) + modFilters[3].processAudioSample(xnR);
			else
				ynR = modFilters[1].processAudioSample(xnR);
			
			outputFrame[0] = ynL*gain*wetMix + xnL*dryMix;
			outputFrame[1] = ynR*gain*wetMix + xnR*dryMix;
			parameters.fcMeter = modFilters[0].getParameters().fcMeter;
			return true;
		}

		return false; // NOT handled
	}


	/** get parameters: note use of custom structure for passing param data */
	/**
	\return FilterModParameters custom data structure
	*/
	FilterModParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param FilterModParameters custom data structure
	*/
	void setParameters(const FilterModParameters& params)
	{
		// --- copy them; note you may choose to ignore certain items
		//     and copy the variables one at a time, or you may test
		//     to see if cook-able variables have changed; if not, then
		//     do not re-cook them as it just wastes CPU
		parameters = params;

		EnvelopeFollowerNewParameters filterParams = modFilters[0].getParameters();

		if (filterParams.fc != parameters.filterFc || filterParams.Q != parameters.filterQ)
		{
			filterParams.fc = parameters.filterFc;
			if (parameters.coupledQ)
			{
				double percent = (parameters.filterFc - 20.0) / 22000.0;
				filterParams.Q = 0.707 + percent * 14.3;	//when f = 20, Q = 0.707 - when f = 20000, Q = 12
			}
			else
				filterParams.Q = parameters.filterQ;
		}

		filterParams.threshold_dB = parameters.threshold_DB;
		filterParams.sensitivity = parameters.sensitivity;
		filterParams.attackTime_mSec = parameters.attackTime_ms;
		filterParams.releaseTime_mSec = parameters.releaseTime_ms;
		filterParams.lfoRate = parameters.lfoRate;
		filterParams.modDown = parameters.modDown;
		filterParams.fcMeter = parameters.fcMeter;
		if (params.detectMode == DetectMode::LFO)
			filterParams.useLFO = true;
		else
			filterParams.useLFO = false;

		int index;
		if (parameters.filterType == FilterType::LPF) index = 11;
		else if (parameters.filterType == FilterType::BPF) index = 5;
		else if (parameters.filterType == FilterType::HPF) index = 4;
		filterParams.filterAlgorithm = filterAlgorithm(index);

		modFilters[0].setParameters(filterParams);
		modFilters[1].setParameters(filterParams);

		if (parameters.engagePF)
		{
			EnvelopeFollowerNewParameters parallelFilterParams = filterParams;
			parallelFilterParams.fc = parameters.filterFc * parameters.parFactor;
			boundValue(parallelFilterParams.fc, 20.0, kMaxFilterFrequency);
			if (parameters.coupledQ)
			{
				parallelFilterParams.Q = 0.707 + ((parallelFilterParams.fc - 20.0) / 22000.0) * 14.3;
			}
			modFilters[2].setParameters(parallelFilterParams);
			modFilters[3].setParameters(parallelFilterParams);
		}
	}

private:
	FilterModParameters parameters; ///< object parameters
	EnvelopeFollowerNew modFilters[4];

	// --- local variables used by this object
	double sampleRate = 0.0;	///< sample rate
};

#endif