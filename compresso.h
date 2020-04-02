#pragma once

#ifndef __Compresso__
#define __Compresso__

#include "fxobjects.h"
#include "pluginstructures.h"
#include "dynamicsprocessornew.h"

/**
\struct CompressoParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the Compresso object.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/

struct CompressoParameters
{
	CompressoParameters() {}

	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	CompressoParameters& operator=(const CompressoParameters& params)	// need this override for collections to work
	{
		// --- it is possible to try to make the object equal to itself
		//     e.g. thisObject = thisObject; so this code catches that
		//     trivial case and just returns this object
		if (this == &params)
			return *this;

		// --- copy from params (argument) INTO our variables
		splitFreq1 = params.splitFreq1;
		splitFreq2 = params.splitFreq2;
		splitFreq3 = params.splitFreq3;

		soloLow = params.soloLow;
		soloLowMid = params.soloLowMid;
		soloHighMid = params.soloHighMid;
		soloHigh = params.soloHigh;

		lowThreshold = params.lowThreshold;
		lowMidThreshold = params.lowMidThreshold;
		highMidThreshold = params.highMidThreshold;
		highThreshold = params.highThreshold;

		lowGain_dB = params.lowGain_dB;
		lowMidGain_dB = params.lowMidGain_dB;
		highMidGain_dB = params.highMidGain_dB;
		highGain_dB = params.highGain_dB;

		lowGRMeter = params.lowGRMeter;
		lowMidGRMeter = params.lowMidGRMeter;
		highMidGRMeter = params.highMidGRMeter;
		highGRMeter = params.highGRMeter;

		lowMeter = params.lowMeter;
		lowMidMeter = params.lowMidMeter;
		highMidMeter = params.highMidMeter;
		highMeter = params.highMeter;

		attack = params.attack;
		release = params.release;

		lowType = params.lowType;
		lowMidType = params.lowMidType;
		highMidType = params.highMidType;
		highType = params.highType;

		auditionDiff = params.auditionDiff;

		// --- MUST be last
		return *this;
	}

	// --- individual parameters
	double splitFreq1 = 0.0;	///< init
	double splitFreq2 = 0.0;
	double splitFreq3 = 0.0;

	bool soloLow = false;
	bool soloLowMid = false;
	bool soloHighMid = false;
	bool soloHigh = false;

	double lowThreshold = 0.0;
	double lowMidThreshold = 0.0;
	double highMidThreshold = 0.0;
	double highThreshold = 0.0;

	double lowGain_dB = 0.0;
	double lowMidGain_dB = 0.0;
	double highMidGain_dB = 0.0;
	double highGain_dB = 0.0;

	float lowGRMeter = 0.f;
	float lowMidGRMeter = 0.f;
	float highMidGRMeter = 0.f;
	float highGRMeter = 0.f;

	float lowMeter = 0.f;
	float lowMidMeter = 0.f;
	float highMidMeter = 0.f;
	float highMeter = 0.f;

	int attack = 0;
	int release = 0;

	int lowType = 0;
	int lowMidType = 0;
	int highMidType = 0;
	int highType = 0;

	int auditionDiff = 0;

};


/**
\class Compresso
\ingroup FX-Objects
\brief
The Compresso object implements ....

Audio I/O:
- Processes mono input to mono output.
- *** Optionally, process frame *** Modify this according to your object functionality

Control I/F:
- Use CompressoParameters structure to get/set object params.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
class Compresso : public IAudioSignalProcessor
{
public:
	Compresso(void) {}	/* C-TOR */
	~Compresso(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- store the sample rate
		sampleRate = _sampleRate;

		// --- do any other per-audio-run inits here
		for (int i = 0; i < 6, i++;)
			splitterFilters[i].reset(_sampleRate);
		for (int i = 0; i < 4; i++)
			bandDynamics[i].reset(_sampleRate);

		DynamicsProcessorNewParameters params[4];
		for (int i = 0; i < 4; i++)
		{
			params[i] = bandDynamics[i].getParameters();
			params[i].enableSidechain = true;
			params[i].attackTime_mSec = 5.0;
			params[i].releaseTime_mSec = 70.0;
			params[i].calculation = dynamicsProcessorType::kCompressor;
			params[i].softKnee = true;
			params[i].kneeWidth_dB = 5.0;
			params[i].ratio = 3.5;
			bandDynamics[i].setParameters(params[i]);
		}

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
	virtual bool canProcessAudioFrame() { return true; }

	/** process audio frame: implement this function if you answer "true" to above query */
	virtual bool processAudioFrame(const float* inputFrame,
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		double xnL = inputFrame[0];
		double xnR = inputFrame[1];

		double ynL = 0.0;

		FilterBankOutput splitFilter_Left1 = splitterFilters[0].processFilterBank(xnL);
		FilterBankOutput splitFilter_Right1 = splitterFilters[1].processFilterBank(xnR);

		FilterBankOutput splitFilter_Left2 = splitterFilters[2].processFilterBank(splitFilter_Left1.HFOut);
		FilterBankOutput splitFilter_Right2 = splitterFilters[3].processFilterBank(splitFilter_Right1.HFOut);

		FilterBankOutput splitFilter_Left3 = splitterFilters[4].processFilterBank(splitFilter_Left2.HFOut);
		FilterBankOutput splitFilter_Right3 = splitterFilters[5].processFilterBank(splitFilter_Right2.HFOut);

		double bandSigL[4];
		bandSigL[0] = splitFilter_Left1.LFOut;
		bandSigL[1] = splitFilter_Left2.LFOut;
		bandSigL[2] = splitFilter_Left3.LFOut;
		bandSigL[3] = splitFilter_Left3.HFOut;
		
		double ynR = 0.0;
		double bandSigR[4];
		double diffSigL[4];
		double diffSigR[4];

		double sidechain[4];
		
		if (inputChannels == 1)
		{
			for (int i = 0; i < 4; i++)
				sidechain[i] = bandSigL[i];
		}
		else
		{
			bandSigR[0] = splitFilter_Right1.LFOut;
			bandSigR[1] = splitFilter_Right2.LFOut;
			bandSigR[2] = splitFilter_Right3.LFOut;
			bandSigR[3] = splitFilter_Right3.HFOut;

			for(int i = 0; i < 4; i++)
				sidechain[i] = 0.5 * bandSigL[i] + 0.5 * bandSigR[i];

		}

		double gain[4];
		for (int i = 0; i < 4; i++)
		{
			bandDynamics[i].processAuxInputAudioSample(sidechain[i]);
			gain[i] = bandDynamics[i].processAudioSample(sidechain[i]);

			if (parameters.auditionDiff) 
			{
				bandSigL[i] = bandSigL[i] * (1.0-bandDynamics[i].getParameters().gainReduction);
				bandSigR[i] = bandSigR[i] * (1.0-bandDynamics[i].getParameters().gainReduction);
			}

			else
			{
				bandSigL[i] = bandSigL[i] * gain[i] * gainCooked[i];
				bandSigR[i] = bandSigR[i] * gain[i] * gainCooked[i];
			}
		}

		if (soloAny)
		{
			if (parameters.soloLow)
				ynL += bandSigL[0];
			if (parameters.soloLowMid)
				ynL += bandSigL[1];
			if (parameters.soloHighMid)
				ynL += bandSigL[2];
			if (parameters.soloHigh)
				ynL += bandSigL[3];
		}

		else ynL = bandSigL[0]
			+ bandSigL[1]
			+ bandSigL[2]
			+ bandSigL[3];


		// mono-mono
		if (inputChannels == 1 && outputChannels == 1)
		{
			parameters.lowMeter = bandSigL[1];
			parameters.lowMidMeter = bandSigL[2];
			parameters.highMidMeter = bandSigL[3];
			parameters.highMeter = bandSigL[4];

			outputFrame[0] = ynL;
			return true;
		}

		//mono-stereo
		if (inputChannels == 1 && outputChannels == 2)
		{
			parameters.lowMeter = bandSigL[1];
			parameters.lowMidMeter = bandSigL[2];
			parameters.highMidMeter = bandSigL[3];
			parameters.highMeter = bandSigL[4];
			
			outputFrame[0] = ynL;
			outputFrame[1] = ynL;
			return true;
		}

		//stereo-stereo
		if (inputChannels == 2 && outputChannels == 2)
		{
			parameters.lowMeter = 0.5 * bandSigL[1] + 0.5 * bandSigR[1];
			parameters.lowMidMeter = 0.5 * bandSigL[2] + 0.5 * bandSigR[2];
			parameters.highMidMeter = 0.5 * bandSigL[3] + 0.5 * bandSigR[3];
			parameters.highMeter = 0.5 * bandSigL[4] + 0.5 * bandSigR[4];
			
			if (soloAny)
			{
				if (parameters.soloLow)
					ynR += bandSigR[0];
				if (parameters.soloLowMid)
					ynR += bandSigR[1];
				if (parameters.soloHighMid)
					ynR += bandSigR[2];
				if (parameters.soloHigh)
					ynR += bandSigR[3];
			}

			else ynR = bandSigR[0]
				+ bandSigR[1]
				+ bandSigR[2]
				+ bandSigR[3];

			outputFrame[0] = ynL;
			outputFrame[1] = ynR;
			return true;
		}
		return false;
	}


	/** get parameters: note use of custom structure for passing param data */
	/**
	\return CompressoParameters custom data structure
	*/
	CompressoParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param CompressoParameters custom data structure
	*/
	void setParameters(const CompressoParameters& params)
	{
		// --- copy them; note you may choose to ignore certain items
		//     and copy the variables one at a time, or you may test
		//     to see if cook-able variables have changed; if not, then
		//     do not re-cook them as it just wastes CPU
		
		if (parameters.lowGain_dB != params.lowGain_dB)
			gainCooked[0] = pow(10.0, parameters.lowGain_dB / 20.0);
		if (parameters.lowMidGain_dB != params.lowMidGain_dB)
			gainCooked[1] = pow(10.0, parameters.lowMidGain_dB / 20.0);
		if (parameters.highMidGain_dB != params.highMidGain_dB)
			gainCooked[2] = pow(10.0, parameters.highMidGain_dB / 20.0);
		if (parameters.highGain_dB != params.highGain_dB)
			gainCooked[3] = pow(10.0, parameters.highGain_dB / 20.0);

		parameters = params;

		DynamicsProcessorNewParameters bandParams[4];
		
		for (int i = 0; i < 4; i++)
			bandParams[i] = bandDynamics[i].getParameters();
		
		bandParams[0].threshold_dB = parameters.lowThreshold;
		bandParams[1].threshold_dB = parameters.lowMidThreshold;
		bandParams[2].threshold_dB = parameters.highMidThreshold;
		bandParams[3].threshold_dB = parameters.highThreshold;

		bandParams[0].attackTime_mSec = (parameters.attack + 1) * 22.0;
		bandParams[1].attackTime_mSec = (parameters.attack + 1) * 12.0;
		bandParams[2].attackTime_mSec = (parameters.attack + 1) * 7.0;
		bandParams[3].attackTime_mSec = (parameters.attack + 1) * 3.0;

		bandParams[0].ratio = 1.7 + parameters.lowType * 2.1;
		bandParams[1].ratio = 1.7 + parameters.lowMidType * 2.1;
		bandParams[2].ratio = 1.7 + parameters.highMidType * 2.1;
		bandParams[3].ratio = 1.7 + parameters.highType * 2.1;

		bandParams[0].kneeWidth_dB = 13 - parameters.lowType * 6;
		bandParams[1].kneeWidth_dB = 13 - parameters.lowMidType * 6;
		bandParams[2].kneeWidth_dB = 13 - parameters.highMidType * 6;
		bandParams[3].kneeWidth_dB = 13 - parameters.highType * 6;

		bandParams[0].releaseTime_mSec = (parameters.release + 1) * 80.0 - 40.0;
		bandParams[1].releaseTime_mSec = (parameters.release + 1) * 60.0 - 30.0;
		bandParams[2].releaseTime_mSec = (parameters.release + 1) * 40.0 - 20.0;
		bandParams[3].releaseTime_mSec = (parameters.release + 1) * 20.0 - 10.0;

		parameters.lowGRMeter = bandParams[0].gainReduction;
		parameters.lowMidGRMeter = bandParams[1].gainReduction;
		parameters.highMidGRMeter = bandParams[2].gainReduction;
		parameters.highGRMeter = bandParams[3].gainReduction;

		for (int i = 0; i < 4; i++)
			bandDynamics[i].setParameters(bandParams[i]);

		bool changedSplitFreq =
			(splitterFilters[0].getParameters().splitFrequency != parameters.splitFreq1) ||
			(splitterFilters[2].getParameters().splitFrequency != parameters.splitFreq2) ||
			(splitterFilters[4].getParameters().splitFrequency != parameters.splitFreq3);

		if (changedSplitFreq)
		{
			LRFilterBankParameters filterParams[6];

			for (int i = 0; i < 6; i++)
				filterParams[i] = splitterFilters[i].getParameters();

			filterParams[0].splitFrequency = parameters.splitFreq1;
			filterParams[1].splitFrequency = parameters.splitFreq1;
			filterParams[2].splitFrequency = parameters.splitFreq2;
			filterParams[3].splitFrequency = parameters.splitFreq2;
			filterParams[4].splitFrequency = parameters.splitFreq3;
			filterParams[5].splitFrequency = parameters.splitFreq3;

			for (int i = 0; i < 6; i++)
				splitterFilters[i].setParameters(filterParams[i]);

		}

		soloAny = (parameters.soloLow || parameters.soloLowMid || parameters.soloHighMid || parameters.soloHigh);

	}

private:
	CompressoParameters parameters; ///< object parameters
	DynamicsProcessorNew bandDynamics[4];
	LRFilterBank splitterFilters[6];

	// --- local variables used by this object
	double sampleRate = 0.0;	///< sample rate
	bool soloAny = false;
	double gainCooked[4] = {1.0, 1.0, 1.0, 1.0};
};

#endif#pragma once
