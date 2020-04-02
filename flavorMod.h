#pragma once

#ifndef __FlavorMod__
#define __FlavorMod__

#include "fxobjects.h"
#include "pluginstructures.h"
#include "superlfo.h"

/**
\struct FlavorModParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the FlavorMod object.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
enum class BandSelect { SUM, LPF, LBPF, HBPF, HPF };
enum class DistType { Creamy, Crusty, Crunchy, Crispy };

struct FlavorModParameters
{
	FlavorModParameters() {}

	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	FlavorModParameters& operator=(const FlavorModParameters& params)	// need this override for collections to work
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

		enableLow = params.enableLow;
		enableLowMid = params.enableLowMid;
		enableHighMid = params.enableHighMid;
		enableHigh = params.enableHigh;
		
		lowSat = params.lowSat;
		lowMidSat = params.lowMidSat;
		highMidSat = params.highMidSat;
		highSat = params.highSat;

		lowAmpdB = params.lowAmpdB;
		lowMidAmpdB = params.lowMidAmpdB;
		highMidAmpdB = params.highMidAmpdB;
		highAmpdB = params.highAmpdB;
		
		dryAmpdB = params.dryAmpdB;
		distType = params.distType;
		distBypass = params.distBypass;
		gainCompOn = params.gainCompOn;
		distRate = params.distRate;
		distDepth = params.distDepth;
		enableLFO = params.enableLFO;

		// --- MUST be last
		return *this;
	}

	// --- individual parameters
	double splitFreq1 = 0.0;	///< init
	double splitFreq2 = 0.0;
	double splitFreq3 = 0.0;

	bool enableLow = true;
	bool enableLowMid = true;
	bool enableHighMid = true;
	bool enableHigh = true;
	
	double lowSat = 0.0;
	double lowMidSat = 0.0;
	double highMidSat = 0.0;
	double highSat = 0.0;

	double lowAmpdB = 0.0;
	double lowMidAmpdB = 0.0;
	double highMidAmpdB = 0.0;
	double highAmpdB = 0.0;

	double distRate = 0.0;
	double distDepth = 0.0;

	double dryAmpdB = 0.0;
	DistType distType = DistType::Crunchy;
	bool distBypass = false;
	bool gainCompOn = true;
	bool enableLFO = false;
};


/**
\class FlavorMod
\ingroup FX-Objects
\brief
The FlavorMod object implements ....

Audio I/O:
- Processes mono input to mono output.
- *** Optionally, process frame *** Modify this according to your object functionality

Control I/F:
- Use FlavorModParameters structure to get/set object params.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
class FlavorMod : public IAudioSignalProcessor
{
public:
	FlavorMod(void) {}	/* C-TOR */
	~FlavorMod(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- store the sample rate
		sampleRate = _sampleRate;

		// --- do any other per-audio-run inits here
		for (int i = 0; i < 6, i++;)
		{
			splitterFilters[i].reset(_sampleRate);
		}
		distLFO.reset(_sampleRate);

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

		double ynL;

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

		//dry is made up of all bands pre-distortion
		double dryL = 0.0;
		for (int i = 0; i < 4; i++)		
			dryL = dryL + bandSigL[i];

		//are we automating the distortion?
		SignalModulatorOutput lfoOutput = distLFO.renderModulatorOutput();
		if (parameters.enableLFO)
		{
			parameters.lowSat = doBipolarModulation(lfoOutput.normalOutput, 1.0, 11.0);
			parameters.lowMidSat = doBipolarModulation(lfoOutput.invertedOutput, 1.0, 11.0);
			parameters.highMidSat = doBipolarModulation(lfoOutput.quadPhaseOutput_neg, 1.0, 11.0);
			parameters.highSat = doBipolarModulation(lfoOutput.quadPhaseOutput_pos, 1.0, 11.0);
		}

		//now distort by band if not bypassed
		if (!parameters.distBypass)
		{
			if (parameters.lowSat > 1.0 && parameters.enableLow)
				bandSigL[0] = performDistortion(bandSigL[0], parameters.lowSat, parameters.distType, parameters.gainCompOn);
			if (parameters.lowMidSat > 1.0 && parameters.enableLowMid)
				bandSigL[1] = performDistortion(bandSigL[1], parameters.lowMidSat, parameters.distType, parameters.gainCompOn);
			if (parameters.highMidSat > 1.0 && parameters.enableHighMid)
				bandSigL[2] = performDistortion(bandSigL[2], parameters.highMidSat, parameters.distType, parameters.gainCompOn);
			if (parameters.highSat > 1.0 && parameters.enableHigh)
				bandSigL[3] = performDistortion(bandSigL[3], parameters.highSat, parameters.distType, parameters.gainCompOn);
		}
		ynL = bandSigL[0] * lowAmp
			+ bandSigL[1] * lowMidAmp
			+ bandSigL[2] * highMidAmp
			+ bandSigL[3] * highAmp
			+ dryL * dryAmp;
		

		// mono-mono
		if (inputChannels == 1 && outputChannels == 1)
		{
			outputFrame[0] = ynL;
			return true;
		}

		//mono-stereo
		if (inputChannels == 1 && outputChannels == 2)
		{
			outputFrame[0] = ynL;
			outputFrame[1] = ynL;
			return true;
		}

		//stereo-stereo
		if (inputChannels == 2 && outputChannels == 2)
		{
			double ynR;
			double bandSigR[4];
			bandSigR[0] = splitFilter_Right1.LFOut;
			bandSigR[1] = splitFilter_Right2.LFOut;
			bandSigR[2] = splitFilter_Right3.LFOut; 
			bandSigR[3] = splitFilter_Right3.HFOut;
			
			//dry is made up of all bands pre-distortion
			double dryR = 0.0;
			for (int i = 0; i < 4; i++)		
				dryR = dryR + bandSigR[i];
			//now distort by band if not bypassed
			if (!parameters.distBypass)
			{
				if (parameters.lowSat > 1.0 && parameters.enableLow)
					bandSigR[0] = performDistortion(bandSigR[0], parameters.lowSat, parameters.distType, parameters.gainCompOn);
				if (parameters.lowMidSat > 1.0 && parameters.enableLowMid)
					bandSigR[1] = performDistortion(bandSigR[1], parameters.lowMidSat, parameters.distType, parameters.gainCompOn);
				if (parameters.highMidSat > 1.0 && parameters.enableHighMid)
					bandSigR[2] = performDistortion(bandSigR[2], parameters.highMidSat, parameters.distType, parameters.gainCompOn);
				if (parameters.highSat > 1.0 && parameters.enableHigh)
					bandSigR[3] = performDistortion(bandSigR[3], parameters.highSat, parameters.distType, parameters.gainCompOn);
			}
			ynR = bandSigR[0] * lowAmp
				+ bandSigR[1] * lowMidAmp
				+ bandSigR[2] * highMidAmp
				+ bandSigR[3] * highAmp
				+ dryR * dryAmp;

  			outputFrame[0] = ynL;
			outputFrame[1] = ynR;
			return true;
		}
		return false;
	}


	/** get parameters: note use of custom structure for passing param data */
	/**
	\return FlavorModParameters custom data structure
	*/
	FlavorModParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param FlavorModParameters custom data structure
	*/
	void setParameters(const FlavorModParameters& params)
	{
		// --- copy them; note you may choose to ignore certain items
		//     and copy the variables one at a time, or you may test
		//     to see if cook-able variables have changed; if not, then
		//     do not re-cook them as it just wastes CPU
		parameters = params;

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

		if (parameters.lowAmpdB == -96.0 || !parameters.enableLow) lowAmp = 0.0;
		else lowAmp = pow(10.0, parameters.lowAmpdB / 20.0);

		if (parameters.lowMidAmpdB == -96.0 || !parameters.enableLowMid) lowMidAmp = 0.0;
		else lowMidAmp = pow(10.0, parameters.lowMidAmpdB / 20.0);

		if (parameters.highMidAmpdB == -96.0 || !parameters.enableHighMid) highMidAmp = 0.0;
		else highMidAmp = pow(10.0, parameters.highMidAmpdB / 20.0);

		if (parameters.highAmpdB == -96.0 || !parameters.enableHigh) highAmp = 0.0;
		else highAmp = pow(10.0, parameters.highAmpdB / 20.0);

		if (parameters.dryAmpdB == -96.0) dryAmp = 0.0;
		else dryAmp = pow(10.0, parameters.dryAmpdB / 20.0);

		//if (parameters.enableLFO)
		//{
			SuperLFOParameters lfoParams = distLFO.getParameters();
			//if (lfoParams.frequency_Hz != parameters.distRate || lfoParams.outputAmplitude != parameters.distRate)
			//{
				lfoParams.frequency_Hz = parameters.distRate;
				lfoParams.outputAmplitude = parameters.distDepth / 100.0;
				lfoParams.waveform = LFOWaveform::kTriangle;
				distLFO.setParameters(lfoParams);
			//}
		//}
			
	}

private:
	FlavorModParameters parameters; ///< object parameters
	LRFilterBank splitterFilters[6];
	SuperLFO distLFO;

	// --- local variables used by this object
	double sampleRate = 0.0;	///< sample rate

	double lowAmp = 0.0;
	double lowMidAmp = 0.0;
	double highMidAmp = 0.0;
	double highAmp = 0.0;

	double dryAmp = 0.0;

	double performDistortion(double x, double k, DistType type, bool comp)		//x is input sample, k ranges from 1 to 11, type selection
	{
		double y = x;
		double compGain;
		if (type == DistType::Creamy)	//mix together arraya and modified sigmoid
		{
			double y1 = 1.5*x*(1.0 - x*x/3.0);			//arraya
			double y2 = 2.0/(1.0 + exp(-6.0*x)) - 1.0;	//sigmoid where k=6
			double k1 = (k - 1.0)/10.0;							//k now ranges from 0 to 1
			y = (1.0 - k1) * y1 + k1 * y2;
			if (comp)
			{
				compGain = 1.0 - (k - 1.0) / 18.0;
				y = y * compGain;
			}
		}

		else if (type == DistType::Crusty)	//arctangent
		{
			y = atan(x * k) / atan(k);
			if (comp)
			{
				compGain = 1.0 - (k - 1.0) / 11.6;
				y = y * compGain;
			}
		}
		
		else if (type == DistType::Crunchy)	// //mix together sigmoid2 and bipolar square root 
		{
			double y1, y2;
			y1 = ((exp(x) - 1.0)*(exp(1.0) + 1.0)) / ((exp(x) + 1.0)*(exp(1.0) - 1.0));

			if (x >= 0.0)
				y2 = sqrt(x);
			else
				y2 = -1.0 * sqrt(-1.0 * x);
			double k1 = (k - 1.0) / 10.0;
			y = (1.0 - k1) * y1 + k1 * y2;
			if (comp)
			{
				compGain = 1.0 - (k - 1.0) / 13.0;
				y = y * compGain;
			}
		}

		else if (type == DistType::Crispy)
		{
			y = tanh(x * k*k) / tanh(k*k);	// overdriven hyperbolic tangent
			if (comp)
			{
				compGain = 1.0 - sqrt((k - 1.0) / 12.0);
				y = y * compGain;
			}
		}
		
		return y;
	}

};

#endif