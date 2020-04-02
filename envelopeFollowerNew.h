#pragma once
/**
\struct EnvelopeFollowerNewParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the EnvelopeFollowerNew object.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/

#include "superlfo.h"

struct EnvelopeFollowerNewParameters
{
	EnvelopeFollowerNewParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	EnvelopeFollowerNewParameters& operator=(const EnvelopeFollowerNewParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		fc = params.fc;
		Q = params.Q;
		attackTime_mSec = params.attackTime_mSec;
		releaseTime_mSec = params.releaseTime_mSec;
		threshold_dB = params.threshold_dB;
		sensitivity = params.sensitivity;
		filterAlgorithm = params.filterAlgorithm;
		useLFO = params.useLFO;
		lfoRate = params.lfoRate;
		modDown = params.modDown;
		fcMeter = params.fcMeter;

		return *this;
	}

	// --- individual parameters
	double fc = 0.0;				///< filter fc
	double Q = 0.707;				///< filter Q
	double attackTime_mSec = 10.0;	///< detector attack time
	double releaseTime_mSec = 10.0;	///< detector release time
	double threshold_dB = 0.0;		///< detector threshold in dB
	double sensitivity = 1.0;		///< detector sensitivity
	double lfoRate = 1.0;
	filterAlgorithm filterAlgorithm = filterAlgorithm::kMMALPF2;
	bool useLFO = false;
	bool modDown = true;
	float fcMeter = 0.f;
};

/**
\class EnvelopeFollowerNew
\ingroup FX-Objects
\brief
The EnvelopeFollowerNew object implements a traditional envelope follower effect modulating a LPR fc value
using the strength of the detected input.

Audio I/O:
- Processes mono input to mono output.

Control I/F:
- Use EnvelopeFollowerNewParameters structure to get/set object params.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class EnvelopeFollowerNew : public IAudioSignalProcessor
{
public:
	EnvelopeFollowerNew() {
		// --- setup the filter
		AudioFilterParameters filterParams;
		filterParams.algorithm = filterAlgorithm::kMMALPF2;
		filterParams.fc = 1000.0;
		filterParams.boostCut_dB = 0.0;
		filterParams.Q = 0.707;
		filter.setParameters(filterParams);

		// --- setup the detector
		AudioDetectorParameters adParams;
		adParams.attackTime_mSec = -1.0;
		adParams.releaseTime_mSec = -1.0;
		adParams.detectMode = TLD_AUDIO_DETECT_MODE_RMS;
		adParams.detect_dB = true;
		adParams.clampToUnityMax = false;
		detector.setParameters(adParams);

		SuperLFOParameters lfoParams;
		lfoParams.frequency_Hz = 1.0;
		lfoParams.mode = LFOMode::kSync;
		lfoParams.waveform = LFOWaveform::kRSH;
		lfoParams.outputAmplitude = 1.0;
		lfo.setParameters(lfoParams);

	}		/* C-TOR */
	~EnvelopeFollowerNew() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		filter.reset(_sampleRate);
		detector.reset(_sampleRate);
		lfo.reset(_sampleRate);

		SuperLFOParameters params = lfo.getParameters();
		params.waveform = LFOWaveform::kRSH;
		lfo.setParameters(params);

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return EnvelopeFollowerNewParameters custom data structure
	*/
	EnvelopeFollowerNewParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param EnvelopeFollowerNewParameters custom data structure
	*/
	void setParameters(const EnvelopeFollowerNewParameters& params)
	{
		AudioFilterParameters filterParams = filter.getParameters();
		AudioDetectorParameters adParams = detector.getParameters();
		SuperLFOParameters lfoParams = lfo.getParameters();

		filterParams.fc = params.fc;
		filterParams.Q = params.Q;
		filterParams.algorithm = params.filterAlgorithm;

		filter.setParameters(filterParams);

		if (params.attackTime_mSec != parameters.attackTime_mSec ||
			params.releaseTime_mSec != parameters.releaseTime_mSec)
		{
			adParams.attackTime_mSec = params.attackTime_mSec;
			adParams.releaseTime_mSec = params.releaseTime_mSec;
			detector.setParameters(adParams);
		}

		lfoParams.frequency_Hz = params.lfoRate;
		lfo.setParameters(lfoParams);

		parameters = params;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the envelope follower to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		AudioFilterParameters filterParams = filter.getParameters();
		
		SignalModulatorOutput lfoOutput = lfo.renderModulatorOutput();	

		if (parameters.useLFO)
		{
			if (filterParams.algorithm == filterAlgorithm::kMMALPF2)
				filterParams.fc = doUnipolarModulationFromMin(lfoOutput.unipolarOutputFromMin, parameters.fc, kMaxFilterFrequency);
			else if (filterParams.algorithm == filterAlgorithm::kBPF2)
				filterParams.fc = doBipolarModulation(lfoOutput.normalOutput, kMinFilterFrequency, kMaxFilterFrequency);
			else if (filterParams.algorithm == filterAlgorithm::kHPF2)
				filterParams.fc = doUnipolarModulationFromMin(lfoOutput.unipolarOutputFromMin, kMinFilterFrequency, parameters.fc);
		}
		else
		{
			double threshValue = pow(10.0, parameters.threshold_dB / 20.0);
			double detect_dB = detector.processAudioSample(xn);
			double detectValue = pow(10.0, detect_dB / 20.0);
			double deltaValue = detectValue - threshValue;
			double modulatorValue = (deltaValue * parameters.sensitivity);
			double min, max, freq, octaves;

			if (filterParams.algorithm == filterAlgorithm::kHPF2)
			{
				min = kMinFilterFrequency;
				max = parameters.fc;
				freq = doUnipolarModulationFromMin(modulatorValue, min, max);
				octaves = freq / min;
				if (!parameters.modDown)
					filterParams.fc = max / octaves;
				else
					filterParams.fc = freq;
			}
			else
			{
				min = parameters.fc;
				max = kMaxFilterFrequency;
				freq = doUnipolarModulationFromMin(modulatorValue, min, max);
				octaves = freq / min;
				if (parameters.modDown)
					filterParams.fc = max / octaves;
				else
					filterParams.fc = freq;
			}
		}

		parameters.fcMeter = filterParams.fc / kMaxFilterFrequency;

		// --- update with new modulated frequency
		filter.setParameters(filterParams);

		// --- perform the filtering operation
		return filter.processAudioSample(xn);
	}

protected:
	EnvelopeFollowerNewParameters parameters; ///< object parameters

	AudioFilter filter;		///< filter to modulate
	AudioDetector detector; ///< detector to track input signal
	SuperLFO lfo;				///< lfo to set cutoff frequency
};