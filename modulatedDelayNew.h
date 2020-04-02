#include <memory>
#include <math.h>
#include "fxobjects.h"
#include "audioDelayStereoFeedback.h"

enum class modDelayNewAlgorithm { kNone, kFlanger, kCircle_Flanger, kChorus, kWide_Chorus, kWow, kFlutter, kVibrato };

/**
\struct ModulatedDelayNewParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the ModulatedDelayNew object.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct ModulatedDelayNewParameters
{
	ModulatedDelayNewParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ModulatedDelayNewParameters& operator=(const ModulatedDelayNewParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		algorithm = params.algorithm;
		trim_Pct = params.trim_Pct;
		lfoDepth_Pct = params.lfoDepth_Pct;
		feedback_Pct = params.feedback_Pct;
		return *this;
	}

	// --- individual parameters
	modDelayNewAlgorithm algorithm = modDelayNewAlgorithm::kNone; ///< mod delay algorithm
	double trim_Pct = 0.0;	///< mod delay LFO rate in Hz
	double lfoDepth_Pct = 0.0;	///< mod delay LFO depth in %
	double feedback_Pct = 0.0;	///< feedback in %
};

/**
\class ModulatedDelayNew
\ingroup FX-Objects
\brief
The ModulatedDelayNew object implements the three basic algorithms: flanger, chorus, vibrato.

Audio I / O :
	-Processes mono input to mono OR stereo output.

Control I / F :
	-Use ModulatedDelayNewParameters structure to get / set object params.

\author Will Pirkle http ://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed.by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ModulatedDelayNew : public IAudioSignalProcessor
{
public:
	ModulatedDelayNew() {
	}		/* C-TOR */
	~ModulatedDelayNew() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- create new buffer, 100mSec long
		delay.reset(_sampleRate);
		monoDelay.reset(_sampleRate);
		delay.createDelayBuffers(_sampleRate, 100.0);
		monoDelay.createDelayBuffers(_sampleRate, 100.0);

		// --- lfo
		lfo.reset(_sampleRate);
		OscillatorParameters params = lfo.getParameters();
		params.waveform = generatorWaveform::kTriangle;
		lfo.setParameters(params);

		return true;
	}

	/** process input sample */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		float input = xn;
		float output = 0.0;
		processAudioFrame(&input, &output, 1, 1);
		return output;
	}

	/** return true: this object can process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay of frames */
	virtual bool processAudioFrame(const float* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- render LFO
		SignalGenData lfoOutput = lfo.renderAudioOutput();

		// --- setup delay modulation
		AudioDelayStereoFeedbackParameters params = delay.getParameters();
		AudioDelayParameters monoParams = monoDelay.getParameters();
		double minDelay_mSec = 0.0;
		double maxDepth_mSec = 0.0;

		OscillatorParameters lfoParams = lfo.getParameters();
		if (parameters.algorithm == modDelayNewAlgorithm::kVibrato || parameters.algorithm == modDelayNewAlgorithm::kWow || parameters.algorithm == modDelayNewAlgorithm::kFlutter)
			lfoParams.waveform = generatorWaveform::kSin;
		else
			lfoParams.waveform = generatorWaveform::kTriangle;

		// --- set delay times, wet/dry and feedback


		if (parameters.algorithm == modDelayNewAlgorithm::kFlanger || parameters.algorithm == modDelayNewAlgorithm::kCircle_Flanger)
		{
			minDelay_mSec = 0.1;
			maxDepth_mSec = 7.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -3.0;
			lfoParams.frequency_Hz = 0.45 + (parameters.trim_Pct / 100.0) * 0.42;
		}

		if (parameters.algorithm == modDelayNewAlgorithm::kChorus || parameters.algorithm == modDelayNewAlgorithm::kWide_Chorus)
		{
			minDelay_mSec = 10.0;
			maxDepth_mSec = 28.0 - (parameters.trim_Pct / 100.0) * 5.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -0.0;
			params.feedbackLeft_Pct = 0.0;
			params.feedbackRight_Pct = 0.0;
			lfoParams.frequency_Hz = 0.27 + (parameters.trim_Pct / 100.0) * 0.26;
		}
		if (parameters.algorithm == modDelayNewAlgorithm::kWow)
		{
			minDelay_mSec = 0.0;
			maxDepth_mSec = 14.0;
			params.wetLevel_dB = 0.0;
			params.dryLevel_dB = -96.0;
			params.feedbackLeft_Pct = 0.0;
			params.feedbackRight_Pct = 0.0;
			lfoParams.frequency_Hz = 0.65 + (parameters.trim_Pct / 100.0) * 0.35;
		}
		if (parameters.algorithm == modDelayNewAlgorithm::kFlutter)
		{
			minDelay_mSec = 0.0;
			maxDepth_mSec = 0.16 - (parameters.trim_Pct / 100.0) * 0.04;
			params.wetLevel_dB = 0.0;
			params.dryLevel_dB = -96.0;
			params.feedbackLeft_Pct = 0.0;
			params.feedbackRight_Pct = 0.0;
			lfoParams.frequency_Hz = 800.0 + (parameters.trim_Pct / 100.0) * 600.0;
		}
		if (parameters.algorithm == modDelayNewAlgorithm::kVibrato)
		{
			minDelay_mSec = 0.0;
			maxDepth_mSec = 7.0 - (parameters.trim_Pct / 100.0) * 4.0;
			params.wetLevel_dB = 0.0;
			params.dryLevel_dB = -96.0;
			params.feedbackLeft_Pct = 0.0;
			params.feedbackRight_Pct = 0.0;
			lfoParams.frequency_Hz = 4.0 + (parameters.trim_Pct / 100.0) * 3.5;
		}

		lfo.setParameters(lfoParams);


		// --- calc modulated delay times
		double depth = parameters.lfoDepth_Pct / 100.0;
		double modulationMin = minDelay_mSec;
		double modulationMax = minDelay_mSec + maxDepth_mSec;

		// --- flanger - unipolar
		if (parameters.algorithm == modDelayNewAlgorithm::kFlanger)
		{
			params.leftDelay_mSec = doUnipolarModulationFromMin(bipolarToUnipolar(depth * lfoOutput.normalOutput),
				modulationMin, modulationMax);
			params.feedbackLeft_Pct = 5.0;
			params.feedbackRight_Pct = 5.0;
		}
		else if (parameters.algorithm == modDelayNewAlgorithm::kWide_Chorus)
		{
			params.leftDelay_mSec = doBipolarModulation(depth * lfoOutput.quadPhaseOutput_pos, modulationMin, modulationMax);
			monoParams.leftDelay_mSec = doBipolarModulation(depth * lfoOutput.normalOutput, modulationMin, modulationMax);
			params.rightDelay_mSec = doBipolarModulation(depth * lfoOutput.quadPhaseOutput_neg, modulationMin, modulationMax);
		}
		else
		{
			params.leftDelay_mSec = doBipolarModulation(depth * lfoOutput.normalOutput, modulationMin, modulationMax);
			params.rightDelay_mSec = params.leftDelay_mSec;
		}

		// --- wide flanger - different left/right feedbacks
		if (parameters.algorithm == modDelayNewAlgorithm::kCircle_Flanger)
		{
			params.feedbackLeft_Pct = doUnipolarModulationFromMax(lfoOutput.quadPhaseOutput_neg, 10.0, 96.0);
			params.feedbackRight_Pct = doUnipolarModulationFromMin(lfoOutput.quadPhaseOutput_pos, 10.0, 96.0);
		}


		// --- modulate the delay
		delay.setParameters(params);

		double monoSample = 0.0;
		monoParams.dryLevel_dB = params.dryLevel_dB;
		monoParams.wetLevel_dB = params.wetLevel_dB;
		monoParams.feedback_Pct = params.feedbackLeft_Pct;

		monoDelay.setParameters(monoParams);
		double monoSampleIn;
		if (inputChannels > 1)
			monoSampleIn = 0.5*inputFrame[0] + 0.5*inputFrame[1];
		else
			monoSampleIn = inputFrame[0];
		monoSample = monoDelay.processAudioSample(monoSampleIn);

		// --- just call the function and pass our info in/out
		bool processed = delay.processAudioFrame(inputFrame, outputFrame, inputChannels, outputChannels);
		if (parameters.algorithm == modDelayNewAlgorithm::kWide_Chorus)
		{
			if (outputChannels > 1)
			{
				outputFrame[0] = outputFrame[0] + 0.5*monoSample;
				outputFrame[1] = outputFrame[1] + 0.5*monoSample;
			}
			else
				outputFrame[0] = outputFrame[0] + monoSample;
		}
		return processed;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ModulatedDelayNewParameters custom data structure
	*/
	ModulatedDelayNewParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ModulatedDelayNewParameters custom data structure
	*/
	void setParameters(ModulatedDelayNewParameters _parameters)
	{
		// --- bulk copy
		parameters = _parameters;
	}

private: 
	ModulatedDelayNewParameters parameters; ///< object parameters
	AudioDelayStereoFeedback delay;	///< the delay to modulate
	AudioDelay monoDelay; // extra for chorus
	LFO lfo;			///< the modulator
};
