#include <memory>
#include <math.h>
#include "fxobjects.h"

/**
\struct AudioDelayStereoFeedbackParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AudioDelayStereoFeedback object.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AudioDelayStereoFeedbackParameters
{
	AudioDelayStereoFeedbackParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AudioDelayStereoFeedbackParameters& operator=(const AudioDelayStereoFeedbackParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		wetLevel_dB = params.wetLevel_dB;
		dryLevel_dB = params.dryLevel_dB;
		feedbackLeft_Pct = params.feedbackLeft_Pct;
		feedbackRight_Pct = params.feedbackRight_Pct;

		leftDelay_mSec = params.leftDelay_mSec;
		rightDelay_mSec = params.rightDelay_mSec;

		return *this;
	}

	// --- individual parameters
	double wetLevel_dB = -3.0;	///< wet output level in dB
	double dryLevel_dB = -3.0;	///< dry output level in dB
	double feedbackLeft_Pct = 0.0;	///< feedback as a % value
	double feedbackRight_Pct = 0.0;	///< feedback as a % value

	double leftDelay_mSec = 0.0;	///< left delay time
	double rightDelay_mSec = 0.0;	///< right delay time
};

/**
\class AudioDelayStereoFeedback
\ingroup FX-Objects
\brief
The AudioDelayStereoFeedback object implements a stereo audio delay with multiple delay algorithms.

Audio I/O:
- Processes mono input to mono output OR stereo output.

Control I/F:
- Use AudioDelayStereoFeedbackParameters structure to get/set object params.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AudioDelayStereoFeedback : public IAudioSignalProcessor
{
public:
	AudioDelayStereoFeedback() {}		/* C-TOR */
	~AudioDelayStereoFeedback() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- if sample rate did not change
		if (sampleRate == _sampleRate)
		{
			// --- just flush buffer and return
			delayBuffer_L.flushBuffer();
			delayBuffer_R.flushBuffer();
			return true;
		}

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffers(_sampleRate, bufferLength_mSec);

		return true;
	}

	/** process MONO audio delay */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- read delay
		double yn = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- create input for delay buffer
		double dn = xn + (parameters.feedbackLeft_Pct / 100.0) * yn;

		// --- write to delay buffer
		delayBuffer_L.writeBuffer(dn);

		// --- form mixture out = dry*xn + wet*yn
		double output = dryMix * xn + wetMix * yn;

		return output;
	}

	/** return true: this object can also process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay in frames */
	virtual bool processAudioFrame(const float* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- if only one output channel, revert to mono operation
		if (outputChannels == 1)
		{
			// --- process left channel only
			outputFrame[0] = processAudioSample(inputFrame[0]);
			return true;
		}

		// --- if we get here we know we have 2 output channels
		//
		// --- pick up inputs
		//
		// --- LEFT channel
		double xnL = inputFrame[0];

		// --- RIGHT channel (duplicate left input if mono-in)
		double xnR = inputChannels > 1 ? inputFrame[1] : xnL;

		// --- read delay LEFT
		double ynL = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- read delay RIGHT
		double ynR = delayBuffer_R.readBuffer(delayInSamples_R);

		// --- create input for delay buffer with LEFT channel info
		double dnL = xnL + (parameters.feedbackLeft_Pct / 100.0) * ynL;

		// --- create input for delay buffer with RIGHT channel info
		double dnR = xnR + (parameters.feedbackRight_Pct / 100.0) * ynR;

		// --- decode

		// --- write to LEFT delay buffer with LEFT channel info
		delayBuffer_L.writeBuffer(dnL);
		// --- write to RIGHT delay buffer with RIGHT channel info
		delayBuffer_R.writeBuffer(dnR);
	

		// --- form mixture out = dry*xn + wet*yn
		double outputL = dryMix * xnL + wetMix * ynL;

		// --- form mixture out = dry*xn + wet*yn
		double outputR = dryMix * xnR + wetMix * ynR;

		// --- set left channel
		outputFrame[0] = outputL;

		// --- set right channel
		outputFrame[1] = outputR;

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AudioDelayStereoFeedbackParameters custom data structure
	*/
	AudioDelayStereoFeedbackParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AudioDelayStereoFeedbackParameters custom data structure
	*/
	void setParameters(AudioDelayStereoFeedbackParameters _parameters)
	{
		// --- check mix in dB for calc
		if (_parameters.dryLevel_dB != parameters.dryLevel_dB)
			dryMix = pow(10.0, _parameters.dryLevel_dB / 20.0);
		if (_parameters.wetLevel_dB != parameters.wetLevel_dB)
			wetMix = pow(10.0, _parameters.wetLevel_dB / 20.0);

		// --- save; rest of updates are cheap on CPU
		parameters = _parameters;

		// --- check update type first:
	
		// --- set left and right delay times
		// --- calculate total delay time in samples + fraction
		double newDelayInSamples_L = parameters.leftDelay_mSec*(samplesPerMSec);
		double newDelayInSamples_R = parameters.rightDelay_mSec*(samplesPerMSec);

		// --- new delay time with fraction
		delayInSamples_L = newDelayInSamples_L;
		delayInSamples_R = newDelayInSamples_R;
	}

	/** creation function */
	void createDelayBuffers(double _sampleRate, double _bufferLength_mSec)
	{
		// --- store for math
		bufferLength_mSec = _bufferLength_mSec;
		sampleRate = _sampleRate;
		samplesPerMSec = sampleRate / 1000.0;

		// --- total buffer length including fractional part
		bufferLength = (unsigned int)(bufferLength_mSec*(samplesPerMSec)) + 1; // +1 for fractional part

																			   // --- create new buffer
		delayBuffer_L.createCircularBuffer(bufferLength);
		delayBuffer_R.createCircularBuffer(bufferLength);
	}

private:
	AudioDelayStereoFeedbackParameters parameters; ///< object parameters

	double sampleRate = 0.0;		///< current sample rate
	double samplesPerMSec = 0.0;	///< samples per millisecond, for easy access calculation
	double delayInSamples_L = 0.0;	///< double includes fractional part
	double delayInSamples_R = 0.0;	///< double includes fractional part
	double bufferLength_mSec = 0.0;	///< buffer length in mSec
	unsigned int bufferLength = 0;	///< buffer length in samples
	double wetMix = 0.707; ///< wet output default = -3dB
	double dryMix = 0.707; ///< dry output default = -3dB

	// --- delay buffer of doubles
	CircularBuffer<double> delayBuffer_L;	///< LEFT delay buffer of doubles
	CircularBuffer<double> delayBuffer_R;	///< RIGHT delay buffer of doubles
};
