#include <memory>
#include <math.h>
#include "fxobjects.h"

/**
\struct AudioDelayNewParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AudioDelayNew object.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AudioDelayNewParameters
{
	AudioDelayNewParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AudioDelayNewParameters& operator=(const AudioDelayNewParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;
		wetLevel_dB = params.wetLevel_dB;
		dryLevel_dB = params.dryLevel_dB;
		feedback_Pct = params.feedback_Pct;
		for (int i = 0; i < 4; i++)
		{
			delay_mSec[i] = params.delay_mSec[i];
			tapGain[i] = params.tapGain[i];
			enable[i] = params.enable[i];
		}
		duckingMode = params.duckingMode;
		duckingSensitivity = params.duckingSensitivity;
		duckingThreshold = params.duckingThreshold;
		attMeter = params.attMeter;
		scMeter = params.scMeter;
		lowCut = params.lowCut;
		return *this;
	}

	// --- individual parameters
	double wetLevel_dB = -3.0;	///< wet output level in dB
	double dryLevel_dB = -3.0;	///< dry output level in dB
	double feedback_Pct;	///< feedback as a % value
	double delay_mSec[4];
	double tapGain[4];
	bool enable[4] = { 0,0,0,0 };
	int duckingMode;
	double duckingSensitivity;
	double duckingThreshold;
	float attMeter = 0.f;
	float scMeter = 0.f;
	int lowCut;
};

/**
\class AudioDelayNew
\ingroup FX-Objects
\brief
The AudioDelayNew object implements a stereo audio delay with multiple delay algorithms.

Audio I/O:
- Processes mono input to mono output OR stereo output.

Control I/F:
- Use AudioDelayNewParameters structure to get/set object params.

\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AudioDelayNew : public IAudioSignalProcessor
{
public:
	AudioDelayNew() {

		AudioDetectorParameters adParams;
		adParams.attackTime_mSec = -1.0;
		adParams.releaseTime_mSec = -1.0;
		adParams.detectMode = TLD_AUDIO_DETECT_MODE_RMS;
		adParams.detect_dB = true;
		adParams.clampToUnityMax = false;
		duckingDetector.setParameters(adParams);
	
	}		/* C-TOR */
	~AudioDelayNew() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- if sample rate did not change
		if (sampleRate == _sampleRate)
		{
			// --- just flush buffer and return
			delayBuffer.flushBuffer();
			return true;
		}

		// --- create new buffer, will store sample rate and length(mSec)
		bufferLength_mSec = 2000.0;
		createDelayBuffers(_sampleRate, bufferLength_mSec);

		feedbackCrusher.reset(_sampleRate);
		BitCrusherParameters feedbackCrusherParams;
		feedbackCrusherParams.quantizedBitDepth = 7.0;
		feedbackCrusher.setParameters(feedbackCrusherParams);

		feedbackFilter.reset(_sampleRate);
		AudioFilterParameters feedbackFilterParams;
		feedbackFilterParams.algorithm = filterAlgorithm::kLPF1;
		feedbackFilterParams.fc = 2200.0;
		feedbackFilterParams.boostCut_dB = 0.0;
		feedbackFilterParams.Q = 0.707;
		feedbackFilter.setParameters(feedbackFilterParams);

		detectHPF.reset(_sampleRate);
		AudioFilterParameters detectHPFParams;
		detectHPFParams.algorithm = filterAlgorithm::kHPF1;
		detectHPFParams.fc = 400.0;
		detectHPFParams.boostCut_dB = 0.0;
		detectHPFParams.Q = 0.707;
		detectHPF.setParameters(detectHPFParams);

		return true;
	}

	/** process MONO audio delay */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- send either sidechain or input audio to detector
		double detectSample, filtered;

		if (parameters.duckingMode == 2 && enableSidechain)
			detectSample = sidechainInputSample;
		else
			detectSample = xn;

		//always process filter, use as detector input if enabled
		filtered = detectHPF.processAudioSample(detectSample);
		if (parameters.lowCut == 1)
			detectSample = filtered * 1.2;

		detect_dB = duckingDetector.processAudioSample(detectSample);

		double detect = pow(10.0, detect_dB / 20.0);
		if (parameters.duckingMode == 2)
			parameters.scMeter = detect;

		if (parameters.duckingMode != 0)								//ducker engaged
		{
			delta_dB = detect_dB - parameters.duckingThreshold;			//dB difference of input above threshold
			if (delta_dB > 0.0)											//positive dB only
			{
				delta_dB = delta_dB *parameters.duckingSensitivity;	//dB difference gets larger when threshold lower
				attenuation = pow(10.0, delta_dB / 20.0);				//linear value above 1.0
				attenuation = 1.0 / attenuation;						//linear value below 1.0
				parameters.attMeter = attenuation;
			}
		}
		else
		{
			attenuation = 1.0;
			parameters.attMeter = 0.0;
		}
			
		

		double feedback;	//only longest delay will feedback
		double yn = 0.0;
		double d[4];
		for (int i = 0; i < 4; i++)
		{
			if (parameters.enable[i])
			{
				d[i] = delayBuffer.readBuffer(delay_InSamples[i]) * tapCooked[i];
				yn += d[i];
				feedback = d[i];
			}
		}

		// --- process feedback, then write into delay buffer
		feedback = feedbackCrusher.processAudioSample(feedback);
		feedback = feedbackFilter.processAudioSample(feedback);

		double dn = xn + (parameters.feedback_Pct / 100.0) * feedback;
		delayBuffer.writeBuffer(dn);

		// --- form mixture out = dry*xn + wet*ducking*yn
		double output = dryMix * xn + wetMix * attenuation * yn;
		return output;	//have amplitudes for each delay (consider having feedback only apply to the LONGEST delay time?)
	}

	/** return false: this object does not process frames */
	virtual bool canProcessAudioFrame() { return false; }

	/** don't implement this */
	virtual bool processAudioFrame(const float* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		return false;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AudioDelayNewParameters custom data structure
	*/
	AudioDelayNewParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AudioDelayNewParameters custom data structure
	*/
	void setParameters(AudioDelayNewParameters _parameters)
	{
		// IDEA: OPTION FOR GRAINIER SOUND - DISABLE LINEAR INTERPOLATION? delayBuffer object has function: void setInterpolate(bool b)

		for (int i = 0; i < 4; i++)
		{
			if(_parameters.tapGain[i] != parameters.tapGain[i])
				tapCooked[i] = pow(10.0, _parameters.tapGain[i] / 20.0);
			//calculate total delay times (samples + fraction)
			delay_InSamples[i] = parameters.delay_mSec[i] * (samplesPerMSec);
		}
		if (_parameters.dryLevel_dB != parameters.dryLevel_dB)
			dryMix = pow(10.0, _parameters.dryLevel_dB / 20.0);
		if (_parameters.wetLevel_dB != parameters.wetLevel_dB)
			wetMix = pow(10.0, _parameters.wetLevel_dB / 20.0);

		AudioDetectorParameters adParams = duckingDetector.getParameters();
		if (adParams.attackTime_mSec != attackTime ||
			adParams.releaseTime_mSec != releaseTime)
		{
			adParams.attackTime_mSec = attackTime;
			adParams.releaseTime_mSec = releaseTime;
			duckingDetector.setParameters(adParams);
		}

		// --- save; rest of updates are cheap on CPU
		parameters = _parameters;
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
		delayBuffer.createCircularBuffer(bufferLength);
	}

	/** enable sidechain input */
	virtual void enableAuxInput(bool enableAuxInput) { enableSidechain = enableAuxInput; }

	/** process the sidechain by saving the value for the upcoming processAudioSample() call */
	virtual double processAuxInputAudioSample(double xn)
	{
		sidechainInputSample = xn;
		//parameters.scMeter = xn;
		return sidechainInputSample;
	}

private:
	AudioDelayNewParameters parameters; ///< object parameters

	double sampleRate;		///< current sample rate
	double samplesPerMSec;	///< samples per millisecond, for easy access calculation
	double delay_InSamples[4];
	double bufferLength_mSec;	///< buffer length in mSec
	unsigned int bufferLength = 0;	///< buffer length in samples
	double wetMix = 0.707; ///< wet output default = -3dB
	double dryMix = 0.707; ///< dry output default = -3dB
	double tapCooked[4];
	double attenuation, attenuation_dB;
	double detect_dB, delta_dB;
	double attackTime = 40.0;
	double releaseTime = 300.0;
	double sidechainInputSample = 0.0;
	bool enableSidechain = false;

	// --- delay buffer of doubles
	CircularBuffer<double> delayBuffer;	///< LEFT delay buffer of doubles
	AudioFilter feedbackFilter;
	BitCrusher feedbackCrusher;
	AudioDetector duckingDetector;
	AudioFilter detectHPF;
};

