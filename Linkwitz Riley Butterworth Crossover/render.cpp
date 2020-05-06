/*
 * 2-way audio crossover
 * using the BeagleBone Black.
 * this code runs on BeagleBone Black with the Bela/BeagleRT environment.
 */

#include <Bela.h>
#include <cmath>
#include <Utilities.h>
#include <math.h>
#include <Scope.h>

// making it easier to access this. 44100 is set in stone anyway
#define SAMPLE_FREQUENCY 44100
// we only need to keep track of current and 2 past outputs and inputs for butterworth filter
#define BW_BUFFER_LENGTH 3
// we need to keep track of the current and past 4 outputs and inputs for linkwitzRiley filter
#define LR_BUFFER_LENGTH 5
// declare the oscscope
Scope scope;

/*float gFrequency = 500;		// Hz
float gDuration = 5;		// Seconds
float gAmplitude = 1;     // Volume control
float gPhase = 0.0;*/ // generator for testing

// boolean indicating whether the linkwitzRiley flag is set
bool isLR = false;

// Global variables for Low Pass Filter coefficients
float gLPFB0 = 0.0;//0.0046;
float gLPFB1 = 0.0;//0.0092;
float gLPFB2 = 0.0;//0.0046;
float gLPFA1 = 0.0;//-1.7991;
float gLPFA2 = 0.0;//0.8175;

// LinkwitzRiley
float gLPFlrB0 = 0.0;
float gLPFlrB1 = 0.0;
float gLPFlrB2 = 0.0;
float gLPFlrB3 = 0.0;
float gLPFlrB4 = 0.0;

float gLPFlrA1 = 0.0;
float gLPFlrA2 = 0.0;
float gLPFlrA3 = 0.0;
float gLPFlrA4 = 0.0;

// Global variables for High Pass Filter coefficients
float gHPFB0 = 0.0;
float gHPFB1 = 0.0;
float gHPFB2 = 0.0;
float gHPFA1 = 0.0;
float gHPFA2 = 0.0;

// LinkwitzRiley
float gHPFlrB0 = 0.0;
float gHPFlrB1 = 0.0;
float gHPFlrB2 = 0.0;
float gHPFlrB3 = 0.0;
float gHPFlrB4 = 0.0;

float gHPFlrA1 = 0.0;
float gHPFlrA2 = 0.0;
float gHPFlrA3 = 0.0;
float gHPFlrA4 = 0.0;

// Butterworth circular buffers
float bwInputBuffer[BW_BUFFER_LENGTH];
float bwLowPassBuffer[BW_BUFFER_LENGTH];
float bwHighPassBuffer[BW_BUFFER_LENGTH];

// LinkwitzRiley circular buffers
float lrInputBuffer[LR_BUFFER_LENGTH];
float lrLowPassBuffer[LR_BUFFER_LENGTH];
float lrHighPassBuffer[LR_BUFFER_LENGTH];

// circular buffer pointer - we only need one
int bufferPointer = 0;

// setup() is called once before the audio rendering starts.
// Use it to perform any initialisation and allocation which is dependent
// on the period size or sample rate.
//
// userData holds an opaque pointer to a data structure that was passed
// in from the call to initAudio().
//
// Return true on success; returning false halts the program.

bool setup(BelaContext *context, void *userData)
{
	float crossoverFrequency;
	float tangent;
	float denominator;
	// struct passed in from main with command line information
	struct userSettings {float frequency; bool linkwitzRiley;};
	// userData contians one such struct so cast userData and dereference it and put it into uD
	struct userSettings uD = *(userSettings*)userData;

	// Retrieve a parameter passed in from the initAudio() call
	if(userData != 0){
		// Checking what the cutoff frequency selected is on the console
		printf("Frequency selected %f", uD.frequency);
		printf("Is it a linkwitzRiley? %d", uD.linkwitzRiley);

		// use our global variable to be able to tell whether it is LR in render() later
		isLR = uD.linkwitzRiley;

		// convert Hz to angular frequency
		crossoverFrequency = *(float *)userData * 2.0 * M_PI;

		// this is the tangent part of the pre-warping of the input frequency
		// The 2/T is not necessary as this cancels out in the bilinear transformation, so we can preemptively just keep track of this
		// N.B., this is refered to as Tau in the report
		tangent = tan(crossoverFrequency / SAMPLE_FREQUENCY / 2.0);
	} else {
		// this is the default behaviour if the user does not input any cutoff frequency on the command line
		crossoverFrequency = 1000.0 * 2.0 * M_PI;
		tangent = tan(crossoverFrequency / SAMPLE_FREQUENCY / 2.0);
	}
	// set up the oscscope for testing
	scope.setup(2, context->audioSampleRate, 1);

	// since the denominator of each coefficient is the same, it is better to save rather than computing every time
	denominator = 1 + sqrt(2) * tangent + pow(tangent, 2);

	// computing all of the coefficients
	gLPFB0 = pow(tangent, 2) / denominator;
	gLPFB1 = (2 * pow(tangent, 2)) / denominator;
	gLPFB2 = gLPFB0; // the first and third coefficients are the same
	gLPFA1 = (2 * pow(tangent, 2) - 2) / denominator;
	gLPFA2 = (1 + pow(tangent, 2) - sqrt(2) * tangent) / denominator;

	gHPFB0 = 1 / denominator;
	gHPFB1 = -2 / denominator;
	gHPFB2 = gHPFB0; // the first and third coefficients are the same
	gHPFA1 = (2 * pow(tangent, 2) - 2) / denominator;
	gHPFA2 = (pow(tangent, 2) - sqrt(2) * tangent + 1) / denominator;

	// crosscheck these against matlab
	printf("LPF B1 is %f, B2 is %f, B3 is %f, A2 is %f, A3 is %f", gLPFB0, gLPFB1, gLPFB2, gLPFA1, gLPFA2);
	printf("HPF B1 is %f, B2 is %f, B3 is %f, A2 is %f, A3 is %f", gHPFB0, gHPFB1, gHPFB2, gHPFA1, gHPFA2);

	// if it is a linkwitzRiley filter, then compute the LR coefficients based on the butterworth ones
	if(isLR){
		gLPFlrB0 = pow(gLPFB0, 2);
		gLPFlrB1 = 2 * gLPFB0 * gLPFB1;
		gLPFlrB2 = 2 * gLPFB0 * gLPFB2 + pow(gLPFB1, 2);
		gLPFlrB3 = 2 * gLPFB1 * gLPFB2;
		gLPFlrB4 = pow(gLPFB2, 2);

		gLPFlrA1 = 2 * gLPFA1;
		gLPFlrA2 = 2 * gLPFA2 + pow(gLPFA1, 2);
		gLPFlrA3 = 2 * gLPFA1 * gLPFA2;
		gLPFlrA4 = pow(gLPFA2, 2);


		gHPFlrB0 = pow(gHPFB0, 2);
		gHPFlrB1 = 2 * gHPFB0 * gHPFB1;
		gHPFlrB2 = 2 * gHPFB0 * gHPFB2 + pow(gHPFB1, 2);
		gHPFlrB3 = 2 * gHPFB1 * gHPFB2;
		gHPFlrB4 = pow(gHPFB2, 2);

		gHPFlrA1 = 2 * gHPFA1;
		gHPFlrA2 = 2 * gHPFA2 + pow(gHPFA1, 2);
		gHPFlrA3 = 2 * gHPFA1 * gHPFA2;
		gHPFlrA4 = pow(gHPFA2, 2);

		// initialise circular buffers for lr filter
		for (int i = 0; i < LR_BUFFER_LENGTH; i++){
			lrInputBuffer[i] = 0.0;
			lrLowPassBuffer[i] = 0.0;
			lrHighPassBuffer[i] = 0.0;
		}
	} else {
		// initialise circular buffers for bw filter
		for (int i = 0; i < BW_BUFFER_LENGTH; i++){
			bwInputBuffer[i] = 0.0;
			bwLowPassBuffer[i] = 0.0;
			bwHighPassBuffer[i] = 0.0;
		}
	}

	return true;
}

// render() is called regularly at the highest priority by the audio engine.
// Input and output are given from the audio hardware and the other
// ADCs and DACs (if available). If only audio is available, numMatrixFrames
// will be 0.

void render(BelaContext *context, void *userData)
{

	for (unsigned int n = 0; n < context->audioFrames; n++){

		// for linkwitzRiley filter
		if (isLR){
				// Input -- Take both channels and mix them into a mono by adding them to the lrInputBuffer
			lrInputBuffer[bufferPointer] = 0.0; // clear the lrInputBuffer
			for (unsigned int channel = 0; channel < context->audioInChannels; channel++){
				// add the inputs from the two channels into a mono mix.
				lrInputBuffer[bufferPointer] += context->audioIn[n * context->audioInChannels + channel];
			}

			// Difference Equation LPF linkwitzRiley
			lrLowPassBuffer[bufferPointer] =
				lrInputBuffer[bufferPointer] * gLPFlrB0 +
				lrInputBuffer[(bufferPointer + 4) % LR_BUFFER_LENGTH] * gLPFlrB1 +
				lrInputBuffer[(bufferPointer + 3) % LR_BUFFER_LENGTH] * gLPFlrB2 +
				lrInputBuffer[(bufferPointer + 2) % LR_BUFFER_LENGTH] * gLPFlrB3 +
				lrInputBuffer[(bufferPointer + 1) % LR_BUFFER_LENGTH] * gLPFlrB4 -
				lrLowPassBuffer[(bufferPointer + 4) % LR_BUFFER_LENGTH] * gLPFlrA1 -
				lrLowPassBuffer[(bufferPointer + 3) % LR_BUFFER_LENGTH] * gLPFlrA2 -
				lrLowPassBuffer[(bufferPointer + 2) % LR_BUFFER_LENGTH] * gLPFlrA3 -
				lrLowPassBuffer[(bufferPointer + 1) % LR_BUFFER_LENGTH] * gLPFlrA4;

			// Difference Equation HPF linkwitzRiley
			lrHighPassBuffer[bufferPointer] =
				lrInputBuffer[bufferPointer] * gHPFlrB0 +
				lrInputBuffer[(bufferPointer + 4) % LR_BUFFER_LENGTH] * gHPFlrB1 +
				lrInputBuffer[(bufferPointer + 3) % LR_BUFFER_LENGTH] * gHPFlrB2 +
				lrInputBuffer[(bufferPointer + 2) % LR_BUFFER_LENGTH] * gHPFlrB3 +
				lrInputBuffer[(bufferPointer + 1) % LR_BUFFER_LENGTH] * gHPFlrB4 -
				lrHighPassBuffer[(bufferPointer + 4) % LR_BUFFER_LENGTH] * gHPFlrA1 -
				lrHighPassBuffer[(bufferPointer + 3) % LR_BUFFER_LENGTH] * gHPFlrA2 -
				lrHighPassBuffer[(bufferPointer + 2) % LR_BUFFER_LENGTH] * gHPFlrA3 -
				lrHighPassBuffer[(bufferPointer + 1) % LR_BUFFER_LENGTH] * gHPFlrA4;

			// test outputs with oscscope
			scope.log(lrLowPassBuffer[bufferPointer], lrHighPassBuffer[bufferPointer]);

			// Output
			for (unsigned int channel = 0; channel < context->audioOutChannels; channel++){
				if (channel == 0)
				{
					// pipe lowpass through L channel
					context->audioOut[n * context->audioOutChannels + channel] = lrLowPassBuffer[bufferPointer];
				} else if (channel == 1)
				{
					// pipe highpass through R channel
					context->audioOut[n * context->audioOutChannels + channel] = lrHighPassBuffer[bufferPointer];
				}

			}
			// increment the bufferPointer wrapping around buffer if need be
			bufferPointer = (bufferPointer + 1) % LR_BUFFER_LENGTH;
		// do the butterworth filter
		} else {

			/*float out = gAmplitude* sin(gPhase);
			gPhase = fmodf(2 * M_PI * gFrequency/context->audioSampleRate + gPhase, 2 * M_PI);*/ //testing

			// Input -- Take both channels and mix them into a mono by adding them to the bwInputBuffer
			bwInputBuffer[bufferPointer] = 0.0; // clear the bwInputBuffer
			for (unsigned int channel = 0; channel < context->audioInChannels; channel++){
				// add the inputs from the two channels into a mono mix.
				bwInputBuffer[bufferPointer] = context->audioIn[n * context->audioInChannels + channel];
			}

			// Difference Equation LPF
			bwLowPassBuffer[bufferPointer] =
				bwInputBuffer[bufferPointer] * gLPFB0 +
				bwInputBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gLPFB1 +
				bwInputBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gLPFB2 -
				bwLowPassBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gLPFA1 -
				bwLowPassBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gLPFA2;

			// Difference Equation HPF
			bwHighPassBuffer[bufferPointer] =
				bwInputBuffer[bufferPointer] * gHPFB0 +
				bwInputBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gHPFB1 +
				bwInputBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gHPFB2 -
				bwHighPassBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gHPFA1 -
				bwHighPassBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gHPFA2;

			// test outputs with oscscope
			scope.log(bwLowPassBuffer[bufferPointer], bwHighPassBuffer[bufferPointer]);

			// Output
			for (unsigned int channel = 0; channel < context->audioOutChannels; channel++){
				if (channel == 0)
				{
					// pipe lowpass through L channel
					context->audioOut[n * context->audioOutChannels + channel] = bwLowPassBuffer[bufferPointer];
				} else if (channel == 1)
				{
					// pipe highpass through R channel
					context->audioOut[n * context->audioOutChannels + channel] = bwHighPassBuffer[bufferPointer];
				}

			}
			// increment the bufferPointer wrapping around buffer if need be
			bufferPointer = (bufferPointer + 1) % BW_BUFFER_LENGTH;
		}

	}

}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext *context, void *userData)
{
	/* TASK:
	 * If you allocate any memory, be sure to release it here.
	 * You may or may not need anything in this function, depending
	 * on your implementation.
	 */
}
