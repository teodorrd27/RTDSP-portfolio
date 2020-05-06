/*
 * drum machine which plays sampled drum sounds in loops.
 *
 * This code runs on BeagleBone Black with the Bela/BeagleRT environment.
 *
 */


#include <Bela.h>
#include <cmath>
#include <Scope.h>
#include "drums.h"
#include "helpers.h"

#define ANALOG_SAMPLE_FREQUENCY 22050
#define BW_BUFFER_LENGTH 3


/* Variables which are given to you: */

/* Drum samples are pre-loaded in these buffers. Length of each
 * buffer is given in gDrumSampleBufferLengths.
 */
extern float *gDrumSampleBuffers[NUMBER_OF_DRUMS];
extern int gDrumSampleBufferLengths[NUMBER_OF_DRUMS];

int gIsPlaying = 0;			/* Whether we should play or not. Implement this in Step 4b. */

/* Read pointer into the current drum sample buffer.
 *
 * TODO (step 3): you will replace this with two arrays, one
 * holding each read pointer, the other saying which buffer
 * each read pointer corresponds to.
 */
int gReadPointers[16];
int gDrumBufferForReadPointer[16];

/* Patterns indicate which drum(s) should play //////on which beat.
 * Each element of gPatterns is an array, whose length is given
 * by gPatternLengths.
 */
extern int *gPatterns[NUMBER_OF_PATTERNS];
extern int gPatternLengths[NUMBER_OF_PATTERNS];

/* These variables indicate which pattern we're playing, and
 * where within the pattern we currently are. Used in Step 4c.
 */
int gCurrentPattern = 3;
int gCurrentIndexInPattern = 0;

/* Triggers from buttons (step 2 etc.). Read these here and
 * do something if they are nonzero (resetting them when done). */
int gButton1;
int gButton2;

int gTriggerButton1;
int gTriggerButton2;

/* This variable holds the interval between events in **milliseconds**
 * To use it (Step 4a), you will need to work out how many samples
 * it corresponds to.
 */
int gEventIntervalMilliseconds = 300;
int gAudioFramesPerAnalogFrame;
int targetSamples;

/* This variable indicates whether samples should be triggered or
 * not. It is used in Step 4b, and should be set in gpio.cpp.
 */
extern int gIsPlaying;

/* This indicates whether we should play the samples backwards.
 */
int gPlaysBackwards = 0;

/* For bonus step only: these variables help implement a fill
 * (temporary pattern) which is triggered by tapping the board.
 */
int gShouldPlayFill = 0;
int gPreviousPattern = 0;
int gPreviousPatternShouldPlay = 0;
/* TODO: Declare any further global variables you need here */

int gSampleInterval;
int gLEDTimer;
/* Accelerometer variables */
float x_calibration;
float y_calibration;
float z_calibration;

float x_axis;
float y_axis;
float z_axis;

float x_gravity = 0;
float y_gravity = 0;
float z_gravity = 1;

int x_trig;
int y_trig;
int z_trig;

int readAccelOnce = 1;
// Schmitt Trigger variables
float schmittTrigMemory[3];
float schmittTrigTolerance;
// Debounce variable
int tapDebounce = 0;
int analogTimer = 0;

// High Pass Butterworth
// Global variables for High Pass Filter coefficients
float bwHighPassBuffer[BW_BUFFER_LENGTH];
float bwInputBuffer[BW_BUFFER_LENGTH];
// circular buffer pointer
int bufferPointer = 0;

float gHPFB0 = 0.0;
float gHPFB1 = 0.0;
float gHPFB2 = 0.0;
float gHPFA1 = 0.0;
float gHPFA2 = 0.0;

// Test the high pass filter
Scope scope;

// setup() is called once before the audio rendering starts.

bool setup(BelaContext *context, void *userData)
{
	/* Step 2: initialise GPIO pins */
	gLEDTimer = -1;
	pinMode(context, 0, P8_07, INPUT);	 // Button 1
	pinMode(context, 0, P8_08, INPUT);	// Button 2
	pinMode(context, 0, P8_09, OUTPUT); // LED
	// initialise gReadPointers to off state -1
	for (int i = 0; i < 16; i++){
		gReadPointers[i] = -1;
		gDrumBufferForReadPointer[i] = -1;
	}
	// we don't need the Schmitt Trigger for this implementation
	schmittTrigTolerance = 0.0;
	gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames; // useful calculation

	/* butterworth implementation */

	// initialise circular buffers
	for (int i = 0; i < BW_BUFFER_LENGTH; i++){
			bwInputBuffer[i] = 0.0;
			bwHighPassBuffer[i] = 0.0;
	}

	float crossoverFrequency = 10000.0 * 2.0 * M_PI;
	float tangent = tan(crossoverFrequency / ANALOG_SAMPLE_FREQUENCY / 2.0);
	float denominator = 1 + sqrt(2) * tangent + pow(tangent, 2);

	gHPFB0 = 1 / denominator;
	gHPFB1 = -2 / denominator;
	gHPFB2 = gHPFB0; // the first and third coefficients are the same
	gHPFA1 = (2 * pow(tangent, 2) - 2) / denominator;
	gHPFA2 = (pow(tangent, 2) - sqrt(2) * tangent + 1) / denominator;

	// monitor butterworth high pass filter output
	scope.setup(2, context->audioSampleRate, 1);
	return true;
}

// render() is called regularly at the highest priority by the audio engine.

void render(BelaContext *context, void *userData)
{
	float out = 0.0; // DANGER: this needs to be a global variable if context->audioFrames == 16 && we need to be able to play more than 16 voices at once
	float analogPot = 0.0; // potentiometer variable

	/* TODO: your audio processing code goes here! */
	for(unsigned int n = 0; n < context->audioFrames; n++){
		out = 0.0;
		checkButtons(context);
		// Do analog frame stuff at half speed
		if (!(n % gAudioFramesPerAnalogFrame)){
			analogPot = analogRead(context, 0, 0);
			analogPot = ((int)(analogPot * 10000)) / 10000.0;
			gEventIntervalMilliseconds = (int)nearbyintf(map(analogPot, 0, 0.831, 50, 1000));
			bool tapOrNot = tapListen((float)((int)(analogRead(context, 0, 4) * 1000)) / 1000);
			if(tapDebounce == 1000){
				if(tapOrNot){
					gShouldPlayFill = 1;
					gPreviousPattern = gCurrentPattern;
					rt_printf("%d gShouldPlayfill", gShouldPlayFill);
					tapDebounce = 0;
				}
			}
			// we only need the following done once every 300 analog frames
			if (analogTimer == 300){
				if(readAccelOnce == 1){
					setupAccel(context);
					rt_printf("%f", bwHighPassBuffer[bufferPointer]);
				}
				calculateAxes(context);
				setPattern();
				analogTimer = 0;
			} // stuff every 300 frames
			if(tapDebounce < 1000){
				tapDebounce++;
			} // debounce housekeeping
			analogTimer++;
		} // analog frame stuff

		// Button trigger on press, not release
		if(gTriggerButton1 == 1){
			gIsPlaying = !gIsPlaying;
			rt_printf("%d gIsPlaying", gIsPlaying);
			gTriggerButton1 = 0;
			gSampleInterval = 0;
		}
		/*if(gTriggerButton2 == 1){
			startPlayingDrum(0);
			startPlayingDrum(6);
			gTriggerButton2 = 0;
		} */ // vestigial button
		if (gSampleInterval == 0){
			targetSamples = (int)nearbyintf(context->audioSampleRate * (float)gEventIntervalMilliseconds * 0.001);
			if (targetSamples > 44100) targetSamples = 44100;
		}
		/** some testing */ // if(gSampleInterval == 0 || gSampleInterval == 22050 || gSampleInterval == 11025 || gSampleInterval == 33075) rt_printf("%d mili\n", targetSamples);
		if((gSampleInterval == targetSamples - 1) && gIsPlaying == 1){
			startNextEvent();
			//rt_printf("%d sample counter\n", gSampleInterval);
			gSampleInterval = -1;
			gLEDTimer = 0;
		}
		flashLED(context);
		// read drums and mix into out
		for(int i = 0; i < 16; i++){
			if(gDrumBufferForReadPointer[i] != -1){

				out += gDrumSampleBuffers[gDrumBufferForReadPointer[i]][gReadPointers[i]];
				if(gPlaysBackwards){
					gReadPointers[i]--;
					if(gReadPointers[i] == 0){
						gReadPointers[i] = -1;
						gDrumBufferForReadPointer[i] = -1;
					}
				} else{
					gReadPointers[i]++;
					if(gReadPointers[i] == gDrumSampleBufferLengths[gDrumBufferForReadPointer[i]]){
						gReadPointers[i] = -1;
						gDrumBufferForReadPointer[i] = -1;
					}
				}
			}
		}
		// put result into output buffer
		for (unsigned int channel = 0; channel < context->audioOutChannels; channel++){
			context->audioOut[n * context->audioOutChannels + channel] = out;
		}

		// tempo housekeeping
		gSampleInterval++;
		if(gSampleInterval == 44100) gSampleInterval = 0;
	}
}

/* Start playing a particular drum sound given by drumIndex.
 */
void startPlayingDrum(int drumIndex) {
	/* TODO in Steps 3a and 3b */
	if(gPlaysBackwards){
		for(int i = 0; i < 16; i++){
			if(gReadPointers[i] == -1){
				gDrumBufferForReadPointer[i] = drumIndex;
				gReadPointers[i] = gDrumSampleBufferLengths[gDrumBufferForReadPointer[i]] - 1;
				break;
			}
		}
	} else {
		for(int i = 0; i < 16; i++){
			if(gReadPointers[i] == -1){
				gReadPointers[i] = 0;
				gDrumBufferForReadPointer[i] = drumIndex;
				break;
			}
		}
	}
}

/* Start playing the next event in the pattern */
void startNextEvent() {
	/* TODO in Step 4 */
	for(int i = 0; i < NUMBER_OF_DRUMS; i++){
		if(eventContainsDrum(gPatterns[gCurrentPattern][gCurrentIndexInPattern], i) == 1){
				startPlayingDrum(i);
		}
	}
	gCurrentIndexInPattern++;
	if(gCurrentIndexInPattern == gPatternLengths[gCurrentPattern]){
		gCurrentIndexInPattern = 0;
		if(gCurrentPattern == FILL_PATTERN){
			gShouldPlayFill = 0; // reset the fill at the end of the previous pattern
			gPreviousPatternShouldPlay = 1;
			rt_printf("gshouldNOTplayfill");
		}
		if(gCurrentPattern == gPreviousPattern){
			gPreviousPatternShouldPlay = 0;
		}
	}
}

/* Returns whether the given event contains the given drum sound */
int eventContainsDrum(int event, int drum) {
	if(event & (1 << drum))
		return 1;
	return 0;
}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext *context, void *userData)
{
	// nothing here
}

// buttons
void checkButtons(BelaContext* context){
	int prevButton1 = gButton1;

	if(digitalRead(context, 0, P8_07) == 0){
		gButton1 = 1;
	} else {
		gButton1 = 0;
	}

	if(prevButton1 < gButton1){
		gTriggerButton1 = 1;
	}

	int prevButton2 = gButton2;

	if(digitalRead(context, 0, P8_08) == 0){
		gButton2 = 1;
	} else {
		gButton2 = 0;
	}

	if(prevButton2 < gButton2){
		gTriggerButton2 = 1;
	}
}
// runs early on in render and sets the calibration axis info
void setupAccel(BelaContext* context){
	x_calibration = map(analogRead(context, 0, 2), 0, 0.831, 0, 1);
	y_calibration = map(analogRead(context, 0, 3), 0, 0.831, 0, 1);
	z_calibration = map(analogRead(context, 0, 4), 0, 0.831, 0, 1);
	rt_printf("%f, %f, %f\n", x_calibration, y_calibration, z_calibration);
	readAccelOnce = 0;
}
// check what pattern to play based on lots of things
void setPattern(){
	if(gCurrentIndexInPattern == 0){
		if(gPreviousPatternShouldPlay){
			rt_printf("%d gprevpatternshouldplay", gPreviousPatternShouldPlay);
			gCurrentPattern = gPreviousPattern;
		} else if(x_trig == 0 && y_trig == 0 && z_trig == 1){
			gCurrentPattern = 0;
			gPlaysBackwards = 0;
		} else if(x_trig == -1 && y_trig == 0 && z_trig == 0){
			gCurrentPattern = 1;
			gPlaysBackwards = 0;
		} else if(x_trig == 1 && y_trig == 0 && z_trig == 0){
			gCurrentPattern = 2;
			gPlaysBackwards = 0;
		} else if(x_trig == 0 && y_trig == -1 && z_trig == 0){
			gCurrentPattern = 3;
			gPlaysBackwards = 0;
		} else if(x_trig == 0 && y_trig == 1 && z_trig == 0){
			gCurrentPattern = 4;
			gPlaysBackwards = 0;
		} else if(x_trig == 0 && y_trig == 0 && z_trig == -1){
			gPlaysBackwards = 1;
		}

		if(gShouldPlayFill){
			gCurrentPattern = FILL_PATTERN;
		}

	}
}
// calculate the orientation of the accelerometer
void calculateAxes(BelaContext* context){
	x_axis = map(analogRead(context, 0, 2), 0, 0.831, 0, 1);
	y_axis = map(analogRead(context, 0, 3), 0, 0.831, 0, 1);
	z_axis = map(analogRead(context, 0, 4), 0, 0.831, 0, 1);

	x_gravity = roundf(map(x_axis - x_calibration, -0.3645, 0.3645, -1.5, 1.5) * 100) / 100;
	y_gravity = roundf(map(y_axis - y_calibration, -0.3645, 0.3645, -1.5, 1.5) * 100) / 100;
	z_gravity = roundf(map(z_axis + (1 - z_calibration), 0.3955, 1.1245, -1.5, 1.5) * 100) / 100;

	/**
		WE DON'T NEED THE SCHMITT TRIGGER, BUT IT'S THERE ANYWAY WITH NO THRESHOLD SET (0.0)
	*/
	// X
	if(schmittTrigMemory[0] > (x_gravity + schmittTrigTolerance) || schmittTrigMemory[0] < (x_gravity - schmittTrigTolerance)){
		if(x_gravity > -0.75 && x_gravity <= 0.75 && (y_trig != 0 || z_trig != 0)){
			x_trig = 0;
		} else if (x_gravity <= -0.75){
			x_trig = -1;
		} else if (x_gravity > 0.75){
			x_trig = 1;
		}
		schmittTrigMemory[0] = x_gravity;
	}
	// Y
	if(schmittTrigMemory[1] > (y_gravity + schmittTrigTolerance) || schmittTrigMemory[1] < (y_gravity - schmittTrigTolerance)){
		if(y_gravity > -0.75 && y_gravity <= 0.75 && (x_trig != 0 || z_trig != 0)){
			y_trig = 0;
		} else if (y_gravity <= -0.75){
			y_trig = -1;
		} else if (y_gravity > 0.75){
			y_trig = 1;
		}
		schmittTrigMemory[1] = y_gravity;
	}
	// Z, this one is a bit funky, it's why we don't need the schmitt trigger
	if(schmittTrigMemory[2] > (z_gravity + schmittTrigTolerance) || schmittTrigMemory[2] < (z_gravity - schmittTrigTolerance)){
		if(z_gravity > -0.75 && z_gravity <= 0.75 && (x_trig != 0 || y_trig != 0)){
			z_trig = 0;
		} else if (z_gravity <= -0.75){
			z_trig = -1;
		} else if (z_gravity > 0.75){
			z_trig = 1;
		}
		schmittTrigMemory[2] = z_gravity;
	}
	//rt_printf("%d x, %d y, %d z\n", x_trig, y_trig, z_trig);
	//rt_printf("X = %f xTrig = %d xAx = %f, Y = %f yTrig = %d yAx = %f, Z = %f zTrig = %d zAx = %f\n", x_gravity, x_trig, x_axis, y_gravity, y_trig, y_axis, z_gravity, z_trig, z_axis);
}

void checkAnalog(BelaContext* context){
	// I used it for testing...?
}

// does what it says on the tin
void flashLED(BelaContext* context){
	if (gLEDTimer == 0){
		digitalWrite(context, 0, P8_09, 1);
	} else if (gLEDTimer > 200) {
		digitalWrite(context, 0, P8_09, 0);
		gLEDTimer = -1;
	}
	if (gLEDTimer != -1){
		gLEDTimer++;
	}
}
// this is the listener for a tap on the z-axis
bool tapListen(float input){
	/** butterowrth stuff in here */
	bwInputBuffer[bufferPointer] = input;
	bwHighPassBuffer[bufferPointer] =
				bwInputBuffer[bufferPointer] * gHPFB0 +
				bwInputBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gHPFB1 +
				bwInputBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gHPFB2 -
				bwHighPassBuffer[(bufferPointer + 2) % BW_BUFFER_LENGTH] * gHPFA1 -
				bwHighPassBuffer[(bufferPointer + 1) % BW_BUFFER_LENGTH] * gHPFA2;
	scope.log(bwHighPassBuffer[bufferPointer]);
	// if there's an insignificant output, return false
	if(bwHighPassBuffer[bufferPointer] < 0.015){
		bufferPointer = (bufferPointer + 1) % BW_BUFFER_LENGTH;
		return false;
	} else if (bwHighPassBuffer[bufferPointer] >= 0.015){
		bufferPointer = (bufferPointer + 1) % BW_BUFFER_LENGTH;
		return true;
	}
	return false;

}
