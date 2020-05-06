/*
 * render.cpp
 */


#include <Bela.h>
#include <ne10/NE10.h>					// neon library
#include <cmath>
#include <Scope.h>
#include "helpers.h"
#include "Intervals.h"
#include "Inversions.h"

#define NUM_CHANNELS 2
#define BUFFER_SIZE 44100		// Notice: BUFFER_SIZE much bigger than gFFTSize
#define FFT_SIZE 8192

Scope scope;

// 4 buttons facilitate passage between states
int gButton1;
int gButton2;
int gButton3;
int gButton4;
// triggers for each of the buttons detect when a button has been pressed and
// an action should be pending
int gTriggerButton1;
int gTriggerButton2;
int gTriggerButton3;
int gTriggerButton4;
// each button facilitates 4 different states. For easy transitioning between
// states, sample dependent timers for each of the buttons are implemented
int gButton1Timer;
int gButton2Timer;
int gButton3Timer;
int gButton4Timer;
// keeps track of how many times (1-4) a button has been pressed within the
// timer's window before the new state comes into effect
int gButton1Clicks;
int gButton2Clicks;
int gButton3Clicks;
int gButton4Clicks;

float analogPot = 0.0; // potentiometer variable
// each mode has 4 states. Hence, although not all in use, there are 256 state
// combinations
int voice_mode;
int position_mode;
int chord_mode;
int inversion_mode;
// stores a pointer to the array returned by the voice swapping functions
double* new_chord;
// arithmetic for diatonic intervals
const int maj_diatonic[] = {1, 3, 4, 6, 8, 10, 11};
/** logical control **/
int diatonic_offset;
int diatonic_offset_trigger;

int inversion_trigger;
int chord_trigger;
int voice_trigger;
int position_trigger;
/*********************/
// keeps track of four part chords for middle part swap
int chord_id = -1;
// handy variable to only read every analog frame occurring twice less often as
// digital frames
int gAudioFramesPerAnalogFrame;
// holds input signal
float gInputBuffer[BUFFER_SIZE];
int gInputBufferPointer = 0;

// hold processed outputs -- 3 voices
float gOutputBuffer[BUFFER_SIZE];
float gOutputBuffer2[BUFFER_SIZE];
float gOutputBuffer3[BUFFER_SIZE];
// circular buffer needs two pointers
int gOutputBufferReadPointer;
int gOutputBufferWritePointer;
// fft circular buffer needs two pointers
int gFFTInputBufferPointer;
int gFFTOutputBufferPointer;
float *gWindowBuffer;

// keep track of samples
int gSampleCount = 0;
// -----------------------------------------------
// Size of STFT window
int gFFTSize = 8192;
// oversampling rate (how many times the windows overlap)
int osamp = 4;
int gHopSize = gFFTSize / osamp; // space between beginning of windows
float gFFTScaleFactor = 0; // scale the output windows for the overlap-add

/** pitch shifter processing variables **/
// keep track of previous window's frequency phases
float gLastPhase[FFT_SIZE];
float gSumPhase[FFT_SIZE];
float gSumPhase2[FFT_SIZE];
float gSumPhase3[FFT_SIZE];
// true frequency arrays
float gAnaFreq[FFT_SIZE];
float gAnaMagn[FFT_SIZE];
// synthesis frequency arrays
float gSynFreq[FFT_SIZE];
float gSynFreq2[FFT_SIZE];
float gSynFreq3[FFT_SIZE];
// synthesis magnitude arrays
float gSynMagn[FFT_SIZE];
float gSynMagn2[FFT_SIZE];
float gSynMagn3[FFT_SIZE];

double real;
double imaginary;
// output magnitude for each voice
double magnitude;
double magnitude2;
double magnitude3;
// output phase for each voice
double phase;
double phase2;
double phase3;
// temporary variable for the phase unwrapping process of each voice
double tmp;
double tmp2;
double tmp3;
// expected frequency calculation between windows
double expct = 2*M_PI*(double)gHopSize/(double)gFFTSize;
double freqPerBin = 22050 / ((double)gFFTSize / 2); // bin width
// unwrapping
int pis;
// pitch shifting for each voice
long pitch_shift_index;
long pitch_shift_index2;
long pitch_shift_index3;

double pitchShift2;
double pitchShift3;
double pitchShift4;
// stores pitch shift for the two voice state
double twoVoiceIntervals[4];
// testing variable to only rt_printf once
int print_counter = 0;
// Sample info

// FFT vars -- one time domain in -- outputs for each voice
ne10_fft_cpx_float32_t* timeDomainIn;
ne10_fft_cpx_float32_t* timeDomainOut;
ne10_fft_cpx_float32_t* timeDomainOut2;
ne10_fft_cpx_float32_t* timeDomainOut3;

ne10_fft_cpx_float32_t* frequencyDomain;
ne10_fft_cpx_float32_t* frequencyDomain2;
ne10_fft_cpx_float32_t* frequencyDomain3;
// configuration
ne10_fft_cfg_float32_t cfg;

// Auxiliary task for calculating FFT -- multithreading
AuxiliaryTask gFFTTask;

void process_fft_background(void *);

// setup() is called once before the audio rendering starts.
// Return true on success; returning false halts the program.

bool setup(BelaContext *context, void *userData)
{
	// buttons setup
	pinMode(context, 0, P8_07, INPUT);	 // Button 1
	pinMode(context, 0, P8_08, INPUT);	// Button 2
	pinMode(context, 0, P8_09, INPUT);	 // Button 3
	pinMode(context, 0, P8_10, INPUT);	// Button 4

	gButton1 = 0;
	gButton2 = 0;
	gButton3 = 0;
	gButton4 = 0;

	gTriggerButton1 = 0;
	gTriggerButton2 = 0;
	gTriggerButton3 = 0;
	gTriggerButton4 = 0;

	gButton1Timer = -1;
	gButton2Timer = -1;
	gButton3Timer = -1;
	gButton4Timer = -1;

	gButton1Clicks = 0;
	gButton2Clicks = 0;
	gButton3Clicks = 0;
	gButton4Clicks = 0;
	// initialise logical triggers and states
	voice_mode = 0;
	position_mode = 0;
	chord_mode = 0;
	inversion_mode = 0;
	// initially there is no pitch shift
	pitchShift2 = 1;
	pitchShift3 = 1;
	pitchShift4 = 1;
	// allocate memory for new_chord
	new_chord = (double*)malloc(3 * sizeof(double));
	// two voice state arrays
	twoVoiceIntervals[0] = Intervals::up[maj_diatonic[1]]; // third up
	twoVoiceIntervals[1] = Intervals::up[maj_diatonic[3]]; // fifth up
	twoVoiceIntervals[2] = Intervals::up[maj_diatonic[4]]; // sixth up
	twoVoiceIntervals[3] = Intervals::up[maj_diatonic[0]]; // second up
	// no diatonic offset to start off with. Just a major third.
	diatonic_offset = 0;
	diatonic_offset_trigger = 0;
	// Testing scope
	scope.setup(1, 44100.0);
	gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames; // useful calculation

	// set up FFT
	gFFTScaleFactor = 1.0f / (float)gFFTSize * 1000;
	rt_printf("scale factor: %f\n", gFFTScaleFactor); // window scaling factor
	// fft arrays initialisation
	timeDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	timeDomainOut = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	timeDomainOut2 = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	timeDomainOut3 = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));

	frequencyDomain = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	frequencyDomain2 = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	frequencyDomain3 = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));

	cfg = ne10_fft_alloc_c2c_float32_neon (gFFTSize);


	// Initialise auxiliary tasks
	if((gFFTTask = Bela_createAuxiliaryTask(&process_fft_background, 90, "fft-calculation")) == 0)
		return false;

	memset(gInputBuffer, 0, BUFFER_SIZE * sizeof(float));
	memset(timeDomainOut, 0, gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	memset(timeDomainOut2, 0, gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	memset(timeDomainOut3, 0, gFFTSize * sizeof (ne10_fft_cpx_float32_t));

	memset(gOutputBuffer, 0, BUFFER_SIZE * sizeof(float));
	memset(gOutputBuffer2, 0, BUFFER_SIZE * sizeof(float));
	memset(gOutputBuffer3, 0, BUFFER_SIZE * sizeof(float));


	// Initialise output pointers
	gOutputBufferReadPointer = 0;
	gOutputBufferWritePointer = gHopSize;

	// Allocate the window buffer based on the FFT size
	gWindowBuffer = (float *)malloc(gFFTSize * sizeof(float));
	if(gWindowBuffer == 0)
		return false;

	// Calculate a Hann window
	for(int n = 0; n < gFFTSize; n++) {
		gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0 * M_PI * n / (float)(gFFTSize - 1)));
	}

	// initialise pitch shifting arrays
	memset(gLastPhase, 0, (FFT_SIZE)*sizeof(float));

	memset(gSumPhase, 0, (FFT_SIZE)*sizeof(float));
	memset(gSumPhase2, 0, (FFT_SIZE)*sizeof(float));
	memset(gSumPhase3, 0, (FFT_SIZE)*sizeof(float));

	memset(gAnaFreq, 0, (FFT_SIZE)*sizeof(float));
	memset(gAnaMagn, 0, (FFT_SIZE)*sizeof(float));

	return true;
}

// This function handles the FFT processing in this example once the buffer has
// been assembled.
void process_fft(float *inBuffer, int inWritePointer, float *outBuffer, float *outBuffer2, float *outBuffer3, int outWritePointer)
{
	// Copy buffer into FFT input
	// earliest sample
	int pointer = (inWritePointer - gFFTSize + BUFFER_SIZE) % BUFFER_SIZE;
	// prepare fft input
	for(int n = 0; n < gFFTSize; n++) {
		timeDomainIn[n].r = (ne10_float32_t) inBuffer[pointer] * gWindowBuffer[n];
		timeDomainIn[n].i = 0;

		// circular buffer pointer
		pointer++;
		if (pointer >= BUFFER_SIZE)
			pointer = 0;
	}

	// Run the FFT
	ne10_fft_c2c_1d_float32_neon (frequencyDomain, timeDomainIn, cfg, 0);
	/** PITCH SHIFTING STAGES **/
	/** ********** Analysis ********* **/
	for(int n = 0; n <= gFFTSize / 2; n++){
		// get real and imaginary values from frequency domain
		real = frequencyDomain[n].r;
		imaginary = frequencyDomain[n].i;

		// compute magnitude and phase
		magnitude = 2 * sqrt(real * real + imaginary * imaginary);
		phase = (float)atan2(imaginary, real);

		// compute phase difference
		tmp = phase - gLastPhase[n];
		gLastPhase[n] = phase; // store phase for future use

		// subtract expected phase difference
		tmp -= (double)n*expct;

		// unwrap delta phase into +/- Pi interval
		pis = tmp/M_PI; // pis --> how many Pi-s in the phase differences' offset?
		if(pis >= 0) pis += pis&1; // if odd add a Pi
		else pis -= pis&1; // if odd subtract a Pi
		tmp -= M_PI * (double)pis; // finally unwrapped into +/- interval

		// get deviation from bin frequency from the +/- Pi
		tmp = osamp*tmp/(2.*M_PI);

		// compute the n-th partials' true frequency
		tmp = (double)n*freqPerBin + tmp*freqPerBin;

		// store magnitude and true frequency in analysis arrdays
		gAnaMagn[n] = magnitude;
		gAnaFreq[n] = tmp;
	}
	/* ********* Processing *********** */
	memset(gSynMagn, 0, (gFFTSize)*sizeof(float));
	memset(gSynFreq, 0, (gFFTSize)*sizeof(float));
	memset(gSynMagn2, 0, (gFFTSize)*sizeof(float));
	memset(gSynFreq2, 0, (gFFTSize)*sizeof(float));
	memset(gSynMagn3, 0, (gFFTSize)*sizeof(float));
	memset(gSynFreq3, 0, (gFFTSize)*sizeof(float));
	// shift up frequencies up to nyquist
	for(int n = 0; n <= gFFTSize / 2; n++){
		pitch_shift_index = n * pitchShift2;
		pitch_shift_index2 = n * pitchShift3;
		pitch_shift_index3 = n * pitchShift4;
		if(pitch_shift_index <= gFFTSize / 2){ // push out frequencies that end up being higher than 22.1KHz
			gSynMagn[pitch_shift_index] = gAnaMagn[n];
			gSynFreq[pitch_shift_index] = gAnaFreq[n] * pitchShift2;
		}
		if(pitch_shift_index2 <= gFFTSize / 2){
			gSynMagn2[pitch_shift_index2] += gAnaMagn[n];
			gSynFreq2[pitch_shift_index2] = gAnaFreq[n] * pitchShift3;
		}
		if(pitch_shift_index3 <= gFFTSize / 2){
			gSynMagn3[pitch_shift_index3] += gAnaMagn[n];
			gSynFreq3[pitch_shift_index3] = gAnaFreq[n] * pitchShift4;
		}
	}

	/* ********* Synthesis of 3 pitch shifted voices ********** */
	for(int n = 0; n <= gFFTSize / 2; n++){ // up to nyquist
		// get magnitude and true frequency from synthesis arrays
		magnitude = gSynMagn[n];
		magnitude2 = gSynMagn2[n];
		magnitude3 = gSynMagn3[n];
		tmp = gSynFreq[n];
		tmp2 = gSynFreq2[n];
		tmp3 = gSynFreq3[n];
		// subtract bin mid frequency
		tmp -= (double)n*freqPerBin;
		tmp2 -= (double)n*freqPerBin;
		tmp3 -= (double)n*freqPerBin;

		// get bin deviation from freq deviation
		tmp /= freqPerBin;
		tmp2 /= freqPerBin;
		tmp3 /= freqPerBin;
		// undo deviation
		tmp = 2.*M_PI*tmp/osamp;
		tmp2 = 2.*M_PI*tmp/osamp;
		tmp3 = 2.*M_PI*tmp/osamp;
		//add the overlap phase advance back in
		tmp += (double)n*expct;
		tmp2 += (double)n*expct;
		tmp3 += (double)n*expct;

		// wrap delta phase to get bin phase
		gSumPhase[n] = fmodf((gSumPhase[n] + tmp), (2 * M_PI));
		gSumPhase2[n] = fmodf((gSumPhase2[n] + tmp2), (2 * M_PI));
		gSumPhase3[n] = fmodf((gSumPhase3[n] + tmp3), (2 * M_PI));

		phase = gSumPhase[n];
		phase2 = gSumPhase2[n];
		phase3 = gSumPhase3[n];

		float cosine = cos(phase);
		float sine = sin(phase);
		float cosine2 = cos(phase2);
		float sine2 = sin(phase2);
		float sine3 = sin(phase3);
		float cosine3 = cos(phase3);
		// get real and imag part and real part back in
		frequencyDomain[n].r = magnitude * cosine;
		frequencyDomain[n].i = magnitude * sine;// * sin(phase);

		frequencyDomain2[n].r = magnitude2 * cosine2;// * cos(phase2);
		frequencyDomain2[n].i = magnitude2 * sine2;// * sin(phase2);

		frequencyDomain3[n].r = magnitude3 * cosine3;// * cos(phase3);
		frequencyDomain3[n].i = magnitude3 * sine3;// * sin(phase3);

	}
	// zero negative frequencies ie from nyquist up
	for (int n = gFFTSize / 2; n < gFFTSize; n++){
		frequencyDomain[n].r = 0;
		frequencyDomain[n].i = 0;

		frequencyDomain2[n].r = 0;
		frequencyDomain2[n].i = 0;

		frequencyDomain3[n].r = 0;
		frequencyDomain3[n].i = 0;
	}

	// Run the inverse FFT on all three synthesised voices
	ne10_fft_c2c_1d_float32_neon (timeDomainOut, frequencyDomain, cfg, 1);
	ne10_fft_c2c_1d_float32_neon (timeDomainOut2, frequencyDomain2, cfg, 1);
	ne10_fft_c2c_1d_float32_neon (timeDomainOut3, frequencyDomain3, cfg, 1);

	// copy time domain out into outputbuffers
	pointer = outWritePointer;
	for(int n=0; n<gFFTSize; n++) {
		outBuffer[pointer] += timeDomainOut[n].r * gFFTScaleFactor;
		outBuffer2[pointer] += timeDomainOut2[n].r * gFFTScaleFactor;
		outBuffer3[pointer] += timeDomainOut3[n].r * gFFTScaleFactor;
		pointer++;
		if(pointer >= BUFFER_SIZE)
			pointer = 0;
	}
}

// Function to process the FFT in a thread at lower priority
void process_fft_background(void *) {
	process_fft(gInputBuffer, gFFTInputBufferPointer, gOutputBuffer, gOutputBuffer2, gOutputBuffer3, gFFTOutputBufferPointer);
}

// render() is called regularly at the highest priority by the audio engine.
void render(BelaContext *context, void *userData)
{
	for(unsigned int n = 0; n < context->audioFrames; n++) {
		checkButtons(context);
		set_diatonic_offset(n, context);
		selectListener();
		voiceLogic();

		// read sample data into input buffer
		gInputBuffer[gInputBufferPointer] = 0.0; // clear out input buffer
		for (unsigned int channel = 0; channel < context->audioInChannels; channel++){
				// add the inputs from the two channels into a mono mix.
				gInputBuffer[gInputBufferPointer] += context->audioIn[n * context->audioInChannels + channel];
		}
		scope.log(gOutputBuffer[gOutputBufferReadPointer]); // bit of testing
		// Copy output buffer to output
		for(int channel = 0; channel < context->audioOutChannels; channel++) {
			// keep own voice output a bit lower since you're singing it as well
			context->audioOut[n * context->audioOutChannels + channel] += gInputBuffer[gInputBufferPointer] /2;
			if(voice_mode > 0){ // 2 voices state // this will be the bass in upper states, and thus should be stronger, hence the x8 multiplier
				context->audioOut[n * context->audioOutChannels + channel] += gOutputBuffer[gOutputBufferReadPointer] * 8;
			}
			if(voice_mode > 1){ // 3 voices state // this will be the alto/ soprano in this and upper state. Should sound clearly.
				context->audioOut[n * context->audioOutChannels + channel] += gOutputBuffer2[gOutputBufferReadPointer] * 4;
			}
			if(voice_mode > 2){ // 4 voices state // this will be the highest voice // Ear always picks out top voice as loudest, hence the mere x2 multiplier
				context->audioOut[n * context->audioOutChannels + channel] += gOutputBuffer3[gOutputBufferReadPointer] * 2;
			}
			// take down the output a bit since there are so many concurrent voices
			context->audioOut[n * context->audioOutChannels + channel] /= 4;
		}

		// clear out output buffers
		gOutputBuffer[gOutputBufferReadPointer] = 0;
		gOutputBuffer2[gOutputBufferReadPointer] = 0;
		gOutputBuffer3[gOutputBufferReadPointer] = 0;

		// advance pointers
		gOutputBufferReadPointer++;
		if(gOutputBufferReadPointer >= BUFFER_SIZE)
			gOutputBufferReadPointer = 0;

		// update write pointer
		gOutputBufferWritePointer++;
		if(gOutputBufferWritePointer >= BUFFER_SIZE)
			gOutputBufferWritePointer = 0;

		gInputBufferPointer++;
		if(gInputBufferPointer >= BUFFER_SIZE)
			gInputBufferPointer = 0;

		// if hop size has ellapsed, execute auxiliary fftTask
		gSampleCount++;
		if(gSampleCount >= gHopSize /* change this condition in part 1 */) {
			gFFTInputBufferPointer = gInputBufferPointer;
			gFFTOutputBufferPointer = gOutputBufferWritePointer;
			Bela_scheduleAuxiliaryTask(gFFTTask);
			// reset sample counter
			gSampleCount = 0;
		}

		/** handle button timers **/
		// Voice state button
		if(gButton1Timer != -1 && gButton1Timer < 44100){
			 gButton1Timer++;
		} else if(gButton1Timer >= 44100) {
			voice_trigger = 1;
			gButton1Timer = -1;
			voice_mode = gButton1Clicks - 1;
			rt_printf("Voice mode, %d", voice_mode);
			gButton1Clicks = 0;
		}
		// Position state button
		if(gButton2Timer != -1 && gButton2Timer < 44100){
			gButton2Timer++;
		} else if (gButton2Timer >= 44100){
			position_trigger = 1;
			gButton2Timer = -1;
			position_mode = gButton2Clicks -1;
			rt_printf("Position mode, %d", position_mode);
			gButton2Clicks = 0;
		}
		// Chord state button
		if(gButton3Timer != -1 && gButton3Timer < 44100){
			gButton3Timer++;
		} else if (gButton3Timer >= 44100){
			chord_trigger = 1;
			gButton3Timer = -1;
			chord_mode = gButton3Clicks -1;
			rt_printf("Chord mode, %d", chord_mode);
			gButton3Clicks = 0;
		}
		// Inversion state button
		if(gButton4Timer != -1 && gButton4Timer < 44100){
			gButton4Timer++;
		} else if (gButton4Timer >= 44100){
			inversion_trigger = 1;
			gButton4Timer = -1;
			inversion_mode = gButton4Clicks -1;
			rt_printf("Inversion mode, %d", inversion_mode);
			gButton4Clicks = 0;
		}
	} // audio frames
} // render

// button handler
void checkButtons(BelaContext* context){
	int prevButton1 = gButton1; // previous 1st button state
	if(digitalRead(context, 0, P8_07) == 0){
		gButton1 = 1;
	} else {
		gButton1 = 0;
	}
	// trigger if 1st button state changed
	if(prevButton1 < gButton1){
		gTriggerButton1 = 1;
	}

	int prevButton2 = gButton2; // previous 2nd button state
	if(digitalRead(context, 0, P8_08) == 0){
		gButton2 = 1;
	} else {
		gButton2 = 0;
	}
	// trigger if 2nd button state changed
	if(prevButton2 < gButton2){
		gTriggerButton2 = 1;
	}

	int prevButton3 = gButton3; // previous 3rd button state
	if(digitalRead(context, 0, P8_09) == 0){
		gButton3 = 1;
	} else {
		gButton3 = 0;
	}
	// trigger if 3rd button state changed
	if(prevButton3 < gButton3){
		gTriggerButton3 = 1;
	}

	int prevButton4 = gButton4; // previous 4th button state
	if(digitalRead(context, 0, P8_10) == 0){
		gButton4 = 1;
	} else {
		gButton4 = 0;
	}
	// trigger if 4th button state changed
	if(prevButton4 < gButton4){
		gTriggerButton4 = 1;
	}
}

/** Handles 2 voice state logic **/
void twoVoiceLogic(){
	// only 2 inversion modes active in 2 voice state
	// potentiometer active in this state
	if(diatonic_offset_trigger || inversion_trigger){
		if(inversion_mode == 0){
			twoVoiceIntervals[0] = Intervals::up[maj_diatonic[1] + diatonic_offset]; // third
			twoVoiceIntervals[1] = Intervals::up[maj_diatonic[3] + diatonic_offset]; // fifth
			twoVoiceIntervals[2] = Intervals::up[maj_diatonic[4] + diatonic_offset]; // sixth
			twoVoiceIntervals[3] = Intervals::up[maj_diatonic[0] + diatonic_offset]; // second
			if(diatonic_offset_trigger == 1)
				diatonic_offset_trigger = 0;
			else inversion_trigger = 0;
		} else if (inversion_mode == 1){
			twoVoiceIntervals[0] = Intervals::down[maj_diatonic[1] + diatonic_offset]; // sixth
			twoVoiceIntervals[1] = Intervals::down[maj_diatonic[3] + diatonic_offset]; // fourth
			twoVoiceIntervals[2] = Intervals::down[maj_diatonic[4] + diatonic_offset]; // third
			twoVoiceIntervals[3] = Intervals::down[maj_diatonic[0] + diatonic_offset]; // seventh
			if(diatonic_offset_trigger) diatonic_offset_trigger = 0;
			else inversion_trigger = 0;
		}
	}
	// check chord_mode and set pitch shifters
	switch(chord_mode){
		case 0: pitchShift2 = twoVoiceIntervals[0];
				break;
		case 1: pitchShift2 = twoVoiceIntervals[1];
				break;
		case 2: pitchShift2 = twoVoiceIntervals[2];
				break;
		case 3: pitchShift2 = twoVoiceIntervals[3];
		//		break;
		default: break;
	}
}

/** handles 3 voice state logic **/
void threeVoiceLogic(){
	// chord button has 4 possible states in 3 voice state
	// inversion button has 4 possible states in 3 voice state
	// potentiometer deactivated
	if(inversion_trigger || chord_trigger || voice_trigger){
		switch (chord_mode){
			case 0: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::maj[0][0];
							pitchShift3 = Inversions::maj[0][1];
							break;
					case 1: pitchShift2 = Inversions::maj[1][0];
							pitchShift3 = Inversions::maj[1][1];
							break;
					case 2: pitchShift2 = Inversions::maj[2][0];
							pitchShift3 = Inversions::maj[2][1];
							break;
					case 3: pitchShift2 = Inversions::maj[3][0];
							pitchShift3 = Inversions::maj[3][1];
							break;
					default: break;
				}
				break;
			}
			case 1: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::min[0][0];
							pitchShift3 = Inversions::min[0][1];
							break;
					case 1: pitchShift2 = Inversions::min[1][0];
							pitchShift3 = Inversions::min[1][1];
							break;
					case 2: pitchShift2 = Inversions::min[2][0];
							pitchShift3 = Inversions::min[2][1];
							break;
					case 3: pitchShift2 = Inversions::min[3][0];
							pitchShift3 = Inversions::min[3][1];
							break;
					default: break;
				}
				break;
			}
			case 2: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::maj_flat_seven[0][0];
							pitchShift3 = Inversions::maj_flat_seven[0][1];
							break;
					case 1: pitchShift2 = Inversions::maj_flat_seven[1][0];
							pitchShift3 = Inversions::maj_flat_seven[1][1];
							break;
					case 2: pitchShift2 = Inversions::maj_flat_seven[2][0];
							pitchShift3 = Inversions::maj_flat_seven[2][1];
							break;
					case 3: pitchShift2 = Inversions::maj_flat_seven[3][0];
							pitchShift3 = Inversions::maj_flat_seven[3][1];
							break;
					default: break;
				}
				break;
			}
			case 3: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::min_flat_seven[0][0];
							pitchShift3 = Inversions::min_flat_seven[0][1];
							break;
					case 1: pitchShift2 = Inversions::min_flat_seven[1][0];
							pitchShift3 = Inversions::min_flat_seven[1][1];
							break;
					case 2: pitchShift2 = Inversions::min_flat_seven[2][0];
							pitchShift3 = Inversions::min_flat_seven[2][1];
							break;
					case 3: pitchShift2 = Inversions::min_flat_seven[3][0];
							pitchShift3 = Inversions::min_flat_seven[3][1];
							break;
					default: break;
				}
				break;
			}
			default: break;
		}
		if(inversion_trigger) inversion_trigger = 0;
		else if (chord_trigger) chord_trigger = 0;
		else voice_trigger = 0;
	}
}

/** handles 4 voice state logic **/
void fourVoiceLogic(){
	// chord button has 4 states in 4 voice state
	// inversion button has 4 possible states in 4 voice state
	// position button has 2 possible states in 4 voice state
	// potentiometer deactivated
	if(inversion_trigger || chord_trigger || voice_trigger){
		rt_printf("Inversion %d, Chord %d ", inversion_mode, chord_mode);
		switch (chord_mode){
			case 0: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::maj_4_part[0][0];
							pitchShift3 = Inversions::maj_4_part[0][1];
							pitchShift4 = Inversions::maj_4_part[0][2];
							chord_id = 0;
							break;
					case 1: pitchShift2 = Inversions::maj_4_part[1][0];
							pitchShift3 = Inversions::maj_4_part[1][1];
							pitchShift4 = Inversions::maj_4_part[1][2];
							chord_id = 1;
							break;
					case 2: pitchShift2 = Inversions::maj_4_part[2][0];
							pitchShift3 = Inversions::maj_4_part[2][1];
							pitchShift4 = Inversions::maj_4_part[2][2];
							chord_id = 2;
							break;
					case 3: pitchShift2 = Inversions::maj_4_part[3][0];
							pitchShift3 = Inversions::maj_4_part[3][1];
							pitchShift4 = Inversions::maj_4_part[3][2];
							chord_id = 3;
							break;
					default: break;
				}
				break;
			}
			case 1: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::min_4_part[0][0];
							pitchShift3 = Inversions::min_4_part[0][1];
							pitchShift4 = Inversions::min_4_part[0][2];
							chord_id = 4;
							break;
					case 1: pitchShift2 = Inversions::min_4_part[1][0];
							pitchShift3 = Inversions::min_4_part[1][1];
							pitchShift4 = Inversions::min_4_part[1][2];
							chord_id = 5;
							break;
					case 2: pitchShift2 = Inversions::min_4_part[2][0];
							pitchShift3 = Inversions::min_4_part[2][1];
							pitchShift4 = Inversions::min_4_part[2][2];
							chord_id = 6;
							break;
					case 3: pitchShift2 = Inversions::min_4_part[3][0];
							pitchShift3 = Inversions::min_4_part[3][1];
							pitchShift4 = Inversions::min_4_part[3][2];
							chord_id = 7;
							break;
					default: break;
				}
				break;
			}
			case 2: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::maj_flat_seven_4_part[0][0];
							pitchShift3 = Inversions::maj_flat_seven_4_part[0][1];
							pitchShift4 = Inversions::maj_flat_seven_4_part[0][2];
							chord_id = 8;
							break;
					case 1: pitchShift2 = Inversions::maj_flat_seven_4_part[1][0];
							pitchShift3 = Inversions::maj_flat_seven_4_part[1][1];
							pitchShift4 = Inversions::maj_flat_seven_4_part[1][2];
							chord_id = 9;
							break;
					case 2: pitchShift2 = Inversions::maj_flat_seven_4_part[2][0];
							pitchShift3 = Inversions::maj_flat_seven_4_part[2][1];
							pitchShift4 = Inversions::maj_flat_seven_4_part[2][2];
							chord_id = 10;
							break;
					case 3: pitchShift2 = Inversions::maj_flat_seven_4_part[3][0];
							pitchShift3 = Inversions::maj_flat_seven_4_part[3][1];
							pitchShift4 = Inversions::maj_flat_seven_4_part[3][2];
							chord_id = 11;
							break;
					default: break;
				}
				break;
			}
			case 3: {
				switch (inversion_mode){
					case 0: pitchShift2 = Inversions::min_flat_seven_4_part[0][0];
							pitchShift3 = Inversions::min_flat_seven_4_part[0][1];
							pitchShift4 = Inversions::min_flat_seven_4_part[0][2];
							chord_id = 12;
							break;
					case 1: pitchShift2 = Inversions::min_flat_seven_4_part[1][0];
							pitchShift3 = Inversions::min_flat_seven_4_part[1][1];
							pitchShift4 = Inversions::min_flat_seven_4_part[1][2];
							chord_id = 13;
							break;
					case 2: pitchShift2 = Inversions::min_flat_seven_4_part[2][0];
							pitchShift3 = Inversions::min_flat_seven_4_part[2][1];
							pitchShift4 = Inversions::min_flat_seven_4_part[2][2];
							chord_id = 14;
							break;
					case 3: pitchShift2 = Inversions::min_flat_seven_4_part[3][0];
							pitchShift3 = Inversions::min_flat_seven_4_part[3][1];
							pitchShift4 = Inversions::min_flat_seven_4_part[3][2];
							chord_id = 15;
							break;
					default: break;
				}
				break;
			}
			default: break;
		} // chord_mode
		if(chord_trigger) chord_trigger = 0;
		else if (inversion_trigger) inversion_trigger = 0;
		else if (voice_trigger) voice_trigger = 0;
		else if (position_trigger) position_trigger = 0;
	} // inversion_chord_voice trigger
	if(position_trigger){
		if(position_mode == 1){
			rt_printf("%d", chord_id);
			switch (chord_id){
				case 0: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 1: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 2: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 3: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;

				case 4: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 5: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 6: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 7: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;

				case 8: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 9: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						break;
				case 10: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;
				case 11: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;

				case 12: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;
				case 13: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;
				case 14: new_chord = Inversions::swap_to_alto(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;
				case 15: new_chord = Inversions::swap_to_tenor(pitchShift2, pitchShift3, pitchShift4);
						pitchShift2 = new_chord[0];
						pitchShift3 = new_chord[1];
						pitchShift4 = new_chord[2];
						 break;

				default: break;
			}
		} else {
			inversion_trigger = 1;
			rt_printf("entered inversion trigger");
		}
		rt_printf("Contents of pitchshifts, %f, %f, %f \n", pitchShift2, pitchShift3, pitchShift4); // testing
		position_trigger = 0;
	} // position_trigger
} // 4voicelogic

/** voice logic router, gets called in render every audio frame **/
void voiceLogic(){
	switch (voice_mode){
		case 1:	{
			twoVoiceLogic();
		}
		case 2: {
			threeVoiceLogic();
		}
		case 3: {
			fourVoiceLogic();
		}
		default: break;
	}
}

/** Potentiometer control for diatonic_offset **/
// active only in 2 voice state
void set_diatonic_offset(int n, BelaContext* context){
	if (!(n % gAudioFramesPerAnalogFrame)){
		analogPot = analogRead(context, 0, 0);
		if(analogPot < 0.277 && diatonic_offset!= -1){
			diatonic_offset = -1;
			diatonic_offset_trigger = 1;
		} else if (analogPot >= 0.277 && analogPot < 0.554 && diatonic_offset != 0){
			diatonic_offset = 0;
			diatonic_offset_trigger = 1;
		} else if (analogPot >= 0.554 && diatonic_offset != 1){
			diatonic_offset = 1;
			diatonic_offset_trigger = 1;
		} // 3 pot states
	} // analog frames
} // function

/** mode selector, state changer with 1 second timer allowing for 256 user friendly state combinations **/
void selectListener(){
	// if button was pressed
	if(gTriggerButton1){
		rt_printf("Change voice\n");
		// start timer
		gButton1Timer = 0;
		// if timer is not running
		if(gButton1Timer == -1){

			gButton1Clicks++;
		// if timer is running
		} else {
			gButton1Clicks++;
			if(gButton1Clicks > 4)
				gButton1Clicks = 1;
		}
		gTriggerButton1 = 0;
	}
	if(gTriggerButton2){
		rt_printf("Change position\n");
		// start timer
		gButton2Timer = 0;
		// if timer is not running
		if(gButton2Timer == -1){

			gButton2Clicks++;
		// if timer is running
		} else {
			gButton2Clicks++;
			if(gButton2Clicks > 4)
				gButton2Clicks = 1;
		}
		gTriggerButton2 = 0;
	}
	if(gTriggerButton3){
		rt_printf("Change chord\n");
		// start timer
		gButton3Timer = 0;
		// if timer is not running
		if(gButton3Timer == -1){

			gButton3Clicks++;
		// if timer is running
		} else {
			gButton3Clicks++;
			if(gButton3Clicks > 4)
				gButton3Clicks = 1;
		}
		gTriggerButton3 = 0;
	}
	if(gTriggerButton4){
		rt_printf("Change inversion\n");
		// start timer
		gButton4Timer = 0;
		// if timer is not running
		if(gButton4Timer == -1){

			gButton4Clicks++;
		// if timer is running
		} else {
			gButton4Clicks++;
			if(gButton4Clicks > 4)
				gButton4Clicks = 1;
		}
		gTriggerButton4 = 0;
	}
}

// cleanup() is called once at the end, after the audio has stopped.

void cleanup(BelaContext *context, void *userData)
{
	NE10_FREE(timeDomainIn);
	NE10_FREE(timeDomainOut);
	NE10_FREE(timeDomainOut2);
	NE10_FREE(timeDomainOut3);
	NE10_FREE(frequencyDomain);
	NE10_FREE(frequencyDomain2);
	NE10_FREE(frequencyDomain3);

	NE10_FREE(cfg);

	free(new_chord); // mallocced for array returned by voice swappers
}
