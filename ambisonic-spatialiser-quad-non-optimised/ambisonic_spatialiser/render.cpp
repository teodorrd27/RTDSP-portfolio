
#include <Bela.h>
#include <Scope.h>
#include <chrono>
#include <iostream>
#include <ne10/NE10.h>


#include "WAVHandler.h"
#include "NeonFFT.h"
#include "BFormat.h"
#include "AmbisonicEncoder.h"
#include "AmbisonicDecoder.h"

/*----------*/
/*----------*/
/* IMU #includes*/
#include <rtdk.h>
#include "Bela_BNO055.h"
/*----------*/
/*----------*/

using namespace std;
using namespace std::chrono;

// fft will contain half this amount of bins
#define HRTF_LENGTH 512
#define BUFFER_SIZE 44100
#define NUMBER_OF_SPEAKERS 4
#define SPEAKER_ARRAY_CONFIG kAmblib_Quad
#define NUM_STREAMS 1
// Bela scope
Scope scope;
// Feeds hrtfs when needed
WAVHandler* hrtf_handler;
WAVHandler* sound_effects;
NeonFFT* neon_fft;

high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
high_resolution_clock::time_point t3;
high_resolution_clock::time_point t4;
high_resolution_clock::time_point t5;
high_resolution_clock::time_point t6;
high_resolution_clock::time_point t7;
high_resolution_clock::time_point t8;



// stereo outputs
float outL;
float outR;

float current_phase = 0.0;
float past_phase = 0.0;
float phase_increment = 0.0;

int gFFTSize = 1024;
int convolution_length = gFFTSize;// + HRTF_LENGTH;
int half_convolution_length = convolution_length / 2;
float osamp = 2;
int gHopSize = gFFTSize / osamp;
int gSampleCount = 0;
int gEffectSampleCount = 0;
float* gWindowBuffer;
float gFFTScaleFactor;

float gSpeakerInputBuffers[NUMBER_OF_SPEAKERS][BUFFER_SIZE];
int gSpeakerInputBufferPointer;

float gSpeakerOutputBuffersL[NUMBER_OF_SPEAKERS][BUFFER_SIZE];
float gSpeakerOutputBuffersR[NUMBER_OF_SPEAKERS][BUFFER_SIZE];

int gSpeakerOutputBufferReadPointer;
int gSpeakerOutputBufferWritePointer;

ne10_fft_cpx_float32_t* timeDomainIn[NUMBER_OF_SPEAKERS];
ne10_fft_cpx_float32_t* frequencyDomain[NUMBER_OF_SPEAKERS];
ne10_fft_cpx_float32_t* frequencyDomainInterimBuffers[NUMBER_OF_SPEAKERS];
ne10_fft_cpx_float32_t* timeDomainOutL[NUMBER_OF_SPEAKERS];
ne10_fft_cpx_float32_t* timeDomainOutR[NUMBER_OF_SPEAKERS];

ne10_float32_t* cpx_hrtf_src_L[NUMBER_OF_SPEAKERS]; // complex interleaved hrtf source buffer left channel hr hi hr hi ...
ne10_float32_t* cpx_freq_src_L[NUMBER_OF_SPEAKERS]; // complex interleaved real time buffer left channel fr fi fi fr ...
ne10_float32_t* cpx_dst_L[NUMBER_OF_SPEAKERS]; // complex destination interleaved buffer left channel a b c d

ne10_float32_t* cpx_hrtf_src_R[NUMBER_OF_SPEAKERS]; // complex interleaved hrtf source buffer right channel hr hi hr hi ...
ne10_float32_t* cpx_freq_src_R[NUMBER_OF_SPEAKERS]; // complex interleaved real time buffer right channel fr fi fi fr ...
ne10_float32_t* cpx_dst_R[NUMBER_OF_SPEAKERS]; // complex destination interleaved buffer right channel a b c d


ne10_fft_cpx_float32_t* hrtf_fd_right[NUMBER_OF_SPEAKERS];//[HRTF_LENGTH];
ne10_fft_cpx_float32_t* hrtf_fd_left[NUMBER_OF_SPEAKERS];//[HRTF_LENGTH];

// configuration
ne10_fft_cfg_float32_t cfg;

int gFFTInputBufferPointer;
int gFFTOutputBufferPointer;






int counter = 0;
float* hrtf_td_right[NUMBER_OF_SPEAKERS];//[HRTF_LENGTH];
float* hrtf_td_left[NUMBER_OF_SPEAKERS];//[HRTF_LENGTH];

float* sound_effects_buffer_td;
ne10_fft_cpx_float32_t* sound_effects_buffer_fd;

float live_input_buffer[BUFFER_SIZE];
int live_input_buffer_write_pointer;
int live_input_buffer_read_pointer;

int sound_effects_buffer_length;
int sound_effects_buffer_pointer;

float* BFormat_live_input_buffer;
int BFormat_buffer_length;
int BFormat_buffer_pointer;



float** ppfSpeakerFeeds = new float*[NUMBER_OF_SPEAKERS];
float** ppfSpeakerFeedsSecondary = new float*[NUMBER_OF_SPEAKERS];
int BF_processing_toggle;
int gSpeakerFeedSampleCounter = 0;

CBFormat live_input_BFormat;
CBFormat tmp_BFormat;
CBFormat hrtf_left_BFormat;
CBFormat hrtf_right_BFormat;
PolarPoint position;
CAmbisonicEncoder myEncoder;
CAmbisonicDecoder myDecoder;



// Auxiliary task for calculating FFT -- multithreading
AuxiliaryTask gFFTTask;
AuxiliaryTask gBFormatTask;
// asynchronous function declarations
void process_fft_background(void *);
void process_BFormat_background(void *);
void process_BFormat();
void process_fft();


/*----------*/
/*----------*/
/*IMU #variables*/

// Change this to change how often the BNO055 IMU is read (in Hz)
int readInterval = 100;

I2C_BNO055 bno; // IMU sensor object
int buttonPin = P8_08; // calibration button pin
int lastButtonValue = 0; // using a pulldown resistor

// Quaternions and Vectors
imu::Quaternion gCal, gCalLeft, gCalRight, gIdleConj = {1, 0, 0, 0};
imu::Quaternion qGravIdle, qGravCal, quat, steering, qRaw;

imu::Vector<3> gRaw;
imu::Vector<3> gGravIdle, gGravCal;
imu::Vector<3> ypr; //yaw pitch and roll angles


int calibrationState = 0; // state machine variable for calibration
int setForward = 0; // flag for setting forward orientation

// variables handling threading
AuxiliaryTask i2cTask;		// Auxiliary task to read I2C
AuxiliaryTask gravityNeutralTask;		// Auxiliary task to read gravity from I2C
AuxiliaryTask gravityDownTask;		// Auxiliary task to read gravity from I2C

int readCount = 0;			// How long until we read again...
int readIntervalSamples = 0; // How many samples between reads

int printThrottle = 0; // used to limit printing frequency
/*----------*/
/*----------*/

// function declarations
void readIMU(void*);
void getNeutralGravity(void*);
void getDownGravity(void*);
void calibrate();
void resetOrientation();
void rotateVectors();
void createVectors();
void loadBFormatBuffer();
void initialiseBuffers();
int initialiseIMU(BelaContext* context);

float gInitialSourcePosition[NUM_STREAMS][3];
int gInitialSourceAzimuth[NUM_STREAMS]={-45};
int gInitialSourceElevation[NUM_STREAMS]={45};

float gUpdatedSourceAzimuth[NUM_STREAMS] = {0.0};
float gUpdatedSourceElevation[NUM_STREAMS] = {0.0};


// function to convert input stream point source locations to 3D vectors


bool setup(BelaContext *context, void *userData)
{
	if(initialiseIMU(context) == -1)
		return false;

	// generate sinewave
	float sinewave[512];
	for(int ni = 0; ni < 512; ni++)
		sinewave[ni] = (float)sin((ni / 128.f) * (M_PI * 2));

	// initialise scope
	scope.setup(8, context->audioSampleRate);

	// all the buffers and pointers
	initialiseBuffers();

	// Initialise auxiliary tasks
	if((gFFTTask = Bela_createAuxiliaryTask(&process_fft_background, 90, "fft-calculation")) == 0)
		return false;
	if((gBFormatTask = Bela_createAuxiliaryTask(&process_BFormat_background, 89, "BFormat-processing")) == 0)
		return false;

	// hrtf_handler
	hrtf_handler = new WAVHandler("hrtfs/quad_setup/quad", NUMBER_OF_SPEAKERS);
	for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		hrtf_handler->spool_into_buffer(hrtf_td_right[i], RIGHT);
		hrtf_handler->spool_into_buffer(hrtf_td_left[i], LEFT);
		for(int j = 0; j < gFFTSize; j++){
			hrtf_td_right[i][j] *= HRTF_LENGTH / 512;
			hrtf_td_left[i][j] *= HRTF_LENGTH / 512;
		}
	}

	/*for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		for(int j = 0; j < gFFTSize; j++){
			hrtf_td_left[i][gFFTSize / 4 * 3 - 1 - j] = hrtf_td_left[i][gFFTSize / 2 - 1 - j];
			hrtf_td_right[i][gFFTSize / 4 * 3 - 1 - j] = hrtf_td_right[i][gFFTSize / 2 - 1 - j];
		}
	}*/
	/*for(int i = gFFTSize / 2 - 20; i < gFFTSize / 2 + 20; i++){
		rt_printf(" %f ", hrtf_td_right[0][i]);
	}
	return false;*/ // debugging

	/*sound_effects = new WAVHandler("sound_effects/drum");
	sound_effects_buffer_length = sound_effects->spool_mono_into_buffer(sound_effects_buffer_td);
	sound_effects_buffer_pointer = 0;*/
	sound_effects_buffer_length = 512;
	sound_effects_buffer_pointer = 0;
	sound_effects_buffer_td = (float*)malloc(sound_effects_buffer_length * sizeof(float));
	for(int i = 0; i < sound_effects_buffer_length; i++){
		sound_effects_buffer_td[i] = sinewave[i];
	}
	rt_printf("there are %d samples in this sound effect buffer", sound_effects_buffer_length);

	BFormat_buffer_length = gFFTSize;
	BFormat_live_input_buffer = (float*)malloc(BFormat_buffer_length * sizeof(float));
	BFormat_buffer_pointer = 0;
	for(int i = 0; i < BFormat_buffer_length; i++){
		BFormat_live_input_buffer[i] = sound_effects_buffer_td[sound_effects_buffer_pointer];
		sound_effects_buffer_pointer++;
		if(sound_effects_buffer_pointer >= sound_effects_buffer_length)
			sound_effects_buffer_pointer = 0;
	}


	// neon_fft_hrtf
	neon_fft = new NeonFFT(convolution_length, 1);
	for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		neon_fft->to_freq_domain(hrtf_fd_right[i], hrtf_td_right[i]);
		neon_fft->to_freq_domain(hrtf_fd_left[i], hrtf_td_left[i]);
		// interleave cpx buffers
		for(int n = 0; n < half_convolution_length; n++){
			cpx_hrtf_src_L[i][n * 4 + 0] = hrtf_fd_left[i][n].r;
			cpx_hrtf_src_L[i][n * 4 + 1] = hrtf_fd_left[i][n].i;
			cpx_hrtf_src_L[i][n * 4 + 2] = hrtf_fd_left[i][n].r;
			cpx_hrtf_src_L[i][n * 4 + 3] = hrtf_fd_left[i][n].i;

			cpx_hrtf_src_R[i][n * 4 + 0] = hrtf_fd_right[i][n].r;
			cpx_hrtf_src_R[i][n * 4 + 1] = hrtf_fd_right[i][n].i;
			cpx_hrtf_src_R[i][n * 4 + 2] = hrtf_fd_right[i][n].r;
			cpx_hrtf_src_R[i][n * 4 + 3] = hrtf_fd_right[i][n].i;
		}
	}
	// test the system to see if everything was interleaved properly
	/*if(cpx_hrtf_src_R[3][1463] != hrtf_fd_right[3][365].i){
		rt_printf("sorry, %f != %f", cpx_hrtf_src_R[6][1463], hrtf_fd_right[6][365]);
		return false;
	}*/
	/*
	for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		neon_fft->to_time_domain(hrtf_td_left_secondary[i], hrtf_fd_left[i]);
		neon_fft->to_time_domain(hrtf_td_right_secondary[i], hrtf_fd_right[i]);
	}*/
	/*for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		for(int j = gFFTSize - 20; j < gFFTSize + 20; j++){
			rt_printf("%f, %f\n", hrtf_td_left[i][j], hrtf_td_right[i][j]);
			///int left = (int) (hrtf_td_left[i][j] * 10000000.f);
			int left_secondary = (int) (hrtf_td_left_secondary[i][j] * 10000000.f);
			int right = (int) (hrtf_td_right[i][j] * 10000000.f);
			int right_secondary = (int) (hrtf_td_right_secondary[i][j] * 10000000.f);

			if(left != left_secondary || right != right_secondary)
				rt_printf("speaker: %d, index: %d, value primary left: %d, value secondary left: %d, value primary right: %d, value secondary right: %d\n", i, j, left, left_secondary, right, right_secondary);
		}
	}*/


	// CBFormat as 1st order 3D,
	live_input_BFormat.Configure(1, true, BFormat_buffer_length);
	// Ambisonic encoder, also 1st order 3D
	myEncoder.Configure(1, true, 0);

	// Set test signal's position in the soundfield

	position.fAzimuth = 0;
	position.fElevation = 0;
	position.fDistance = 5;
	myEncoder.SetPosition(position);
	myEncoder.Refresh();

	// Encode test signal into BFormat buffer
	myEncoder.Process(BFormat_live_input_buffer, BFormat_buffer_length, &live_input_BFormat);

	// Ambisonic decoder, also 1st order 3D
	myDecoder.Configure(1, true, SPEAKER_ARRAY_CONFIG, NUMBER_OF_SPEAKERS);

	// Allocate buffers for speaker feeds

	for(int niSpeaker = 0; niSpeaker < NUMBER_OF_SPEAKERS; niSpeaker++){
	    ppfSpeakerFeeds[niSpeaker] = new float[BFormat_buffer_length];
	    ppfSpeakerFeedsSecondary[niSpeaker] = new float[BFormat_buffer_length];
	}
	// Decode to get the speaker feeds
	myDecoder.Process(&live_input_BFormat, BFormat_buffer_length, ppfSpeakerFeeds);
	BF_processing_toggle = 0;

	return true;
}
// FFT function
void process_fft(){//float* inBuffer, int inWritePointer, float *outBufferL, float *outBufferR, int outWritePointer, ne10_fft_cpx_float32_t* timeDomainInBle, ne10_fft_cpx_float32_t* frequencyDomainBle, ne10_fft_cpx_float32_t* frequencyDomainInterimBufferble, ne10_fft_cpx_float32_t* timeDomainOutLble, ne10_fft_cpx_float32_t* timeDomainOutRble, ne10_fft_cpx_float32_t* hrtf_fd_leftble, ne10_fft_cpx_float32_t* hrtf_fd_rightble){
	int pointer;
	int tmpInputPointer = (gFFTInputBufferPointer - gFFTSize + BUFFER_SIZE) % BUFFER_SIZE;
	int tmpOutputPointer = gFFTOutputBufferPointer;
	//t5 = high_resolution_clock::now();


	for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){

		pointer = tmpInputPointer;

		//rt_printf("%d pointer fft \n", inWritePointer);
		/*for(int i = 0; i < gFFTSize; i++){
			timeDomainIn[i].r = 0.0;
			timeDomainIn[i].i = 0.0;
		}*/
		for(int n = 0; n < gFFTSize; n++) {
			timeDomainIn[i][n].r = (ne10_float32_t) gSpeakerInputBuffers[i][pointer] * gWindowBuffer[n];
			timeDomainIn[i][n].i = 0.0;
			//scope.log(timeDomainIn[n].r);
			// circular buffer pointer
			pointer++;
			if (pointer >= BUFFER_SIZE)
				pointer = 0;
		}

		// Run the FFT

		ne10_fft_c2c_1d_float32_neon (frequencyDomain[i], timeDomainIn[i], cfg, 0);

		// convolve with hrtfs
		for(int n = 0; n < half_convolution_length; n++){
			/*tmpMulCst[0] = hrtf_fd_left[i][n].r;
			tmpMulCst[1] = hrtf_fd_left[i][n].i;
			tmpMulCst[2] = tmpMulCst[0];
			tmpMulCst[3] = tmpMulCst[1];
			tmpMulSrc[0] = frequencyDomain[i][n].r;
			tmpMulSrc[1] = frequencyDomain[i][n].i;
			tmpMulSrc[2] = tmpMulSrc[1];
			tmpMulSrc[3] = tmpMulSrc[0];
			ne10_mul_float_neon (tmpMulDst, tmpMulSrc, tmpMulCst, 4);
			frequencyDomainInterimBuffers[i][n].r = tmpMulDst[0] - tmpMulDst[1];
			frequencyDomainInterimBuffers[i][n].i = tmpMulDst[2] + tmpMulDst[3];*/

			frequencyDomainInterimBuffers[i][n].r = frequencyDomain[i][n].r * hrtf_fd_left[i][n].r - frequencyDomain[i][n].i * hrtf_fd_left[i][n].i;
			frequencyDomainInterimBuffers[i][n].i = frequencyDomain[i][n].i * hrtf_fd_left[i][n].r + frequencyDomain[i][n].r * hrtf_fd_left[i][n].i;
			/*int index = n * 4;
			cpx_freq_src_L[i][index + 0] = frequencyDomain[i][n].r;
			cpx_freq_src_L[i][index + 1] = frequencyDomain[i][n].i;
			cpx_freq_src_L[i][index + 2] = frequencyDomain[i][n].i;
			cpx_freq_src_L[i][index + 3] = frequencyDomain[i][n].r;*/
		}

		//ne10_mul_float_neon(cpx_dst_L[i], cpx_hrtf_src_L[i], cpx_freq_src_L[i], half_convolution_length);

		//ne10_mul_float_neon(cpx_dst_R[i], cpx_hrtf_src_R[i], cpx_freq_src_L[i], half_convolution_length);


		/*for(int n = 0; n < half_convolution_length; n++){
			int index = n * 4;
			frequencyDomainInterimBuffers[i][n].r = cpx_dst_L[i][index + 0] - cpx_dst_L[i][index + 1];
			frequencyDomainInterimBuffers[i][n].i = cpx_dst_L[i][index + 2] + cpx_dst_L[i][index + 3];
		}*/




		// run IFFT for left channel
		ne10_fft_c2c_1d_float32_neon (timeDomainOutL[i], frequencyDomainInterimBuffers[i], cfg, 1);

		/*for(int n = 0; n < half_convolution_length; n++){
			int index = n * 4;
			frequencyDomainInterimBuffers[i][n].r = cpx_dst_R[i][index + 0] - cpx_dst_R[i][index + 1];
			frequencyDomainInterimBuffers[i][n].i = cpx_dst_R[i][index + 2] + cpx_dst_R[i][index + 3];
		}*/
		for(int n = 0; n < half_convolution_length; n++){
			/*tmpMulCst[0] = hrtf_fd_right[i][n].r;
			tmpMulCst[1] = hrtf_fd_right[i][n].i;
			tmpMulCst[2] = tmpMulCst[0];
			tmpMulCst[3] = tmpMulCst[1];
			tmpMulSrc[0] = frequencyDomain[i][n].r;
			tmpMulSrc[1] = frequencyDomain[i][n].i;
			tmpMulSrc[2] = tmpMulSrc[1];
			tmpMulSrc[3] = tmpMulSrc[0];
			ne10_mul_float_neon (tmpMulDst, tmpMulSrc, tmpMulCst, 4);
			//frequencyDomainInterimBuffers[i][n].r = tmpMulDst[0] - tmpMulDst[1];
			//frequencyDomainInterimBuffers[i][n].i = tmpMulDst[2] + tmpMulDst[3];*/
			frequencyDomainInterimBuffers[i][n].r = frequencyDomain[i][n].r * hrtf_fd_right[i][n].r - frequencyDomain[i][n].i * hrtf_fd_right[i][n].i;
			frequencyDomainInterimBuffers[i][n].i = frequencyDomain[i][n].i * hrtf_fd_right[i][n].r + frequencyDomain[i][n].r * hrtf_fd_right[i][n].i;
		}

		// run IFFT for right channel
		//t5 = high_resolution_clock::now();

		ne10_fft_c2c_1d_float32_neon (timeDomainOutR[i], frequencyDomainInterimBuffers[i], cfg, 1);

		// write samples into main buffer
		pointer = tmpOutputPointer;
		//rt_printf("outpointer %d\n", outWritePointer);

		for(int n = 0; n < convolution_length; n++){
			gSpeakerOutputBuffersL[i][pointer] += timeDomainOutL[i][n].r * gFFTScaleFactor;
			gSpeakerOutputBuffersR[i][pointer] += timeDomainOutR[i][n].r * gFFTScaleFactor;
			//scope.log(outBufferL[pointer]);
			//scope.log(timeDomainOutR[n].r, timeDomainOutR[n].r);
			pointer++;
			if(pointer >= BUFFER_SIZE)
				pointer = 0;
		}


	}
	//t6 = high_resolution_clock::now();

	//rt_printf("subsection lasted %d microseconds %d buffer ptr\n", (int)duration_cast<microseconds>(t6 - t5).count(), gFFTInputBufferPointer);
}


void process_BFormat(){
	//rt_printf("triggered bformat");
	//t1 = high_resolution_clock::now(); // test time
	//rt_printf("azimuth is %f", position.fAzimuth);
	//rt_printf("azimuth is %f", gUpdatedSourceAzimuth[0]);

	float whatever = (float)((int)(gUpdatedSourceAzimuth[0] * 10000000) / 10000000.0);
	position.fAzimuth = whatever;//gUpdatedSourceAzimuth[0];
	whatever = (float)((int)(gUpdatedSourceElevation[0] * 10000000) / 10000000.0);
	/*if(position.fAzimuth >= 6.28)
		position.fAzimuth -= 6.28;*/
	position.fElevation = whatever;
	position.fDistance = 5;
	myEncoder.SetPosition(position);
	myEncoder.Refresh();
	loadBFormatBuffer();

	// Encode test signal into BFormat buffer

	myEncoder.Process(BFormat_live_input_buffer, BFormat_buffer_length, &live_input_BFormat);


	// toggle between speaker feed buffers so that bformat doesn't interfere with the reading into the fft
	if(BF_processing_toggle){
		//t7 = high_resolution_clock::now();
		myDecoder.Process(&live_input_BFormat, BFormat_buffer_length, ppfSpeakerFeeds);
		//t8 = high_resolution_clock::now();
		//rt_printf("bformat buffer lasted %d milliseconds \n", (int)duration_cast<milliseconds>(t8 - t7).count());
		//BF_processing_toggle = !BF_processing_toggle;
	}else if(!BF_processing_toggle){
		myDecoder.Process(&live_input_BFormat, BFormat_buffer_length, ppfSpeakerFeedsSecondary);
		//BF_processing_toggle = !BF_processing_toggle;
	}


	//t2 = high_resolution_clock::now(); // test time
	//rt_printf("bformat task lasted %d ", (int)duration_cast<milliseconds>(t2 - t1).count()); //approx 57 milliseconds to complete inside fft, 108-135 outside as thread
}

void process_fft_background(void *) {
	t3 = high_resolution_clock::now();
	//for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
	process_fft();//gSpeakerInputBuffers[i], gFFTInputBufferPointer, gSpeakerOutputBuffersL[i], gSpeakerOutputBuffersR[i], gFFTOutputBufferPointer, timeDomainIn[i], frequencyDomain[i], frequencyDomainInterimBuffers[i], timeDomainOutL[i], timeDomainOutR[i], hrtf_fd_left[i], hrtf_fd_right[i]);
	t4 = high_resolution_clock::now();
	//rt_printf("fft task lasted %d milliseconds \n", (int)duration_cast<milliseconds>(t4 - t3).count());
	rt_printf("buffer counter is : %d, duration %d milliseconds \n", gFFTInputBufferPointer, (int)duration_cast<milliseconds>(t4 - t3).count());
}

void process_BFormat_background(void *) {
	t1 = high_resolution_clock::now(); // test time
	rotateVectors();
	process_BFormat();
	t2 = high_resolution_clock::now(); // test time
	rt_printf("bformat task lasted %d milliseconds \n", (int)duration_cast<milliseconds>(t2 - t1).count()); //approx 57 milliseconds to complete inside fft, 108-135 outside as thread
}


void render(BelaContext *context, void *userData)
{

	for(unsigned int n = 0; n < context->audioFrames; n++){

		/*----------*/
	    /*----------*/
	    /*IMU #setup routine*/

	    // this schedules the imu sensor readings
	    if(++readCount >= readIntervalSamples) {
	      readCount = 0;
	      Bela_scheduleAuxiliaryTask(i2cTask);
	    }

	    // print IMU values, but not every sample
	    printThrottle++;
	    if(printThrottle >= 44100){
	      //rt_printf("Tracker Value: %d %d %d \n",gVBAPTracking[0],gVBAPTracking[1],gVBAPTracking[2]); //print horizontal head-track value
	      //rt_printf("%f %f %f\n", ypr[0], ypr[1], ypr[2]);
	      //rt_printf("Positions Update: %d %d\n",gVBAPUpdatePositions[0],gVBAPUpdatePositions[9]); //print horizontal head-track value
	      imu::Vector<3> qForward = gIdleConj.toEuler();
	      printThrottle = 0;
	    }

	    //read the value of the button
	    /*int buttonValue = digitalRead(context, 0, buttonPin);

	    // if button wasn't pressed before and is pressed now
	    if( buttonValue != lastButtonValue && buttonValue == 1 ){
	      // then run calibration to set looking forward (gGravIdle)
	      // and looking down (gGravCal)
	      switch(calibrationState) {
	      case 0: // first time button was pressed
	        setForward = 1;
	        // run task to get gravity values when sensor in neutral position
	        Bela_scheduleAuxiliaryTask(gravityNeutralTask);
	        calibrationState = 1;	// progress calibration state
	        break;
	      case 1: // second time button was pressed
	        // run task to get gravity values when sensor 'looking down' (for head-tracking)
	        Bela_scheduleAuxiliaryTask(gravityDownTask);
	        calibrationState = 0; // reset calibration state for next time
	        break;
	      }
	    }
	    lastButtonValue = buttonValue;*/
	    /*----------*/
	    /*----------*/


		// clear out buffers
		for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
			gSpeakerInputBuffers[i][gSpeakerInputBufferPointer] = 0.0;



		// take input
		for(unsigned int channel = 0; channel < context->audioInChannels; channel++){
			live_input_buffer[live_input_buffer_write_pointer] = context->audioIn[n * context->audioInChannels + channel];
		}
		//live_input_buffer[live_input_buffer_write_pointer] /= 8;

		//scope.log(gSpeakerInputBuffers[0][tmp_pointer], gSpeakerInputBuffers[1][tmp_pointer], gSpeakerInputBuffers[2][tmp_pointer], gSpeakerInputBuffers[3][tmp_pointer]);

		if(!BF_processing_toggle){
			for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
				gSpeakerInputBuffers[i][gSpeakerInputBufferPointer] += ppfSpeakerFeeds[i][gSpeakerFeedSampleCounter];//live_input_buffer[live_input_buffer_write_pointer];//
			//scope.log(ppfSpeakerFeeds[0][gSpeakerFeedSampleCounter], ppfSpeakerFeeds[1][gSpeakerFeedSampleCounter], ppfSpeakerFeeds[2][gSpeakerFeedSampleCounter], ppfSpeakerFeeds[3][gSpeakerFeedSampleCounter]);

		} else if (BF_processing_toggle){
			for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
				gSpeakerInputBuffers[i][gSpeakerInputBufferPointer] += ppfSpeakerFeedsSecondary[i][gSpeakerFeedSampleCounter]; //live_input_buffer[live_input_buffer_write_pointer];//
			//scope.log(ppfSpeakerFeedsSecondary[0][gSpeakerFeedSampleCounter], ppfSpeakerFeedsSecondary[1][gSpeakerFeedSampleCounter], ppfSpeakerFeedsSecondary[2][gSpeakerFeedSampleCounter], ppfSpeakerFeedsSecondary[3][gSpeakerFeedSampleCounter]);
		}


		for (unsigned int channel = 0; channel < context->audioOutChannels; channel++){
			if(channel == 0){
				for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
					context->audioOut[n * context->audioOutChannels + channel] += gSpeakerOutputBuffersL[i][gSpeakerOutputBufferReadPointer];
				//scope.log(gSpeakerOutputBuffersL[0][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[1][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[2][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[3][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[4][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[5][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[6][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersL[7][gSpeakerOutputBufferReadPointer]);
				outL = context->audioOut[n * context->audioOutChannels + channel];
				context->audioOut[n * context->audioOutChannels + channel] /= 4;

			}
			else{
				for(int i = 0; i < NUMBER_OF_SPEAKERS; i++)
					context->audioOut[n * context->audioOutChannels + channel] += gSpeakerOutputBuffersR[i][gSpeakerOutputBufferReadPointer];
				//scope.log(gSpeakerOutputBuffersR[0][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[1][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[2][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[3][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[4][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[5][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[6][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[7][gSpeakerOutputBufferReadPointer]);
				outR = context->audioOut[n * context->audioOutChannels + channel];

				context->audioOut[n * context->audioOutChannels + channel] /= 4;
			}
		}
		//scope.log(gSpeakerOutputBuffersR[0][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[1][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[2][gSpeakerOutputBufferReadPointer], gSpeakerOutputBuffersR[3][gSpeakerOutputBufferReadPointer]);

		// clear out output buffers
		for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
			gSpeakerOutputBuffersL[i][gSpeakerOutputBufferReadPointer] = 0;
			gSpeakerOutputBuffersR[i][gSpeakerOutputBufferReadPointer] = 0;
		}

		// advance output buffer READ pointers
		gSpeakerOutputBufferReadPointer++;
		if(gSpeakerOutputBufferReadPointer >= BUFFER_SIZE)
			gSpeakerOutputBufferReadPointer = 0;
		// advance output buffer WRITE pointers
		gSpeakerOutputBufferWritePointer++;
		if(gSpeakerOutputBufferWritePointer >= BUFFER_SIZE)
			gSpeakerOutputBufferWritePointer = 0;
		// advance input buffer pointers
		gSpeakerInputBufferPointer++;
		if(gSpeakerInputBufferPointer >= BUFFER_SIZE)
			gSpeakerInputBufferPointer = 0;

		//if(gSpeakerFeedSampleCounter == 0)
		//	Bela_scheduleAuxiliaryTask(gBFormatTask);
		gSpeakerFeedSampleCounter++;
		if(gSpeakerFeedSampleCounter >= BFormat_buffer_length){
			//rt_printf("buffer pointer %d\n", gSpeakerInputBufferPointer);
			Bela_scheduleAuxiliaryTask(gBFormatTask);
			BF_processing_toggle = !BF_processing_toggle;
			gSpeakerFeedSampleCounter = 0;
		}


		//scope.log(gSpeakerInputBuffers[0][gSpeakerInputBufferPointer], gSpeakerInputBuffers[1][gSpeakerInputBufferPointer], gSpeakerInputBuffers[2][gSpeakerInputBufferPointer], gSpeakerInputBuffers[3][gSpeakerInputBufferPointer]);
		gSampleCount++;
		if(gSampleCount >= gHopSize){
			//rt_printf("buffer pointer %d \n", gSpeakerInputBufferPointer);
			gFFTInputBufferPointer = gSpeakerInputBufferPointer;
			gFFTOutputBufferPointer = gSpeakerOutputBufferWritePointer;
			Bela_scheduleAuxiliaryTask(gFFTTask);
			gSampleCount = 0;
		}


		live_input_buffer_write_pointer++;
		if(live_input_buffer_write_pointer >= BUFFER_SIZE)
			live_input_buffer_write_pointer = 0;

	}
}

// Auxiliary task to read from the I2C board
void readIMU(void*)
{
	// get calibration status
	uint8_t sys, gyro, accel, mag;
	bno.getCalibration(&sys, &gyro, &accel, &mag);
	// status of 3 means fully calibrated
	//rt_printf("CALIBRATION STATUSES\n");
	//rt_printf("System: %d   Gyro: %d Accel: %d  Mag: %d\n", sys, gyro, accel, mag);

	// quaternion data routine from MrHeadTracker
  	imu::Quaternion qRaw = bno.getQuat(); //get sensor raw quaternion data

  	if( setForward ) {
  		gIdleConj = qRaw.conjugate(); // sets what is looking forward
  		setForward = 0; // reset flag so only happens once
  	}

  	steering = gIdleConj * qRaw; // calculate relative rotation data
  	quat = gCalLeft * steering; // transform it to calibrated coordinate system
  	quat = quat * gCalRight;

  	ypr = quat.toEuler(); // transform from quaternion to Euler
}

// Auxiliary task to read from the I2C board
void getNeutralGravity(void*) {
	// read in gravity value
  	imu::Vector<3> gravity = bno.getVector(I2C_BNO055::VECTOR_GRAVITY);
  	gravity = gravity.scale(-1);
  	gravity.normalize();
  	gGravIdle = gravity;
}

// Auxiliary task to read from the I2C board
void getDownGravity(void*) {
	// read in gravity value
  	imu::Vector<3> gravity = bno.getVector(I2C_BNO055::VECTOR_GRAVITY);
  	gravity = gravity.scale(-1);
  	gravity.normalize();
  	gGravCal = gravity;
  	// run calibration routine as we should have both gravity values
  	calibrate();
}

// calibration of coordinate system from MrHeadTracker
// see http://www.aes.org/e-lib/browse.cfm?elib=18567 for full paper
// describing algorithm
void calibrate() {
  	imu::Vector<3> g, gravCalTemp, x, y, z;
  	g = gGravIdle; // looking forward in neutral position

  	z = g.scale(-1);
  	z.normalize();

  	gravCalTemp = gGravCal; // looking down
  	y = gravCalTemp.cross(g);
  	y.normalize();

  	x = y.cross(z);
  	x.normalize();

  	imu::Matrix<3> rot;
  	rot.cell(0, 0) = x.x();
  	rot.cell(1, 0) = x.y();
  	rot.cell(2, 0) = x.z();
  	rot.cell(0, 1) = y.x();
  	rot.cell(1, 1) = y.y();
  	rot.cell(2, 1) = y.z();
  	rot.cell(0, 2) = z.x();
  	rot.cell(1, 2) = z.y();
  	rot.cell(2, 2) = z.z();

  	gCal.fromMatrix(rot);

  	resetOrientation();
}

// from MrHeadTracker
// resets values used for looking forward
void resetOrientation() {
  	gCalLeft = gCal.conjugate();
  	gCalRight = gCal;
}

void cleanup(BelaContext *context, void *userData)
{
	delete hrtf_handler;
	// De-allocate speaker feed buffers
	for(int niSpeaker = 0; niSpeaker < NUMBER_OF_SPEAKERS; niSpeaker++){
			delete [] ppfSpeakerFeeds[niSpeaker];
			delete [] ppfSpeakerFeedsSecondary[niSpeaker];
			NE10_FREE(timeDomainIn[niSpeaker]);
			NE10_FREE(frequencyDomain[niSpeaker]);
			NE10_FREE(frequencyDomainInterimBuffers[niSpeaker]);
			NE10_FREE(timeDomainOutL[niSpeaker]);
			NE10_FREE(timeDomainOutR[niSpeaker]);
	}
	delete [] ppfSpeakerFeeds;
	delete [] ppfSpeakerFeedsSecondary;
	delete neon_fft;
	free(sound_effects_buffer_td);
	free(BFormat_live_input_buffer);
	free(gWindowBuffer);
}

void loadBFormatBuffer(){
	/*for(int i = 0; i < BFormat_buffer_length; i++){
		BFormat_buffer[i] = sound_effects_buffer_td[sound_effects_buffer_pointer];
		sound_effects_buffer_pointer++;
		if(sound_effects_buffer_pointer >= sound_effects_buffer_length)
			sound_effects_buffer_pointer = 0;
	}*/
	//t7 = high_resolution_clock::now();
	for(int i = 0; i < BFormat_buffer_length; i++){
		BFormat_live_input_buffer[i] = live_input_buffer[live_input_buffer_read_pointer];
		//scope.log(BFormat_buffer[i]);
		live_input_buffer_read_pointer++;
		if(live_input_buffer_read_pointer >= BUFFER_SIZE)
			live_input_buffer_read_pointer = 0;
	}
	//t8 = high_resolution_clock::now();
	//rt_printf("loading bformat buffer lasted %d milliseconds \n", (int)duration_cast<milliseconds>(t8 - t7).count());

}
void initialiseBuffers(){

	// Allocate the window buffer based on the FFT size
	gWindowBuffer = (float *)malloc(gFFTSize * sizeof(float));
	// Calculate a Hann window
	for(int n = 0; n < gFFTSize; n++) {
		gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0f * M_PI * n / (float)(gFFTSize - 1)));
	}
	gFFTScaleFactor = 1.0 / (float)gFFTSize * 1000;

	// speaker amount dependent buffers
	for(int i = 0; i < NUMBER_OF_SPEAKERS; i++){
		// buffer size dependent buffers
		for(int j = 0; j < BUFFER_SIZE; j++){
			gSpeakerInputBuffers[i][j] = 0;
			gSpeakerOutputBuffersL[i][j] = 0;
			gSpeakerOutputBuffersR[i][j] = 0;
		}
		/*for(int j = 0; j < HRTF_LENGTH; j++){
			hrtf_td_right[i][j] = 0;
			hrtf_td_left[i][j] = 0;
			hrtf_fd_right[i][j].r = 0;
			hrtf_fd_right[i][j].i = 0;
			hrtf_fd_left[i][j].r = 0;
			hrtf_fd_left[i][j].i = 0;
		}*/

		timeDomainIn[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		frequencyDomain[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		frequencyDomainInterimBuffers[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		timeDomainOutL[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		timeDomainOutR[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		hrtf_fd_left[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		hrtf_fd_right[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC((convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		hrtf_td_left[i] = (float*)NE10_MALLOC((convolution_length) * sizeof(float));
		hrtf_td_right[i] = (float*)NE10_MALLOC((convolution_length) * sizeof(float));

		memset(timeDomainIn[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(frequencyDomain[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(frequencyDomainInterimBuffers[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(timeDomainOutL[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(timeDomainOutR[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(hrtf_fd_left[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(hrtf_fd_right[i], 0, (convolution_length) * sizeof(ne10_fft_cpx_float32_t));
		memset(hrtf_td_left[i], 0, (convolution_length) * sizeof(float));
		memset(hrtf_td_right[i], 0, (convolution_length) * sizeof(float));

		// complex neon multiply buffers
		cpx_hrtf_src_L[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));
		cpx_freq_src_L[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));
		cpx_dst_L[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));

		cpx_hrtf_src_R[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));
		cpx_freq_src_R[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));
		cpx_dst_R[i] = (ne10_float32_t*)NE10_MALLOC(4 * half_convolution_length * sizeof(ne10_float32_t));

		memset(cpx_hrtf_src_L[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));
		memset(cpx_freq_src_L[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));
		memset(cpx_dst_L[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));

		memset(cpx_hrtf_src_R[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));
		memset(cpx_freq_src_R[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));
		memset(cpx_dst_R[i], 0, 4 * half_convolution_length * sizeof(ne10_float32_t));
	}

	// buffer size length dependent speaker number independent buffers
	for(int i = 0; i < BUFFER_SIZE; i++)
		live_input_buffer[i] = 0;

	gSpeakerInputBufferPointer = 0;
	gSpeakerOutputBufferReadPointer = 0;
	gSpeakerOutputBufferWritePointer = gHopSize;
	live_input_buffer_write_pointer = 0;
	live_input_buffer_read_pointer = 0;
	cfg = ne10_fft_alloc_c2c_float32_neon (convolution_length);
}

int initialiseIMU(BelaContext* context){
	if(!bno.begin()) {
		rt_printf("Error initialising BNO055\n");
		return -1;
	}

	//rt_printf("Initialised BNO055\n");

	// use external crystal for better accuracy
	bno.setExtCrystalUse(true);

	// get the system status of the sensor to make sure everything is ok
	uint8_t sysStatus, selfTest, sysError;
	bno.getSystemStatus(&sysStatus, &selfTest, &sysError);
	rt_printf("System Status: %d (0 is Idle)   Self Test: %d (15 is all good)   System Error: %d (0 is no error)\n", sysStatus, selfTest, sysError);

	// set sensor reading in a separate thread
	// so it doesn't interfere with the audio processing
	i2cTask = Bela_createAuxiliaryTask(&readIMU, 5, "bela-bno");
	readIntervalSamples = context->audioSampleRate / readInterval;

	gravityNeutralTask = Bela_createAuxiliaryTask(&getNeutralGravity, 5, "bela-neu-gravity");
	gravityDownTask = Bela_createAuxiliaryTask(&getDownGravity, 5, "bela-down-gravity");

	// set up button pin
	pinMode(context, 0, buttonPin, INPUT);

	createVectors();
	return 0;
}


// create IMU vectors
void createVectors(){
	float aziRad = 0.0;
	float eleRad = 0.0;
	for(int i = 0; i < NUM_STREAMS; i++){
		aziRad = gInitialSourceAzimuth[i] * M_PI / 180;
		eleRad = gInitialSourceElevation[i] * M_PI / 180;
		gInitialSourcePosition[i][0] = cos(eleRad) * cos(aziRad);
		gInitialSourcePosition[i][1] = cos(eleRad) * sin(aziRad);
		gInitialSourcePosition[i][2] = sin(eleRad);
		rt_printf("\nSource %d – X: %f\t Y: %f\t Z: %f\n", i, \
		gInitialSourcePosition[i][0],gInitialSourcePosition[i][1],gInitialSourcePosition[i][2]);
	}
}

// calculating IMU rotation
void rotateVectors(){
	float yawRot[3] = {0};
	float yawSin = sin(ypr[0]);
	float yawCos = cos(ypr[0]);
	float pitchRot[3]={0};
	float pitchSin = sin(-ypr[1]);
	float pitchCos = cos(-ypr[1]);
	float rollRot[3]={0};
	float rollSin = sin(ypr[2]);
	float rollCos = cos(ypr[2]);
	for(int i=0; i < NUM_STREAMS;i++){
		yawRot[0] = yawCos * gInitialSourcePosition[i][0] + -yawSin * gInitialSourcePosition[i][1];
		yawRot[1] = yawSin * gInitialSourcePosition[i][0] + yawCos * gInitialSourcePosition[i][1];
		yawRot[2] = gInitialSourcePosition[i][2];

		pitchRot[0] = pitchCos*yawRot[0] + pitchSin*yawRot[2];
		pitchRot[1] = yawRot[1];
		pitchRot[2] = -pitchSin*yawRot[0] + pitchCos*yawRot[2];

		rollRot[0] = pitchRot[0];
		rollRot[1] = rollCos * pitchRot[1] + -rollSin*pitchRot[2];
		rollRot[2] = rollSin *pitchRot[1] + rollCos*pitchRot[2];
		gUpdatedSourceAzimuth[i] = atan2(rollRot[1], rollRot[0]);
		gUpdatedSourceElevation[i] = asin(rollRot[2] / (sqrt(pow(rollRot[0],2) + pow(rollRot[1],2) + pow(rollRot[2], 2))));
	}
	//rt_printf("Azimuth %f - Elevation %f\n",gUpdatedSourceAzimuth[0], gUpdatedSourceElevation[0]);
}
