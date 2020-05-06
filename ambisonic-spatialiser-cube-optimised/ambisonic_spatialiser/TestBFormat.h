#include <ne10/NE10.h>
#include <iostream>
#include <fstream>

#include "WAVHandler.h"

#include "BFormat.h"
#include "AmbisonicEncoder.h"

#ifndef TESTBFORMAT_H
#define TESTBFORMAT_H

#define NUMBER_OF_HRTFS 187

class TestBFormat{

  public:
    CBFormat impulse_bformat;
    PolarPoint position;
    CAmbisonicEncoder myEncoder;

    int convolution_length;
    int half_convolution_length;

    float* impulse_signal_td;

    float* wxyz_impulse_td[4];
    ne10_fft_cpx_float32_t* wxyz_left_hrtf_fd[4];
    ne10_fft_cpx_float32_t* wxyz_right_hrtf_fd[4];

    ne10_fft_cpx_float32_t* timeDomainIn[4];
    ne10_fft_cpx_float32_t* frequencyDomain[4];
    ne10_fft_cpx_float32_t* frequencyDomainInterimBuffers[4];
    ne10_fft_cpx_float32_t* timeDomainOutL[4];
    ne10_fft_cpx_float32_t* timeDomainOutR[4];

    float* outputL[NUMBER_OF_HRTFS];
    float* outputR[NUMBER_OF_HRTFS];


    ne10_fft_cfg_float32_t cfg;



    TestBFormat(int convolution_length, ne10_fft_cpx_float32_t* wxyz_left_fd[4], ne10_fft_cpx_float32_t* wxyz_right_fd[4]);



  private:
    WAVHandler* white_noise;

    void initialise_buffers();
    void process_BFormat();
    void export_csv();
    ~TestBFormat();

};

#endif /** TESTBFORMAT_H **/
