#include "WAVHandler.h"
#include "NeonFFT.h"
#include <string>
#include <ne10/NE10.h>
#include <iostream>
#include <fstream>

#ifndef TESTHRTF_H
#define TESTHRTF_H

class TestHRTF{

  public:
    float** hrtf_td_left;
    float** hrtf_td_right;
    ne10_fft_cpx_float32_t** hrtf_fd_left;
    ne10_fft_cpx_float32_t** hrtf_fd_right;

    float** imp_res_td_left;
    float** imp_res_td_right;
    ne10_fft_cpx_float32_t* tmp_fd_convolution_output_left;
    ne10_fft_cpx_float32_t* tmp_fd_convolution_output_right;

    float* impulse_signal_td;
    ne10_fft_cpx_float32_t* impulse_signal_fd;

    TestHRTF(string containing_folder, int number_of_hrirs, int convolution_length);
    ~TestHRTF();

  private:
    WAVHandler* hrtfs;
    WAVHandler* white_noise;
    int convolution_length;
    int number_of_hrirs;
    NeonFFT* neon_fft;

    void initialise_buffers();
    void hrtfs_impulse_convolve();
    void td_to_fd();
    void export_csv();

};

#endif /** TESTHRTF_H **/
