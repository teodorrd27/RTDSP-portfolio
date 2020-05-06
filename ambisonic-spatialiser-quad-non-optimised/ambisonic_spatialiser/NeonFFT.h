/** NeonFFT.h **/
#include <ne10/NE10.h>
#include <math.h>

#ifndef NEONFFT_H
#define NEONFFT_H

#define FFT_SIZE 1024

class NeonFFT{
  private:
    int fftSize;
    int overlap;
    int hopSize;
    float scaleFactor;
    float* windowBuffer;

    ne10_fft_cfg_float32_t config;

    void calculate_hann_window();

  public:
    ne10_fft_cpx_float32_t* timeDomainIn;
    ne10_fft_cpx_float32_t* timeDomainOut;
    ne10_fft_cpx_float32_t* frequencyDomainIn;
    ne10_fft_cpx_float32_t* frequencyDomainOut;
    // constructor
    NeonFFT(int fft_size, int overlap);

    // destructor
    ~NeonFFT();

    void to_freq_domain(ne10_fft_cpx_float32_t* outBuffer, float* samples);

    void to_time_domain(float* outBuffer, ne10_fft_cpx_float32_t* samples);

    void convolve_with_hrtf(ne10_fft_cpx_float32_t* hrtf, ne10_fft_cpx_float32_t* inBuffer);
};
#endif // NEONFFT_H
