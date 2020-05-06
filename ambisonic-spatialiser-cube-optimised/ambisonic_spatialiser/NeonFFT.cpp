/** NeonFFT.cpp **/
#include "NeonFFT.h"

void NeonFFT::calculate_hann_window(){
  for(int i = 0; i < fftSize; i++){
    windowBuffer[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (float)(fftSize - 1)));
  }
}

NeonFFT::NeonFFT(int fft_size, int overlap){
  this->overlap = overlap;
  fftSize = fft_size;
  hopSize = fftSize / this->overlap;

  timeDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC(fftSize * sizeof(ne10_fft_cpx_float32_t));
  memset(timeDomainIn, 0, fftSize * sizeof(ne10_fft_cpx_float32_t));

  timeDomainOut = (ne10_fft_cpx_float32_t*) NE10_MALLOC(fftSize * sizeof(ne10_fft_cpx_float32_t));
  memset(timeDomainOut, 0, fftSize * sizeof(ne10_fft_cpx_float32_t));

  frequencyDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC(fftSize * sizeof(ne10_fft_cpx_float32_t));
  memset(frequencyDomainIn, 0, fftSize * sizeof(ne10_fft_cpx_float32_t));

  frequencyDomainOut = (ne10_fft_cpx_float32_t*) NE10_MALLOC(fftSize * sizeof(ne10_fft_cpx_float32_t));
  memset(frequencyDomainOut, 0, fftSize * sizeof(ne10_fft_cpx_float32_t));

  config = ne10_fft_alloc_c2c_float32_neon(fftSize);

  windowBuffer = (float*) malloc(fftSize * sizeof(float));
  calculate_hann_window();
}

    // destructor
NeonFFT::~NeonFFT(){
  NE10_FREE(timeDomainIn);
  NE10_FREE(timeDomainOut);
  NE10_FREE(frequencyDomainIn);
  NE10_FREE(frequencyDomainOut);
  NE10_FREE(config);
}

void NeonFFT::to_freq_domain(ne10_fft_cpx_float32_t* outBuffer, float* samples){
  for(int i = 0; i < fftSize; i++){
    timeDomainIn[i].r = samples[i];
    timeDomainIn[i].i = 0;
  }
  ne10_fft_c2c_1d_float32_neon(frequencyDomainOut, timeDomainIn, config, 0);
  for(int i = 0; i < fftSize; i++){
    outBuffer[i].r = frequencyDomainOut[i].r;
    outBuffer[i].i = frequencyDomainOut[i].i;
  }
}

void NeonFFT::to_time_domain(float* outBuffer, ne10_fft_cpx_float32_t* samples){
  for(int i = 0; i < fftSize; i++){
    frequencyDomainIn[i].r = samples[i].r;
    frequencyDomainIn[i].i = samples[i].i;
  }
  ne10_fft_c2c_1d_float32_neon(timeDomainOut, frequencyDomainIn, config, 1);
  for(int i = 0; i < fftSize; i++){
    outBuffer[i] = timeDomainOut[i].r;
  }
}

void NeonFFT::convolve_with_hrtf(ne10_fft_cpx_float32_t* hrtf, ne10_fft_cpx_float32_t* inBuffer){}
/*  for(int i = 0; i < fftSize; i++){
    frequencyDomain[i].r = inBuffer[i].r * hrtf[i].r - inBuffer[i].i * hrtf[i].i;
    frequencyDomain[i].i = inBuffer[i].i * hrtf[i].r + inBuffer[i].r * hrtf[i].i;
  }
  // IFFT
  ne10_fft_c2c_1d_float32_neon(inBuffer, frequencyDomain, config, 1);
}*/
