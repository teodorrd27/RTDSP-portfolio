#include "TestBFormat.h"

TestBFormat::TestBFormat(int convolution_length, ne10_fft_cpx_float32_t* wxyz_left_fd[4], ne10_fft_cpx_float32_t* wxyz_right_fd[4]){
  this->convolution_length = convolution_length;
  half_convolution_length = convolution_length / 2;
  //white_noise = new WAVHandler("white_noise/white_noise", 1);

  initialise_buffers();

  for(int i = 0; i < 4; i++){
    memcpy(wxyz_left_hrtf_fd[i], wxyz_left_fd[i], convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memcpy(wxyz_right_hrtf_fd[i], wxyz_right_fd[i], convolution_length * sizeof(ne10_fft_cpx_float32_t));
  }
  //white_noise->spool_mono_into_buffer(impulse_signal_td);

  impulse_bformat.Configure(1, true, convolution_length);
  myEncoder.Configure(1, true, 0);
  position.fDistance = 5;

  process_BFormat();
  export_csv();
}


void TestBFormat::initialise_buffers(){
  impulse_signal_td = (float*)malloc(convolution_length * sizeof(float));
  memset(impulse_signal_td, 0, convolution_length * sizeof(float));
  for(int i = 0; i < 1024; i++){
    impulse_signal_td[i] = sin(2 * M_PI * 464.4 * i / 44100);
  }// impulse signal is 1 at location 0, 0 otherwise
  //impulse_signal_td[0] = 1.0;


  for(int i = 0; i < 4; i++){
    wxyz_left_hrtf_fd[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(wxyz_left_hrtf_fd[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));
    wxyz_right_hrtf_fd[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(wxyz_right_hrtf_fd[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    wxyz_impulse_td[i] = (float*)malloc(convolution_length * sizeof(float));
    memset(wxyz_impulse_td[i], 0, convolution_length * sizeof(float));

    timeDomainIn[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(timeDomainIn[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    frequencyDomain[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(frequencyDomain[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    frequencyDomainInterimBuffers[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(frequencyDomainInterimBuffers[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    timeDomainOutL[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(timeDomainOutL[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    timeDomainOutR[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(timeDomainOutR[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));
  }

  for(int i = 0; i < NUMBER_OF_HRTFS; i++){
    outputL[i] = (float*)malloc(convolution_length * sizeof(float));
    memset(outputL[i], 0, convolution_length * sizeof(float));

    outputR[i] = (float*)malloc(convolution_length * sizeof(float));
    memset(outputR[i], 0, convolution_length * sizeof(float));
  }

  cfg = ne10_fft_alloc_c2c_float32_neon (convolution_length);

}

void TestBFormat::process_BFormat(){
  int pattern[24] = {10, 7, 8, 7, 9, 7, 8, 7, 9, 7, 8, 7, 9, 7, 8, 7, 9, 7, 8, 7, 9, 7, 8, 7};
  int azimuths[NUMBER_OF_HRTFS];
  int elevations[NUMBER_OF_HRTFS];
  int index = 0;
  int current_pattern;
  int pattern_index = 0;
  for(int azim = 0; azim >= -345; azim -= 15){
    current_pattern = pattern[pattern_index];
    //printf("current_pattern %d", current_pattern);
    if(current_pattern == 10){
        for(int elev = -45; elev <= 90; elev += 15){
          azimuths[index] = azim;
          elevations[index] = elev;
          printf("index: %d, azim: %d, elev: %d", index, azim, elev);
          index++;
        }
    } else if(current_pattern == 9){
        for(int elev = -45; elev <= 75; elev += 15){
          azimuths[index] = azim;
          elevations[index] = elev;
          printf("index: %d, azim: %d, elev: %d", index, azim, elev);
          index++;
        }
    } else if(current_pattern == 8){
        for(int elev = -45; elev <= 60; elev += 15){
          azimuths[index] = azim;
          elevations[index] = elev;
          printf("index: %d, azim: %d, elev: %d", index, azim, elev);
          index++;
        }
    } else if(current_pattern == 7){
        for(int elev = -45; elev <= 45; elev += 15){
          azimuths[index] = azim;
          elevations[index] = elev;
          printf("index: %d, azim: %d, elev: %d", index, azim, elev);
          index++;
        }
    }
    pattern_index++;
  }
  for(int pos = 0; pos < NUMBER_OF_HRTFS; pos++){
    //printf("Position: azim %d, elev %d\n", azimuths[pos], elevations[pos]);

    position.fAzimuth = DegreesToRadians(azimuths[pos]);
    position.fElevation = DegreesToRadians(elevations[pos]);
    //position.fAzimuth = DegreesToRadians(azimuths[pos]);
    myEncoder.SetPosition(position);
    myEncoder.Refresh();

    myEncoder.Process(impulse_signal_td, convolution_length, &impulse_bformat);
    for(int i = 0; i < 4; i++){
      impulse_bformat.ExtractStream(wxyz_impulse_td[i], i, convolution_length);
    }
    /*for(int n = 0; n < convolution_length; n++)
      wxyz_impulse_td[0][n] *= 0.70710678118;*/
    for(int i = 0; i < 4; i++){
      for(int n = 0; n < convolution_length; n++){
        timeDomainIn[i][n].r = (ne10_float32_t) wxyz_impulse_td[i][n];
        timeDomainIn[i][n].i = 0.0;
      }
      // perform FFT
      ne10_fft_c2c_1d_float32_neon (frequencyDomain[i], timeDomainIn[i], cfg, 0);
      // convolution left
      for(int n = 0; n < convolution_length; n++){
        frequencyDomainInterimBuffers[i][n].r = frequencyDomain[i][n].r * wxyz_left_hrtf_fd[i][n].r - frequencyDomain[i][n].i * wxyz_left_hrtf_fd[i][n].i;
  			frequencyDomainInterimBuffers[i][n].i = frequencyDomain[i][n].i * wxyz_left_hrtf_fd[i][n].r + frequencyDomain[i][n].r * wxyz_left_hrtf_fd[i][n].i;
      }
      // perform IFFT left
      ne10_fft_c2c_1d_float32_neon (timeDomainOutL[i], frequencyDomainInterimBuffers[i], cfg, 1);
      // convolution right
      for(int n = 0; n < convolution_length; n++){
        frequencyDomainInterimBuffers[i][n].r = frequencyDomain[i][n].r * wxyz_right_hrtf_fd[i][n].r - frequencyDomain[i][n].i * wxyz_right_hrtf_fd[i][n].i;
  			frequencyDomainInterimBuffers[i][n].i = frequencyDomain[i][n].i * wxyz_right_hrtf_fd[i][n].r + frequencyDomain[i][n].r * wxyz_right_hrtf_fd[i][n].i;
      }
      // perform IFFT right
      ne10_fft_c2c_1d_float32_neon (timeDomainOutR[i], frequencyDomainInterimBuffers[i], cfg, 1);

      for(int n = 0; n < convolution_length; n++){
        outputL[pos][n] += timeDomainOutL[i][n].r;
        outputR[pos][n] += timeDomainOutR[i][n].r;
      }
    }
    for(int n = 0; n < convolution_length; n++){
      outputL[pos][n] /= 8;
      outputR[pos][n] /= 8;
    }
  }
}

void TestBFormat::export_csv(){

  std::ofstream outL("tests/testLamb_IRC1013_impulse.csv");
  std::ofstream outR("tests/testRamb_IRC1013_impulse.csv");

  for(int i = 0; i < convolution_length; i++){
  	outL << (float)i << ',';
  	outR << (float)i << ',';
  }
  outL << '\n';
  outR << '\n';
  for(int i = 0; i < NUMBER_OF_HRTFS; i++){
    for(int j = 0; j < convolution_length; j++){
      outL << (float)outputL[i][j] << ',';
    }
    outL << '\n';
  }
  for(int i = 0; i < NUMBER_OF_HRTFS; i++){
    for(int j = 0; j < convolution_length; j++){
      outR << (float)outputR[i][j] << ',';
    }
    outR << '\n';
  }
}


TestBFormat::~TestBFormat(){
  free(impulse_signal_td);
  for(int i = 0; i < 4; i++){
    NE10_FREE(wxyz_left_hrtf_fd[i]);

    free(wxyz_impulse_td[i]);

    NE10_FREE(timeDomainIn[i]);
    NE10_FREE(frequencyDomain[i]);
    NE10_FREE(frequencyDomainInterimBuffers[i]);
    NE10_FREE(timeDomainOutL[i]);
    NE10_FREE(timeDomainOutR[i]);
  }
  for(int i = 0; i < NUMBER_OF_HRTFS; i++){
    free(outputL[i]);
    free(outputR[i]);
  }
}
