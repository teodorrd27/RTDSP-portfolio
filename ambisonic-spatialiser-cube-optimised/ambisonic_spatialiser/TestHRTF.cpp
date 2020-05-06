#include "TestHRTF.h"

TestHRTF::TestHRTF(string containing_folder, int number_of_hrirs, int convolution_length){
    this->number_of_hrirs = number_of_hrirs;
    this->convolution_length = convolution_length;
    string filepath = "hrtfs/" + containing_folder + "/hrtf";
    hrtfs = new WAVHandler(filepath, number_of_hrirs);
    //white_noise = new WAVHandler("white_noise/white_noise", 1);
    neon_fft = new NeonFFT(convolution_length, 1);

    // allocate memory for all buffers in the class
    initialise_buffers();

    for(int i = 0; i < number_of_hrirs; i++){
      hrtfs->spool_into_buffer(hrtf_td_left[i], LEFT);
      hrtfs->spool_into_buffer(hrtf_td_right[i], RIGHT);
    }
    //white_noise->spool_mono_into_buffer(impulse_signal_td);
    // use the neon fft to convert all time domain buffers in this class to frequency domain
    td_to_fd();

    // convolution
    hrtfs_impulse_convolve();

    // export results to csv
    export_csv();


}

void TestHRTF::initialise_buffers(){
  hrtf_td_left = (float**)malloc(number_of_hrirs * sizeof(float*));
  hrtf_td_right = (float**)malloc(number_of_hrirs * sizeof(float*));
  hrtf_fd_left = (ne10_fft_cpx_float32_t**)NE10_MALLOC(number_of_hrirs * sizeof(ne10_fft_cpx_float32_t*));
  hrtf_fd_right = (ne10_fft_cpx_float32_t**)NE10_MALLOC(number_of_hrirs * sizeof(ne10_fft_cpx_float32_t*));

  imp_res_td_left = (float**)malloc(number_of_hrirs * sizeof(float*));
  imp_res_td_right = (float**)malloc(number_of_hrirs * sizeof(float*));

  for(int i = 0; i < number_of_hrirs; i++){
    hrtf_td_left[i] = (float*)malloc(convolution_length * sizeof(float));
    hrtf_td_right[i] = (float*)malloc(convolution_length * sizeof(float));
    hrtf_fd_left[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
    hrtf_fd_right[i] = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));

    imp_res_td_left[i] = (float*)malloc(convolution_length * sizeof(float));
    imp_res_td_right[i] = (float*)malloc(convolution_length * sizeof(float));

    memset(hrtf_td_left[i], 0, convolution_length * sizeof(float));
    memset(hrtf_td_right[i], 0, convolution_length * sizeof(float));
    memset(hrtf_fd_left[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));
    memset(hrtf_fd_right[i], 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

    memset(imp_res_td_left[i], 0, convolution_length * sizeof(float));
    memset(imp_res_td_right[i], 0, convolution_length * sizeof(float));

  }
  impulse_signal_td = (float*)malloc(convolution_length * sizeof(float));
  impulse_signal_fd = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
  memset(impulse_signal_td, 0, convolution_length * sizeof(float));
  /*for(int i = 0; i < 1024; i++){
    impulse_signal_td[i] = sin(2 * M_PI * 464.4 * i / 44100);
  }*/
  impulse_signal_td[0] = 1.0;
  memset(impulse_signal_fd, 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));

  tmp_fd_convolution_output_left = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
  memset(tmp_fd_convolution_output_left, 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));
  tmp_fd_convolution_output_right = (ne10_fft_cpx_float32_t*)NE10_MALLOC(convolution_length * sizeof(ne10_fft_cpx_float32_t));
  memset(tmp_fd_convolution_output_right, 0, convolution_length * sizeof(ne10_fft_cpx_float32_t));
}

void TestHRTF::td_to_fd(){
  neon_fft->to_freq_domain(impulse_signal_fd, impulse_signal_td);
  for(int i = 0; i < number_of_hrirs; i++){
    neon_fft->to_freq_domain(hrtf_fd_left[i], hrtf_td_left[i]);
    neon_fft->to_freq_domain(hrtf_fd_right[i], hrtf_td_right[i]);
  }
}

void TestHRTF::hrtfs_impulse_convolve(){
  for(int i = 0; i < number_of_hrirs; i++){
    for(int j = 0; j < convolution_length; j++){
      tmp_fd_convolution_output_left[j].r = hrtf_fd_left[i][j].r * impulse_signal_fd[j].r - hrtf_fd_left[i][j].i * impulse_signal_fd[j].i;
      tmp_fd_convolution_output_left[j].i = hrtf_fd_left[i][j].i * impulse_signal_fd[j].r + hrtf_fd_left[i][j].r * impulse_signal_fd[j].i;
      tmp_fd_convolution_output_right[j].r = hrtf_fd_right[i][j].r * impulse_signal_fd[j].r - hrtf_fd_right[i][j].i * impulse_signal_fd[j].i;
      tmp_fd_convolution_output_right[j].i = hrtf_fd_right[i][j].i * impulse_signal_fd[j].r + hrtf_fd_right[i][j].r * impulse_signal_fd[j].i;
    }
    neon_fft->to_time_domain(imp_res_td_left[i], tmp_fd_convolution_output_left);
    neon_fft->to_time_domain(imp_res_td_right[i], tmp_fd_convolution_output_right);
  }
}

void TestHRTF::export_csv(){

  std::ofstream outL("tests/testLhrtf_IRC1013_impulse.csv");
  std::ofstream outR("tests/testRhrtf_IRC1013_impulse.csv");

  for(int i = 0; i < convolution_length; i++){
  	outL << (float)i << ',';
  	outR << (float)i << ',';
  }
  outL << '\n';
  outR << '\n';
  for(int i = 0; i < number_of_hrirs; i++){
    for(int j = 0; j < convolution_length; j++){
      outL << (float)imp_res_td_left[i][j] << ',';
    }
    outL << '\n';
  }
  for(int i = 0; i < number_of_hrirs; i++){
    for(int j = 0; j < convolution_length; j++){
      outR << (float)imp_res_td_right[i][j] << ',';
    }
    outR << '\n';
  }
}

TestHRTF::~TestHRTF(){
  for(int i = 0; i < number_of_hrirs; i++){
    free(hrtf_td_left[i]);
    free(hrtf_td_right[i]);

    free(imp_res_td_left[i]);
    free(imp_res_td_right[i]);

    NE10_FREE(hrtf_fd_left[i]);
    NE10_FREE(hrtf_fd_right[i]);
  }
  free(hrtf_td_left);
  free(hrtf_td_right);

  free(imp_res_td_left);
  free(imp_res_td_right);

  NE10_FREE(hrtf_fd_left);
  NE10_FREE(hrtf_fd_right);
  free(impulse_signal_td);
  NE10_FREE(impulse_signal_fd);
  NE10_FREE(tmp_fd_convolution_output_left);
  NE10_FREE(tmp_fd_convolution_output_right);

  delete hrtfs;
  delete neon_fft;
}
