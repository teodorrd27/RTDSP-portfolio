/***** WAVHandler.cpp *****/
#include "WAVHandler.h"

/** this class handles the WAV samples by calling the wav_loader function **/

    void WAVHandler::deinterleave_buffers(){
      for(int i = 0; i < number_of_wavs; i++){
        left_wavs[i] = (float*)malloc(wav_buffer_lengths[i] / 2 * sizeof(float));
        right_wavs[i] = (float*)malloc(wav_buffer_lengths[i] / 2 * sizeof(float));

        for(int j = 0; j < wav_length; j++){
          left_wavs[i][j] = wav_sample_buffers[i][j * 2];
          right_wavs[i][j] = wav_sample_buffers[i][j * 2 + 1];
        }
      }
    }

    // returns 0 when the wav is single channel and 1 if it is stereo
    int WAVHandler::initWAVS(){
      if(WavLoader::load_wav(wavPrefix, wav_sample_buffers, wav_buffer_lengths, number_of_wavs, &wav_length) == -1){
        return 0;
      }
      else{
      	
      	rt_printf("%d wav length is ", wav_length);
        return 1;
      }
    }

    // constructor
    WAVHandler::WAVHandler(string wav_prefix, int number_of_wavs){
      this->number_of_wavs = number_of_wavs;
      wav_buffer_lengths = (int*)malloc(number_of_wavs * sizeof(int));
      wav_sample_buffers = (float**)malloc(number_of_wavs * sizeof(float*));
      left_wavs = (float**)malloc(number_of_wavs * sizeof(float*));
      right_wavs = (float**)malloc(number_of_wavs * sizeof(float*));

      wavPrefix = wav_prefix;
      leftWavBufferPointer = 0;
      leftWavBufferIndex = 0;
      rightWavBufferPointer = 0;
      rightWavBufferIndex = 0;
      // if channel is stereo, deinterleave
      if(initWAVS())
        deinterleave_buffers();
      rt_printf("WAVHandler says: HELLO!");
    }
    // destructor
    WAVHandler::~WAVHandler(){
      WavLoader::unload_wav(wav_sample_buffers, number_of_wavs);
      WavLoader::unload_wav(left_wavs, number_of_wavs);
      WavLoader::unload_wav(right_wavs, number_of_wavs);
      free(wav_buffer_lengths);
      free(left_wavs);
      free(right_wavs);
      free(wav_sample_buffers);
      if(wav_sample_buffers[0] == NULL && left_wavs[0] == NULL && right_wavs[0] == NULL){
      	rt_printf("WAVHandler says: BYE!");
      }
    }

    float WAVHandler::feed_next_left_float(){
      leftWavBufferPointer++;
      if(leftWavBufferPointer == wav_length){
        leftWavBufferPointer = 0;
        leftWavBufferIndex++;
        if(leftWavBufferIndex == number_of_wavs)
          leftWavBufferIndex = 0;
      }
      //rt_printf(" %f ", left_wavs[leftWavBufferIndex][leftWavBufferPointer]);
      return left_wavs[leftWavBufferIndex][leftWavBufferPointer];
    }

    float WAVHandler::feed_next_right_float(){
      rightWavBufferPointer++;
      if(rightWavBufferPointer == wav_length){
        rightWavBufferPointer = 0;
        rightWavBufferIndex++;
        if(rightWavBufferIndex == number_of_wavs)
          rightWavBufferIndex = 0;
      }
      return right_wavs[rightWavBufferIndex][rightWavBufferPointer];
    }

    int WAVHandler::spool_into_buffer(float* buffer, int left_right){
      if(left_right == LEFT){
        for(int i = 0; i < wav_length; i++){
          buffer[i] = feed_next_left_float();
          //rt_printf(" %f ", buffer[i]);
        }
      } else if(left_right == RIGHT){
        for(int i = 0; i < wav_length; i++){
          buffer[i] = feed_next_right_float();
        }
      } else if(left_right == MONO){
        buffer = (float*) malloc(wav_length * sizeof(float));
        for(int i = 0; i < wav_length; i++){
          buffer[i] = wav_sample_buffers[0][i];
        }
      }
      return wav_length;
    }

    int WAVHandler::spool_mono_into_buffer(float* &buffer){
        buffer = (float*) malloc(wav_length * sizeof(float));
        for(int i = 0; i < wav_length; i++){
          buffer[i] = wav_sample_buffers[0][i];
        }
        return wav_length;
    }
