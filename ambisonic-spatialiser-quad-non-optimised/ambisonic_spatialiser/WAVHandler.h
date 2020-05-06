/***** WAVHandler.h *****
  ** Implements an easy way to read in wavs from wav files (#TODO and sofa files)
  ** Uses circular buffers that are managed internally to keep track of the current
  ** samples to be fed to the caller.
*/
#include <string>

#include "wav_loader.h"

using namespace std;

#ifndef WAVHANDLER_H
#define WAVHANDLER_H

#define NUMBER_OF_WAVS 4
#define LEFT 1
#define RIGHT 2
#define MONO 0

/** this class handles the wav samples by calling the wav_loader function **/
class WAVHandler{

  private:
    string wavPrefix;

    int wav_length;
    int number_of_wavs;
    int leftWavBufferPointer;
    int rightWavBufferPointer;
    int leftWavBufferIndex;
    int rightWavBufferIndex;

    int* wav_buffer_lengths;

    float** left_wavs;
    float** right_wavs;
    float** wav_sample_buffers;

    void deinterleave_buffers();

    int initWAVS();

  public:
    // constructor
    WAVHandler(string wav_prefix, int number_of_wavs = 1);
    // destructor
    ~WAVHandler();

    /**
      *  use these to feed in next input that is currently on the wavhandler
      *  They can be used to request the next sample in the left and right
      *  circular buffers for each wav.
      */
    float feed_next_left_float();
    float feed_next_right_float();
    int spool_into_buffer(float* buffer, int left_right);
    int spool_mono_into_buffer(float* &buffer);
};


#endif // WAVHANDLER_H
