/***** wav_loader.h *****
  ** Uses the sndfile.h library to read in wavs when sample_buffers can be provided
  ** in order to store.
*/
#include <sndfile.h>
#include <cstdlib>
#include <string>

using namespace std;

#ifndef WAV_LOADER_H
#define WAV_LOADER_H


class WavLoader{
  public:
    static int load_wav(string wavPrefix, float* sample_buffers[], int buffer_lengths[], int how_many_wavs, int* out_wave_length);
    static void unload_wav(float* sample_buffers[], int index);
  private:
    WavLoader(){};
};


#endif // LOADER_H
