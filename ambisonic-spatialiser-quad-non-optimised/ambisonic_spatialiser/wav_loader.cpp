/***** wav_loader.cpp *****/
#include "wav_loader.h"


int WavLoader::load_wav(string wavPrefix, float* sample_buffers[], int buffer_lengths[], int how_many_wavs, int* out_wave_length) {
	/* Load drums from WAV files */
	SF_INFO sfinfo;
	SNDFILE *sndfile ;
	char filename[32];

	for(int i = 0; i < how_many_wavs; i++) {
		snprintf(filename, 32, "./%s%d.wav", wavPrefix.c_str(), i);

		if (!(sndfile = sf_open (filename, SFM_READ, &sfinfo))) {
			printf("Couldn't open file %s\n", filename);

			/* Free already loaded sounds */
			for(int j = 0; j < i; j++)
				free(sample_buffers[j]);
			return 1;
		}
		if(sfinfo.channels > 2){
			printf("File %s has more than 2 channels\n", filename);
			for(int j = 0; j < i; j++)
				free(sample_buffers[j]);
			return 1;
		}

		buffer_lengths[i] = sfinfo.frames * sfinfo.channels;
		sample_buffers[i] = (float*)malloc(buffer_lengths[i] * sizeof(float));
		if(sample_buffers[i] == NULL) {
			printf("Error: couldn't allocate buffer for %s\n", filename);

			/* Free already loaded sounds */
			for(int j = 0; j < i; j++)
				free(sample_buffers[j]);
			return 1;
		}

		int readcount = sf_readf_float(sndfile, sample_buffers[i], sfinfo.frames);

		/* Pad with zeros in case we couldn't read whole file */
		for(int k = readcount * sfinfo.channels; k < buffer_lengths[i]; k++)
			sample_buffers[i][k] = 0;

		sf_close(sndfile);
	}
	rt_printf("Channels in this wav", sfinfo.channels);
	*out_wave_length = sfinfo.frames;
	if(sfinfo.channels == 1)
		return -1;
	return 0;
}

void WavLoader::unload_wav(float* sample_buffers[], int index){
	for(int i = 0; i < index; i++){
		free(sample_buffers[i]);
		sample_buffers[i] = NULL;
	}
}
