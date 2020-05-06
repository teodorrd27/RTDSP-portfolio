/** OperationalSampleBuffer.h **/
#include <stdlib.h>

#ifndef SAMPLEBUFFERS_H
#define SAMPLEBUFFERS_H

class OperationalSampleBuffer{

  private:

  public:
    int bufferPointer;
    int bufferLength;
    bool stereo;
    float* buffer;

    OperationalSampleBuffer(int buffer_length);

    int spool_into_operational_buffer(float* inBuffer, int inBufferLength, bool stereo);
};

#endif // SAMPLEBUFFERS_H
