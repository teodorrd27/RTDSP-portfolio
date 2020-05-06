/** OperationalSampleBuffer.cpp **/
#include <Bela.h>
#include "OperationalSampleBuffer.h"

OperationalSampleBuffer::OperationalSampleBuffer(int buffer_length){
  bufferLength = buffer_length;
  buffer = (float*) malloc(bufferLength * sizeof(float));
  for(int i = 0; i < bufferLength; i++)
    buffer[i] = 0;
  bufferPointer = 0;
  stereo = false;
  for(int i = 0; i < bufferLength; i++){
    buffer[i] = 0;
  }
}

int OperationalSampleBuffer::spool_into_operational_buffer(float* inBuffer, int inBufferLength, bool stereo){
  if(this->stereo != stereo){
    rt_printf("The operational buffer already contains samples of %s, buffer pointer not advanced", stereo ? "stereo" : "mono");
    return bufferPointer;
  }
  for(int i = 0; i < inBufferLength; i++){
    buffer[bufferPointer] = inBuffer[i];
    bufferPointer++;
    if(bufferPointer >= bufferLength)
      bufferPointer = 0;
  }
  return bufferPointer;
}
