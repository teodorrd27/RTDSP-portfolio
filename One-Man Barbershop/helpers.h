/***** helpers.h *****/

#include <Bela.h>

#ifndef HELPERS_H
#define HELPERS_H
	// checks state of buttons
	void checkButtons(BelaContext* context);
	// listens for changes in state
	void selectListener();
	// directs control to function with adecquate logic for x amount of voices
	void voiceLogic();
	// in the 2 voice state, this uses the potentiometer to tune the upper voice
	// up or down a semitone in the event that the singer moves between notes of
	// the major scale
	void set_diatonic_offset(int n, BelaContext* context);
#endif // helpers_h
