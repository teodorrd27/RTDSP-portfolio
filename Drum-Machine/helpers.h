/***** button.h *****/
#ifndef HELPERS_H
#define HELPERS_H

	void checkButtons(BelaContext* context);
	void flashLED(BelaContext* context);
	void calculateAxes(BelaContext* context);
	void checkAnalog(BelaContext* context);
	void setupAccel(BelaContext* context);
	void setPattern();
	void bwHPF(float input);
	bool tapListen(float input);
#endif