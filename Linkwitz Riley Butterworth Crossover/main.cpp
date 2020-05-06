/*
 * First assignment for ECS732 RTDSP, to implement a 2-way audio crossover
 * using the BeagleBone Black.
 *
 * This code runs on the BeagleBone Black with the Bela/BeagleRT environment.
 */

#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <libgen.h>
#include <signal.h>
#include <getopt.h>
#include <Bela.h>

using namespace std;

// Handle Ctrl-C by requesting that the audio rendering stop
void interrupt_handler(int var)
{
	gShouldStop = true;
}

// Print usage information
void usage(const char * processName)
{
	cerr << "Usage: " << processName << " [options]" << endl;

	Bela_usage();

	cerr << "   --help [-h]:                Print this menu\n";
}

int main(int argc, char *argv[])
{
	BelaInitSettings settings;	// Standard audio settings
	//float frequency = 1000.0;	// Frequency of crossover
	struct userSettings { float frequency; bool linkwitzRiley;};
	// this struct will be passed in to Bela_initAudio()
	struct userSettings uS = {1000.0, false};

	struct option customOptions[] =
	{
		{"help", 0, NULL, 'h'},
		{"frequency", 1, NULL, 'f'},
		{"linkwitzRiley", 0, NULL, 'l'},
		{NULL, 0, NULL, 0}
	};

	// Set default settings
	Bela_defaultSettings(&settings);
	settings.setup = setup;
	settings.render = render;
	settings.cleanup = cleanup;

	// Parse command-line arguments
	while (1) {
		int c;
		if ((c = Bela_getopt_long(argc, argv, "hf:l", customOptions, &settings)) < 0){
				break;
		}
		switch (c) {
		case 'h':
				usage(basename(argv[0]));
				exit(0);
        case 'f':
        		uS.frequency = atof(optarg);
        		if(uS.frequency < 20.0)
        			uS.frequency = 20.0;
        		if(uS.frequency > 5000.0)
        			uS.frequency = 5000.0;
        		break;
        case 'l':
        		uS.linkwitzRiley = true;
        		break;
		case '?':
		default:
				usage(basename(argv[0]));
				exit(1);
		}
	}

	// Initialise the PRU audio device
	if(Bela_initAudio(&settings, &uS) != 0) {
		cout << "Error: unable to initialise audio" << endl;
		return -1;
	}



	// Start the audio device running
	if(Bela_startAudio()) {
		cout << "Error: unable to start real-time audio" << endl;
		return -1;
	}

	// Set up interrupt handler to catch Control-C
	signal(SIGINT, interrupt_handler);
	signal(SIGTERM, interrupt_handler);

	// Run until told to stop
	while(!gShouldStop) {
		usleep(100000);
	}

	// Stop the audio device
	Bela_stopAudio();

	// Clean up any resources allocated for audio
	Bela_cleanupAudio();

	// All done!
	return 0;
}
