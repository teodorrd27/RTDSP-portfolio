/***** Intervals.cpp *****/
#include "Intervals.h"
		// implementation of intervals class
		const double Intervals::octave_d = (double)1 / 2;
		const double Intervals::maj7th_d = (double)8 / 15;
		const double Intervals::min7th_d = (double)5 / 9;
		const double Intervals::maj6th_d = (double)3 / 5;
		const double Intervals::min6th_d = (double)5 / 8;
		const double Intervals::per5th_d = (double)2 / 3;
		const double Intervals::aug4th_d = (double)32 / 45;
		const double Intervals::per4th_d = (double)3 / 4;
		const double Intervals::maj3rd_d = (double)4 / 5;
		const double Intervals::min3rd_d = (double)5 / 6;
		const double Intervals::maj2nd_d = (double)8 / 9;
		const double Intervals::min2nd_d = (double)15 / 16;

		const double Intervals::min2nd_u = (double)16 / 15;
		const double Intervals::maj2nd_u = (double)9 / 8;
		const double Intervals::min3rd_u = (double)6 / 5;
		const double Intervals::maj3rd_u = (double)5 / 4;
		const double Intervals::per4th_u = (double)4 / 3;
		const double Intervals::aug4th_u = (double)45 / 32;
		const double Intervals::per5th_u = (double)3 / 2;
		const double Intervals::min6th_u = (double)8 / 5;
		const double Intervals::maj6th_u = (double)5 / 3;
		const double Intervals::min7th_u = (double)9 / 5;
		const double Intervals::maj7th_u = (double)15 / 8;
		const double Intervals::octave_u = 2.;
		// array holding intervals allows for arithmetic operations
		const double Intervals::up[] = {Intervals::min2nd_u, Intervals::maj2nd_u, Intervals::min3rd_u, Intervals::maj3rd_u, Intervals::per4th_u, Intervals::aug4th_u, Intervals::per5th_u, Intervals::min6th_u, Intervals::maj6th_u, Intervals::min7th_u, Intervals::maj7th_u, Intervals::octave_u};
		const double Intervals::down[] = {Intervals::maj7th_d, Intervals::min7th_d, Intervals::maj6th_d, Intervals::min6th_d, Intervals::per5th_d, Intervals::aug4th_d, Intervals::per4th_d, Intervals::maj3rd_d, Intervals::min3rd_d, Intervals::maj2nd_d, Intervals::min2nd_d, Intervals::octave_d};
