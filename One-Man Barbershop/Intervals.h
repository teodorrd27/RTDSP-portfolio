/***** Intervals.h *****/


#ifndef INTERVALS_H_
#define INTERVALS_H_
// harmonic intervals
class Intervals {

	public:
		static const double octave_d;
		static const double maj7th_d;
		static const double min7th_d;
		static const double maj6th_d;
		static const double min6th_d;
		static const double per5th_d;
		static const double aug4th_d;
		static const double per4th_d;
		static const double maj3rd_d;
		static const double min3rd_d;
		static const double maj2nd_d;
		static const double min2nd_d;

		static const double min2nd_u;
		static const double maj2nd_u;
		static const double min3rd_u;
		static const double maj3rd_u;
		static const double per4th_u;
		static const double aug4th_u;
		static const double per5th_u;
		static const double min6th_u;
		static const double maj6th_u;
		static const double min7th_u;
		static const double maj7th_u;
		static const double octave_u;

		// this array provides support for diatonic offseting
		static const double up[12];
		static const double down[12];
};


#endif //INTERVALS_H_
