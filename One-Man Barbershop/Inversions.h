/***** Inversions.h *****/

#ifndef INVERSIONS_H_
#define INVERSIONS_H_
// Inversions class serves the three part and four part states
class Inversions {

	public:
		// three part chords
		static double maj[4][2];
		static double min[4][2];
		static double maj_flat_seven[4][2];
		static double min_flat_seven[4][2];
		// four part chords
		static double maj_4_part[4][3];
		static double min_4_part[4][3];
		static double maj_flat_seven_4_part[4][3];
		static double min_flat_seven_4_part[4][3];
		// in the four part state, the below function enables the singer (tenor/alto)
		// to move between middle parts
		static double* swap_to_alto(double bass, double alto, double soprano);
		static double* swap_to_tenor(double bass, double tenor, double soprano);
};

#endif // inversions
