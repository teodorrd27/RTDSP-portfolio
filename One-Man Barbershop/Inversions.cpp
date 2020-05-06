/***** Inversions.cpp *****/
#include "Inversions.h"
#include "Intervals.h"
#include <Bela.h>

/** PLEASE REFER TO REPORT FOR DESIGN CHOICES REGARDING THE BELOW CHORDS, INVERSIONS, AND HARMONIC VOICING **/

/** ******************************* Major chords ************************************ **/

double Inversions::maj[4][2] = {{Intervals::per5th_d, Intervals::maj6th_u}, //root // you are 5th
								{Intervals::min6th_d, Intervals::per5th_u}, // 1st inv // you are root
								{Intervals::maj6th_d, Intervals::min6th_u}, // 2nd inv // you are 3rd
								{Intervals::maj3rd_d, Intervals::min3rd_u}} // tight root // you are 3rd
							   ;
/** ********************************************************************************* **/

/** ******************************* Minor chords ************************************ **/

double Inversions::min[4][2] = {{Intervals::per5th_d, Intervals::min6th_u}, // root // you are 5th
								{Intervals::maj6th_d, Intervals::per5th_u}, // 1st inv // you are root
								{Intervals::min6th_d, Intervals::maj6th_u}, // 2nd inv // you are 3rd
								{Intervals::min3rd_d, Intervals::maj3rd_u}} // tight root // you are 3rd
							   ;
/** ********************************************************************************* **/


/** **************************** Major flat 7th chords ****************************** **/

double Inversions::maj_flat_seven[4][2] =  {{Intervals::min7th_d, Intervals::aug4th_u}, // root // you are 7th
											{Intervals::aug4th_d, Intervals::maj6th_u}, // 1st inv // you are 7th
											{Intervals::maj6th_d, Intervals::aug4th_u}, // 2nd inv // you are 3rd
											{Intervals::aug4th_d, Intervals::min6th_u}} // 3rd inv // you are 3rd
										   ;
/** ********************************************************************************* **/


/** **************************** Minor flat 7th chords ****************************** **/

double Inversions::min_flat_seven[4][2] =  {{Intervals::min7th_d, Intervals::per4th_u}, // root // you are 7th
											{Intervals::per5th_d, Intervals::maj2nd_u}, // 1st inv // you are 7th
											{Intervals::per4th_d, Intervals::min7th_u}, // 2nd inv // you are root
											{Intervals::maj6th_d, Intervals::per4th_u}} // 3rd inv // you are 5th
										   ;
/** ********************************************************************************* **/


/// indexes: UP: min2nd 0 --> octave 11
///			 DOWN: maj7th 0 --> min2nd 10 --> octave 11
/** ******************************* Major chords 4 part ************************************ **/

double Inversions::maj_4_part[4][3] = {{Intervals::per5th_d, /*5th*/ Intervals::per4th_u ,Intervals::maj6th_u}, //root // you are 5th // double the root
								{Intervals::min6th_d, Intervals::per4th_d, /*root*/ Intervals::per5th_u}, // 1st inv // you are root // double the fifth // alto
								{Intervals::maj6th_d, Intervals::maj3rd_d, /*3rd*/ Intervals::min6th_u}, // 2nd inv // you are 3rd // double root // alto
								{Intervals::maj3rd_d, /*3rd*/ Intervals::min3rd_u, Intervals::min6th_u}} // tight root // you are 3rd
							   ;
/** ********************************************************************************* **/

/** ******************************* Minor chords 4 part ************************************ **/

double Inversions::min_4_part[4][3] = {{Intervals::per5th_d, /*5th*/ Intervals::per4th_u, Intervals::min6th_u}, // root // you are 5th // double the root
								{Intervals::maj6th_d, Intervals::per4th_d, /*root*/ Intervals::per5th_u}, // 1st inv // you are root // double the fifth // alto
								{Intervals::min6th_d, Intervals::min3rd_d, /*3rd*/ Intervals::maj6th_u}, // 2nd inv // you are 3rd // double root // alto
								{Intervals::min3rd_d, /*3rd*/ Intervals::maj3rd_u, Intervals::maj6th_u}} // tight root // you are 3rd
							   ;
/** ********************************************************************************* **/


/** **************************** Major flat 7th chords 4 part ****************************** **/

double Inversions::maj_flat_seven_4_part[4][3] =  {{Intervals::min7th_d, Intervals::min3rd_d, /*7th*/ Intervals::aug4th_u}, // root // you are 7th // alto
											{Intervals::aug4th_d, /*7th*/ Intervals::maj2nd_u, Intervals::maj6th_u}, // 1st inv // you are 7th
											{Intervals::maj6th_d, Intervals::maj3rd_d, /*3rd*/ Intervals::aug4th_u}, // 2nd inv // you are 3rd // alto
											{Intervals::aug4th_d, /*3rd*/ Intervals::min3rd_u, Intervals::min6th_u}} // 3rd inv // you are 3rd
										   ;
/** ********************************************************************************* **/


/** **************************** Minor flat 7th chords 4 part ****************************** **/

double Inversions::min_flat_seven_4_part[4][3] =  {{Intervals::min7th_d, Intervals::min3rd_d, /*7th*/ Intervals::per4th_u}, // root // you are 7th // alto
											{Intervals::per5th_d, Intervals::min3rd_d, /*7th*/ Intervals::maj2nd_u}, // 1st inv // you are 7th //alto
											{Intervals::per4th_d, /*root*/ Intervals::min3rd_u, Intervals::min7th_u}, // 2nd inv // you are root //
											{Intervals::maj6th_d, Intervals::maj3rd_d, /*5th*/ Intervals::per4th_u}} // 3rd inv // you are 5th // alto
										   ;
/** ********************************************************************************* **/

/** The below functions return a pointer to an array of 3 doubles. The results
are securely stored in allocated space within the render file **/
// In the four part state, if the singer is a tenor, he swaps positions with the
// alto.
double* Inversions::swap_to_alto(double bass, double alto, double soprano){
	// since alto has to become 1, divide everything by alto.
	bass = bass / alto;
	double tenor = 1.0 / alto; // the machine will be singing a tenor
	soprano = soprano / alto;
	static double new_chord[] = {bass, tenor, soprano};
	return new_chord;
}
// Vice versa
double* Inversions::swap_to_tenor(double bass, double tenor, double soprano){
	// since tenor has to become 1, divide everything by tenor.
	bass = bass / tenor;
	double alto = 1.0 / tenor; // the machine will be singing an alto
	soprano = soprano / tenor;
	static double new_chord[] = {bass, alto, soprano};
	return new_chord;
}
