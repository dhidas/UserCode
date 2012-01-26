#ifndef GUARD_dead_noisy_h
#define GUARD_dead_noisy_h


#include <iostream>
#include <math.h>
#include "TH2F.h"

void         Find_Dead_And_Noisy(TH2F * map_pixels, TH2F & map_maskedOff, TH2F & map_dead);
void         Find_Dead_And_Noisy(TH2F * map_pixels, TH2F & map_maskedOff, TH2F & map_dead, TH2F & map_abs_dead);
void         Fill_loc(TH2F & hist, float z);
inline bool  is_fiducial(int ir, int ic, int ir_max, int ic_max);
float        get_stdev(float * list, float mean);
inline float get_efficiency(float pixel, float local_mean);
float        get_mean(float *list);
int          get_neighbors(TH2F* map_pixels, int ir, int ic, TH2F & Bad1, float *neighbors);
int          get_neighbors(TH2F* map_pixels, int ir, int ic, TH2F & Bad1, TH2F & Bad2, float *neighbors);
bool         is_very_dead(float pixel, float global_mean, float global_stdev, float global_efficiency);
inline bool  is_abs_dead(float pixel, float global_mean, float global_stdev, float global_efficiency);
inline bool  is_dead(float pixel, float local_mean, float local_stdev, float local_efficiency);
inline bool  is_very_noisy(float pixel, float global_mean, float global_stdev, float global_efficiency);
inline bool  is_noisy(float pixel, float local_mean, float local_stdev, float local_efficiency);
//end prototypes

#endif
