//
// File: AHV_Model_types.h
//
// Code generated for Simulink model 'AHV_Model'.
//
// Model version                  : 1.5458
// Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
// C/C++ source code generated on : Tue Mar  8 16:53:12 2022
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Windows64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef RTW_HEADER_AHV_Model_types_h_
#define RTW_HEADER_AHV_Model_types_h_
#include "rtwtypes.h"

// Model Code Variants
#ifndef DEFINED_TYPEDEF_FOR_Wave_OceanCurrent_
#define DEFINED_TYPEDEF_FOR_Wave_OceanCurrent_

typedef struct {
  real_T spectrum_type;
  real_T hs;
  real_T omega_peak;
  real_T psi_mean;
  real_T gamma;
  real_T spread;
  real_T depth;
  real_T nfreq;
  real_T ndir;
  real_T energylim;
  real_T freq_cutoff;
  real_T dir_cutoff;
  real_T rand_freq;
  real_T rand_dir;
  real_T Current_direction;
  real_T Current_speed;
} Wave_OceanCurrent;

#endif

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T

#ifndef typedef_emxArray_real_T_AHV_Model_T
#define typedef_emxArray_real_T_AHV_Model_T

typedef struct emxArray_real_T emxArray_real_T_AHV_Model_T;

#endif                                 //typedef_emxArray_real_T_AHV_Model_T
#endif                                 // RTW_HEADER_AHV_Model_types_h_

//
// File trailer for generated code.
//
// [EOF]
//
