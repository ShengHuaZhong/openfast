//
// File: AHV_Model_private.h
//
// Code generated for Simulink model 'AHV_Model'.
//
// Model version                  : 1.4920
// Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
// C/C++ source code generated on : Mon Jan 24 13:06:33 2022
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Windows64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef RTW_HEADER_AHV_Model_private_h_
#define RTW_HEADER_AHV_Model_private_h_
#include "rtwtypes.h"

// Private macros used by the generated code to access rtModel
#ifndef rtmSetFirstInitCond
# define rtmSetFirstInitCond(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
# define rtmIsFirstInitCond(rtm)       ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T rt_roundd_snf(real_T u);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern void rt_mldivide_U1d6x6_U2d_4sw8yi9v(const real_T u0[36], const real_T
  u1[6], real_T y[6]);
extern uint32_T plook_u32d_evencka(real_T u, real_T bp0, real_T bpSpace,
  uint32_T maxIndex);
extern uint32_T plook_lincp(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, uint32_T *prevIndex);
extern real_T intrp3d_l_pw(const uint32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[]);
extern uint32_T linsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex);
extern uint32_T plook_evenc(real_T u, real_T bp0, real_T bpSpace, uint32_T
  maxIndex, real_T *fraction);

// private model entry point functions
extern void AHV_Model_derivatives();

#endif                                 // RTW_HEADER_AHV_Model_private_h_

//
// File trailer for generated code.
//
// [EOF]
//
