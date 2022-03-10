//
// File: AHV_Model.cpp
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
#include "AHV_Model.h"
#include "AHV_Model_private.h"

// Named constants for Chart: '<S66>/Chart'
const uint8_T AHV_Model_IN_hold = 1U;
const uint8_T AHV_Model_IN_unhold = 2U;

// Named constants for Chart: '<S94>/Chart1'
const uint8_T AHV_Model_IN_NO_ACTIVE_CHILD = 0U;
const uint8_T AHV_Model_IN_change_angel = 1U;
const uint8_T AHV_Model_IN_hold_angel = 2U;
const uint8_T AHV_Model_IN_local = 1U;
const uint8_T AHV_Model_IN_unlocal = 2U;

// Named constants for Chart: '<S94>/Chart3'
const uint8_T AHV_Model_IN_NO_ACTIVE_CHILD_j = 0U;
const uint8_T AHV_Model_IN_change_angel_n = 1U;
const uint8_T AHV_Model_IN_hold_angel_l = 2U;
const uint8_T AHV_Model_IN_local_f = 1U;
const uint8_T AHV_Model_IN_unlocal_h = 2U;

// Named constants for Chart: '<S190>/Chart1'
const uint8_T AHV_Model_IN_NO_ACTIVE_CHILD_jn = 0U;
const uint8_T AHV_Model_IN_change_angel_g = 1U;
const uint8_T AHV_Model_IN_hold_angel_d = 2U;
const uint8_T AHV_Model_IN_local_j = 1U;
const uint8_T AHV_Model_IN_unlocal_a = 2U;

// Exported block signals
real_T spectrum_type;                  // '<Root>/spectrum_type'
real_T hs;                             // '<Root>/hs'
real_T omega_peak;                     // '<Root>/omega_peak'
real_T psi_mean;                       // '<Root>/psi_mean'
real_T gamma_value;                          // '<Root>/gamma_value
real_T spread;                         // '<Root>/spread'
real_T depth;                          // '<Root>/depth'
real_T nfreq;                          // '<Root>/nfreq'
real_T ndir;                           // '<Root>/ndir'
real_T energylim;                      // '<Root>/energylim'
real_T freq_cutoff;                    // '<Root>/freq_cutoff'
real_T dir_cutoff;                     // '<Root>/dir_cutoff'
real_T rand_freq;                      // '<Root>/rand_freq'
real_T rand_dir;                       // '<Root>/rand_dir'
real_T Current_direction;              // '<Root>/Current_direction'
real_T Current_speed;                  // '<Root>/Current_speed'
real_T Vessel_init1[6];                // '<Root>/Vessel_init1'
boolean_T Drving_Mode1;                   // '<Root>/Drving_Mode1'
real_T heading_mode1;                  // '<Root>/heading_mode1'
boolean_T Hold_Position1;              // '<Root>/Hold_Position1'
real_T Rudde_angle1;                   // '<Root>/Rudde_ angle1'
real_T Thruster_percentage1;           // '<Root>/Thruster_percentage1'
real_T heading_angle_ref1;             // '<Root>/heading_angle_ref1'
real_T Vessel_X_Ref1;                  // '<Root>/Vessel_X_Ref1'
real_T Vessel_Y_Ref1;                  // '<Root>/Vessel_Y_Ref1'
real_T traget_speed1;                  // '<Root>/traget_speed1'
real_T tau_cable1[3];                  // '<Root>/tau_cable1'
real_T ahv_fairlead1[3];               // '<Root>/ahv_fairlead1'
real_T Vessel_init2[6];                // '<Root>/Vessel_init2'
real_T Drving_Mode2;                   // '<Root>/Drving_Mode2'
real_T heading_mode2;                  // '<Root>/heading_mode2'
boolean_T Hold_Position2;              // '<Root>/Hold_Position2'
real_T Rudde_angle2;                   // '<Root>/Rudde_ angle2'
real_T Thruster_percentage2;           // '<Root>/Thruster_percentage2'
real_T heading_angle_ref2;             // '<Root>/heading_angle_ref2'
real_T Vessel_X_Ref2;                  // '<Root>/Vessel_X_Ref2'
real_T Vessel_Y_Ref2;                  // '<Root>/Vessel_Y_Ref2'
real_T traget_speed2;                  // '<Root>/traget_speed2'
real_T tau_cable2[3];                  // '<Root>/tau_cable2'
real_T ahv_fairlead2[3];               // '<Root>/ahv_fairlead2'
real_T tau_cable3[3];                  // '<Root>/tau_cable3'
real_T ahv_fairlead3[3];               // '<Root>/ahv_fairlead3'
real_T tau_cable4[3];                  // '<Root>/tau_cable4'
real_T ahv_fairlead4[3];               // '<Root>/ahv_fairlead4'
real_T Vessel_init3[6];                // '<Root>/Vessel_init3'
real_T Drving_Mode3;                   // '<Root>/Drving_Mode3'
real_T heading_mode3;                  // '<Root>/heading_mode3'
boolean_T Hold_Position3;              // '<Root>/Hold_Position3'
real_T Rudde_angle3;                   // '<Root>/Rudde_ angle3'
real_T Thruster_percentage3;           // '<Root>/Thruster_percentage3'
real_T heading_angle_ref3;             // '<Root>/heading_angle_ref3'
real_T Vessel_X_Ref3;                  // '<Root>/Vessel_X_Ref3'
real_T Vessel_Y_Ref3;                  // '<Root>/Vessel_Y_Ref3'
real_T traget_speed3;                  // '<Root>/traget_speed3'
real_T Vessel_init4[6];                // '<Root>/Vessel_init4'
real_T Drving_Mode4;                   // '<Root>/Drving_Mode4'
real_T heading_mode4;                  // '<Root>/heading_mode4'
boolean_T Hold_Position4;              // '<Root>/Hold_Position4'
real_T Rudde_angle4;                   // '<Root>/Rudde_ angle4'
real_T Thruster_percentage4;           // '<Root>/Thruster_percentage4'
real_T heading_angle_ref4;             // '<Root>/heading_angle_ref4'
real_T Vessel_X_Ref4;                  // '<Root>/Vessel_X_Ref4'
real_T Vessel_Y_Ref4;                  // '<Root>/Vessel_Y_Ref4'
real_T traget_speed4;                  // '<Root>/traget_speed4'
real_T eta_AHV1[6];                    // '<Root>/eta_AHV1'
real_T nu1[6];                         // '<Root>/nu1'
real_T AHV_speed1;                     // '<Root>/AHV_speed1'
real_T eta_AHV2[6];                    // '<Root>/eta_AHV2'
real_T nu2[6];                         // '<Root>/nu2'
real_T AHV_speed2;                     // '<Root>/AHV_speed2'
real_T eta_AHV3[6];                    // '<Root>/eta_AHV3'
real_T nu3[6];                         // '<Root>/nu3'
real_T AHV_speed3;                     // '<Root>/AHV_speed3'
real_T eta_AHV4[6];                    // '<Root>/eta_AHV4'
real_T nu4[6];                         // '<Root>/nu4'
real_T AHV_speed4;                     // '<Root>/AHV_speed4'
uint32_T plook_u32d_evencka(real_T u, real_T bp0, real_T bpSpace, uint32_T
  maxIndex)
{
  uint32_T bpIndex;
  real_T fbpIndex;

  // Prelookup - Index only
  // Index Search method: 'even'
  // Extrapolation method: 'Clip'
  // Use previous index: 'off'
  // Use last breakpoint for index at or above upper limit: 'on'
  // Remove protection against out-of-range input in generated code: 'off'

  if (u <= bp0) {
    bpIndex = 0U;
  } else {
    fbpIndex = (u - bp0) * (1.0 / bpSpace);
    if (fbpIndex < maxIndex) {
      bpIndex = static_cast<uint32_T>(fbpIndex);
    } else {
      bpIndex = maxIndex;
    }
  }

  return bpIndex;
}

uint32_T plook_lincp(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                     *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  // Prelookup - Index and Fraction
  // Index Search method: 'linear'
  // Extrapolation method: 'Clip'
  // Use previous index: 'on'
  // Use last breakpoint for index at or above upper limit: 'off'
  // Remove protection against out-of-range input in generated code: 'off'

  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = linsearch_u32d(u, bp, *prevIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex - 1U;
    *fraction = 1.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp3d_l_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                    table[], const uint32_T stride[])
{
  real_T yL_2d;
  uint32_T offset_2d;
  real_T yL_1d;
  uint32_T offset_0d;

  // Column-major Interpolation 3-D
  // Interpolation method: 'Linear point-slope'
  // Use last breakpoint for index at or above upper limit: 'off'
  // Overflow mode: 'portable wrapping'

  offset_2d = (bpIndex[2U] * stride[2U] + bpIndex[1U] * stride[1U]) + bpIndex[0U];
  yL_1d = (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] +
    table[offset_2d];
  offset_0d = offset_2d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) - yL_1d) * frac[1U] + yL_1d;
  offset_2d += stride[2U];
  yL_1d = (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] +
    table[offset_2d];
  offset_0d = offset_2d + stride[1U];
  return (((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
             table[offset_0d]) - yL_1d) * frac[1U] + yL_1d) - yL_2d) * frac[2U]
    + yL_2d;
}

uint32_T linsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex)
{
  uint32_T bpIndex;

  // Linear Search
  for (bpIndex = startIndex; u < bp[bpIndex]; bpIndex--) {
  }

  while (u >= bp[bpIndex + 1U]) {
    bpIndex++;
  }

  return bpIndex;
}

uint32_T plook_evenc(real_T u, real_T bp0, real_T bpSpace, uint32_T maxIndex,
                     real_T *fraction)
{
  uint32_T bpIndex;
  real_T invSpc;
  real_T fbpIndex;

  // Prelookup - Index and Fraction
  // Index Search method: 'even'
  // Extrapolation method: 'Clip'
  // Use previous index: 'off'
  // Use last breakpoint for index at or above upper limit: 'off'
  // Remove protection against out-of-range input in generated code: 'off'

  if (u <= bp0) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else {
    invSpc = 1.0 / bpSpace;
    fbpIndex = (u - bp0) * invSpc;
    if (fbpIndex < maxIndex) {
      bpIndex = static_cast<uint32_T>(fbpIndex);
      *fraction = (u - (static_cast<real_T>(bpIndex) * bpSpace + bp0)) * invSpc;
    } else {
      bpIndex = maxIndex - 1U;
      *fraction = 1.0;
    }
  }

  return bpIndex;
}

//
// This function updates continuous states using the ODE4 fixed-step
// solver algorithm
//
void AH_Model_v1ModelClass::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = static_cast<ODE4_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 472;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  AHV_Model_derivatives();

  // f1 = f(t + (h/2), y + (h/2)*f0)
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  this->step();
  AHV_Model_derivatives();

  // f2 = f(t + (h/2), y + (h/2)*f1)
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  this->step();
  AHV_Model_derivatives();

  // f3 = f(t + h, y + h*f2)
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  this->step();
  AHV_Model_derivatives();

  // tnew = t + h
  // ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3)
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

//
// System reset for iterator system:
//    '<S21>/Cross-flow drag trapezoidal integration'
//    '<S117>/Cross-flow drag trapezoidal integration'
//    '<S213>/Cross-flow drag trapezoidal integration'
//    '<S309>/Cross-flow drag trapezoidal integration'
//
void AH_Model_v1ModelClass::Crossflowdragtrapezoidali_Reset(real_T
  *memory2_PreviousInput, real_T *memory1_PreviousInput)
{
  // InitializeConditions for Memory: '<S25>/memory2'
  *memory2_PreviousInput = 0.0;

  // InitializeConditions for Memory: '<S25>/memory1'
  *memory1_PreviousInput = 0.0;
}

//
// Output and update for iterator system:
//    '<S21>/Cross-flow drag trapezoidal integration'
//    '<S117>/Cross-flow drag trapezoidal integration'
//    '<S213>/Cross-flow drag trapezoidal integration'
//    '<S309>/Cross-flow drag trapezoidal integration'
//
void AH_Model_v1ModelClass::Crossflowdragtrapezoidalintegra(real_T rtu_N, real_T
  rtu_dx, real_T rtu_v_r, real_T rtu_r, real_T *rty_sum1, real_T *rty_sum2,
  real_T rtp_Lpp)
{
  int32_T tmp;
  int32_T i;
  real_T rtb_x_lh;
  real_T rtb_memory1;
  real_T memory2_PreviousInput;
  real_T memory1_PreviousInput;

  // Outputs for Iterator SubSystem: '<S21>/Cross-flow drag trapezoidal integration' incorporates:
  //   ForIterator: '<S25>/For Iterator'

  Crossflowdragtrapezoidali_Reset(&memory2_PreviousInput, &memory1_PreviousInput);
  if (rtu_N < 2.147483648E+9) {
    if (rtu_N >= -2.147483648E+9) {
      tmp = static_cast<int32_T>(rtu_N);
    } else {
      tmp = MIN_int32_T;
    }
  } else {
    tmp = MAX_int32_T;
  }

  if (tmp > 2147483646) {
    tmp = 2147483646;
  } else {
    if (tmp < 0) {
      tmp = 0;
    }
  }

  for (i = 1; i <= tmp; i++) {
    // Sum: '<S25>/Sum3' incorporates:
    //   Constant: '<S25>/Lpp//2'
    //   Constant: '<S25>/constant'
    //   Product: '<S25>/Product1'
    //   Sum: '<S25>/Sum4'

    rtb_x_lh = (static_cast<real_T>(i) - 1.0) * rtu_dx - rtp_Lpp / 2.0;

    // Sum: '<S25>/Sum1' incorporates:
    //   Product: '<S25>/x * r'

    rtb_memory1 = rtu_r * rtb_x_lh + rtu_v_r;

    // If: '<S25>/If i=1 or i=N' incorporates:
    //   Abs: '<S25>/Abs'
    //   Inport: '<S27>/In1'
    //   Product: '<S25>/Product'
    //   Product: '<S25>/Product3'

    if ((i == 1) || (i == rtu_N)) {
      // Outputs for IfAction SubSystem: '<S25>/multiply with 0.5' incorporates:
      //   ActionPort: '<S26>/Action Port'

      // Gain: '<S26>/Gain' incorporates:
      //   Abs: '<S25>/Abs'
      //   Product: '<S25>/Product'
      //   Product: '<S25>/Product3'

      rtb_memory1 = 0.5 * (std::abs(rtb_memory1) * rtb_memory1 * rtu_dx);

      // End of Outputs for SubSystem: '<S25>/multiply with 0.5'
    } else {
      // Outputs for IfAction SubSystem: '<S25>/multiply with 1' incorporates:
      //   ActionPort: '<S27>/Action Port'

      rtb_memory1 = std::abs(rtb_memory1) * rtb_memory1 * rtu_dx;

      // End of Outputs for SubSystem: '<S25>/multiply with 1'
    }

    // End of If: '<S25>/If i=1 or i=N'

    // Sum: '<S25>/Sum2' incorporates:
    //   Memory: '<S25>/memory2'
    //   Product: '<S25>/Product4'

    *rty_sum2 = rtb_memory1 * rtb_x_lh + memory2_PreviousInput;

    // Sum: '<S25>/Sum' incorporates:
    //   Memory: '<S25>/memory1'

    *rty_sum1 = memory1_PreviousInput + rtb_memory1;

    // Update for Memory: '<S25>/memory2'
    memory2_PreviousInput = *rty_sum2;

    // Update for Memory: '<S25>/memory1'
    memory1_PreviousInput = *rty_sum1;
  }

  // End of Outputs for SubSystem: '<S21>/Cross-flow drag trapezoidal integration' 
}

//
// Output and update for atomic system:
//    '<S66>/Chart'
//    '<S162>/Chart'
//    '<S258>/Chart'
//    '<S354>/Chart'
//
void AH_Model_v1ModelClass::AHV_Model_Chart(boolean_T rtu_hold, real_T rtu_x_ref,
  real_T rtu_y_ref, real_T rtu_x_hold, real_T rtu_y_hold, real_T rtu_yaw_ref,
  real_T rtu_yaw_hold, real_T *rty_x_ref_rel, real_T *rty_y_ref_rel, real_T
  *rty_yaw_ref_rel, DW_Chart_AHV_Model_T *localDW)
{
  // Chart: '<S66>/Chart'
  if (localDW->is_active_c129_AHV_Model == 0U) {
    localDW->is_active_c129_AHV_Model = 1U;
    if (rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_hold;
      localDW->x_local = rtu_x_hold;
      localDW->y_local = rtu_y_hold;
      localDW->yaw_local = rtu_yaw_hold;
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    } else {
      localDW->is_c129_AHV_Model = AHV_Model_IN_unhold;
      *rty_yaw_ref_rel = rtu_yaw_ref;
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    }
  } else if (localDW->is_c129_AHV_Model == AHV_Model_IN_hold) {
    if (!rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_unhold;
      *rty_yaw_ref_rel = rtu_yaw_ref;
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    } else {
      *rty_yaw_ref_rel = localDW->yaw_local;
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    }
  } else {
    // case IN_unhold:
    if (rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_hold;
      localDW->x_local = rtu_x_hold;
      localDW->y_local = rtu_y_hold;
      localDW->yaw_local = rtu_yaw_hold;
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    } else {
      *rty_yaw_ref_rel = rtu_yaw_ref;
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    }
  }

  // End of Chart: '<S66>/Chart'
}

//
// Output and update for action system:
//    '<S78>/If Action Subsystem'
//    '<S174>/If Action Subsystem'
//    '<S270>/If Action Subsystem'
//    '<S366>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem(real_T rtu_In1, real_T
  rtu_In1_b, real_T rtu_In1_l, real_T rty_Kp[3], const real_T rtp_Kp_A[9])
{
  int32_T i;

  // Gain: '<S83>/Kp_u' incorporates:
  //   SignalConversion generated from: '<S83>/Kp_u'

  for (i = 0; i < 3; i++) {
    rty_Kp[i] = 0.0;
    rty_Kp[i] += rtp_Kp_A[i] * rtu_In1;
    rty_Kp[i] += rtp_Kp_A[i + 3] * rtu_In1_b;
    rty_Kp[i] += rtp_Kp_A[i + 6] * rtu_In1_l;
  }

  // End of Gain: '<S83>/Kp_u'
}

//
// Output and update for action system:
//    '<S78>/If Action Subsystem1'
//    '<S174>/If Action Subsystem1'
//    '<S270>/If Action Subsystem1'
//    '<S366>/If Action Subsystem1'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem1(real_T rtu_In1, real_T
  rtu_In1_i, real_T rtu_In1_k, real_T rty_Kp[3], const real_T rtp_Kp[9])
{
  int32_T i;

  // Gain: '<S84>/Kp_u' incorporates:
  //   SignalConversion generated from: '<S84>/Kp_u'

  for (i = 0; i < 3; i++) {
    rty_Kp[i] = 0.0;
    rty_Kp[i] += rtp_Kp[i] * rtu_In1;
    rty_Kp[i] += rtp_Kp[i + 3] * rtu_In1_i;
    rty_Kp[i] += rtp_Kp[i + 6] * rtu_In1_k;
  }

  // End of Gain: '<S84>/Kp_u'
}

//
// Output and update for action system:
//    '<S85>/If Action Subsystem'
//    '<S86>/If Action Subsystem'
//    '<S181>/If Action Subsystem'
//    '<S182>/If Action Subsystem'
//    '<S277>/If Action Subsystem'
//    '<S278>/If Action Subsystem'
//    '<S373>/If Action Subsystem'
//    '<S374>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem_f(real_T rtu_Increment,
  real_T rtu_main_T, real_T *rty_Out1)
{
  // Switch: '<S87>/Switch' incorporates:
  //   Gain: '<S87>/Gain'

  if (rtu_main_T > 0.0) {
    *rty_Out1 = rtu_Increment;
  } else {
    *rty_Out1 = -rtu_Increment;
  }

  // End of Switch: '<S87>/Switch'
}

//
// System initialize for atomic system:
//    '<S94>/Chart1'
//    '<S94>/Chart2'
//    '<S190>/Chart2'
//    '<S286>/Chart2'
//    '<S382>/Chart2'
//
void AH_Model_v1ModelClass::AHV_Model_Chart1_Init(real_T rtu_Rudder_angle,
  real_T rtu_eta_yaw_deg, real_T *rty_heading_deg, DW_Chart1_AHV_Model_T
  *localDW)
{
  localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD;
  localDW->yaw_local_deg = 0.0;
  *rty_heading_deg = 0.0;

  // Chart: '<S94>/Chart1'
  if (rtu_Rudder_angle == 0.0) {
    // 不改变航向
    localDW->is_c1_AHV_Model = AHV_Model_IN_hold_angel;
    *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    if (localDW->local == 0.0) {
      localDW->is_hold_angel = AHV_Model_IN_local;
      localDW->yaw_local_deg = rtu_eta_yaw_deg;
    } else {
      localDW->is_hold_angel = AHV_Model_IN_unlocal;
    }
  } else {
    localDW->is_c1_AHV_Model = AHV_Model_IN_change_angel;
    *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
  }

  // End of Chart: '<S94>/Chart1'
}

//
// Output and update for atomic system:
//    '<S94>/Chart1'
//    '<S94>/Chart2'
//    '<S190>/Chart2'
//    '<S286>/Chart2'
//    '<S382>/Chart2'
//
void AH_Model_v1ModelClass::AHV_Model_Chart1(real_T rtu_Rudder_angle, real_T
  rtu_eta_yaw_deg, real_T *rty_heading_deg, DW_Chart1_AHV_Model_T *localDW)
{
  *rty_heading_deg = 0.0;
  if (localDW->is_c1_AHV_Model == AHV_Model_IN_change_angel) {
    if (rtu_Rudder_angle == 0.0) {
      // 不改变航向
      localDW->is_c1_AHV_Model = AHV_Model_IN_hold_angel;
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
      if (localDW->local == 0.0) {
        localDW->is_hold_angel = AHV_Model_IN_local;
        localDW->yaw_local_deg = rtu_eta_yaw_deg;
      } else {
        localDW->is_hold_angel = AHV_Model_IN_unlocal;
      }
    } else {
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    }
  } else {
    // case IN_hold_angel:
    if (rtu_Rudder_angle != 0.0) {
      // 改变航向
      if (localDW->is_hold_angel == AHV_Model_IN_local) {
        localDW->yaw_local_deg = rtu_eta_yaw_deg;
        localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD;
      } else {
        localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD;
      }

      localDW->is_c1_AHV_Model = AHV_Model_IN_change_angel;
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    } else {
      *rty_heading_deg = localDW->yaw_local_deg;
      if (localDW->is_hold_angel == AHV_Model_IN_local) {
        if (localDW->local == 0.0) {
          localDW->yaw_local_deg = rtu_eta_yaw_deg;
          localDW->is_hold_angel = AHV_Model_IN_unlocal;
        } else {
          localDW->yaw_local_deg = rtu_eta_yaw_deg;
        }
      } else {
        // case IN_unlocal:
      }
    }
  }
}

//
// System initialize for atomic system:
//    '<S94>/Chart3'
//    '<S190>/Chart3'
//    '<S286>/Chart3'
//    '<S382>/Chart3'
//
void AH_Model_v1ModelClass::AHV_Model_Chart3_Init(real_T rtu_eta_yaw_deg, real_T
  rtu_Rudder_angle, real_T *rty_heading_deg, DW_Chart3_AHV_Model_T *localDW)
{
  localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_j;
  localDW->yaw_local_deg = 0.0;
  *rty_heading_deg = 0.0;

  // Chart: '<S94>/Chart3'
  if (rtu_Rudder_angle == 0.0) {
    // 不改变航向
    localDW->is_c3_AHV_Model = AHV_Model_IN_hold_angel_l;
    *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    if (localDW->local == 0.0) {
      localDW->is_hold_angel = AHV_Model_IN_local_f;
      localDW->yaw_local_deg = rtu_eta_yaw_deg;
    } else {
      localDW->is_hold_angel = AHV_Model_IN_unlocal_h;
    }
  } else {
    localDW->is_c3_AHV_Model = AHV_Model_IN_change_angel_n;
    *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
  }

  // End of Chart: '<S94>/Chart3'
}

//
// Output and update for atomic system:
//    '<S94>/Chart3'
//    '<S190>/Chart3'
//    '<S286>/Chart3'
//    '<S382>/Chart3'
//
void AH_Model_v1ModelClass::AHV_Model_Chart3(real_T rtu_eta_yaw_deg, real_T
  rtu_Rudder_angle, real_T *rty_heading_deg, DW_Chart3_AHV_Model_T *localDW)
{
  *rty_heading_deg = 0.0;
  if (localDW->is_c3_AHV_Model == AHV_Model_IN_change_angel_n) {
    if (rtu_Rudder_angle == 0.0) {
      // 不改变航向
      localDW->is_c3_AHV_Model = AHV_Model_IN_hold_angel_l;
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
      if (localDW->local == 0.0) {
        localDW->is_hold_angel = AHV_Model_IN_local_f;
        localDW->yaw_local_deg = rtu_eta_yaw_deg;
      } else {
        localDW->is_hold_angel = AHV_Model_IN_unlocal_h;
      }
    } else {
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    }
  } else {
    // case IN_hold_angel:
    if (rtu_Rudder_angle != 0.0) {
      // 改变航向
      if (localDW->is_hold_angel == AHV_Model_IN_local_f) {
        localDW->yaw_local_deg = rtu_eta_yaw_deg;
        localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_j;
      } else {
        localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_j;
      }

      localDW->is_c3_AHV_Model = AHV_Model_IN_change_angel_n;
      *rty_heading_deg = rtu_eta_yaw_deg + rtu_Rudder_angle;
    } else {
      *rty_heading_deg = localDW->yaw_local_deg;
      if (localDW->is_hold_angel == AHV_Model_IN_local_f) {
        if (localDW->local == 0.0) {
          localDW->yaw_local_deg = rtu_eta_yaw_deg;
          localDW->is_hold_angel = AHV_Model_IN_unlocal_h;
        } else {
          localDW->yaw_local_deg = rtu_eta_yaw_deg;
        }
      } else {
        // case IN_unlocal:
      }
    }
  }
}

//
// Output and update for action system:
//    '<S8>/If Action Subsystem1'
//    '<S104>/If Action Subsystem1'
//    '<S200>/If Action Subsystem1'
//    '<S296>/If Action Subsystem1'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem1_k(real_T
  rtu_heading_mode, real_T rtu_Head_angel, real_T rtu_target_surge, real_T
  rtu_main_add, real_T rtu_fu_add, real_T rtu_target_sway, real_T
  *rty_Target_yaw, real_T *rty_Target_x, real_T *rty_surge_add, real_T
  *rty_sway_add, real_T *rty_Target_y, real_T *rty_heading_mode_dp, real_T
  *rty_SurgeForce, real_T *rty_SwayForce, real_T *rty_YawMoment)
{
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // SignalConversion generated from: '<S92>/SurgeForce' incorporates:
    //   Constant: '<S92>/Surge Force'

    *rty_SurgeForce = 1.0;

    // SignalConversion generated from: '<S92>/SwayForce' incorporates:
    //   Constant: '<S92>/Sway Force'

    *rty_SwayForce = 1.0;

    // SignalConversion generated from: '<S92>/YawMoment' incorporates:
    //   Constant: '<S92>/Yaw Moment'

    *rty_YawMoment = 1.0;
  }

  // Inport: '<S92>/Head_angel'
  *rty_Target_yaw = rtu_Head_angel;

  // Inport: '<S92>/target_surge'
  *rty_Target_x = rtu_target_surge;

  // Inport: '<S92>/main_add'
  *rty_surge_add = rtu_main_add;

  // Inport: '<S92>/fu_add'
  *rty_sway_add = rtu_fu_add;

  // Inport: '<S92>/target_sway'
  *rty_Target_y = rtu_target_sway;

  // Inport: '<S92>/heading_mode'
  *rty_heading_mode_dp = rtu_heading_mode;
}

//
// Output and update for action system:
//    '<S93>/If Action Subsystem2'
//    '<S189>/If Action Subsystem2'
//    '<S285>/If Action Subsystem2'
//    '<S381>/If Action Subsystem2'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem2(boolean_T
  *rty_KeepPosition)
{
  if (1) {
    // SignalConversion generated from: '<S99>/KeepPosition' incorporates:
    //   Constant: '<S99>/Constant'

    *rty_KeepPosition = true;
  }
}

//
// System initialize for action system:
//    '<S104>/If Action Subsystem'
//    '<S200>/If Action Subsystem'
//    '<S296>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Mo_IfActionSubsystem_j_Init(const real_T
  rtu_eta[6], B_IfActionSubsystem_AHV_Mod_d_T *localB,
  DW_IfActionSubsystem_AHV_M_ci_T *localDW)
{
  // InitializeConditions for RateLimiter: '<S190>/Rate Limiter1'
  localDW->LastMajorTime = (rtInf);

  // InitializeConditions for RateLimiter: '<S190>/Rate Limiter2'
  localDW->LastMajorTime_l = (rtInf);

  // SystemInitialize for Chart: '<S190>/Chart2'
  AHV_Model_Chart1_Init(localB->RateLimiter1, localB->Gain4,
                        &localB->heading_deg_n, &localDW->sf_Chart2);

  // SystemInitialize for Chart: '<S190>/Chart3'
  AHV_Model_Chart3_Init(rtu_eta[0], localB->Product, &localB->heading_deg,
                        &localDW->sf_Chart3);

  // SystemInitialize for Chart: '<S190>/Chart1'
  localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_jn;
  localDW->yaw_local_deg = 0.0;

  // Chart: '<S190>/Chart1'
  if (localB->Product1 == 0.0) {
    // 不改变航向
    localDW->is_c9_AHV_Model = AHV_Model_IN_hold_angel_d;
    if (localDW->local == 0.0) {
      localDW->is_hold_angel = AHV_Model_IN_local_j;
      localDW->yaw_local_deg = rtu_eta[1];
    } else {
      localDW->is_hold_angel = AHV_Model_IN_unlocal_a;
    }
  } else {
    localDW->is_c9_AHV_Model = AHV_Model_IN_change_angel_g;
  }

  // End of Chart: '<S190>/Chart1'
}

//
// Outputs for action system:
//    '<S104>/If Action Subsystem'
//    '<S200>/If Action Subsystem'
//    '<S296>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem_c(const real_T rtu_eta[6],
  real_T rtu_Rudder_angle, real_T rtu_Thruster_percentage1, real_T
  rtu_heading_mode, boolean_T rtu_HoldPosition, real_T *rty_Target_yaw_headdeg,
  real_T *rty_Target_x_head, real_T *rty_surge_add_head, real_T
  *rty_sway_add_head, real_T *rty_Target_y_head, real_T *rty_heading_mode_head,
  real_T *rty_SurgeForce, real_T *rty_SwayForce, real_T *rty_YawMoment,
  B_IfActionSubsystem_AHV_Mod_d_T *localB, DW_IfActionSubsystem_AHV_M_ci_T
  *localDW)
{
  real_T riseValLimit;
  real_T rtb_Floor_c;
  real_T rtb_Gain5_c;

  // RateLimiter: '<S190>/Rate Limiter1'
  if (localDW->LastMajorTime == (rtInf)) {
    localB->RateLimiter1 = rtu_Rudder_angle;
  } else {
    rtb_Floor_c = (&AHV_Model_M)->Timing.t[0] - localDW->LastMajorTime;
    riseValLimit = rtb_Floor_c * 0.1;
    rtb_Gain5_c = rtu_Rudder_angle - localDW->PrevY;
    if (rtb_Gain5_c > riseValLimit) {
      localB->RateLimiter1 = localDW->PrevY + riseValLimit;
    } else {
      rtb_Floor_c *= -0.1;
      if (rtb_Gain5_c < rtb_Floor_c) {
        localB->RateLimiter1 = localDW->PrevY + rtb_Floor_c;
      } else {
        localB->RateLimiter1 = rtu_Rudder_angle;
      }
    }
  }

  // End of RateLimiter: '<S190>/Rate Limiter1'

  // Gain: '<S190>/Gain4'
  localB->Gain4 = 57.324840764331206 * rtu_eta[5];
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // Chart: '<S190>/Chart2'
    AHV_Model_Chart1(localB->RateLimiter1, localB->Gain4, &localB->heading_deg_n,
                     &localDW->sf_Chart2);
  }

  // RateLimiter: '<S190>/Rate Limiter2'
  if (localDW->LastMajorTime_l == (rtInf)) {
    localB->RateLimiter2 = localB->heading_deg_n;
  } else {
    rtb_Floor_c = (&AHV_Model_M)->Timing.t[0] - localDW->LastMajorTime_l;
    riseValLimit = rtb_Floor_c * 0.1;
    rtb_Gain5_c = localB->heading_deg_n - localDW->PrevY_c;
    if (rtb_Gain5_c > riseValLimit) {
      localB->RateLimiter2 = localDW->PrevY_c + riseValLimit;
    } else {
      rtb_Floor_c *= -0.1;
      if (rtb_Gain5_c < rtb_Floor_c) {
        localB->RateLimiter2 = localDW->PrevY_c + rtb_Floor_c;
      } else {
        localB->RateLimiter2 = localB->heading_deg_n;
      }
    }
  }

  // End of RateLimiter: '<S190>/Rate Limiter2'

  // Gain: '<S190>/Gain2'
  rtb_Gain5_c = 0.017444444444444446 * localB->RateLimiter2;

  // Gain: '<S190>/Gain3'
  rtb_Floor_c = 50.0 * rtu_Thruster_percentage1;

  // Product: '<S190>/Product' incorporates:
  //   Trigonometry: '<S190>/x'

  localB->Product = std::cos(rtb_Gain5_c) * rtb_Floor_c;
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // Chart: '<S190>/Chart3'
    AHV_Model_Chart3(rtu_eta[0], localB->Product, &localB->heading_deg,
                     &localDW->sf_Chart3);

    // SignalConversion generated from: '<S187>/Target_x_head'
    *rty_Target_x_head = localB->heading_deg;
  }

  // Product: '<S190>/Product1' incorporates:
  //   Trigonometry: '<S190>/y'

  localB->Product1 = rtb_Floor_c * std::sin(rtb_Gain5_c);
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // Chart: '<S190>/Chart1'
    if (localDW->is_c9_AHV_Model == 1) {
      if (localB->Product1 == 0.0) {
        // 不改变航向
        localDW->is_c9_AHV_Model = AHV_Model_IN_hold_angel_d;

        // SignalConversion generated from: '<S187>/Target_y_head'
        *rty_Target_y_head = rtu_eta[1] + localB->Product1;
        if (localDW->local == 0.0) {
          localDW->is_hold_angel = AHV_Model_IN_local_j;
          localDW->yaw_local_deg = rtu_eta[1];
        } else {
          localDW->is_hold_angel = AHV_Model_IN_unlocal_a;
        }
      } else {
        // SignalConversion generated from: '<S187>/Target_y_head'
        *rty_Target_y_head = rtu_eta[1] + localB->Product1;
      }
    } else {
      // case IN_hold_angel:
      if (localB->Product1 != 0.0) {
        // 改变航向
        if (localDW->is_hold_angel == 1) {
          localDW->yaw_local_deg = rtu_eta[1];
          localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_jn;
        } else {
          localDW->is_hold_angel = AHV_Model_IN_NO_ACTIVE_CHILD_jn;
        }

        localDW->is_c9_AHV_Model = AHV_Model_IN_change_angel_g;

        // SignalConversion generated from: '<S187>/Target_y_head'
        *rty_Target_y_head = rtu_eta[1] + localB->Product1;
      } else {
        // SignalConversion generated from: '<S187>/Target_y_head'
        *rty_Target_y_head = localDW->yaw_local_deg;
        if (localDW->is_hold_angel == 1) {
          if (localDW->local == 0.0) {
            localDW->yaw_local_deg = rtu_eta[1];
            localDW->is_hold_angel = AHV_Model_IN_unlocal_a;
          } else {
            localDW->yaw_local_deg = rtu_eta[1];
          }
        } else {
          // case IN_unlocal:
        }
      }
    }

    // End of Chart: '<S190>/Chart1'
  }

  // Switch: '<S190>/Switch1' incorporates:
  //   Constant: '<S190>/Constant'
  //   Constant: '<S190>/Constant1'
  //   Switch: '<S190>/Switch'

  if (rtu_HoldPosition) {
    *rty_sway_add_head = 10.0;
    *rty_surge_add_head = 1.0;
  } else {
    *rty_sway_add_head = localB->Product1;
    *rty_surge_add_head = localB->Product;
  }

  // End of Switch: '<S190>/Switch1'

  // Gain: '<S190>/Gain5'
  rtb_Gain5_c *= 57.324840764331206;

  // Product: '<S194>/Divide' incorporates:
  //   Constant: '<S194>/Constant3'

  rtb_Floor_c = rtb_Gain5_c / 360.0;

  // Rounding: '<S194>/Floor'
  if (rtb_Floor_c < 0.0) {
    rtb_Floor_c = std::ceil(rtb_Floor_c);
  } else {
    rtb_Floor_c = std::floor(rtb_Floor_c);
  }

  // End of Rounding: '<S194>/Floor'

  // Sum: '<S194>/Sum2' incorporates:
  //   Constant: '<S194>/Constant3'
  //   Product: '<S194>/Product2'

  *rty_Target_yaw_headdeg = rtb_Gain5_c - rtb_Floor_c * 360.0;
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // SignalConversion generated from: '<S187>/SurgeForce' incorporates:
    //   Constant: '<S187>/Surge Force'

    *rty_SurgeForce = 1.0;

    // SignalConversion generated from: '<S187>/SwayForce' incorporates:
    //   Constant: '<S187>/Sway Force'

    *rty_SwayForce = 1.0;

    // SignalConversion generated from: '<S187>/YawMoment' incorporates:
    //   Constant: '<S187>/Yaw Moment'

    *rty_YawMoment = 1.0;
  }

  // Gain: '<S187>/Gain'
  *rty_heading_mode_head = 0.0 * rtu_heading_mode;
}

//
// Update for action system:
//    '<S104>/If Action Subsystem'
//    '<S200>/If Action Subsystem'
//    '<S296>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Mo_IfActionSubsystem_Update
  (B_IfActionSubsystem_AHV_Mod_d_T *localB, DW_IfActionSubsystem_AHV_M_ci_T
   *localDW)
{
  // Update for RateLimiter: '<S190>/Rate Limiter1' incorporates:
  //   RateLimiter: '<S190>/Rate Limiter2'

  localDW->PrevY = localB->RateLimiter1;
  localDW->LastMajorTime = (&AHV_Model_M)->Timing.t[0];

  // Update for RateLimiter: '<S190>/Rate Limiter2'
  localDW->PrevY_c = localB->RateLimiter2;
  localDW->LastMajorTime_l = localDW->LastMajorTime;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(static_cast<real_T>(u0_0), static_cast<real_T>(u1_0));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

void rt_mldivide_U1d6x6_U2d_4sw8yi9v(const real_T u0[36], const real_T u1[6],
  real_T y[6])
{
  real_T A[36];
  int8_T ipiv[6];
  int32_T jj;
  int32_T c;
  int32_T c_0;
  int32_T kAcol;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T iy;
  int32_T jy;
  int32_T j;
  int32_T ijA;
  std::memcpy(&A[0], &u0[0], 36U * sizeof(real_T));
  for (iy = 0; iy < 6; iy++) {
    ipiv[iy] = static_cast<int8_T>(iy + 1);
  }

  for (kAcol = 0; kAcol < 5; kAcol++) {
    c = kAcol * 7 + 2;
    jj = kAcol * 7;
    c_0 = 6 - kAcol;
    iy = 1;
    ix = c - 2;
    smax = std::abs(A[jj]);
    for (jy = 2; jy <= c_0; jy++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        iy = jy;
        smax = s;
      }
    }

    if (A[(c + iy) - 3] != 0.0) {
      if (iy - 1 != 0) {
        iy += kAcol;
        ipiv[kAcol] = static_cast<int8_T>(iy);
        ix = kAcol;
        iy--;
        for (jy = 0; jy < 6; jy++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 6;
          iy += 6;
        }
      }

      iy = c - kAcol;
      for (ix = c; ix <= iy + 4; ix++) {
        A[ix - 1] /= A[jj];
      }
    }

    c_0 = 4 - kAcol;
    jy = jj + 6;
    for (j = 0; j <= c_0; j++) {
      if (A[jy] != 0.0) {
        smax = -A[jy];
        ix = c - 1;
        iy = jj - kAcol;
        for (ijA = jj + 8; ijA <= iy + 12; ijA++) {
          A[ijA - 1] += A[ix] * smax;
          ix++;
        }
      }

      jy += 6;
      jj += 6;
    }
  }

  for (iy = 0; iy < 6; iy++) {
    y[iy] = u1[iy];
  }

  for (ix = 0; ix < 5; ix++) {
    if (ix + 1 != ipiv[ix]) {
      smax = y[ix];
      kAcol = ipiv[ix] - 1;
      y[ix] = y[kAcol];
      y[kAcol] = smax;
    }
  }

  for (jy = 0; jy < 6; jy++) {
    kAcol = 6 * jy - 1;
    if (y[jy] != 0.0) {
      for (ix = jy + 2; ix < 7; ix++) {
        y[ix - 1] -= A[ix + kAcol] * y[jy];
      }
    }
  }

  for (jy = 5; jy >= 0; jy--) {
    kAcol = 6 * jy;
    if (y[jy] != 0.0) {
      y[jy] /= A[jy + kAcol];
      iy = jy - 1;
      for (ix = 0; ix <= iy; ix++) {
        y[ix] -= A[ix + kAcol] * y[jy];
      }
    }
  }
}

// Model step function
void AH_Model_v1ModelClass::step()
{
  // local block i/o variables
  real_T rtb_Merge3;
  real_T rtb_Merge4;
  real_T rtb_Merge3_k;
  real_T rtb_Merge4_j;
  real_T rtb_Merge3_j;
  real_T rtb_Merge4_l;
  real_T rtb_Merge3_c;
  real_T rtb_Merge4_g;
  real_T rtb_nu[3];
  real_T rtb_Gain2_g[3];
  real_T rtb_v[3];
  real_T rtb_k_u[3];
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // set solver stop time
    rtsiSetSolverStopTime(&(&AHV_Model_M)->solverInfo,(((&AHV_Model_M)
      ->Timing.clockTick0+1)*(&AHV_Model_M)->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if (rtmIsMinorTimeStep((&AHV_Model_M))) {
    (&AHV_Model_M)->Timing.t[0] = rtsiGetT(&(&AHV_Model_M)->solverInfo);
  }

  {
    real_T rtb_Product1_lf;
    real_T rtb_Integrator_m[6];
    real_T rtb_Product2_n;
    real_T rtb_Product1_e3;
    real_T rtb_rxpi_c;
    real_T rtb_Product2_a;
    real_T rtb_Product1_c;
    real_T rtb_rxpi_n;
    real_T rtb_Sum_j;
    real_T rtb_y_dot;
    real_T rtb_sintheta[3];
    real_T rtb_costheta[3];
    real_T rtb_TransferFcn;
    real_T rtb_Gain1_n[3];
    real_T rtb_Switch2;
    real_T rtb_TmpSignalConversionAtProduc[9];
    real_T rtb_TmpSignalConversionAtProd_i[9];
    real_T rtb_TmpSignalConversionAtProd_o[9];
    real_T rtb_TmpSignalConversionAtProd_n[9];
    real_T rtb_TmpSignalConversionAtProd_f[9];
    real_T rtb_TmpSignalConversionAtPro_fx[9];
    real_T rtb_TmpSignalConversionAtProd_j[9];
    real_T rtb_TmpSignalConversionAtPro_ni[9];
    real_T rtb_rxpi_h;
    real_T rtb_Product_g;
    real_T rtb_Product1_lq;
    real_T rtb_Product2_b;
    real_T rtb_Integrator_f[6];
    real_T rtb_Sum[6];
    real_T SwayFailure;
    real_T rtb_tau_WF[6];
    real_T rtb_tau_WD[6];
    real_T rtb_Row2_e;
    real_T rtb_tau_WF_em[6];
    real_T rtb_tau_WD_n[6];
    real_T rtb_tau_WF_o[6];
    real_T rtb_tau_WD_p[6];
    real_T rtb_Switch3;
    real_T psi_mean_h;
    int32_T i;
    real_T tmp;
    real_T tmp_0;
    real_T tmp_1;
    real_T tmp_2;
    real_T tmp_3;
    real_T rtb_costheta_0[6];
    real_T tmp_4[4];
    real_T tmp_5[4];
    real_T tmp_6[4];
    real_T tmp_7[4];
    int32_T i_0;
    real_T tmp_8;
    real_T rtb_costheta_1;
    real_T rtb_TmpSignalConversionAtProd_9;
    real_T rtb_costheta_tmp;
    real_T rtb_costheta_tmp_0;
    real_T rtb_sintheta_tmp;
    real_T rtb_costheta_tmp_1;
    real_T rtb_sintheta_tmp_0;
    real_T rtb_costheta_tmp_2;
    real_T rtb_sintheta_tmp_1;
    real_T rtb_costheta_tmp_3;
    real_T rtb_TmpSignalConversionAtProd_0;
    real_T rtb_costheta_tmp_4;
    real_T rtb_sintheta_tmp_2;
    real_T rtb_costheta_tmp_5;
    real_T rtb_costheta_tmp_6;
    real_T rtb_costheta_tmp_7;
    real_T rtb_costheta_tmp_8;
    real_T rtb_costheta_tmp_9;
    int32_T tmp_9;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S10>/Subsystem3'
      // Inport: '<Root>/spectrum_type' incorporates:
      //   Inport: '<Root>/depth'
      //   Inport: '<Root>/dir_cutoff'
      //   Inport: '<Root>/energylim'
      //   Inport: '<Root>/freq_cutoff'
      //   Inport: '<Root>/gamma_value'
      //   Inport: '<Root>/hs'
      //   Inport: '<Root>/ndir'
      //   Inport: '<Root>/nfreq'
      //   Inport: '<Root>/omega_peak'
      //   Inport: '<Root>/psi_mean'
      //   Inport: '<Root>/rand_dir'
      //   Inport: '<Root>/rand_freq'
      //   Inport: '<Root>/spread'

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread, depth,
                nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq,
                rand_dir, AHV_Model_B.Zeta_a_c, AHV_Model_B.Omega_d,
                AHV_Model_B.Phase_a, AHV_Model_B.Wavenum_b, AHV_Model_B.Psi_g,
                &psi_mean_h, &AHV_Model_B.Subsystem3, &AHV_Model_DW.Subsystem3,
                20.0, 10.0);

      // End of Outputs for SubSystem: '<S10>/Subsystem3'
    }

    // Integrator: '<S15>/Integrator1' incorporates:
    //   Inport: '<Root>/Vessel_init1'

    for (i = 0; i < 6; i++) {
      if (AHV_Model_DW.Integrator1_IWORK != 0) {
        AHV_Model_X.Integrator1_CSTATE[i] = Vessel_init1[i];
      }

      AHV_Model_B.Integrator1[i] = AHV_Model_X.Integrator1_CSTATE[i];
    }

    // End of Integrator: '<S15>/Integrator1'

    // Outputs for Atomic SubSystem: '<S10>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_B.Integrator1, AHV_Model_B.Psi_g,
                   AHV_Model_B.Wavenum_b, AHV_Model_B.Omega_d,
                   AHV_Model_B.Phase_a, AHV_Model_B.Zeta_a_c, rtb_tau_WF_o,
                   rtb_tau_WD_p, &AHV_Model_B.WaveloadsU0,
                   &AHV_Model_ConstB.WaveloadsU0, &AHV_Model_DW.WaveloadsU0);

    // End of Outputs for SubSystem: '<S10>/Wave loads (U=0)'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S106>/Subsystem3'
      // Inport: '<Root>/spectrum_type' incorporates:
      //   Inport: '<Root>/depth'
      //   Inport: '<Root>/dir_cutoff'
      //   Inport: '<Root>/energylim'
      //   Inport: '<Root>/freq_cutoff'
      //   Inport: '<Root>/gamma_value'
      //   Inport: '<Root>/hs'
      //   Inport: '<Root>/ndir'
      //   Inport: '<Root>/nfreq'
      //   Inport: '<Root>/omega_peak'
      //   Inport: '<Root>/psi_mean'
      //   Inport: '<Root>/rand_dir'
      //   Inport: '<Root>/rand_freq'
      //   Inport: '<Root>/spread'

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread, depth,
                nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq,
                rand_dir, AHV_Model_B.Zeta_a_l, AHV_Model_B.Omega_f,
                AHV_Model_B.Phase_i, AHV_Model_B.Wavenum_do, AHV_Model_B.Psi_b,
                &psi_mean_h, &AHV_Model_B.Subsystem3_d,
                &AHV_Model_DW.Subsystem3_d, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S106>/Subsystem3'
    }

    // Integrator: '<S111>/Integrator1' incorporates:
    //   Inport: '<Root>/Vessel_init2'

    for (i = 0; i < 6; i++) {
      if (AHV_Model_DW.Integrator1_IWORK_o != 0) {
        AHV_Model_X.Integrator1_CSTATE_n[i] = Vessel_init2[i];
      }

      AHV_Model_B.Integrator1_b[i] = AHV_Model_X.Integrator1_CSTATE_n[i];
    }

    // End of Integrator: '<S111>/Integrator1'

    // Outputs for Atomic SubSystem: '<S106>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_B.Integrator1_b, AHV_Model_B.Psi_b,
                   AHV_Model_B.Wavenum_do, AHV_Model_B.Omega_f,
                   AHV_Model_B.Phase_i, AHV_Model_B.Zeta_a_l, rtb_tau_WF_em,
                   rtb_tau_WD_n, &AHV_Model_B.WaveloadsU0_h,
                   &AHV_Model_ConstB.WaveloadsU0_h, &AHV_Model_DW.WaveloadsU0_h);

    // End of Outputs for SubSystem: '<S106>/Wave loads (U=0)'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S202>/Subsystem3'
      // Inport: '<Root>/spectrum_type' incorporates:
      //   Inport: '<Root>/depth'
      //   Inport: '<Root>/dir_cutoff'
      //   Inport: '<Root>/energylim'
      //   Inport: '<Root>/freq_cutoff'
      //   Inport: '<Root>/gamma_value'
      //   Inport: '<Root>/hs'
      //   Inport: '<Root>/ndir'
      //   Inport: '<Root>/nfreq'
      //   Inport: '<Root>/omega_peak'
      //   Inport: '<Root>/psi_mean'
      //   Inport: '<Root>/rand_dir'
      //   Inport: '<Root>/rand_freq'
      //   Inport: '<Root>/spread'

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread, depth,
                nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq,
                rand_dir, AHV_Model_B.Zeta_a_f, AHV_Model_B.Omega_c,
                AHV_Model_B.Phase_f, AHV_Model_B.Wavenum_d, AHV_Model_B.Psi_k,
                &psi_mean_h, &AHV_Model_B.Subsystem3_h,
                &AHV_Model_DW.Subsystem3_h, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S202>/Subsystem3'
    }

    // Integrator: '<S207>/Integrator1' incorporates:
    //   Inport: '<Root>/Vessel_init3'

    for (i = 0; i < 6; i++) {
      if (AHV_Model_DW.Integrator1_IWORK_b != 0) {
        AHV_Model_X.Integrator1_CSTATE_h[i] = Vessel_init3[i];
      }

      AHV_Model_B.Integrator1_p[i] = AHV_Model_X.Integrator1_CSTATE_h[i];
    }

    // End of Integrator: '<S207>/Integrator1'

    // Outputs for Atomic SubSystem: '<S202>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_B.Integrator1_p, AHV_Model_B.Psi_k,
                   AHV_Model_B.Wavenum_d, AHV_Model_B.Omega_c,
                   AHV_Model_B.Phase_f, AHV_Model_B.Zeta_a_f, rtb_tau_WF,
                   rtb_tau_WD, &AHV_Model_B.WaveloadsU0_p,
                   &AHV_Model_ConstB.WaveloadsU0_p, &AHV_Model_DW.WaveloadsU0_p);

    // End of Outputs for SubSystem: '<S202>/Wave loads (U=0)'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S298>/Subsystem3'
      // Inport: '<Root>/spectrum_type' incorporates:
      //   Inport: '<Root>/depth'
      //   Inport: '<Root>/dir_cutoff'
      //   Inport: '<Root>/energylim'
      //   Inport: '<Root>/freq_cutoff'
      //   Inport: '<Root>/gamma_value'
      //   Inport: '<Root>/hs'
      //   Inport: '<Root>/ndir'
      //   Inport: '<Root>/nfreq'
      //   Inport: '<Root>/omega_peak'
      //   Inport: '<Root>/psi_mean'
      //   Inport: '<Root>/rand_dir'
      //   Inport: '<Root>/rand_freq'
      //   Inport: '<Root>/spread'

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread, depth,
                nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq,
                rand_dir, AHV_Model_B.Zeta_a, AHV_Model_B.Omega,
                AHV_Model_B.Phase, AHV_Model_B.Wavenum, AHV_Model_B.Psi,
                &psi_mean_h, &AHV_Model_B.Subsystem3_o,
                &AHV_Model_DW.Subsystem3_o, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S298>/Subsystem3'
    }

    // Integrator: '<S303>/Integrator1' incorporates:
    //   Inport: '<Root>/Vessel_init4'

    for (i = 0; i < 6; i++) {
      if (AHV_Model_DW.Integrator1_IWORK_c != 0) {
        AHV_Model_X.Integrator1_CSTATE_hn[i] = Vessel_init4[i];
      }

      AHV_Model_B.Integrator1_n[i] = AHV_Model_X.Integrator1_CSTATE_hn[i];
    }

    // End of Integrator: '<S303>/Integrator1'

    // Outputs for Atomic SubSystem: '<S298>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_B.Integrator1_n, AHV_Model_B.Psi,
                   AHV_Model_B.Wavenum, AHV_Model_B.Omega, AHV_Model_B.Phase,
                   AHV_Model_B.Zeta_a, rtb_Integrator_f, rtb_Sum,
                   &AHV_Model_B.WaveloadsU0_j, &AHV_Model_ConstB.WaveloadsU0_j,
                   &AHV_Model_DW.WaveloadsU0_j);

    // End of Outputs for SubSystem: '<S298>/Wave loads (U=0)'

    // Trigonometry: '<S23>/sin(theta)' incorporates:
    //   Fcn: '<S24>/T21 '
    //   Trigonometry: '<S63>/sin(theta)'

    rtb_y_dot = std::sin(AHV_Model_B.Integrator1[3]);

    // Trigonometry: '<S23>/cos(theta)' incorporates:
    //   Fcn: '<S24>/T23'
    //   Fcn: '<S24>/T31 '
    //   Trigonometry: '<S63>/cos(theta)'

    rtb_costheta_tmp = std::cos(AHV_Model_B.Integrator1[3]);
    rtb_costheta_tmp_0 = std::cos(AHV_Model_B.Integrator1[4]);

    // Trigonometry: '<S23>/sin(theta)' incorporates:
    //   Trigonometry: '<S63>/sin(theta)'

    rtb_TransferFcn = std::sin(AHV_Model_B.Integrator1[4]);

    // Trigonometry: '<S23>/cos(theta)' incorporates:
    //   Fcn: '<S70>/Fcn'
    //   Fcn: '<S70>/Fcn1'
    //   Fcn: '<S71>/Row1'
    //   Fcn: '<S71>/Row2'
    //   Fcn: '<S72>/Row1'
    //   Fcn: '<S72>/Row2'
    //   Trigonometry: '<S63>/cos(theta)'

    psi_mean_h = std::cos(AHV_Model_B.Integrator1[5]);

    // Trigonometry: '<S23>/sin(theta)' incorporates:
    //   Fcn: '<S70>/Fcn'
    //   Fcn: '<S70>/Fcn1'
    //   Fcn: '<S71>/Row1'
    //   Fcn: '<S71>/Row2'
    //   Fcn: '<S72>/Row1'
    //   Trigonometry: '<S63>/sin(theta)'

    rtb_sintheta_tmp = std::sin(AHV_Model_B.Integrator1[5]);

    // SignalConversion generated from: '<S20>/Product' incorporates:
    //   Fcn: '<S23>/R11'
    //   Fcn: '<S23>/R12'
    //   Fcn: '<S23>/R21 '
    //   Fcn: '<S23>/R31 '
    //   Trigonometry: '<S23>/cos(theta)'
    //   Trigonometry: '<S23>/sin(theta)'

    rtb_TmpSignalConversionAtProduc[0] = psi_mean_h * rtb_costheta_tmp_0;
    rtb_TmpSignalConversionAtProduc[1] = rtb_sintheta_tmp * rtb_costheta_tmp_0;
    rtb_TmpSignalConversionAtProduc[2] = -rtb_TransferFcn;
    rtb_TmpSignalConversionAtProduc[3] = psi_mean_h * rtb_TransferFcn *
      rtb_y_dot - rtb_sintheta_tmp * rtb_costheta_tmp;

    // Fcn: '<S23>/R22' incorporates:
    //   Fcn: '<S23>/R23'
    //   Trigonometry: '<S23>/sin(theta)'

    rtb_TmpSignalConversionAtProd_9 = rtb_sintheta_tmp * rtb_TransferFcn;

    // SignalConversion generated from: '<S20>/Product' incorporates:
    //   Fcn: '<S23>/R13'
    //   Fcn: '<S23>/R22'
    //   Fcn: '<S23>/R23'
    //   Fcn: '<S23>/R32'
    //   Fcn: '<S23>/R33'
    //   Trigonometry: '<S23>/cos(theta)'
    //   Trigonometry: '<S23>/sin(theta)'

    rtb_TmpSignalConversionAtProduc[4] = rtb_TmpSignalConversionAtProd_9 *
      rtb_y_dot + psi_mean_h * rtb_costheta_tmp;
    rtb_TmpSignalConversionAtProduc[5] = rtb_costheta_tmp_0 * rtb_y_dot;
    rtb_TmpSignalConversionAtProduc[6] = rtb_costheta_tmp * rtb_TransferFcn *
      psi_mean_h + rtb_sintheta_tmp * rtb_y_dot;
    rtb_TmpSignalConversionAtProduc[7] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp - psi_mean_h * rtb_y_dot;
    rtb_TmpSignalConversionAtProduc[8] = rtb_costheta_tmp_0 * rtb_costheta_tmp;

    // SignalConversion generated from: '<S11>/Product1' incorporates:
    //   Fcn: '<S63>/R11'

    rtb_TmpSignalConversionAtProd_i[0] = psi_mean_h * rtb_costheta_tmp_0;

    // Fcn: '<S63>/R21 ' incorporates:
    //   Fcn: '<S63>/R31 '

    rtb_TmpSignalConversionAtProd_9 = psi_mean_h * rtb_TransferFcn;

    // SignalConversion generated from: '<S11>/Product1' incorporates:
    //   Fcn: '<S63>/R12'
    //   Fcn: '<S63>/R21 '
    //   Fcn: '<S63>/R31 '

    rtb_TmpSignalConversionAtProd_i[1] = rtb_TmpSignalConversionAtProd_9 *
      rtb_y_dot + -rtb_sintheta_tmp * rtb_costheta_tmp;
    rtb_TmpSignalConversionAtProd_i[2] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp + rtb_sintheta_tmp * rtb_y_dot;
    rtb_TmpSignalConversionAtProd_i[3] = rtb_sintheta_tmp * rtb_costheta_tmp_0;

    // Fcn: '<S63>/R22' incorporates:
    //   Fcn: '<S63>/R32'

    rtb_TmpSignalConversionAtProd_9 = rtb_sintheta_tmp * rtb_TransferFcn;

    // SignalConversion generated from: '<S11>/Product1' incorporates:
    //   Fcn: '<S63>/R13'
    //   Fcn: '<S63>/R22'
    //   Fcn: '<S63>/R23'
    //   Fcn: '<S63>/R32'

    rtb_TmpSignalConversionAtProd_i[4] = rtb_TmpSignalConversionAtProd_9 *
      rtb_y_dot + psi_mean_h * rtb_costheta_tmp;
    rtb_TmpSignalConversionAtProd_i[5] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp + -psi_mean_h * rtb_y_dot;
    rtb_TmpSignalConversionAtProd_i[6] = -rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[7] = rtb_costheta_tmp_0 * rtb_y_dot;

    // Trigonometry: '<S119>/sin(theta)' incorporates:
    //   Fcn: '<S120>/T21 '
    //   Trigonometry: '<S159>/sin(theta)'

    rtb_TransferFcn = std::sin(AHV_Model_B.Integrator1_b[3]);

    // Trigonometry: '<S119>/cos(theta)' incorporates:
    //   Fcn: '<S120>/T31 '
    //   Trigonometry: '<S159>/cos(theta)'

    rtb_costheta_tmp_1 = std::cos(AHV_Model_B.Integrator1_b[3]);

    // Trigonometry: '<S119>/sin(theta)' incorporates:
    //   Trigonometry: '<S159>/sin(theta)'

    rtb_sintheta_tmp_0 = std::sin(AHV_Model_B.Integrator1_b[4]);

    // Trigonometry: '<S119>/cos(theta)' incorporates:
    //   Fcn: '<S120>/T23'
    //   Trigonometry: '<S159>/cos(theta)'

    rtb_costheta_tmp_2 = std::cos(AHV_Model_B.Integrator1_b[4]);

    // Trigonometry: '<S119>/sin(theta)' incorporates:
    //   Fcn: '<S166>/Fcn'
    //   Fcn: '<S166>/Fcn1'
    //   Fcn: '<S167>/Row1'
    //   Fcn: '<S167>/Row2'
    //   Trigonometry: '<S159>/sin(theta)'

    rtb_sintheta_tmp_1 = std::sin(AHV_Model_B.Integrator1_b[5]);

    // Trigonometry: '<S119>/cos(theta)' incorporates:
    //   Fcn: '<S166>/Fcn'
    //   Fcn: '<S166>/Fcn1'
    //   Fcn: '<S167>/Row1'
    //   Fcn: '<S167>/Row2'
    //   Fcn: '<S168>/Row2'
    //   Trigonometry: '<S159>/cos(theta)'

    rtb_costheta_tmp_3 = std::cos(AHV_Model_B.Integrator1_b[5]);

    // SignalConversion generated from: '<S11>/Product1' incorporates:
    //   Fcn: '<S63>/R33'

    rtb_TmpSignalConversionAtProd_i[8] = rtb_costheta_tmp_0 * rtb_costheta_tmp;

    // SignalConversion generated from: '<S116>/Product' incorporates:
    //   Fcn: '<S119>/R11'
    //   Fcn: '<S119>/R12'
    //   Fcn: '<S119>/R21 '
    //   Fcn: '<S119>/R31 '
    //   Trigonometry: '<S119>/cos(theta)'
    //   Trigonometry: '<S119>/sin(theta)'

    rtb_TmpSignalConversionAtProd_o[0] = rtb_costheta_tmp_3 * rtb_costheta_tmp_2;
    rtb_TmpSignalConversionAtProd_o[1] = rtb_sintheta_tmp_1 * rtb_costheta_tmp_2;
    rtb_TmpSignalConversionAtProd_o[2] = -rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_o[3] = rtb_costheta_tmp_3 * rtb_sintheta_tmp_0
      * rtb_TransferFcn - rtb_sintheta_tmp_1 * rtb_costheta_tmp_1;

    // Fcn: '<S119>/R22' incorporates:
    //   Fcn: '<S119>/R23'
    //   Trigonometry: '<S119>/sin(theta)'

    rtb_TmpSignalConversionAtProd_9 = rtb_sintheta_tmp_1 * rtb_sintheta_tmp_0;

    // SignalConversion generated from: '<S116>/Product' incorporates:
    //   Fcn: '<S119>/R13'
    //   Fcn: '<S119>/R22'
    //   Fcn: '<S119>/R23'
    //   Fcn: '<S119>/R32'
    //   Fcn: '<S119>/R33'
    //   Trigonometry: '<S119>/cos(theta)'
    //   Trigonometry: '<S119>/sin(theta)'

    rtb_TmpSignalConversionAtProd_o[4] = rtb_TmpSignalConversionAtProd_9 *
      rtb_TransferFcn + rtb_costheta_tmp_3 * rtb_costheta_tmp_1;
    rtb_TmpSignalConversionAtProd_o[5] = rtb_costheta_tmp_2 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_o[6] = rtb_costheta_tmp_1 * rtb_sintheta_tmp_0
      * rtb_costheta_tmp_3 + rtb_sintheta_tmp_1 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_o[7] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_1 - rtb_costheta_tmp_3 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_o[8] = rtb_costheta_tmp_2 * rtb_costheta_tmp_1;

    // SignalConversion generated from: '<S107>/Product1' incorporates:
    //   Fcn: '<S159>/R11'

    rtb_TmpSignalConversionAtProd_n[0] = rtb_costheta_tmp_3 * rtb_costheta_tmp_2;

    // Fcn: '<S159>/R21 ' incorporates:
    //   Fcn: '<S159>/R31 '

    rtb_TmpSignalConversionAtProd_9 = rtb_costheta_tmp_3 * rtb_sintheta_tmp_0;

    // SignalConversion generated from: '<S107>/Product1' incorporates:
    //   Fcn: '<S159>/R12'
    //   Fcn: '<S159>/R21 '
    //   Fcn: '<S159>/R31 '

    rtb_TmpSignalConversionAtProd_n[1] = rtb_TmpSignalConversionAtProd_9 *
      rtb_TransferFcn + -rtb_sintheta_tmp_1 * rtb_costheta_tmp_1;
    rtb_TmpSignalConversionAtProd_n[2] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_1 + rtb_sintheta_tmp_1 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_n[3] = rtb_sintheta_tmp_1 * rtb_costheta_tmp_2;

    // Fcn: '<S159>/R22' incorporates:
    //   Fcn: '<S159>/R32'

    rtb_TmpSignalConversionAtProd_9 = rtb_sintheta_tmp_1 * rtb_sintheta_tmp_0;

    // SignalConversion generated from: '<S107>/Product1' incorporates:
    //   Fcn: '<S159>/R13'
    //   Fcn: '<S159>/R22'
    //   Fcn: '<S159>/R23'
    //   Fcn: '<S159>/R32'

    rtb_TmpSignalConversionAtProd_n[4] = rtb_TmpSignalConversionAtProd_9 *
      rtb_TransferFcn + rtb_costheta_tmp_3 * rtb_costheta_tmp_1;
    rtb_TmpSignalConversionAtProd_n[5] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_1 + -rtb_costheta_tmp_3 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_n[6] = -rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_n[7] = rtb_costheta_tmp_2 * rtb_TransferFcn;

    // Trigonometry: '<S215>/sin(theta)' incorporates:
    //   Fcn: '<S216>/T21 '
    //   Trigonometry: '<S255>/sin(theta)'

    rtb_sintheta_tmp_0 = std::sin(AHV_Model_B.Integrator1_p[3]);

    // Trigonometry: '<S215>/cos(theta)' incorporates:
    //   Fcn: '<S216>/T31 '
    //   Trigonometry: '<S255>/cos(theta)'

    rtb_costheta_tmp_4 = std::cos(AHV_Model_B.Integrator1_p[3]);

    // Trigonometry: '<S215>/sin(theta)' incorporates:
    //   Trigonometry: '<S255>/sin(theta)'

    rtb_sintheta_tmp_2 = std::sin(AHV_Model_B.Integrator1_p[4]);

    // Trigonometry: '<S215>/cos(theta)' incorporates:
    //   Fcn: '<S216>/T23'
    //   Trigonometry: '<S255>/cos(theta)'

    rtb_costheta_tmp_5 = std::cos(AHV_Model_B.Integrator1_p[4]);

    // Trigonometry: '<S215>/sin(theta)' incorporates:
    //   Fcn: '<S262>/Fcn'
    //   Fcn: '<S262>/Fcn1'
    //   Fcn: '<S263>/Row1'
    //   Fcn: '<S263>/Row2'
    //   Fcn: '<S264>/Row1'
    //   Trigonometry: '<S255>/sin(theta)'

    rtb_TmpSignalConversionAtProd_9 = std::sin(AHV_Model_B.Integrator1_p[5]);

    // Trigonometry: '<S215>/cos(theta)' incorporates:
    //   Fcn: '<S262>/Fcn'
    //   Fcn: '<S262>/Fcn1'
    //   Fcn: '<S263>/Row1'
    //   Fcn: '<S263>/Row2'
    //   Fcn: '<S264>/Row1'
    //   Fcn: '<S264>/Row2'
    //   Trigonometry: '<S255>/cos(theta)'

    rtb_costheta_tmp_6 = std::cos(AHV_Model_B.Integrator1_p[5]);

    // SignalConversion generated from: '<S107>/Product1' incorporates:
    //   Fcn: '<S159>/R33'

    rtb_TmpSignalConversionAtProd_n[8] = rtb_costheta_tmp_2 * rtb_costheta_tmp_1;

    // SignalConversion generated from: '<S212>/Product' incorporates:
    //   Fcn: '<S215>/R11'
    //   Fcn: '<S215>/R12'
    //   Fcn: '<S215>/R21 '
    //   Fcn: '<S215>/R31 '
    //   Trigonometry: '<S215>/cos(theta)'
    //   Trigonometry: '<S215>/sin(theta)'

    rtb_TmpSignalConversionAtProd_f[0] = rtb_costheta_tmp_6 * rtb_costheta_tmp_5;
    rtb_TmpSignalConversionAtProd_f[1] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_5;
    rtb_TmpSignalConversionAtProd_f[2] = -rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtProd_f[3] = rtb_costheta_tmp_6 * rtb_sintheta_tmp_2
      * rtb_sintheta_tmp_0 - rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_4;

    // Fcn: '<S215>/R22' incorporates:
    //   Fcn: '<S215>/R23'
    //   Trigonometry: '<S215>/sin(theta)'

    rtb_TmpSignalConversionAtProd_0 = rtb_TmpSignalConversionAtProd_9 *
      rtb_sintheta_tmp_2;

    // SignalConversion generated from: '<S212>/Product' incorporates:
    //   Fcn: '<S215>/R13'
    //   Fcn: '<S215>/R22'
    //   Fcn: '<S215>/R23'
    //   Fcn: '<S215>/R32'
    //   Fcn: '<S215>/R33'
    //   Trigonometry: '<S215>/cos(theta)'
    //   Trigonometry: '<S215>/sin(theta)'

    rtb_TmpSignalConversionAtProd_f[4] = rtb_TmpSignalConversionAtProd_0 *
      rtb_sintheta_tmp_0 + rtb_costheta_tmp_6 * rtb_costheta_tmp_4;
    rtb_TmpSignalConversionAtProd_f[5] = rtb_costheta_tmp_5 * rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_f[6] = rtb_costheta_tmp_4 * rtb_sintheta_tmp_2
      * rtb_costheta_tmp_6 + rtb_TmpSignalConversionAtProd_9 *
      rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_f[7] = rtb_TmpSignalConversionAtProd_0 *
      rtb_costheta_tmp_4 - rtb_costheta_tmp_6 * rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_f[8] = rtb_costheta_tmp_5 * rtb_costheta_tmp_4;

    // SignalConversion generated from: '<S203>/Product1' incorporates:
    //   Fcn: '<S255>/R11'

    rtb_TmpSignalConversionAtPro_fx[0] = rtb_costheta_tmp_6 * rtb_costheta_tmp_5;

    // Fcn: '<S255>/R21 ' incorporates:
    //   Fcn: '<S255>/R31 '

    rtb_TmpSignalConversionAtProd_0 = rtb_costheta_tmp_6 * rtb_sintheta_tmp_2;

    // SignalConversion generated from: '<S203>/Product1' incorporates:
    //   Fcn: '<S255>/R12'
    //   Fcn: '<S255>/R21 '
    //   Fcn: '<S255>/R31 '

    rtb_TmpSignalConversionAtPro_fx[1] = rtb_TmpSignalConversionAtProd_0 *
      rtb_sintheta_tmp_0 + -rtb_TmpSignalConversionAtProd_9 * rtb_costheta_tmp_4;
    rtb_TmpSignalConversionAtPro_fx[2] = rtb_TmpSignalConversionAtProd_0 *
      rtb_costheta_tmp_4 + rtb_TmpSignalConversionAtProd_9 * rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtPro_fx[3] = rtb_TmpSignalConversionAtProd_9 *
      rtb_costheta_tmp_5;

    // Fcn: '<S255>/R22' incorporates:
    //   Fcn: '<S255>/R32'

    rtb_TmpSignalConversionAtProd_0 = rtb_TmpSignalConversionAtProd_9 *
      rtb_sintheta_tmp_2;

    // SignalConversion generated from: '<S203>/Product1' incorporates:
    //   Fcn: '<S255>/R13'
    //   Fcn: '<S255>/R22'
    //   Fcn: '<S255>/R23'
    //   Fcn: '<S255>/R32'

    rtb_TmpSignalConversionAtPro_fx[4] = rtb_TmpSignalConversionAtProd_0 *
      rtb_sintheta_tmp_0 + rtb_costheta_tmp_6 * rtb_costheta_tmp_4;
    rtb_TmpSignalConversionAtPro_fx[5] = rtb_TmpSignalConversionAtProd_0 *
      rtb_costheta_tmp_4 + -rtb_costheta_tmp_6 * rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtPro_fx[6] = -rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtPro_fx[7] = rtb_costheta_tmp_5 * rtb_sintheta_tmp_0;

    // Trigonometry: '<S311>/sin(theta)' incorporates:
    //   Fcn: '<S312>/T21 '
    //   Trigonometry: '<S351>/sin(theta)'

    rtb_sintheta_tmp_2 = std::sin(AHV_Model_B.Integrator1_n[3]);

    // Trigonometry: '<S311>/cos(theta)' incorporates:
    //   Fcn: '<S312>/T31 '
    //   Trigonometry: '<S351>/cos(theta)'

    rtb_costheta_tmp_7 = std::cos(AHV_Model_B.Integrator1_n[3]);

    // Trigonometry: '<S311>/sin(theta)' incorporates:
    //   Trigonometry: '<S351>/sin(theta)'

    rtb_rxpi_h = std::sin(AHV_Model_B.Integrator1_n[4]);

    // Trigonometry: '<S311>/cos(theta)' incorporates:
    //   Fcn: '<S312>/T23'
    //   Trigonometry: '<S351>/cos(theta)'

    rtb_costheta_tmp_8 = std::cos(AHV_Model_B.Integrator1_n[4]);

    // Trigonometry: '<S311>/sin(theta)' incorporates:
    //   Fcn: '<S358>/Fcn'
    //   Fcn: '<S358>/Fcn1'
    //   Fcn: '<S359>/Row1'
    //   Fcn: '<S359>/Row2'
    //   Fcn: '<S360>/Row1'
    //   Trigonometry: '<S351>/sin(theta)'

    rtb_TmpSignalConversionAtProd_0 = std::sin(AHV_Model_B.Integrator1_n[5]);

    // Trigonometry: '<S311>/cos(theta)' incorporates:
    //   Fcn: '<S358>/Fcn'
    //   Fcn: '<S358>/Fcn1'
    //   Fcn: '<S359>/Row1'
    //   Fcn: '<S359>/Row2'
    //   Fcn: '<S360>/Row1'
    //   Fcn: '<S360>/Row2'
    //   Trigonometry: '<S351>/cos(theta)'

    rtb_costheta_tmp_9 = std::cos(AHV_Model_B.Integrator1_n[5]);

    // SignalConversion generated from: '<S203>/Product1' incorporates:
    //   Fcn: '<S255>/R33'

    rtb_TmpSignalConversionAtPro_fx[8] = rtb_costheta_tmp_5 * rtb_costheta_tmp_4;

    // SignalConversion generated from: '<S308>/Product' incorporates:
    //   Fcn: '<S311>/R11'
    //   Fcn: '<S311>/R12'
    //   Fcn: '<S311>/R21 '
    //   Fcn: '<S311>/R31 '
    //   Trigonometry: '<S311>/cos(theta)'
    //   Trigonometry: '<S311>/sin(theta)'

    rtb_TmpSignalConversionAtProd_j[0] = rtb_costheta_tmp_9 * rtb_costheta_tmp_8;
    rtb_TmpSignalConversionAtProd_j[1] = rtb_TmpSignalConversionAtProd_0 *
      rtb_costheta_tmp_8;
    rtb_TmpSignalConversionAtProd_j[2] = -rtb_rxpi_h;
    rtb_TmpSignalConversionAtProd_j[3] = rtb_costheta_tmp_9 * rtb_rxpi_h *
      rtb_sintheta_tmp_2 - rtb_TmpSignalConversionAtProd_0 * rtb_costheta_tmp_7;

    // Fcn: '<S311>/R22' incorporates:
    //   Fcn: '<S311>/R23'
    //   Trigonometry: '<S311>/sin(theta)'

    rtb_Product1_lq = rtb_TmpSignalConversionAtProd_0 * rtb_rxpi_h;

    // SignalConversion generated from: '<S308>/Product' incorporates:
    //   Fcn: '<S311>/R13'
    //   Fcn: '<S311>/R22'
    //   Fcn: '<S311>/R23'
    //   Fcn: '<S311>/R32'
    //   Fcn: '<S311>/R33'
    //   Trigonometry: '<S311>/cos(theta)'
    //   Trigonometry: '<S311>/sin(theta)'

    rtb_TmpSignalConversionAtProd_j[4] = rtb_Product1_lq * rtb_sintheta_tmp_2 +
      rtb_costheta_tmp_9 * rtb_costheta_tmp_7;
    rtb_TmpSignalConversionAtProd_j[5] = rtb_costheta_tmp_8 * rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtProd_j[6] = rtb_costheta_tmp_7 * rtb_rxpi_h *
      rtb_costheta_tmp_9 + rtb_TmpSignalConversionAtProd_0 * rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtProd_j[7] = rtb_Product1_lq * rtb_costheta_tmp_7 -
      rtb_costheta_tmp_9 * rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtProd_j[8] = rtb_costheta_tmp_8 * rtb_costheta_tmp_7;

    // Fcn: '<S351>/R31 '
    rtb_Product_g = rtb_costheta_tmp_9 * rtb_rxpi_h * rtb_costheta_tmp_7 +
      rtb_TmpSignalConversionAtProd_0 * rtb_sintheta_tmp_2;

    // Fcn: '<S351>/R12'
    rtb_Product1_lq = rtb_TmpSignalConversionAtProd_0 * rtb_costheta_tmp_8;

    // Fcn: '<S351>/R22'
    rtb_Product2_b = rtb_TmpSignalConversionAtProd_0 * rtb_rxpi_h *
      rtb_sintheta_tmp_2 + rtb_costheta_tmp_9 * rtb_costheta_tmp_7;

    // SignalConversion generated from: '<S299>/Product1' incorporates:
    //   Fcn: '<S351>/R11'
    //   Fcn: '<S351>/R13'
    //   Fcn: '<S351>/R21 '
    //   Fcn: '<S351>/R23'
    //   Fcn: '<S351>/R32'
    //   Fcn: '<S351>/R33'

    rtb_TmpSignalConversionAtPro_ni[0] = rtb_costheta_tmp_9 * rtb_costheta_tmp_8;
    rtb_TmpSignalConversionAtPro_ni[1] = rtb_costheta_tmp_9 * rtb_rxpi_h *
      rtb_sintheta_tmp_2 + -rtb_TmpSignalConversionAtProd_0 * rtb_costheta_tmp_7;
    rtb_TmpSignalConversionAtPro_ni[2] = rtb_Product_g;
    rtb_TmpSignalConversionAtPro_ni[3] = rtb_Product1_lq;
    rtb_TmpSignalConversionAtPro_ni[4] = rtb_Product2_b;
    rtb_TmpSignalConversionAtPro_ni[5] = rtb_TmpSignalConversionAtProd_0 *
      rtb_rxpi_h * rtb_costheta_tmp_7 + -rtb_costheta_tmp_9 * rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtPro_ni[6] = -rtb_rxpi_h;
    rtb_TmpSignalConversionAtPro_ni[7] = rtb_costheta_tmp_8 * rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtPro_ni[8] = rtb_costheta_tmp_8 * rtb_costheta_tmp_7;

    // Gain: '<S302>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_Product1_lf = 0.017453292519943295 * Current_direction;

    // Product: '<S302>/Product' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S302>/cos'

    rtb_rxpi_h = std::cos(rtb_Product1_lf) * Current_speed;

    // Product: '<S302>/Product1' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S302>/sin'

    rtb_Product1_lf = std::sin(rtb_Product1_lf) * Current_speed;

    // Sum: '<S303>/Sum6' incorporates:
    //   Gain: '<S298>/Current on//off '
    //   Integrator: '<S303>/Integrator'

    AHV_Model_B.nu_r[0] = AHV_Model_X.Integrator_CSTATE[0] - rtb_rxpi_h;
    AHV_Model_B.nu_r[1] = AHV_Model_X.Integrator_CSTATE[1] - rtb_Product1_lf;
    AHV_Model_B.nu_r[2] = AHV_Model_X.Integrator_CSTATE[2];
    AHV_Model_B.nu_r[3] = AHV_Model_X.Integrator_CSTATE[3];
    AHV_Model_B.nu_r[4] = AHV_Model_X.Integrator_CSTATE[4];
    AHV_Model_B.nu_r[5] = AHV_Model_X.Integrator_CSTATE[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S309>/Cross-flow drag trapezoidal integration' 
      // Constant: '<S309>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_f,
        AHV_Model_B.nu_r[1], AHV_Model_B.nu_r[5], &AHV_Model_B.Sum_p,
        &AHV_Model_B.Sum2, 82.800003);

      // End of Outputs for SubSystem: '<S309>/Cross-flow drag trapezoidal integration' 

      // Product: '<S309>/dx1' incorporates:
      //   Constant: '<S309>/2D drag coefficient '
      //   Constant: '<S309>/Transversal area//Lpp'
      //   Constant: '<S309>/rho'
      //   Gain: '<S309>/Gain1'

      AHV_Model_B.dx1 = -0.5 * AHV_Model_B.Sum_p * 1025.0 * 5.4 *
        0.69954741847648816;

      // Product: '<S309>/dx2' incorporates:
      //   Constant: '<S309>/2D drag coefficient '
      //   Constant: '<S309>/Transversal area//Lpp'
      //   Constant: '<S309>/rho'
      //   Gain: '<S309>/Gain'

      AHV_Model_B.dx2 = -0.5 * AHV_Model_B.Sum2 * 1025.0 * 5.4 *
        0.69954741847648816;
    }

    // If: '<S285>/If1' incorporates:
    //   Inport: '<Root>/Drving_Mode3'
    //   Inport: '<Root>/Hold_Position3'
    //   Inport: '<Root>/Thruster_percentage3'
    //   Inport: '<S292>/HoldPosition'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If1_ActiveSubsystem = static_cast<int8_T>((!(Drving_Mode3 ==
        0.0)) || (!(Thruster_percentage3 == 0.0)));
    }

    switch (AHV_Model_DW.If1_ActiveSubsystem) {
     case 0:
      // Outputs for IfAction SubSystem: '<S285>/If Action Subsystem2' incorporates:
      //   ActionPort: '<S291>/Action Port'

      AHV_Model_IfActionSubsystem2(&AHV_Model_B.Merge_pf);

      // End of Outputs for SubSystem: '<S285>/If Action Subsystem2'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S285>/If Action Subsystem3' incorporates:
      //   ActionPort: '<S292>/Action Port'

      AHV_Model_B.Merge_pf = Hold_Position3;

      // End of Outputs for SubSystem: '<S285>/If Action Subsystem3'
      break;
    }

    // End of If: '<S285>/If1'

    // If: '<S189>/If1' incorporates:
    //   Inport: '<Root>/Drving_Mode2'
    //   Inport: '<Root>/Hold_Position2'
    //   Inport: '<Root>/Thruster_percentage2'
    //   Inport: '<S196>/HoldPosition'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If1_ActiveSubsystem_p = static_cast<int8_T>((!(Drving_Mode2 ==
        0.0)) || (!(Thruster_percentage2 == 0.0)));
    }

    switch (AHV_Model_DW.If1_ActiveSubsystem_p) {
     case 0:
      // Outputs for IfAction SubSystem: '<S189>/If Action Subsystem2' incorporates:
      //   ActionPort: '<S195>/Action Port'

      AHV_Model_IfActionSubsystem2(&AHV_Model_B.Merge_l);

      // End of Outputs for SubSystem: '<S189>/If Action Subsystem2'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S189>/If Action Subsystem3' incorporates:
      //   ActionPort: '<S196>/Action Port'

      AHV_Model_B.Merge_l = Hold_Position2;

      // End of Outputs for SubSystem: '<S189>/If Action Subsystem3'
      break;
    }

    // End of If: '<S189>/If1'

    // If: '<S93>/If1' incorporates:
    //   Inport: '<Root>/Drving_Mode1'
    //   Inport: '<Root>/Hold_Position1'
    //   Inport: '<Root>/Thruster_percentage1'
    //   Inport: '<S100>/HoldPosition'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If1_ActiveSubsystem_g = static_cast<int8_T>((!(Drving_Mode1 ==
        0.0)) || (!(Thruster_percentage1 == 0.0)));
    }

    switch (AHV_Model_DW.If1_ActiveSubsystem_g) {
     case 0:
      // Outputs for IfAction SubSystem: '<S93>/If Action Subsystem2' incorporates:
      //   ActionPort: '<S99>/Action Port'

      AHV_Model_IfActionSubsystem2(&AHV_Model_B.Merge_g);

      // End of Outputs for SubSystem: '<S93>/If Action Subsystem2'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S93>/If Action Subsystem3' incorporates:
      //   ActionPort: '<S100>/Action Port'

      AHV_Model_B.Merge_g = Hold_Position1;

      // End of Outputs for SubSystem: '<S93>/If Action Subsystem3'
      break;
    }

    // End of If: '<S93>/If1'

    // Switch: '<S1>/Switch3' incorporates:
    //   Constant: '<S1>/Constant20'
    //   Gain: '<S1>/Gain'
    //   Inport: '<Root>/Hold_Position1'
    //   Inport: '<Root>/traget_speed1'
    //   Switch: '<S1>/Switch1'

    if (Hold_Position1) {
      rtb_Switch3 = 3.0;
    } else {
      if (traget_speed1 >= 0.2) {
        // Switch: '<S1>/Switch1' incorporates:
        //   Inport: '<Root>/traget_speed1'

        tmp_8 = traget_speed1;
      } else {
        // Switch: '<S1>/Switch1' incorporates:
        //   Constant: '<S1>/traget speed limit'

        tmp_8 = 0.2;
      }

      rtb_Switch3 = 5.0 * tmp_8;
    }

    // End of Switch: '<S1>/Switch3'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Switch: '<S9>/Switch2' incorporates:
      //   Abs: '<S9>/Abs'
      //   Constant: '<S9>/Constant18'
      //   Product: '<S9>/Divide'
      //   Sum: '<S9>/Add'
      //   UnitDelay: '<S9>/Unit Delay'

      if (AHV_Model_DW.UnitDelay_DSTATE[1] + AHV_Model_DW.UnitDelay_DSTATE[0] !=
          0.0) {
        rtb_Switch2 = std::abs(AHV_Model_DW.UnitDelay_DSTATE[1] /
          AHV_Model_DW.UnitDelay_DSTATE[0]);
      } else {
        rtb_Switch2 = 1.0;
      }

      // End of Switch: '<S9>/Switch2'

      // Switch: '<S9>/Switch' incorporates:
      //   Constant: '<S9>/Constant21'

      if (rtb_Switch2 > 1.0) {
        AHV_Model_B.Switch = 1.0;
      } else {
        AHV_Model_B.Switch = rtb_Switch2;
      }

      // End of Switch: '<S9>/Switch'
    }

    // Saturate: '<S1>/Sway Failure' incorporates:
    //   Inport: '<Root>/Rudde_ angle1'

    if (Rudde_angle1 > 35.0) {
      SwayFailure = 35.0;
    } else if (Rudde_angle1 < -35.0) {
      SwayFailure = -35.0;
    } else {
      SwayFailure = Rudde_angle1;
    }

    // End of Saturate: '<S1>/Sway Failure'

    // If: '<S8>/If' incorporates:
    //   Inport: '<Root>/Drving_Mode1'
    //   Inport: '<Root>/Vessel_X_Ref1'
    //   Inport: '<Root>/Vessel_Y_Ref1'
    //   Inport: '<Root>/heading_angle_ref1'
    //   Inport: '<Root>/heading_mode1'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem = static_cast<int8_T>(!(Drving_Mode1 ==
        0.0));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem) {
     case 0:
      // Outputs for IfAction SubSystem: '<S8>/If Action Subsystem' incorporates:
      //   ActionPort: '<S91>/Action Port'

      // RateLimiter: '<S94>/Rate Limiter1'
      if (AHV_Model_DW.LastMajorTime_e == (rtInf)) {
        AHV_Model_B.RateLimiter1 = SwayFailure;
      } else {
        rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_e;
        rtb_Product1_lf = rtb_Switch3 * 0.1;
        rtb_Switch2 = SwayFailure - AHV_Model_DW.PrevY_a;
        if (rtb_Switch2 > rtb_Product1_lf) {
          AHV_Model_B.RateLimiter1 = AHV_Model_DW.PrevY_a + rtb_Product1_lf;
        } else {
          rtb_Switch3 *= -0.1;
          if (rtb_Switch2 < rtb_Switch3) {
            AHV_Model_B.RateLimiter1 = AHV_Model_DW.PrevY_a + rtb_Switch3;
          } else {
            AHV_Model_B.RateLimiter1 = SwayFailure;
          }
        }
      }

      // End of RateLimiter: '<S94>/Rate Limiter1'

      // Gain: '<S94>/Gain4'
      AHV_Model_B.Gain4 = 57.324840764331206 * AHV_Model_B.Integrator1[5];
      if (rtmIsMajorTimeStep((&AHV_Model_M))) {
        // Chart: '<S94>/Chart2'
        AHV_Model_Chart1(AHV_Model_B.RateLimiter1, AHV_Model_B.Gain4,
                         &AHV_Model_B.heading_deg_g, &AHV_Model_DW.sf_Chart2);
      }

      // RateLimiter: '<S94>/Rate Limiter2'
      if (AHV_Model_DW.LastMajorTime_h == (rtInf)) {
        AHV_Model_B.RateLimiter2 = AHV_Model_B.heading_deg_g;
      } else {
        rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_h;
        rtb_Product1_lf = rtb_Switch3 * 0.1;
        rtb_Switch2 = AHV_Model_B.heading_deg_g - AHV_Model_DW.PrevY_ic;
        if (rtb_Switch2 > rtb_Product1_lf) {
          AHV_Model_B.RateLimiter2 = AHV_Model_DW.PrevY_ic + rtb_Product1_lf;
        } else {
          rtb_Switch3 *= -0.1;
          if (rtb_Switch2 < rtb_Switch3) {
            AHV_Model_B.RateLimiter2 = AHV_Model_DW.PrevY_ic + rtb_Switch3;
          } else {
            AHV_Model_B.RateLimiter2 = AHV_Model_B.heading_deg_g;
          }
        }
      }

      // End of RateLimiter: '<S94>/Rate Limiter2'

      // Gain: '<S94>/Gain2'
      rtb_rxpi_h = 0.017444444444444446 * AHV_Model_B.RateLimiter2;

      // Gain: '<S94>/Gain3' incorporates:
      //   Inport: '<Root>/Thruster_percentage1'

      rtb_Product_g = 50.0 * Thruster_percentage1;

      // Product: '<S94>/Product' incorporates:
      //   Trigonometry: '<S94>/x'

      AHV_Model_B.Product = std::cos(rtb_rxpi_h) * rtb_Product_g;
      if (rtmIsMajorTimeStep((&AHV_Model_M))) {
        // SignalConversion generated from: '<S91>/Target_x_head' incorporates:
        //   Chart: '<S94>/Chart3'

        AHV_Model_Chart3(AHV_Model_B.Integrator1[0], AHV_Model_B.Product,
                         &AHV_Model_B.Merge2, &AHV_Model_DW.sf_Chart3);
      }

      // Product: '<S94>/Product1' incorporates:
      //   Trigonometry: '<S94>/y'

      AHV_Model_B.Product1 = rtb_Product_g * std::sin(rtb_rxpi_h);
      if (rtmIsMajorTimeStep((&AHV_Model_M))) {
        // SignalConversion generated from: '<S91>/Target_y_head' incorporates:
        //   Chart: '<S94>/Chart1'

        AHV_Model_Chart1(AHV_Model_B.Product1, AHV_Model_B.Integrator1[1],
                         &AHV_Model_B.Merge5, &AHV_Model_DW.sf_Chart1);
      }

      // Switch: '<S94>/Switch1' incorporates:
      //   Constant: '<S94>/Constant'
      //   Constant: '<S94>/Constant1'
      //   Switch: '<S94>/Switch'

      if (AHV_Model_B.Merge_g) {
        rtb_Merge4_g = 10.0;
        rtb_Merge3_c = 1.0;
      } else {
        rtb_Merge4_g = AHV_Model_B.Product1;
        rtb_Merge3_c = AHV_Model_B.Product;
      }

      // End of Switch: '<S94>/Switch1'

      // Gain: '<S94>/Gain5'
      rtb_rxpi_h *= 57.324840764331206;

      // Product: '<S98>/Divide' incorporates:
      //   Constant: '<S98>/Constant3'

      rtb_Product_g = rtb_rxpi_h / 360.0;

      // Rounding: '<S98>/Floor'
      if (rtb_Product_g < 0.0) {
        rtb_Product_g = std::ceil(rtb_Product_g);
      } else {
        rtb_Product_g = std::floor(rtb_Product_g);
      }

      // End of Rounding: '<S98>/Floor'

      // Sum: '<S98>/Sum2' incorporates:
      //   Constant: '<S98>/Constant3'
      //   Product: '<S98>/Product2'

      rtb_Product_g = rtb_rxpi_h - rtb_Product_g * 360.0;
      if (rtmIsMajorTimeStep((&AHV_Model_M))) {
        // SignalConversion generated from: '<S91>/SurgeForce' incorporates:
        //   Constant: '<S91>/Surge Force'

        AHV_Model_B.Merge7_e = 1.0;

        // SignalConversion generated from: '<S91>/SwayForce' incorporates:
        //   Constant: '<S91>/Sway Force'

        AHV_Model_B.Merge8_mo = 1.0;

        // SignalConversion generated from: '<S91>/YawMoment' incorporates:
        //   Constant: '<S91>/Yaw Moment'

        AHV_Model_B.Merge9_p = 1.0;
      }

      // Gain: '<S91>/Gain' incorporates:
      //   Inport: '<Root>/heading_mode1'

      rtb_rxpi_h = 0.0 * heading_mode1;

      // End of Outputs for SubSystem: '<S8>/If Action Subsystem'
      break;

     case 1:
      // Switch: '<S1>/Switch4' incorporates:
      //   Constant: '<S1>/Constant22'
      //   Inport: '<Root>/Hold_Position1'

      if (Hold_Position1) {
        tmp_8 = 3.0;
      } else {
        tmp_8 = AHV_Model_B.Switch;
      }

      // End of Switch: '<S1>/Switch4'

      // Outputs for IfAction SubSystem: '<S8>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S92>/Action Port'

      AHV_Model_IfActionSubsystem1_k(heading_mode1, heading_angle_ref1,
        Vessel_X_Ref1, rtb_Switch3, tmp_8, Vessel_Y_Ref1, &rtb_Product_g,
        &AHV_Model_B.Merge2, &rtb_Merge3_c, &rtb_Merge4_g, &AHV_Model_B.Merge5,
        &rtb_rxpi_h, &AHV_Model_B.Merge7_e, &AHV_Model_B.Merge8_mo,
        &AHV_Model_B.Merge9_p);

      // End of Outputs for SubSystem: '<S8>/If Action Subsystem1'
      break;
    }

    // End of If: '<S8>/If'

    // Switch: '<S66>/Switch' incorporates:
    //   Gain: '<S69>/Gain'
    //   Sum: '<S64>/Sum3'
    //   Sum: '<S64>/Sum4'
    //   Trigonometry: '<S64>/Trigonometric Function'

    if (rtb_rxpi_h > 0.5) {
      rtb_Product_g = 57.295779513082323 * rt_atan2d_snf
        (AHV_Model_B.Integrator1[1], AHV_Model_B.Integrator1[0]);
    }

    // End of Switch: '<S66>/Switch'

    // RateLimiter: '<S66>/Rate Limiter'
    if (AHV_Model_DW.LastMajorTime == (rtInf)) {
      AHV_Model_B.RateLimiter = rtb_Product_g;
    } else {
      rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime;
      rtb_Product1_lf = rtb_Switch3 * 10.0;
      rtb_Switch2 = rtb_Product_g - AHV_Model_DW.PrevY;
      if (rtb_Switch2 > rtb_Product1_lf) {
        AHV_Model_B.RateLimiter = AHV_Model_DW.PrevY + rtb_Product1_lf;
      } else {
        rtb_Switch3 *= -10.0;
        if (rtb_Switch2 < rtb_Switch3) {
          AHV_Model_B.RateLimiter = AHV_Model_DW.PrevY + rtb_Switch3;
        } else {
          AHV_Model_B.RateLimiter = rtb_Product_g;
        }
      }
    }

    // End of RateLimiter: '<S66>/Rate Limiter'

    // Gain: '<S76>/Gain1'
    rtb_rxpi_h = 0.017453292519943295 * AHV_Model_B.RateLimiter;

    // Saturate: '<S77>/x_Saturation'
    if (rtb_rxpi_h > 1.0E+10) {
      rtb_rxpi_h = 1.0E+10;
    } else {
      if (rtb_rxpi_h < -1.0E+10) {
        rtb_rxpi_h = -1.0E+10;
      }
    }

    // End of Saturate: '<S77>/x_Saturation'

    // Signum: '<S77>/x_Sign'
    if (rtb_rxpi_h < 0.0) {
      rtb_Product_g = -1.0;
    } else if (rtb_rxpi_h > 0.0) {
      rtb_Product_g = 1.0;
    } else if (rtb_rxpi_h == 0.0) {
      rtb_Product_g = 0.0;
    } else {
      rtb_Product_g = (rtNaN);
    }

    // End of Signum: '<S77>/x_Sign'

    // Gain: '<S77>/pi'
    rtb_Product_g *= 3.1415926535897931;

    // Sum: '<S77>/Sum1'
    rtb_rxpi_h += rtb_Product_g;

    // Math: '<S77>/Math Function' incorporates:
    //   Constant: '<S77>/Constant'

    rtb_rxpi_h = rt_remd_snf(rtb_rxpi_h, 6.2831853071795862);
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S66>/Chart' incorporates:
      //   Sum: '<S77>/Sum'

      AHV_Model_Chart(AHV_Model_B.Merge_g, AHV_Model_B.Merge2,
                      AHV_Model_B.Merge5, AHV_Model_B.Integrator1[0],
                      AHV_Model_B.Integrator1[1], rtb_rxpi_h - rtb_Product_g,
                      AHV_Model_B.Integrator1[5], &AHV_Model_B.x_ref_rel_m,
                      &AHV_Model_B.y_ref_rel_d, &AHV_Model_B.yaw_ref_rel_c,
                      &AHV_Model_DW.sf_Chart);

      // Switch: '<S105>/Switch2' incorporates:
      //   Abs: '<S105>/Abs'
      //   Constant: '<S105>/Constant18'
      //   Product: '<S105>/Divide'
      //   Sum: '<S105>/Add'
      //   UnitDelay: '<S105>/Unit Delay'

      if (AHV_Model_DW.UnitDelay_DSTATE_h[1] + AHV_Model_DW.UnitDelay_DSTATE_h[0]
          != 0.0) {
        rtb_Switch2 = std::abs(AHV_Model_DW.UnitDelay_DSTATE_h[1] /
          AHV_Model_DW.UnitDelay_DSTATE_h[0]);
      } else {
        rtb_Switch2 = 1.0;
      }

      // End of Switch: '<S105>/Switch2'

      // Switch: '<S105>/Switch' incorporates:
      //   Constant: '<S105>/Constant21'

      if (rtb_Switch2 > 1.0) {
        AHV_Model_B.Switch_a = 1.0;
      } else {
        AHV_Model_B.Switch_a = rtb_Switch2;
      }

      // End of Switch: '<S105>/Switch'
    }

    // Switch: '<S2>/Switch3' incorporates:
    //   Constant: '<S2>/Constant20'
    //   Gain: '<S2>/Gain'
    //   Inport: '<Root>/Hold_Position2'
    //   Inport: '<Root>/traget_speed2'
    //   Switch: '<S2>/Switch1'

    if (Hold_Position2) {
      rtb_Switch3 = 3.0;
    } else {
      if (traget_speed2 >= 0.2) {
        // Switch: '<S2>/Switch1' incorporates:
        //   Inport: '<Root>/traget_speed2'

        tmp_8 = traget_speed2;
      } else {
        // Switch: '<S2>/Switch1' incorporates:
        //   Constant: '<S2>/traget speed limit'

        tmp_8 = 0.2;
      }

      rtb_Switch3 = 5.0 * tmp_8;
    }

    // End of Switch: '<S2>/Switch3'

    // If: '<S104>/If' incorporates:
    //   Inport: '<Root>/Drving_Mode2'
    //   Inport: '<Root>/Thruster_percentage2'
    //   Inport: '<Root>/Vessel_X_Ref2'
    //   Inport: '<Root>/Vessel_Y_Ref2'
    //   Inport: '<Root>/heading_angle_ref2'
    //   Inport: '<Root>/heading_mode2'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_l = static_cast<int8_T>(!(Drving_Mode2 ==
        0.0));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_l) {
     case 0:
      // Saturate: '<S2>/Sway Failure' incorporates:
      //   Inport: '<Root>/Rudde_ angle2'

      if (Rudde_angle2 > 35.0) {
        tmp_8 = 35.0;
      } else if (Rudde_angle2 < -35.0) {
        tmp_8 = -35.0;
      } else {
        tmp_8 = Rudde_angle2;
      }

      // End of Saturate: '<S2>/Sway Failure'

      // Outputs for IfAction SubSystem: '<S104>/If Action Subsystem' incorporates:
      //   ActionPort: '<S187>/Action Port'

      AHV_Model_IfActionSubsystem_c(AHV_Model_B.Integrator1_b, tmp_8,
        Thruster_percentage2, heading_mode2, AHV_Model_B.Merge_l, &rtb_Product_g,
        &AHV_Model_B.Merge2_i, &rtb_Merge3_j, &rtb_Merge4_l,
        &AHV_Model_B.Merge5_g, &rtb_rxpi_h, &AHV_Model_B.Merge7_m,
        &AHV_Model_B.Merge8_m, &AHV_Model_B.Merge9_l,
        &AHV_Model_B.IfActionSubsystem_c, &AHV_Model_DW.IfActionSubsystem_c);

      // End of Outputs for SubSystem: '<S104>/If Action Subsystem'
      break;

     case 1:
      // Switch: '<S2>/Switch4' incorporates:
      //   Constant: '<S2>/Constant22'
      //   Inport: '<Root>/Hold_Position2'

      if (Hold_Position2) {
        tmp_8 = 3.0;
      } else {
        tmp_8 = AHV_Model_B.Switch_a;
      }

      // End of Switch: '<S2>/Switch4'

      // Outputs for IfAction SubSystem: '<S104>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S188>/Action Port'

      AHV_Model_IfActionSubsystem1_k(heading_mode2, heading_angle_ref2,
        Vessel_X_Ref2, rtb_Switch3, tmp_8, Vessel_Y_Ref2, &rtb_Product_g,
        &AHV_Model_B.Merge2_i, &rtb_Merge3_j, &rtb_Merge4_l,
        &AHV_Model_B.Merge5_g, &rtb_rxpi_h, &AHV_Model_B.Merge7_m,
        &AHV_Model_B.Merge8_m, &AHV_Model_B.Merge9_l);

      // End of Outputs for SubSystem: '<S104>/If Action Subsystem1'
      break;
    }

    // End of If: '<S104>/If'

    // Switch: '<S162>/Switch' incorporates:
    //   Gain: '<S165>/Gain'
    //   Sum: '<S160>/Sum3'
    //   Sum: '<S160>/Sum4'
    //   Trigonometry: '<S160>/Trigonometric Function'

    if (rtb_rxpi_h > 0.5) {
      rtb_Product_g = 57.295779513082323 * rt_atan2d_snf
        (AHV_Model_B.Integrator1_b[1], AHV_Model_B.Integrator1_b[0]);
    }

    // End of Switch: '<S162>/Switch'

    // RateLimiter: '<S162>/Rate Limiter'
    if (AHV_Model_DW.LastMajorTime_g == (rtInf)) {
      AHV_Model_B.RateLimiter_b = rtb_Product_g;
    } else {
      rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_g;
      rtb_Product1_lf = rtb_Switch3 * 10.0;
      rtb_Switch2 = rtb_Product_g - AHV_Model_DW.PrevY_p;
      if (rtb_Switch2 > rtb_Product1_lf) {
        AHV_Model_B.RateLimiter_b = AHV_Model_DW.PrevY_p + rtb_Product1_lf;
      } else {
        rtb_Switch3 *= -10.0;
        if (rtb_Switch2 < rtb_Switch3) {
          AHV_Model_B.RateLimiter_b = AHV_Model_DW.PrevY_p + rtb_Switch3;
        } else {
          AHV_Model_B.RateLimiter_b = rtb_Product_g;
        }
      }
    }

    // End of RateLimiter: '<S162>/Rate Limiter'

    // Gain: '<S172>/Gain1'
    rtb_rxpi_h = 0.017453292519943295 * AHV_Model_B.RateLimiter_b;

    // Saturate: '<S173>/x_Saturation'
    if (rtb_rxpi_h > 1.0E+10) {
      rtb_rxpi_h = 1.0E+10;
    } else {
      if (rtb_rxpi_h < -1.0E+10) {
        rtb_rxpi_h = -1.0E+10;
      }
    }

    // End of Saturate: '<S173>/x_Saturation'

    // Signum: '<S173>/x_Sign'
    if (rtb_rxpi_h < 0.0) {
      rtb_Product_g = -1.0;
    } else if (rtb_rxpi_h > 0.0) {
      rtb_Product_g = 1.0;
    } else if (rtb_rxpi_h == 0.0) {
      rtb_Product_g = 0.0;
    } else {
      rtb_Product_g = (rtNaN);
    }

    // End of Signum: '<S173>/x_Sign'

    // Gain: '<S173>/pi'
    rtb_Product_g *= 3.1415926535897931;

    // Sum: '<S173>/Sum1'
    rtb_rxpi_h += rtb_Product_g;

    // Math: '<S173>/Math Function' incorporates:
    //   Constant: '<S173>/Constant'

    rtb_rxpi_h = rt_remd_snf(rtb_rxpi_h, 6.2831853071795862);
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S162>/Chart' incorporates:
      //   Sum: '<S173>/Sum'

      AHV_Model_Chart(AHV_Model_B.Merge_l, AHV_Model_B.Merge2_i,
                      AHV_Model_B.Merge5_g, AHV_Model_B.Integrator1_b[0],
                      AHV_Model_B.Integrator1_b[1], rtb_rxpi_h - rtb_Product_g,
                      AHV_Model_B.Integrator1_b[5], &AHV_Model_B.x_ref_rel_l,
                      &AHV_Model_B.y_ref_rel_g, &AHV_Model_B.yaw_ref_rel_g,
                      &AHV_Model_DW.sf_Chart_e);

      // Switch: '<S201>/Switch2' incorporates:
      //   Abs: '<S201>/Abs'
      //   Constant: '<S201>/Constant18'
      //   Product: '<S201>/Divide'
      //   Sum: '<S201>/Add'
      //   UnitDelay: '<S201>/Unit Delay'

      if (AHV_Model_DW.UnitDelay_DSTATE_c[1] + AHV_Model_DW.UnitDelay_DSTATE_c[0]
          != 0.0) {
        rtb_Switch2 = std::abs(AHV_Model_DW.UnitDelay_DSTATE_c[1] /
          AHV_Model_DW.UnitDelay_DSTATE_c[0]);
      } else {
        rtb_Switch2 = 1.0;
      }

      // End of Switch: '<S201>/Switch2'

      // Switch: '<S201>/Switch' incorporates:
      //   Constant: '<S201>/Constant21'

      if (rtb_Switch2 > 1.0) {
        AHV_Model_B.Switch_b = 1.0;
      } else {
        AHV_Model_B.Switch_b = rtb_Switch2;
      }

      // End of Switch: '<S201>/Switch'
    }

    // Switch: '<S3>/Switch3' incorporates:
    //   Constant: '<S3>/Constant20'
    //   Gain: '<S3>/Gain'
    //   Inport: '<Root>/Hold_Position3'
    //   Inport: '<Root>/traget_speed3'
    //   Switch: '<S3>/Switch1'

    if (Hold_Position3) {
      rtb_Switch3 = 3.0;
    } else {
      if (traget_speed3 >= 0.2) {
        // Switch: '<S3>/Switch1' incorporates:
        //   Inport: '<Root>/traget_speed3'

        tmp_8 = traget_speed3;
      } else {
        // Switch: '<S3>/Switch1' incorporates:
        //   Constant: '<S3>/traget speed limit'

        tmp_8 = 0.2;
      }

      rtb_Switch3 = 5.0 * tmp_8;
    }

    // End of Switch: '<S3>/Switch3'

    // If: '<S200>/If' incorporates:
    //   Inport: '<Root>/Drving_Mode3'
    //   Inport: '<Root>/Thruster_percentage3'
    //   Inport: '<Root>/Vessel_X_Ref3'
    //   Inport: '<Root>/Vessel_Y_Ref3'
    //   Inport: '<Root>/heading_angle_ref3'
    //   Inport: '<Root>/heading_mode3'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_a = static_cast<int8_T>(!(Drving_Mode3 ==
        0.0));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_a) {
     case 0:
      // Saturate: '<S3>/Sway Failure' incorporates:
      //   Inport: '<Root>/Rudde_ angle3'

      if (Rudde_angle3 > 35.0) {
        tmp_8 = 35.0;
      } else if (Rudde_angle3 < -35.0) {
        tmp_8 = -35.0;
      } else {
        tmp_8 = Rudde_angle3;
      }

      // End of Saturate: '<S3>/Sway Failure'

      // Outputs for IfAction SubSystem: '<S200>/If Action Subsystem' incorporates:
      //   ActionPort: '<S283>/Action Port'

      AHV_Model_IfActionSubsystem_c(AHV_Model_B.Integrator1_p, tmp_8,
        Thruster_percentage3, heading_mode3, AHV_Model_B.Merge_pf,
        &rtb_Product_g, &AHV_Model_B.Merge2_h, &rtb_Merge3_k, &rtb_Merge4_j,
        &AHV_Model_B.Merge5_l, &rtb_rxpi_h, &AHV_Model_B.Merge7_h,
        &AHV_Model_B.Merge8_j, &AHV_Model_B.Merge9_e,
        &AHV_Model_B.IfActionSubsystem_e, &AHV_Model_DW.IfActionSubsystem_e);

      // End of Outputs for SubSystem: '<S200>/If Action Subsystem'
      break;

     case 1:
      // Switch: '<S3>/Switch4' incorporates:
      //   Constant: '<S3>/Constant22'
      //   Inport: '<Root>/Hold_Position3'

      if (Hold_Position3) {
        tmp_8 = 3.0;
      } else {
        tmp_8 = AHV_Model_B.Switch_b;
      }

      // End of Switch: '<S3>/Switch4'

      // Outputs for IfAction SubSystem: '<S200>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S284>/Action Port'

      AHV_Model_IfActionSubsystem1_k(heading_mode3, heading_angle_ref3,
        Vessel_X_Ref3, rtb_Switch3, tmp_8, Vessel_Y_Ref3, &rtb_Product_g,
        &AHV_Model_B.Merge2_h, &rtb_Merge3_k, &rtb_Merge4_j,
        &AHV_Model_B.Merge5_l, &rtb_rxpi_h, &AHV_Model_B.Merge7_h,
        &AHV_Model_B.Merge8_j, &AHV_Model_B.Merge9_e);

      // End of Outputs for SubSystem: '<S200>/If Action Subsystem1'
      break;
    }

    // End of If: '<S200>/If'

    // Switch: '<S258>/Switch' incorporates:
    //   Gain: '<S261>/Gain'
    //   Sum: '<S256>/Sum3'
    //   Sum: '<S256>/Sum4'
    //   Trigonometry: '<S256>/Trigonometric Function'

    if (rtb_rxpi_h > 0.5) {
      rtb_Product_g = 57.295779513082323 * rt_atan2d_snf
        (AHV_Model_B.Integrator1_p[1], AHV_Model_B.Integrator1_p[0]);
    }

    // End of Switch: '<S258>/Switch'

    // RateLimiter: '<S258>/Rate Limiter'
    if (AHV_Model_DW.LastMajorTime_a == (rtInf)) {
      AHV_Model_B.RateLimiter_c = rtb_Product_g;
    } else {
      rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_a;
      rtb_Product1_lf = rtb_Switch3 * 10.0;
      rtb_Switch2 = rtb_Product_g - AHV_Model_DW.PrevY_i;
      if (rtb_Switch2 > rtb_Product1_lf) {
        AHV_Model_B.RateLimiter_c = AHV_Model_DW.PrevY_i + rtb_Product1_lf;
      } else {
        rtb_Switch3 *= -10.0;
        if (rtb_Switch2 < rtb_Switch3) {
          AHV_Model_B.RateLimiter_c = AHV_Model_DW.PrevY_i + rtb_Switch3;
        } else {
          AHV_Model_B.RateLimiter_c = rtb_Product_g;
        }
      }
    }

    // End of RateLimiter: '<S258>/Rate Limiter'

    // Gain: '<S268>/Gain1'
    rtb_rxpi_h = 0.017453292519943295 * AHV_Model_B.RateLimiter_c;

    // Saturate: '<S269>/x_Saturation'
    if (rtb_rxpi_h > 1.0E+10) {
      rtb_rxpi_h = 1.0E+10;
    } else {
      if (rtb_rxpi_h < -1.0E+10) {
        rtb_rxpi_h = -1.0E+10;
      }
    }

    // End of Saturate: '<S269>/x_Saturation'

    // Signum: '<S269>/x_Sign'
    if (rtb_rxpi_h < 0.0) {
      rtb_Product_g = -1.0;
    } else if (rtb_rxpi_h > 0.0) {
      rtb_Product_g = 1.0;
    } else if (rtb_rxpi_h == 0.0) {
      rtb_Product_g = 0.0;
    } else {
      rtb_Product_g = (rtNaN);
    }

    // End of Signum: '<S269>/x_Sign'

    // Gain: '<S269>/pi'
    rtb_Product_g *= 3.1415926535897931;

    // Sum: '<S269>/Sum1'
    rtb_rxpi_h += rtb_Product_g;

    // Math: '<S269>/Math Function' incorporates:
    //   Constant: '<S269>/Constant'

    rtb_rxpi_h = rt_remd_snf(rtb_rxpi_h, 6.2831853071795862);
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S258>/Chart' incorporates:
      //   Sum: '<S269>/Sum'

      AHV_Model_Chart(AHV_Model_B.Merge_pf, AHV_Model_B.Merge2_h,
                      AHV_Model_B.Merge5_l, AHV_Model_B.Integrator1_p[0],
                      AHV_Model_B.Integrator1_p[1], rtb_rxpi_h - rtb_Product_g,
                      AHV_Model_B.Integrator1_p[5], &AHV_Model_B.x_ref_rel_p,
                      &AHV_Model_B.y_ref_rel_c, &AHV_Model_B.yaw_ref_rel_m,
                      &AHV_Model_DW.sf_Chart_g);

      // Switch: '<S297>/Switch2' incorporates:
      //   Abs: '<S297>/Abs'
      //   Constant: '<S297>/Constant18'
      //   Product: '<S297>/Divide'
      //   Sum: '<S297>/Add'
      //   UnitDelay: '<S297>/Unit Delay'

      if (AHV_Model_DW.UnitDelay_DSTATE_a[1] + AHV_Model_DW.UnitDelay_DSTATE_a[0]
          != 0.0) {
        rtb_Switch2 = std::abs(AHV_Model_DW.UnitDelay_DSTATE_a[1] /
          AHV_Model_DW.UnitDelay_DSTATE_a[0]);
      } else {
        rtb_Switch2 = 1.0;
      }

      // End of Switch: '<S297>/Switch2'

      // Switch: '<S297>/Switch' incorporates:
      //   Constant: '<S297>/Constant21'

      if (rtb_Switch2 > 1.0) {
        AHV_Model_B.Switch_p = 1.0;
      } else {
        AHV_Model_B.Switch_p = rtb_Switch2;
      }

      // End of Switch: '<S297>/Switch'
    }

    // Switch: '<S4>/Switch3' incorporates:
    //   Constant: '<S4>/Constant20'
    //   Gain: '<S4>/Gain'
    //   Inport: '<Root>/Hold_Position4'
    //   Inport: '<Root>/traget_speed4'
    //   Switch: '<S4>/Switch1'

    if (Hold_Position4) {
      rtb_rxpi_h = 3.0;
    } else {
      if (traget_speed4 >= 0.2) {
        // Switch: '<S4>/Switch1' incorporates:
        //   Inport: '<Root>/traget_speed4'

        tmp_8 = traget_speed4;
      } else {
        // Switch: '<S4>/Switch1' incorporates:
        //   Constant: '<S4>/traget speed limit'

        tmp_8 = 0.2;
      }

      rtb_rxpi_h = 5.0 * tmp_8;
    }

    // End of Switch: '<S4>/Switch3'

    // If: '<S381>/If1' incorporates:
    //   Inport: '<Root>/Drving_Mode4'
    //   Inport: '<Root>/Hold_Position4'
    //   Inport: '<Root>/Thruster_percentage4'
    //   Inport: '<S388>/HoldPosition'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If1_ActiveSubsystem_l = static_cast<int8_T>((!(Drving_Mode4 ==
        0.0)) || (!(Thruster_percentage4 == 0.0)));
    }

    switch (AHV_Model_DW.If1_ActiveSubsystem_l) {
     case 0:
      // Outputs for IfAction SubSystem: '<S381>/If Action Subsystem2' incorporates:
      //   ActionPort: '<S387>/Action Port'

      AHV_Model_IfActionSubsystem2(&AHV_Model_B.Merge_ik);

      // End of Outputs for SubSystem: '<S381>/If Action Subsystem2'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S381>/If Action Subsystem3' incorporates:
      //   ActionPort: '<S388>/Action Port'

      AHV_Model_B.Merge_ik = Hold_Position4;

      // End of Outputs for SubSystem: '<S381>/If Action Subsystem3'
      break;
    }

    // End of If: '<S381>/If1'

    // If: '<S296>/If' incorporates:
    //   Inport: '<Root>/Drving_Mode4'
    //   Inport: '<Root>/Thruster_percentage4'
    //   Inport: '<Root>/Vessel_X_Ref4'
    //   Inport: '<Root>/Vessel_Y_Ref4'
    //   Inport: '<Root>/heading_angle_ref4'
    //   Inport: '<Root>/heading_mode4'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_k = static_cast<int8_T>(!(Drving_Mode4 ==
        0.0));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_k) {
     case 0:
      // Saturate: '<S4>/Sway Failure' incorporates:
      //   Inport: '<Root>/Rudde_ angle4'

      if (Rudde_angle4 > 35.0) {
        tmp_8 = 35.0;
      } else if (Rudde_angle4 < -35.0) {
        tmp_8 = -35.0;
      } else {
        tmp_8 = Rudde_angle4;
      }

      // End of Saturate: '<S4>/Sway Failure'

      // Outputs for IfAction SubSystem: '<S296>/If Action Subsystem' incorporates:
      //   ActionPort: '<S379>/Action Port'

      AHV_Model_IfActionSubsystem_c(AHV_Model_B.Integrator1_n, tmp_8,
        Thruster_percentage4, heading_mode4, AHV_Model_B.Merge_ik,
        &rtb_Product2_b, &AHV_Model_B.Merge2_j, &rtb_Merge3, &rtb_Merge4,
        &AHV_Model_B.Merge5_o, &rtb_Product1_lq, &AHV_Model_B.Merge7,
        &AHV_Model_B.Merge8, &AHV_Model_B.Merge9,
        &AHV_Model_B.IfActionSubsystem_eo, &AHV_Model_DW.IfActionSubsystem_eo);

      // End of Outputs for SubSystem: '<S296>/If Action Subsystem'
      break;

     case 1:
      // Switch: '<S4>/Switch4' incorporates:
      //   Constant: '<S4>/Constant22'
      //   Inport: '<Root>/Hold_Position4'

      if (Hold_Position4) {
        tmp_8 = 3.0;
      } else {
        tmp_8 = AHV_Model_B.Switch_p;
      }

      // End of Switch: '<S4>/Switch4'

      // Outputs for IfAction SubSystem: '<S296>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S380>/Action Port'

      AHV_Model_IfActionSubsystem1_k(heading_mode4, heading_angle_ref4,
        Vessel_X_Ref4, rtb_rxpi_h, tmp_8, Vessel_Y_Ref4, &rtb_Product2_b,
        &AHV_Model_B.Merge2_j, &rtb_Merge3, &rtb_Merge4, &AHV_Model_B.Merge5_o,
        &rtb_Product1_lq, &AHV_Model_B.Merge7, &AHV_Model_B.Merge8,
        &AHV_Model_B.Merge9);

      // End of Outputs for SubSystem: '<S296>/If Action Subsystem1'
      break;
    }

    // End of If: '<S296>/If'

    // Integrator: '<S353>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init4'
    //   SignalConversion generated from: '<S353>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK != 0) {
      AHV_Model_X.Integrator3_CSTATE[0] = Vessel_init4[0];
      AHV_Model_X.Integrator3_CSTATE[1] = Vessel_init4[1];
      AHV_Model_X.Integrator3_CSTATE[2] = Vessel_init4[5];
    }

    // Saturate: '<S361>/x_Saturation' incorporates:
    //   Integrator: '<S353>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE[2] > 1.0E+10) {
      rtb_rxpi_h = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE[2] < -1.0E+10) {
      rtb_rxpi_h = -1.0E+10;
    } else {
      rtb_rxpi_h = AHV_Model_X.Integrator3_CSTATE[2];
    }

    // End of Saturate: '<S361>/x_Saturation'

    // Signum: '<S361>/x_Sign'
    if (rtb_rxpi_h < 0.0) {
      rtb_Product_g = -1.0;
    } else if (rtb_rxpi_h > 0.0) {
      rtb_Product_g = 1.0;
    } else if (rtb_rxpi_h == 0.0) {
      rtb_Product_g = 0.0;
    } else {
      rtb_Product_g = (rtNaN);
    }

    // End of Signum: '<S361>/x_Sign'

    // Gain: '<S361>/pi'
    rtb_Product_g *= 3.1415926535897931;

    // Sum: '<S361>/Sum1'
    rtb_rxpi_h += rtb_Product_g;

    // Math: '<S361>/Math Function' incorporates:
    //   Constant: '<S361>/Constant'

    rtb_rxpi_h = rt_remd_snf(rtb_rxpi_h, 6.2831853071795862);

    // Sum: '<S361>/Sum'
    rtb_rxpi_h -= rtb_Product_g;

    // Switch: '<S354>/Switch' incorporates:
    //   Gain: '<S357>/Gain'
    //   Sum: '<S352>/Sum3'
    //   Sum: '<S352>/Sum4'
    //   Trigonometry: '<S352>/Trigonometric Function'

    if (rtb_Product1_lq > 0.5) {
      rtb_Product2_b = 57.295779513082323 * rt_atan2d_snf
        (AHV_Model_B.Integrator1_n[1], AHV_Model_B.Integrator1_n[0]);
    }

    // End of Switch: '<S354>/Switch'

    // RateLimiter: '<S354>/Rate Limiter'
    if (AHV_Model_DW.LastMajorTime_c == (rtInf)) {
      AHV_Model_B.RateLimiter_m = rtb_Product2_b;
    } else {
      rtb_Switch3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_c;
      rtb_Product1_lf = rtb_Switch3 * 10.0;
      rtb_Switch2 = rtb_Product2_b - AHV_Model_DW.PrevY_o;
      if (rtb_Switch2 > rtb_Product1_lf) {
        AHV_Model_B.RateLimiter_m = AHV_Model_DW.PrevY_o + rtb_Product1_lf;
      } else {
        rtb_Switch3 *= -10.0;
        if (rtb_Switch2 < rtb_Switch3) {
          AHV_Model_B.RateLimiter_m = AHV_Model_DW.PrevY_o + rtb_Switch3;
        } else {
          AHV_Model_B.RateLimiter_m = rtb_Product2_b;
        }
      }
    }

    // End of RateLimiter: '<S354>/Rate Limiter'

    // Gain: '<S364>/Gain1'
    rtb_Product_g = 0.017453292519943295 * AHV_Model_B.RateLimiter_m;

    // Saturate: '<S365>/x_Saturation'
    if (rtb_Product_g > 1.0E+10) {
      rtb_Product_g = 1.0E+10;
    } else {
      if (rtb_Product_g < -1.0E+10) {
        rtb_Product_g = -1.0E+10;
      }
    }

    // End of Saturate: '<S365>/x_Saturation'

    // Signum: '<S365>/x_Sign'
    if (rtb_Product_g < 0.0) {
      rtb_Product1_lq = -1.0;
    } else if (rtb_Product_g > 0.0) {
      rtb_Product1_lq = 1.0;
    } else if (rtb_Product_g == 0.0) {
      rtb_Product1_lq = 0.0;
    } else {
      rtb_Product1_lq = (rtNaN);
    }

    // End of Signum: '<S365>/x_Sign'

    // Gain: '<S365>/pi'
    rtb_Product1_lq *= 3.1415926535897931;

    // Sum: '<S365>/Sum1'
    rtb_Product_g += rtb_Product1_lq;

    // Math: '<S365>/Math Function' incorporates:
    //   Constant: '<S365>/Constant'

    rtb_Product_g = rt_remd_snf(rtb_Product_g, 6.2831853071795862);
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S354>/Chart' incorporates:
      //   Sum: '<S365>/Sum'

      AHV_Model_Chart(AHV_Model_B.Merge_ik, AHV_Model_B.Merge2_j,
                      AHV_Model_B.Merge5_o, AHV_Model_B.Integrator1_n[0],
                      AHV_Model_B.Integrator1_n[1], rtb_Product_g -
                      rtb_Product1_lq, AHV_Model_B.Integrator1_n[5],
                      &AHV_Model_B.x_ref_rel, &AHV_Model_B.y_ref_rel,
                      &AHV_Model_B.yaw_ref_rel, &AHV_Model_DW.sf_Chart_c);
    }

    // Sum: '<S355>/Sum2' incorporates:
    //   Integrator: '<S353>/Integrator3'

    AHV_Model_B.regulation_error[0] = AHV_Model_B.x_ref_rel -
      AHV_Model_X.Integrator3_CSTATE[0];
    AHV_Model_B.regulation_error[1] = AHV_Model_B.y_ref_rel -
      AHV_Model_X.Integrator3_CSTATE[1];
    AHV_Model_B.regulation_error[2] = AHV_Model_B.yaw_ref_rel - rtb_rxpi_h;

    // Signum: '<S370>/x_Sign' incorporates:
    //   Fcn: '<S355>/yaw_angle'

    if (rtb_rxpi_h < 0.0) {
      rtb_Product_g = -1.0;
    } else if (rtb_rxpi_h > 0.0) {
      rtb_Product_g = 1.0;
    } else if (rtb_rxpi_h == 0.0) {
      rtb_Product_g = 0.0;
    } else {
      rtb_Product_g = (rtNaN);
    }

    // End of Signum: '<S370>/x_Sign'

    // Gain: '<S370>/pi'
    rtb_Product1_lq = 3.1415926535897931 * rtb_Product_g;

    // Sum: '<S370>/Sum' incorporates:
    //   Constant: '<S370>/Constant'
    //   Fcn: '<S355>/yaw_angle'
    //   Math: '<S370>/Math Function'
    //   Sum: '<S370>/Sum1'

    rtb_Product_g = rt_remd_snf(rtb_rxpi_h + rtb_Product1_lq, 6.2831853071795862)
      - rtb_Product1_lq;

    // Saturate: '<S369>/x_Saturation'
    if (AHV_Model_B.regulation_error[2] > 1.0E+10) {
      rtb_Product1_lq = 1.0E+10;
    } else if (AHV_Model_B.regulation_error[2] < -1.0E+10) {
      rtb_Product1_lq = -1.0E+10;
    } else {
      rtb_Product1_lq = AHV_Model_B.regulation_error[2];
    }

    // End of Saturate: '<S369>/x_Saturation'

    // Signum: '<S369>/x_Sign'
    if (rtb_Product1_lq < 0.0) {
      rtb_Product2_b = -1.0;
    } else if (rtb_Product1_lq > 0.0) {
      rtb_Product2_b = 1.0;
    } else if (rtb_Product1_lq == 0.0) {
      rtb_Product2_b = 0.0;
    } else {
      rtb_Product2_b = (rtNaN);
    }

    // End of Signum: '<S369>/x_Sign'

    // Gain: '<S369>/pi'
    rtb_Product2_b *= 3.1415926535897931;

    // Sum: '<S369>/Sum1'
    rtb_Product1_lq += rtb_Product2_b;

    // Math: '<S369>/Math Function' incorporates:
    //   Constant: '<S369>/Constant'

    rtb_Product1_lq = rt_remd_snf(rtb_Product1_lq, 6.2831853071795862);

    // Sum: '<S369>/Sum'
    AHV_Model_B.rxpi_j5 = rtb_Product1_lq - rtb_Product2_b;

    // Fcn: '<S368>/Row3'
    AHV_Model_B.Row3 = AHV_Model_B.rxpi_j5;

    // Fcn: '<S368>/Row1' incorporates:
    //   Fcn: '<S368>/Row2'

    rtb_Product1_lq = std::sin(rtb_Product_g);
    rtb_Product2_b = std::cos(rtb_Product_g);
    rtb_Product_g = rtb_Product2_b * AHV_Model_B.regulation_error[0] +
      rtb_Product1_lq * AHV_Model_B.regulation_error[1];

    // If: '<S373>/If' incorporates:
    //   Abs: '<S373>/Abs'
    //   Inport: '<S376>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_i = static_cast<int8_T>(!(std::abs
        (rtb_Product_g) <= rtb_Merge3));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_i) {
     case 0:
      // Outputs for IfAction SubSystem: '<S373>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S376>/Action Port'

      AHV_Model_B.Merge = rtb_Product_g;

      // End of Outputs for SubSystem: '<S373>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S373>/If Action Subsystem' incorporates:
      //   ActionPort: '<S375>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge3, rtb_Product_g,
        &AHV_Model_B.Merge);

      // End of Outputs for SubSystem: '<S373>/If Action Subsystem'
      break;
    }

    // End of If: '<S373>/If'

    // Fcn: '<S368>/Row2'
    rtb_Product1_lq = -rtb_Product1_lq * AHV_Model_B.regulation_error[0] +
      rtb_Product2_b * AHV_Model_B.regulation_error[1];

    // If: '<S374>/If' incorporates:
    //   Abs: '<S374>/Abs'
    //   Inport: '<S378>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_j = static_cast<int8_T>(!(std::abs
        (rtb_Product1_lq) <= rtb_Merge4));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_j) {
     case 0:
      // Outputs for IfAction SubSystem: '<S374>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S378>/Action Port'

      AHV_Model_B.Merge_e = rtb_Product1_lq;

      // End of Outputs for SubSystem: '<S374>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S374>/If Action Subsystem' incorporates:
      //   ActionPort: '<S377>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge4, rtb_Product1_lq,
        &AHV_Model_B.Merge_e);

      // End of Outputs for SubSystem: '<S374>/If Action Subsystem'
      break;
    }

    // End of If: '<S374>/If'

    // If: '<S366>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_d = static_cast<int8_T>
        (!AHV_Model_B.Merge_ik);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_d) {
     case 0:
      // Outputs for IfAction SubSystem: '<S366>/If Action Subsystem' incorporates:
      //   ActionPort: '<S371>/Action Port'

      AHV_Model_IfActionSubsystem(AHV_Model_B.Merge, AHV_Model_B.Merge_e,
        AHV_Model_B.Row3, rtb_nu, AHV_Model_ConstP.pooled2);

      // End of Outputs for SubSystem: '<S366>/If Action Subsystem'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S366>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S372>/Action Port'

      AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge, AHV_Model_B.Merge_e,
        AHV_Model_B.Row3, rtb_nu, AHV_Model_ConstP.pooled1);

      // End of Outputs for SubSystem: '<S366>/If Action Subsystem1'
      break;
    }

    // End of If: '<S366>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S355>/Ki' incorporates:
      //   DiscreteIntegrator: '<S355>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki[i] = 0.0;
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[0];
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[1];
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[2];
      }

      // End of Gain: '<S355>/Ki'
    }

    // Sum: '<S355>/Sum3'
    rtb_sintheta[0] = rtb_nu[0] + AHV_Model_B.Ki[0];

    // Integrator: '<S353>/Integrator4'
    rtb_nu[0] = AHV_Model_X.Integrator4_CSTATE[0];

    // Sum: '<S355>/Sum3'
    rtb_sintheta[1] = rtb_nu[1] + AHV_Model_B.Ki[1];

    // Integrator: '<S353>/Integrator4'
    rtb_nu[1] = AHV_Model_X.Integrator4_CSTATE[1];

    // Sum: '<S355>/Sum3'
    rtb_sintheta[2] = rtb_nu[2] + AHV_Model_B.Ki[2];

    // Integrator: '<S353>/Integrator4'
    rtb_nu[2] = AHV_Model_X.Integrator4_CSTATE[2];

    // Sum: '<S355>/Sum1' incorporates:
    //   Gain: '<S355>/Kd'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = rtb_sintheta[i] - ((AHV_Model_ConstP.pooled70[i + 3] *
        rtb_nu[1] + AHV_Model_ConstP.pooled70[i] * rtb_nu[0]) +
        AHV_Model_ConstP.pooled70[i + 6] * rtb_nu[2]);
    }

    // End of Sum: '<S355>/Sum1'

    // Saturate: '<S356>/Surge Force Saturation'
    if (rtb_costheta[0] > 2.0E+6) {
      rtb_costheta_1 = 2.0E+6;
    } else if (rtb_costheta[0] < -2.0E+6) {
      rtb_costheta_1 = -2.0E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[0];
    }

    // End of Saturate: '<S356>/Surge Force Saturation'

    // Product: '<S356>/Product'
    rtb_Product_g = AHV_Model_B.Merge7 * rtb_costheta_1;

    // Saturate: '<S356>/Sway Force Saturation'
    if (rtb_costheta[1] > 1.5E+6) {
      rtb_costheta_1 = 1.5E+6;
    } else if (rtb_costheta[1] < -1.5E+6) {
      rtb_costheta_1 = -1.5E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[1];
    }

    // End of Saturate: '<S356>/Sway Force Saturation'

    // Product: '<S356>/Product1'
    rtb_Product1_lq = AHV_Model_B.Merge8 * rtb_costheta_1;

    // Saturate: '<S356>/Yaw Moment Saturation'
    if (rtb_costheta[2] > 2.0E+7) {
      rtb_costheta_1 = 2.0E+7;
    } else if (rtb_costheta[2] < -2.0E+7) {
      rtb_costheta_1 = -2.0E+7;
    } else {
      rtb_costheta_1 = rtb_costheta[2];
    }

    // End of Saturate: '<S356>/Yaw Moment Saturation'

    // Product: '<S356>/Product2'
    rtb_Product2_b = AHV_Model_B.Merge9 * rtb_costheta_1;

    // Product: '<S303>/G*eta' incorporates:
    //   Constant: '<S303>/Spring stiffness'

    for (i = 0; i < 6; i++) {
      rtb_Integrator_m[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Integrator_m[i] += AHV_Model_ConstP.pooled55[6 * i_0 + i] *
          AHV_Model_B.Integrator1_n[i_0];
      }
    }

    // End of Product: '<S303>/G*eta'

    // Sum: '<S303>/Sum1' incorporates:
    //   Constant: '<S300>/Constant'

    tmp_8 = rtb_Integrator_m[0];
    tmp = rtb_Integrator_m[1];
    tmp_0 = rtb_Integrator_m[2];
    tmp_1 = rtb_Integrator_m[3];
    tmp_2 = rtb_Integrator_m[4];
    tmp_3 = rtb_Integrator_m[5];
    rtb_Integrator_m[0] = rtb_Product_g - tmp_8;
    rtb_Integrator_m[1] = rtb_Product1_lq - tmp;
    rtb_Integrator_m[2] = 0.0 - tmp_0;
    rtb_Integrator_m[3] = 0.0 - tmp_1;
    rtb_Integrator_m[4] = 0.0 - tmp_2;
    rtb_Integrator_m[5] = rtb_Product2_b - tmp_3;
    for (i = 0; i < 3; i++) {
      // Product: '<S299>/Product1' incorporates:
      //   Inport: '<Root>/tau_cable4'

      rtb_costheta_1 = rtb_TmpSignalConversionAtPro_ni[i + 6] * tau_cable4[2] +
        (rtb_TmpSignalConversionAtPro_ni[i + 3] * tau_cable4[1] +
         rtb_TmpSignalConversionAtPro_ni[i] * tau_cable4[0]);

      // Gain: '<S298>/tau_cableGain1'
      rtb_costheta_0[i] = rtb_costheta_1;

      // Product: '<S299>/Product1'
      rtb_costheta[i] = rtb_costheta_1;
    }

    // Gain: '<S298>/tau_cableGain1' incorporates:
    //   Inport: '<Root>/ahv_fairlead4'
    //   Product: '<S350>/ae'
    //   Product: '<S350>/af'
    //   Product: '<S350>/bd'
    //   Product: '<S350>/bf'
    //   Product: '<S350>/cd'
    //   Product: '<S350>/ce'
    //   Sum: '<S350>/Sum'
    //   Sum: '<S350>/Sum1'
    //   Sum: '<S350>/Sum2'

    rtb_costheta_0[3] = ahv_fairlead4[1] * rtb_costheta[2] - ahv_fairlead4[2] *
      rtb_costheta[1];
    rtb_costheta_0[4] = ahv_fairlead4[2] * rtb_costheta[0] - ahv_fairlead4[0] *
      rtb_costheta[2];
    rtb_costheta_0[5] = ahv_fairlead4[0] * rtb_costheta[1] - ahv_fairlead4[1] *
      rtb_costheta[0];

    // Sum: '<S303>/Sum7' incorporates:
    //   Gain: '<S298>/tau_cableGain1'
    //   Sum: '<S298>/Sum'
    //   Sum: '<S303>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Integrator_m[i] = ((rtb_Integrator_f[i] + rtb_Sum[i]) +
        rtb_Integrator_m[i]) + AHV_Model_ConstP.pooled77[i] * rtb_costheta_0[i];
    }

    // End of Sum: '<S303>/Sum7'

    // Sum: '<S303>/Sum5' incorporates:
    //   Abs: '<S309>/Abs1'
    //   Gain: '<S309>/Gain3'
    //   Product: '<S309>/Product7'
    //   Product: '<S309>/Product8'

    tmp_8 = rtb_Integrator_m[0];
    tmp = rtb_Integrator_m[1];
    tmp_0 = rtb_Integrator_m[2];
    tmp_1 = rtb_Integrator_m[3];
    tmp_2 = rtb_Integrator_m[4];
    tmp_3 = rtb_Integrator_m[5];
    rtb_Integrator_m[0] = -(std::abs(AHV_Model_B.nu_r[0]) * AHV_Model_B.nu_r[0] *
      AHV_Model_ConstB.Product5_hj) + tmp_8;
    rtb_Integrator_m[1] = tmp + AHV_Model_B.dx1;
    rtb_Integrator_m[2] = tmp_0;
    rtb_Integrator_m[3] = tmp_1;
    rtb_Integrator_m[4] = tmp_2;
    rtb_Integrator_m[5] = tmp_3 + AHV_Model_B.dx2;

    // Sum: '<S303>/Sum4' incorporates:
    //   Constant: '<S303>/damping'
    //   Product: '<S303>/D*eta '

    for (i = 0; i < 6; i++) {
      tmp_8 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        tmp_8 += AHV_Model_ConstP.pooled58[6 * i_0 + i] * AHV_Model_B.nu_r[i_0];
      }

      rtb_Integrator_m[i] -= tmp_8;
    }

    // End of Sum: '<S303>/Sum4'

    // Product: '<S322>/Product2' incorporates:
    //   Integrator: '<S322>/Integrator'

    rtb_Switch3 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Switch3 += AHV_Model_ConstB.MathFunction_e[i] *
        AHV_Model_X.Integrator_CSTATE_d[i];
    }

    // Sum: '<S303>/Sum' incorporates:
    //   Constant: '<S322>/B44_inf'
    //   Constant: '<S322>/D44'
    //   Product: '<S322>/Product2'
    //   Product: '<S322>/Product3'
    //   Product: '<S322>/Product4'
    //   StateSpace: '<S310>/Dp(1,1)'
    //   StateSpace: '<S310>/Dp(2,2)'
    //   StateSpace: '<S310>/Dp(2,4)'
    //   StateSpace: '<S310>/Dp(2,6)'
    //   StateSpace: '<S310>/Dp(3,3)'
    //   StateSpace: '<S310>/Dp(3,5)'
    //   StateSpace: '<S310>/Dp(4,2)'
    //   StateSpace: '<S310>/Dp(4,6)'
    //   StateSpace: '<S310>/Dp(5,3)'
    //   StateSpace: '<S310>/Dp(5,5)'
    //   StateSpace: '<S310>/Dp(6,2)'
    //   StateSpace: '<S310>/Dp(6,4)'
    //   StateSpace: '<S310>/Dp(6,6)'
    //   Sum: '<S310>/Sum'
    //   Sum: '<S310>/Sum1'
    //   Sum: '<S310>/Sum3'
    //   Sum: '<S310>/Sum4'
    //   Sum: '<S310>/Sum5'
    //   Sum: '<S322>/Sum1'
    //   Sum: '<S322>/Sum2'

    tmp_8 = rtb_Integrator_m[0];
    tmp = rtb_Integrator_m[1];
    tmp_0 = rtb_Integrator_m[2];
    tmp_1 = rtb_Integrator_m[3];
    tmp_2 = rtb_Integrator_m[4];
    tmp_3 = rtb_Integrator_m[5];
    rtb_Integrator_m[0] = tmp_8 - (((((-192192.528024618 *
      AHV_Model_X.Dp11_CSTATE[0] + -4350.6445868200408 *
      AHV_Model_X.Dp11_CSTATE[1]) + -26196.671588656758 *
      AHV_Model_X.Dp11_CSTATE[2]) + -69550.629589629229 *
      AHV_Model_X.Dp11_CSTATE[3]) + -6667.9011063139569 *
      AHV_Model_X.Dp11_CSTATE[4]) + 28360.700724836512 * AHV_Model_B.nu_r[0]);
    rtb_Integrator_m[1] = tmp - (((((((-1.9219252802461751E+6 *
      AHV_Model_X.Dp22_CSTATE[0] + -43506.445868266834 *
      AHV_Model_X.Dp22_CSTATE[1]) + -261966.71588657616 *
      AHV_Model_X.Dp22_CSTATE[2]) + -695506.29589629406 *
      AHV_Model_X.Dp22_CSTATE[3]) + -66679.0110631567 * AHV_Model_X.Dp22_CSTATE
      [4]) + 283607.00724836387 * AHV_Model_B.nu_r[1]) + (((((163377.73376110662
      * AHV_Model_X.Dp24_CSTATE[0] + -2.24222697524951E+6 *
      AHV_Model_X.Dp24_CSTATE[1]) + 782719.9224317756 * AHV_Model_X.Dp24_CSTATE
      [2]) + 906439.02406890341 * AHV_Model_X.Dp24_CSTATE[3]) +
      48702.05056677026 * AHV_Model_X.Dp24_CSTATE[4]) + 342290.49224698928 *
      AHV_Model_B.nu_r[3])) + (((((-377217.71116477274 *
      AHV_Model_X.Dp26_CSTATE[0] + 9.7759770791847613E+6 *
      AHV_Model_X.Dp26_CSTATE[1]) + 1.2719330237464791E+6 *
      AHV_Model_X.Dp26_CSTATE[2]) + 3.67278239005485E+6 *
      AHV_Model_X.Dp26_CSTATE[3]) + -523575.51382008579 *
      AHV_Model_X.Dp26_CSTATE[4]) + 1.4324880870688609E+6 * AHV_Model_B.nu_r[5]));
    rtb_Integrator_m[2] = tmp_0 - ((((((-2.8230572248621574E+6 *
      AHV_Model_X.Dp33_CSTATE[0] + 2971.0437052436932 * AHV_Model_X.Dp33_CSTATE
      [1]) + -112975.07052200216 * AHV_Model_X.Dp33_CSTATE[2]) +
      -459867.38057091518 * AHV_Model_X.Dp33_CSTATE[3]) + 35755.997807249405 *
      AHV_Model_X.Dp33_CSTATE[4]) + 425986.43085235963 * AHV_Model_B.nu_r[2]) +
      (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE[0] +
           -189169.67892356269 * AHV_Model_X.Dp35_CSTATE[1]) +
          1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE[2]) +
         1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE[3]) +
        2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE[4]) +
       6.2116453321941122E+6 * AHV_Model_B.nu_r[4]));
    rtb_Integrator_m[3] = tmp_1 - (((((((163377.73376110662 *
      AHV_Model_X.Dp42_CSTATE[0] + -2.24222697524951E+6 *
      AHV_Model_X.Dp42_CSTATE[1]) + 782719.9224317756 * AHV_Model_X.Dp42_CSTATE
      [2]) + 906439.02406890341 * AHV_Model_X.Dp42_CSTATE[3]) +
      48702.05056677026 * AHV_Model_X.Dp42_CSTATE[4]) + 342290.49224698928 *
      AHV_Model_B.nu_r[1]) + ((2.3074903324953854E+7 * AHV_Model_B.nu_r[3] +
      rtb_Switch3) + 4.97124674352808E+7 * AHV_Model_B.nu_r[3])) +
      (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE[0] +
           2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE[1]) +
          -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE[2]) +
         1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE[3]) +
        -70143.531785937739 * AHV_Model_X.Dp46_CSTATE[4]) +
       -8.0522169282609783E+6 * AHV_Model_B.nu_r[5]));
    rtb_Integrator_m[4] = tmp_2 - ((((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp53_CSTATE[0] + -189169.67892356269 *
      AHV_Model_X.Dp53_CSTATE[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp53_CSTATE[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp53_CSTATE[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp53_CSTATE[4]) + 6.2116453321941122E+6 * AHV_Model_B.nu_r[2])
      + (((((-1.3510265271416564E+7 * AHV_Model_X.Dp55_CSTATE[0] +
             2.1371969925573766E+9 * AHV_Model_X.Dp55_CSTATE[1]) +
            9.8710379288756877E+7 * AHV_Model_X.Dp55_CSTATE[2]) +
           6.1487736935078049E+8 * AHV_Model_X.Dp55_CSTATE[3]) +
          5.7220438213245332E+8 * AHV_Model_X.Dp55_CSTATE[4]) +
         3.0813057877131975E+8 * AHV_Model_B.nu_r[4]));
    rtb_Integrator_m[5] = tmp_3 - (((((((-377217.71116477274 *
      AHV_Model_X.Dp62_CSTATE[0] + 9.7759770791847613E+6 *
      AHV_Model_X.Dp62_CSTATE[1]) + 1.2719330237464791E+6 *
      AHV_Model_X.Dp62_CSTATE[2]) + 3.67278239005485E+6 *
      AHV_Model_X.Dp62_CSTATE[3]) + -523575.51382008579 *
      AHV_Model_X.Dp62_CSTATE[4]) + 1.4324880870688609E+6 * AHV_Model_B.nu_r[1])
      + (((((-1.5746048421550707E+6 * AHV_Model_X.Dp64_CSTATE[0] +
             2.9884305406230576E+7 * AHV_Model_X.Dp64_CSTATE[1]) +
            -4.3953040127794016E+6 * AHV_Model_X.Dp64_CSTATE[2]) +
           1.205872403100216E+7 * AHV_Model_X.Dp64_CSTATE[3]) +
          -70143.531785937739 * AHV_Model_X.Dp64_CSTATE[4]) +
         -8.0522169282609783E+6 * AHV_Model_B.nu_r[3])) +
      (((((-8.07757267677547E+8 * AHV_Model_X.Dp66_CSTATE[0] +
           -1.0798841515040772E+7 * AHV_Model_X.Dp66_CSTATE[1]) +
          -1.27285250509213E+8 * AHV_Model_X.Dp66_CSTATE[2]) +
         -3.0027175358715254E+8 * AHV_Model_X.Dp66_CSTATE[3]) +
        1.7969070654437996E+7 * AHV_Model_X.Dp66_CSTATE[4]) +
       1.1943182502908507E+8 * AHV_Model_B.nu_r[5]));

    // Product: '<S303>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_a, rtb_Integrator_m,
      AHV_Model_B.Minvtau);

    // Outport: '<Root>/nu4' incorporates:
    //   Integrator: '<S303>/Integrator'

    for (i = 0; i < 6; i++) {
      nu4[i] = AHV_Model_X.Integrator_CSTATE[i];
    }

    // End of Outport: '<Root>/nu4'

    // Gain: '<S206>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_Product2_n = 0.017453292519943295 * Current_direction;

    // Product: '<S206>/Product' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S206>/cos'

    rtb_Switch2 = std::cos(rtb_Product2_n) * Current_speed;

    // Product: '<S206>/Product1' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S206>/sin'

    rtb_Product2_n = std::sin(rtb_Product2_n) * Current_speed;

    // Sum: '<S207>/Sum6' incorporates:
    //   Gain: '<S202>/Current on//off '
    //   Integrator: '<S207>/Integrator'

    AHV_Model_B.nu_r_m[0] = AHV_Model_X.Integrator_CSTATE_h[0] - rtb_Switch2;
    AHV_Model_B.nu_r_m[1] = AHV_Model_X.Integrator_CSTATE_h[1] - rtb_Product2_n;
    AHV_Model_B.nu_r_m[2] = AHV_Model_X.Integrator_CSTATE_h[2];
    AHV_Model_B.nu_r_m[3] = AHV_Model_X.Integrator_CSTATE_h[3];
    AHV_Model_B.nu_r_m[4] = AHV_Model_X.Integrator_CSTATE_h[4];
    AHV_Model_B.nu_r_m[5] = AHV_Model_X.Integrator_CSTATE_h[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S213>/Cross-flow drag trapezoidal integration' 
      // Constant: '<S213>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_g,
        AHV_Model_B.nu_r_m[1], AHV_Model_B.nu_r_m[5], &AHV_Model_B.Sum_ps,
        &AHV_Model_B.Sum2_p, 82.800003);

      // End of Outputs for SubSystem: '<S213>/Cross-flow drag trapezoidal integration' 

      // Product: '<S213>/dx1' incorporates:
      //   Constant: '<S213>/2D drag coefficient '
      //   Constant: '<S213>/Transversal area//Lpp'
      //   Constant: '<S213>/rho'
      //   Gain: '<S213>/Gain1'

      AHV_Model_B.dx1_a = -0.5 * AHV_Model_B.Sum_ps * 1025.0 * 5.4 *
        0.69954741847648816;

      // Product: '<S213>/dx2' incorporates:
      //   Constant: '<S213>/2D drag coefficient '
      //   Constant: '<S213>/Transversal area//Lpp'
      //   Constant: '<S213>/rho'
      //   Gain: '<S213>/Gain'

      AHV_Model_B.dx2_c = -0.5 * AHV_Model_B.Sum2_p * 1025.0 * 5.4 *
        0.69954741847648816;
    }

    // Integrator: '<S257>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init3'
    //   SignalConversion generated from: '<S257>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_i != 0) {
      AHV_Model_X.Integrator3_CSTATE_c[0] = Vessel_init3[0];
      AHV_Model_X.Integrator3_CSTATE_c[1] = Vessel_init3[1];
      AHV_Model_X.Integrator3_CSTATE_c[2] = Vessel_init3[5];
    }

    // Saturate: '<S265>/x_Saturation' incorporates:
    //   Integrator: '<S257>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_c[2] > 1.0E+10) {
      rtb_Switch2 = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_c[2] < -1.0E+10) {
      rtb_Switch2 = -1.0E+10;
    } else {
      rtb_Switch2 = AHV_Model_X.Integrator3_CSTATE_c[2];
    }

    // End of Saturate: '<S265>/x_Saturation'

    // Signum: '<S265>/x_Sign'
    if (rtb_Switch2 < 0.0) {
      rtb_Switch3 = -1.0;
    } else if (rtb_Switch2 > 0.0) {
      rtb_Switch3 = 1.0;
    } else if (rtb_Switch2 == 0.0) {
      rtb_Switch3 = 0.0;
    } else {
      rtb_Switch3 = (rtNaN);
    }

    // End of Signum: '<S265>/x_Sign'

    // Gain: '<S265>/pi'
    rtb_Switch3 *= 3.1415926535897931;

    // Sum: '<S265>/Sum1'
    rtb_Switch2 += rtb_Switch3;

    // Math: '<S265>/Math Function' incorporates:
    //   Constant: '<S265>/Constant'

    rtb_Switch2 = rt_remd_snf(rtb_Switch2, 6.2831853071795862);

    // Sum: '<S265>/Sum'
    rtb_Switch2 -= rtb_Switch3;

    // Sum: '<S259>/Sum2' incorporates:
    //   Integrator: '<S257>/Integrator3'

    AHV_Model_B.regulation_error_d[0] = AHV_Model_B.x_ref_rel_p -
      AHV_Model_X.Integrator3_CSTATE_c[0];
    AHV_Model_B.regulation_error_d[1] = AHV_Model_B.y_ref_rel_c -
      AHV_Model_X.Integrator3_CSTATE_c[1];
    AHV_Model_B.regulation_error_d[2] = AHV_Model_B.yaw_ref_rel_m - rtb_Switch2;

    // Signum: '<S274>/x_Sign' incorporates:
    //   Fcn: '<S259>/yaw_angle'

    if (rtb_Switch2 < 0.0) {
      rtb_Switch3 = -1.0;
    } else if (rtb_Switch2 > 0.0) {
      rtb_Switch3 = 1.0;
    } else if (rtb_Switch2 == 0.0) {
      rtb_Switch3 = 0.0;
    } else {
      rtb_Switch3 = (rtNaN);
    }

    // End of Signum: '<S274>/x_Sign'

    // Gain: '<S274>/pi'
    rtb_Product1_lf = 3.1415926535897931 * rtb_Switch3;

    // Sum: '<S274>/Sum' incorporates:
    //   Constant: '<S274>/Constant'
    //   Fcn: '<S259>/yaw_angle'
    //   Math: '<S274>/Math Function'
    //   Sum: '<S274>/Sum1'

    rtb_Switch3 = rt_remd_snf(rtb_Switch2 + rtb_Product1_lf, 6.2831853071795862)
      - rtb_Product1_lf;

    // Saturate: '<S273>/x_Saturation'
    if (AHV_Model_B.regulation_error_d[2] > 1.0E+10) {
      rtb_Product1_lf = 1.0E+10;
    } else if (AHV_Model_B.regulation_error_d[2] < -1.0E+10) {
      rtb_Product1_lf = -1.0E+10;
    } else {
      rtb_Product1_lf = AHV_Model_B.regulation_error_d[2];
    }

    // End of Saturate: '<S273>/x_Saturation'

    // Signum: '<S273>/x_Sign'
    if (rtb_Product1_lf < 0.0) {
      SwayFailure = -1.0;
    } else if (rtb_Product1_lf > 0.0) {
      SwayFailure = 1.0;
    } else if (rtb_Product1_lf == 0.0) {
      SwayFailure = 0.0;
    } else {
      SwayFailure = (rtNaN);
    }

    // End of Signum: '<S273>/x_Sign'

    // Gain: '<S273>/pi'
    SwayFailure *= 3.1415926535897931;

    // Sum: '<S273>/Sum1'
    rtb_Product1_lf += SwayFailure;

    // Math: '<S273>/Math Function' incorporates:
    //   Constant: '<S273>/Constant'

    rtb_Product1_lf = rt_remd_snf(rtb_Product1_lf, 6.2831853071795862);

    // Sum: '<S273>/Sum'
    AHV_Model_B.rxpi_k = rtb_Product1_lf - SwayFailure;

    // Fcn: '<S272>/Row3'
    AHV_Model_B.Row3_j = AHV_Model_B.rxpi_k;

    // Fcn: '<S272>/Row1' incorporates:
    //   Fcn: '<S272>/Row2'

    rtb_Product1_lf = std::sin(rtb_Switch3);
    rtb_Switch3 = std::cos(rtb_Switch3);
    SwayFailure = rtb_Switch3 * AHV_Model_B.regulation_error_d[0] +
      rtb_Product1_lf * AHV_Model_B.regulation_error_d[1];

    // If: '<S277>/If' incorporates:
    //   Abs: '<S277>/Abs'
    //   Inport: '<S280>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_h = static_cast<int8_T>(!(std::abs
        (SwayFailure) <= rtb_Merge3_k));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_h) {
     case 0:
      // Outputs for IfAction SubSystem: '<S277>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S280>/Action Port'

      AHV_Model_B.Merge_o = SwayFailure;

      // End of Outputs for SubSystem: '<S277>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S277>/If Action Subsystem' incorporates:
      //   ActionPort: '<S279>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge3_k, SwayFailure,
        &AHV_Model_B.Merge_o);

      // End of Outputs for SubSystem: '<S277>/If Action Subsystem'
      break;
    }

    // End of If: '<S277>/If'

    // Fcn: '<S272>/Row2'
    rtb_Switch3 = -rtb_Product1_lf * AHV_Model_B.regulation_error_d[0] +
      rtb_Switch3 * AHV_Model_B.regulation_error_d[1];

    // If: '<S278>/If' incorporates:
    //   Abs: '<S278>/Abs'
    //   Inport: '<S282>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_m = static_cast<int8_T>(!(std::abs
        (rtb_Switch3) <= rtb_Merge4_j));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_m) {
     case 0:
      // Outputs for IfAction SubSystem: '<S278>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S282>/Action Port'

      AHV_Model_B.Merge_p = rtb_Switch3;

      // End of Outputs for SubSystem: '<S278>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S278>/If Action Subsystem' incorporates:
      //   ActionPort: '<S281>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge4_j, rtb_Switch3,
        &AHV_Model_B.Merge_p);

      // End of Outputs for SubSystem: '<S278>/If Action Subsystem'
      break;
    }

    // End of If: '<S278>/If'

    // If: '<S270>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_c = static_cast<int8_T>
        (!AHV_Model_B.Merge_pf);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_c) {
     case 0:
      // Outputs for IfAction SubSystem: '<S270>/If Action Subsystem' incorporates:
      //   ActionPort: '<S275>/Action Port'

      AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_o, AHV_Model_B.Merge_p,
        AHV_Model_B.Row3_j, rtb_Gain2_g, AHV_Model_ConstP.pooled2);

      // End of Outputs for SubSystem: '<S270>/If Action Subsystem'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S270>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S276>/Action Port'

      AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_o, AHV_Model_B.Merge_p,
        AHV_Model_B.Row3_j, rtb_Gain2_g, AHV_Model_ConstP.pooled1);

      // End of Outputs for SubSystem: '<S270>/If Action Subsystem1'
      break;
    }

    // End of If: '<S270>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S259>/Ki' incorporates:
      //   DiscreteIntegrator: '<S259>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki_p[i] = 0.0;
        AHV_Model_B.Ki_p[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_d[0];
        AHV_Model_B.Ki_p[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_d[1];
        AHV_Model_B.Ki_p[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_d[2];
      }

      // End of Gain: '<S259>/Ki'
    }

    // Sum: '<S259>/Sum3'
    rtb_sintheta[0] = rtb_Gain2_g[0] + AHV_Model_B.Ki_p[0];

    // Integrator: '<S257>/Integrator4'
    rtb_Gain2_g[0] = AHV_Model_X.Integrator4_CSTATE_f[0];

    // Sum: '<S259>/Sum3'
    rtb_sintheta[1] = rtb_Gain2_g[1] + AHV_Model_B.Ki_p[1];

    // Integrator: '<S257>/Integrator4'
    rtb_Gain2_g[1] = AHV_Model_X.Integrator4_CSTATE_f[1];

    // Sum: '<S259>/Sum3'
    rtb_sintheta[2] = rtb_Gain2_g[2] + AHV_Model_B.Ki_p[2];

    // Integrator: '<S257>/Integrator4'
    rtb_Gain2_g[2] = AHV_Model_X.Integrator4_CSTATE_f[2];

    // Sum: '<S259>/Sum1' incorporates:
    //   Gain: '<S259>/Kd'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = rtb_sintheta[i] - ((AHV_Model_ConstP.pooled70[i + 3] *
        rtb_Gain2_g[1] + AHV_Model_ConstP.pooled70[i] * rtb_Gain2_g[0]) +
        AHV_Model_ConstP.pooled70[i + 6] * rtb_Gain2_g[2]);
    }

    // End of Sum: '<S259>/Sum1'

    // Saturate: '<S260>/Surge Force Saturation'
    if (rtb_costheta[0] > 2.0E+6) {
      rtb_costheta_1 = 2.0E+6;
    } else if (rtb_costheta[0] < -2.0E+6) {
      rtb_costheta_1 = -2.0E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[0];
    }

    // End of Saturate: '<S260>/Surge Force Saturation'

    // Product: '<S260>/Product'
    rtb_Switch3 = AHV_Model_B.Merge7_h * rtb_costheta_1;

    // Saturate: '<S260>/Sway Force Saturation'
    if (rtb_costheta[1] > 1.5E+6) {
      rtb_costheta_1 = 1.5E+6;
    } else if (rtb_costheta[1] < -1.5E+6) {
      rtb_costheta_1 = -1.5E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[1];
    }

    // End of Saturate: '<S260>/Sway Force Saturation'

    // Product: '<S260>/Product1'
    rtb_Product1_lf = AHV_Model_B.Merge8_j * rtb_costheta_1;

    // Saturate: '<S260>/Yaw Moment Saturation'
    if (rtb_costheta[2] > 2.0E+7) {
      rtb_costheta_1 = 2.0E+7;
    } else if (rtb_costheta[2] < -2.0E+7) {
      rtb_costheta_1 = -2.0E+7;
    } else {
      rtb_costheta_1 = rtb_costheta[2];
    }

    // End of Saturate: '<S260>/Yaw Moment Saturation'

    // Product: '<S260>/Product2'
    SwayFailure = AHV_Model_B.Merge9_e * rtb_costheta_1;

    // Product: '<S207>/G*eta' incorporates:
    //   Constant: '<S207>/Spring stiffness'

    for (i = 0; i < 6; i++) {
      rtb_Sum[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Sum[i] += AHV_Model_ConstP.pooled55[6 * i_0 + i] *
          AHV_Model_B.Integrator1_p[i_0];
      }
    }

    // End of Product: '<S207>/G*eta'

    // Sum: '<S207>/Sum1' incorporates:
    //   Constant: '<S204>/Constant'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = rtb_Switch3 - tmp_8;
    rtb_Sum[1] = rtb_Product1_lf - tmp;
    rtb_Sum[2] = 0.0 - tmp_0;
    rtb_Sum[3] = 0.0 - tmp_1;
    rtb_Sum[4] = 0.0 - tmp_2;
    rtb_Sum[5] = SwayFailure - tmp_3;
    for (i = 0; i < 3; i++) {
      // Product: '<S203>/Product1' incorporates:
      //   Inport: '<Root>/tau_cable3'

      rtb_costheta_1 = rtb_TmpSignalConversionAtPro_fx[i + 6] * tau_cable3[2] +
        (rtb_TmpSignalConversionAtPro_fx[i + 3] * tau_cable3[1] +
         rtb_TmpSignalConversionAtPro_fx[i] * tau_cable3[0]);

      // Gain: '<S202>/tau_cableGain1'
      rtb_costheta_0[i] = rtb_costheta_1;

      // Product: '<S203>/Product1'
      rtb_costheta[i] = rtb_costheta_1;
    }

    // Gain: '<S202>/tau_cableGain1' incorporates:
    //   Inport: '<Root>/ahv_fairlead3'
    //   Product: '<S254>/ae'
    //   Product: '<S254>/af'
    //   Product: '<S254>/bd'
    //   Product: '<S254>/bf'
    //   Product: '<S254>/cd'
    //   Product: '<S254>/ce'
    //   Sum: '<S254>/Sum'
    //   Sum: '<S254>/Sum1'
    //   Sum: '<S254>/Sum2'

    rtb_costheta_0[3] = ahv_fairlead3[1] * rtb_costheta[2] - ahv_fairlead3[2] *
      rtb_costheta[1];
    rtb_costheta_0[4] = ahv_fairlead3[2] * rtb_costheta[0] - ahv_fairlead3[0] *
      rtb_costheta[2];
    rtb_costheta_0[5] = ahv_fairlead3[0] * rtb_costheta[1] - ahv_fairlead3[1] *
      rtb_costheta[0];

    // Sum: '<S207>/Sum7' incorporates:
    //   Gain: '<S202>/tau_cableGain1'
    //   Sum: '<S202>/Sum'
    //   Sum: '<S207>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Sum[i] = ((rtb_tau_WF[i] + rtb_tau_WD[i]) + rtb_Sum[i]) +
        AHV_Model_ConstP.pooled77[i] * rtb_costheta_0[i];
    }

    // End of Sum: '<S207>/Sum7'

    // Sum: '<S207>/Sum5' incorporates:
    //   Abs: '<S213>/Abs1'
    //   Gain: '<S213>/Gain3'
    //   Product: '<S213>/Product7'
    //   Product: '<S213>/Product8'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = -(std::abs(AHV_Model_B.nu_r_m[0]) * AHV_Model_B.nu_r_m[0] *
                   AHV_Model_ConstB.Product5_f) + tmp_8;
    rtb_Sum[1] = tmp + AHV_Model_B.dx1_a;
    rtb_Sum[2] = tmp_0;
    rtb_Sum[3] = tmp_1;
    rtb_Sum[4] = tmp_2;
    rtb_Sum[5] = tmp_3 + AHV_Model_B.dx2_c;

    // Sum: '<S207>/Sum4' incorporates:
    //   Constant: '<S207>/damping'
    //   Product: '<S207>/D*eta '

    for (i = 0; i < 6; i++) {
      tmp_8 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        tmp_8 += AHV_Model_ConstP.pooled58[6 * i_0 + i] * AHV_Model_B.nu_r_m[i_0];
      }

      rtb_Sum[i] -= tmp_8;
    }

    // End of Sum: '<S207>/Sum4'

    // Product: '<S226>/Product2' incorporates:
    //   Integrator: '<S226>/Integrator'

    rtb_Row2_e = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Row2_e += AHV_Model_ConstB.MathFunction_n[i] *
        AHV_Model_X.Integrator_CSTATE_c[i];
    }

    // Sum: '<S207>/Sum' incorporates:
    //   Constant: '<S226>/B44_inf'
    //   Constant: '<S226>/D44'
    //   Product: '<S226>/Product2'
    //   Product: '<S226>/Product3'
    //   Product: '<S226>/Product4'
    //   StateSpace: '<S214>/Dp(1,1)'
    //   StateSpace: '<S214>/Dp(2,2)'
    //   StateSpace: '<S214>/Dp(2,4)'
    //   StateSpace: '<S214>/Dp(2,6)'
    //   StateSpace: '<S214>/Dp(3,3)'
    //   StateSpace: '<S214>/Dp(3,5)'
    //   StateSpace: '<S214>/Dp(4,2)'
    //   StateSpace: '<S214>/Dp(4,6)'
    //   StateSpace: '<S214>/Dp(5,3)'
    //   StateSpace: '<S214>/Dp(5,5)'
    //   StateSpace: '<S214>/Dp(6,2)'
    //   StateSpace: '<S214>/Dp(6,4)'
    //   StateSpace: '<S214>/Dp(6,6)'
    //   Sum: '<S214>/Sum'
    //   Sum: '<S214>/Sum1'
    //   Sum: '<S214>/Sum3'
    //   Sum: '<S214>/Sum4'
    //   Sum: '<S214>/Sum5'
    //   Sum: '<S226>/Sum1'
    //   Sum: '<S226>/Sum2'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = tmp_8 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE_l[0] +
      -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE_l[1]) + -26196.671588656758 *
      AHV_Model_X.Dp11_CSTATE_l[2]) + -69550.629589629229 *
      AHV_Model_X.Dp11_CSTATE_l[3]) + -6667.9011063139569 *
      AHV_Model_X.Dp11_CSTATE_l[4]) + 28360.700724836512 * AHV_Model_B.nu_r_m[0]);
    rtb_Sum[1] = tmp - (((((((-1.9219252802461751E+6 *
      AHV_Model_X.Dp22_CSTATE_j[0] + -43506.445868266834 *
      AHV_Model_X.Dp22_CSTATE_j[1]) + -261966.71588657616 *
      AHV_Model_X.Dp22_CSTATE_j[2]) + -695506.29589629406 *
      AHV_Model_X.Dp22_CSTATE_j[3]) + -66679.0110631567 *
      AHV_Model_X.Dp22_CSTATE_j[4]) + 283607.00724836387 * AHV_Model_B.nu_r_m[1])
                         + (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_b[0]
      + -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_b[1]) + 782719.9224317756
      * AHV_Model_X.Dp24_CSTATE_b[2]) + 906439.02406890341 *
      AHV_Model_X.Dp24_CSTATE_b[3]) + 48702.05056677026 *
      AHV_Model_X.Dp24_CSTATE_b[4]) + 342290.49224698928 * AHV_Model_B.nu_r_m[3]))
                        + (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE_n[0]
      + 9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE_n[1]) +
      1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE_n[2]) +
      3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE_n[3]) + -523575.51382008579 *
      AHV_Model_X.Dp26_CSTATE_n[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_m[5]));
    rtb_Sum[2] = tmp_0 - ((((((-2.8230572248621574E+6 *
      AHV_Model_X.Dp33_CSTATE_m[0] + 2971.0437052436932 *
      AHV_Model_X.Dp33_CSTATE_m[1]) + -112975.07052200216 *
      AHV_Model_X.Dp33_CSTATE_m[2]) + -459867.38057091518 *
      AHV_Model_X.Dp33_CSTATE_m[3]) + 35755.997807249405 *
      AHV_Model_X.Dp33_CSTATE_m[4]) + 425986.43085235963 * AHV_Model_B.nu_r_m[2])
                          + (((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp35_CSTATE_f[0] + -189169.67892356269 *
      AHV_Model_X.Dp35_CSTATE_f[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp35_CSTATE_f[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp35_CSTATE_f[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp35_CSTATE_f[4]) + 6.2116453321941122E+6 *
      AHV_Model_B.nu_r_m[4]));
    rtb_Sum[3] = tmp_1 - (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE_b[0]
      + -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE_b[1]) + 782719.9224317756
      * AHV_Model_X.Dp42_CSTATE_b[2]) + 906439.02406890341 *
      AHV_Model_X.Dp42_CSTATE_b[3]) + 48702.05056677026 *
      AHV_Model_X.Dp42_CSTATE_b[4]) + 342290.49224698928 * AHV_Model_B.nu_r_m[1])
      + ((2.3074903324953854E+7 * AHV_Model_B.nu_r_m[3] + rtb_Row2_e) +
         4.97124674352808E+7 * AHV_Model_B.nu_r_m[3])) +
                          (((((-1.5746048421550707E+6 *
      AHV_Model_X.Dp46_CSTATE_o[0] + 2.9884305406230576E+7 *
      AHV_Model_X.Dp46_CSTATE_o[1]) + -4.3953040127794016E+6 *
      AHV_Model_X.Dp46_CSTATE_o[2]) + 1.205872403100216E+7 *
      AHV_Model_X.Dp46_CSTATE_o[3]) + -70143.531785937739 *
      AHV_Model_X.Dp46_CSTATE_o[4]) + -8.0522169282609783E+6 *
      AHV_Model_B.nu_r_m[5]));
    rtb_Sum[4] = tmp_2 - ((((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp53_CSTATE_i[0] + -189169.67892356269 *
      AHV_Model_X.Dp53_CSTATE_i[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp53_CSTATE_i[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp53_CSTATE_i[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp53_CSTATE_i[4]) + 6.2116453321941122E+6 *
      AHV_Model_B.nu_r_m[2]) + (((((-1.3510265271416564E+7 *
      AHV_Model_X.Dp55_CSTATE_k[0] + 2.1371969925573766E+9 *
      AHV_Model_X.Dp55_CSTATE_k[1]) + 9.8710379288756877E+7 *
      AHV_Model_X.Dp55_CSTATE_k[2]) + 6.1487736935078049E+8 *
      AHV_Model_X.Dp55_CSTATE_k[3]) + 5.7220438213245332E+8 *
      AHV_Model_X.Dp55_CSTATE_k[4]) + 3.0813057877131975E+8 *
      AHV_Model_B.nu_r_m[4]));
    rtb_Sum[5] = tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE_g
      [0] + 9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE_g[1]) +
      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE_g[2]) +
      3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE_g[3]) + -523575.51382008579 *
      AHV_Model_X.Dp62_CSTATE_g[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_m[1]) + (((((-1.5746048421550707E+6 *
      AHV_Model_X.Dp64_CSTATE_p[0] + 2.9884305406230576E+7 *
      AHV_Model_X.Dp64_CSTATE_p[1]) + -4.3953040127794016E+6 *
      AHV_Model_X.Dp64_CSTATE_p[2]) + 1.205872403100216E+7 *
      AHV_Model_X.Dp64_CSTATE_p[3]) + -70143.531785937739 *
      AHV_Model_X.Dp64_CSTATE_p[4]) + -8.0522169282609783E+6 *
      AHV_Model_B.nu_r_m[3])) + (((((-8.07757267677547E+8 *
      AHV_Model_X.Dp66_CSTATE_f[0] + -1.0798841515040772E+7 *
      AHV_Model_X.Dp66_CSTATE_f[1]) + -1.27285250509213E+8 *
      AHV_Model_X.Dp66_CSTATE_f[2]) + -3.0027175358715254E+8 *
      AHV_Model_X.Dp66_CSTATE_f[3]) + 1.7969070654437996E+7 *
      AHV_Model_X.Dp66_CSTATE_f[4]) + 1.1943182502908507E+8 *
      AHV_Model_B.nu_r_m[5]));

    // Product: '<S207>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_g, rtb_Sum,
      AHV_Model_B.Minvtau_n);

    // Outport: '<Root>/nu3' incorporates:
    //   Integrator: '<S207>/Integrator'

    for (i = 0; i < 6; i++) {
      nu3[i] = AHV_Model_X.Integrator_CSTATE_h[i];
    }

    // End of Outport: '<Root>/nu3'

    // Gain: '<S110>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_Product2_a = 0.017453292519943295 * Current_direction;

    // Product: '<S110>/Product' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S110>/cos'

    rtb_rxpi_c = std::cos(rtb_Product2_a) * Current_speed;

    // Product: '<S110>/Product1' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S110>/sin'

    rtb_Product2_a = std::sin(rtb_Product2_a) * Current_speed;

    // Sum: '<S111>/Sum6' incorporates:
    //   Gain: '<S106>/Current on//off '
    //   Integrator: '<S111>/Integrator'

    AHV_Model_B.nu_r_j[0] = AHV_Model_X.Integrator_CSTATE_o[0] - rtb_rxpi_c;
    AHV_Model_B.nu_r_j[1] = AHV_Model_X.Integrator_CSTATE_o[1] - rtb_Product2_a;
    AHV_Model_B.nu_r_j[2] = AHV_Model_X.Integrator_CSTATE_o[2];
    AHV_Model_B.nu_r_j[3] = AHV_Model_X.Integrator_CSTATE_o[3];
    AHV_Model_B.nu_r_j[4] = AHV_Model_X.Integrator_CSTATE_o[4];
    AHV_Model_B.nu_r_j[5] = AHV_Model_X.Integrator_CSTATE_o[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S117>/Cross-flow drag trapezoidal integration' 
      // Constant: '<S117>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_o,
        AHV_Model_B.nu_r_j[1], AHV_Model_B.nu_r_j[5], &AHV_Model_B.Sum_nu,
        &AHV_Model_B.Sum2_h, 82.800003);

      // End of Outputs for SubSystem: '<S117>/Cross-flow drag trapezoidal integration' 

      // Product: '<S117>/dx1' incorporates:
      //   Constant: '<S117>/2D drag coefficient '
      //   Constant: '<S117>/Transversal area//Lpp'
      //   Constant: '<S117>/rho'
      //   Gain: '<S117>/Gain1'

      AHV_Model_B.dx1_a0 = -0.5 * AHV_Model_B.Sum_nu * 1025.0 * 5.4 *
        0.69954741847648816;

      // Product: '<S117>/dx2' incorporates:
      //   Constant: '<S117>/2D drag coefficient '
      //   Constant: '<S117>/Transversal area//Lpp'
      //   Constant: '<S117>/rho'
      //   Gain: '<S117>/Gain'

      AHV_Model_B.dx2_f = -0.5 * AHV_Model_B.Sum2_h * 1025.0 * 5.4 *
        0.69954741847648816;
    }

    // Integrator: '<S161>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init2'
    //   SignalConversion generated from: '<S161>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_c != 0) {
      AHV_Model_X.Integrator3_CSTATE_l[0] = Vessel_init2[0];
      AHV_Model_X.Integrator3_CSTATE_l[1] = Vessel_init2[1];
      AHV_Model_X.Integrator3_CSTATE_l[2] = Vessel_init2[5];
    }

    // Saturate: '<S169>/x_Saturation' incorporates:
    //   Integrator: '<S161>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_l[2] > 1.0E+10) {
      rtb_rxpi_c = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_l[2] < -1.0E+10) {
      rtb_rxpi_c = -1.0E+10;
    } else {
      rtb_rxpi_c = AHV_Model_X.Integrator3_CSTATE_l[2];
    }

    // End of Saturate: '<S169>/x_Saturation'

    // Signum: '<S169>/x_Sign'
    if (rtb_rxpi_c < 0.0) {
      rtb_Row2_e = -1.0;
    } else if (rtb_rxpi_c > 0.0) {
      rtb_Row2_e = 1.0;
    } else if (rtb_rxpi_c == 0.0) {
      rtb_Row2_e = 0.0;
    } else {
      rtb_Row2_e = (rtNaN);
    }

    // End of Signum: '<S169>/x_Sign'

    // Gain: '<S169>/pi'
    rtb_Row2_e *= 3.1415926535897931;

    // Sum: '<S169>/Sum1'
    rtb_rxpi_c += rtb_Row2_e;

    // Math: '<S169>/Math Function' incorporates:
    //   Constant: '<S169>/Constant'

    rtb_rxpi_c = rt_remd_snf(rtb_rxpi_c, 6.2831853071795862);

    // Sum: '<S169>/Sum'
    rtb_rxpi_c -= rtb_Row2_e;

    // Sum: '<S163>/Sum2' incorporates:
    //   Integrator: '<S161>/Integrator3'

    AHV_Model_B.regulation_error_k[0] = AHV_Model_B.x_ref_rel_l -
      AHV_Model_X.Integrator3_CSTATE_l[0];
    AHV_Model_B.regulation_error_k[1] = AHV_Model_B.y_ref_rel_g -
      AHV_Model_X.Integrator3_CSTATE_l[1];
    AHV_Model_B.regulation_error_k[2] = AHV_Model_B.yaw_ref_rel_g - rtb_rxpi_c;

    // Signum: '<S178>/x_Sign' incorporates:
    //   Fcn: '<S163>/yaw_angle'

    if (rtb_rxpi_c < 0.0) {
      rtb_Row2_e = -1.0;
    } else if (rtb_rxpi_c > 0.0) {
      rtb_Row2_e = 1.0;
    } else if (rtb_rxpi_c == 0.0) {
      rtb_Row2_e = 0.0;
    } else {
      rtb_Row2_e = (rtNaN);
    }

    // End of Signum: '<S178>/x_Sign'

    // Gain: '<S178>/pi'
    rtb_Product1_e3 = 3.1415926535897931 * rtb_Row2_e;

    // Sum: '<S178>/Sum' incorporates:
    //   Constant: '<S178>/Constant'
    //   Fcn: '<S163>/yaw_angle'
    //   Math: '<S178>/Math Function'
    //   Sum: '<S178>/Sum1'

    rtb_Row2_e = rt_remd_snf(rtb_rxpi_c + rtb_Product1_e3, 6.2831853071795862) -
      rtb_Product1_e3;

    // Saturate: '<S177>/x_Saturation'
    if (AHV_Model_B.regulation_error_k[2] > 1.0E+10) {
      rtb_Product1_e3 = 1.0E+10;
    } else if (AHV_Model_B.regulation_error_k[2] < -1.0E+10) {
      rtb_Product1_e3 = -1.0E+10;
    } else {
      rtb_Product1_e3 = AHV_Model_B.regulation_error_k[2];
    }

    // End of Saturate: '<S177>/x_Saturation'

    // Signum: '<S177>/x_Sign'
    if (rtb_Product1_e3 < 0.0) {
      rtb_Product2_n = -1.0;
    } else if (rtb_Product1_e3 > 0.0) {
      rtb_Product2_n = 1.0;
    } else if (rtb_Product1_e3 == 0.0) {
      rtb_Product2_n = 0.0;
    } else {
      rtb_Product2_n = (rtNaN);
    }

    // End of Signum: '<S177>/x_Sign'

    // Gain: '<S177>/pi'
    rtb_Product2_n *= 3.1415926535897931;

    // Sum: '<S177>/Sum1'
    rtb_Product1_e3 += rtb_Product2_n;

    // Math: '<S177>/Math Function' incorporates:
    //   Constant: '<S177>/Constant'

    rtb_Product1_e3 = rt_remd_snf(rtb_Product1_e3, 6.2831853071795862);

    // Sum: '<S177>/Sum'
    AHV_Model_B.rxpi_dyz = rtb_Product1_e3 - rtb_Product2_n;

    // Fcn: '<S176>/Row3'
    AHV_Model_B.Row3_h = AHV_Model_B.rxpi_dyz;

    // Fcn: '<S176>/Row1' incorporates:
    //   Fcn: '<S176>/Row2'

    rtb_Product1_e3 = std::sin(rtb_Row2_e);
    rtb_Row2_e = std::cos(rtb_Row2_e);
    rtb_Product2_n = rtb_Row2_e * AHV_Model_B.regulation_error_k[0] +
      rtb_Product1_e3 * AHV_Model_B.regulation_error_k[1];

    // If: '<S181>/If' incorporates:
    //   Abs: '<S181>/Abs'
    //   Inport: '<S184>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_e = static_cast<int8_T>(!(std::abs
        (rtb_Product2_n) <= rtb_Merge3_j));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_e) {
     case 0:
      // Outputs for IfAction SubSystem: '<S181>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S184>/Action Port'

      AHV_Model_B.Merge_j = rtb_Product2_n;

      // End of Outputs for SubSystem: '<S181>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S181>/If Action Subsystem' incorporates:
      //   ActionPort: '<S183>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge3_j, rtb_Product2_n,
        &AHV_Model_B.Merge_j);

      // End of Outputs for SubSystem: '<S181>/If Action Subsystem'
      break;
    }

    // End of If: '<S181>/If'

    // Fcn: '<S176>/Row2'
    rtb_Row2_e = -rtb_Product1_e3 * AHV_Model_B.regulation_error_k[0] +
      rtb_Row2_e * AHV_Model_B.regulation_error_k[1];

    // If: '<S182>/If' incorporates:
    //   Abs: '<S182>/Abs'
    //   Inport: '<S186>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_g = static_cast<int8_T>(!(std::abs
        (rtb_Row2_e) <= rtb_Merge4_l));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_g) {
     case 0:
      // Outputs for IfAction SubSystem: '<S182>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S186>/Action Port'

      AHV_Model_B.Merge_i = rtb_Row2_e;

      // End of Outputs for SubSystem: '<S182>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S182>/If Action Subsystem' incorporates:
      //   ActionPort: '<S185>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge4_l, rtb_Row2_e,
        &AHV_Model_B.Merge_i);

      // End of Outputs for SubSystem: '<S182>/If Action Subsystem'
      break;
    }

    // End of If: '<S182>/If'

    // If: '<S174>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_dx = static_cast<int8_T>
        (!AHV_Model_B.Merge_l);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_dx) {
     case 0:
      // Outputs for IfAction SubSystem: '<S174>/If Action Subsystem' incorporates:
      //   ActionPort: '<S179>/Action Port'

      AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_j, AHV_Model_B.Merge_i,
        AHV_Model_B.Row3_h, rtb_k_u, AHV_Model_ConstP.pooled2);

      // End of Outputs for SubSystem: '<S174>/If Action Subsystem'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S174>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S180>/Action Port'

      AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_j, AHV_Model_B.Merge_i,
        AHV_Model_B.Row3_h, rtb_k_u, AHV_Model_ConstP.pooled1);

      // End of Outputs for SubSystem: '<S174>/If Action Subsystem1'
      break;
    }

    // End of If: '<S174>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S163>/Ki' incorporates:
      //   DiscreteIntegrator: '<S163>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki_k[i] = 0.0;
        AHV_Model_B.Ki_k[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_p[0];
        AHV_Model_B.Ki_k[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_p[1];
        AHV_Model_B.Ki_k[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_p[2];
      }

      // End of Gain: '<S163>/Ki'
    }

    // Sum: '<S163>/Sum3'
    rtb_sintheta[0] = rtb_k_u[0] + AHV_Model_B.Ki_k[0];

    // Integrator: '<S161>/Integrator4'
    rtb_k_u[0] = AHV_Model_X.Integrator4_CSTATE_p[0];

    // Sum: '<S163>/Sum3'
    rtb_sintheta[1] = rtb_k_u[1] + AHV_Model_B.Ki_k[1];

    // Integrator: '<S161>/Integrator4'
    rtb_k_u[1] = AHV_Model_X.Integrator4_CSTATE_p[1];

    // Sum: '<S163>/Sum3'
    rtb_sintheta[2] = rtb_k_u[2] + AHV_Model_B.Ki_k[2];

    // Integrator: '<S161>/Integrator4'
    rtb_k_u[2] = AHV_Model_X.Integrator4_CSTATE_p[2];

    // Sum: '<S163>/Sum1' incorporates:
    //   Gain: '<S163>/Kd'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = rtb_sintheta[i] - ((AHV_Model_ConstP.pooled70[i + 3] *
        rtb_k_u[1] + AHV_Model_ConstP.pooled70[i] * rtb_k_u[0]) +
        AHV_Model_ConstP.pooled70[i + 6] * rtb_k_u[2]);
    }

    // End of Sum: '<S163>/Sum1'

    // Saturate: '<S164>/Surge Force Saturation'
    if (rtb_costheta[0] > 2.0E+6) {
      rtb_costheta_1 = 2.0E+6;
    } else if (rtb_costheta[0] < -2.0E+6) {
      rtb_costheta_1 = -2.0E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[0];
    }

    // End of Saturate: '<S164>/Surge Force Saturation'

    // Product: '<S164>/Product'
    rtb_Row2_e = AHV_Model_B.Merge7_m * rtb_costheta_1;

    // Saturate: '<S164>/Sway Force Saturation'
    if (rtb_costheta[1] > 1.5E+6) {
      rtb_costheta_1 = 1.5E+6;
    } else if (rtb_costheta[1] < -1.5E+6) {
      rtb_costheta_1 = -1.5E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[1];
    }

    // End of Saturate: '<S164>/Sway Force Saturation'

    // Product: '<S164>/Product1'
    rtb_Product1_e3 = AHV_Model_B.Merge8_m * rtb_costheta_1;

    // Saturate: '<S164>/Yaw Moment Saturation'
    if (rtb_costheta[2] > 2.0E+7) {
      rtb_costheta_1 = 2.0E+7;
    } else if (rtb_costheta[2] < -2.0E+7) {
      rtb_costheta_1 = -2.0E+7;
    } else {
      rtb_costheta_1 = rtb_costheta[2];
    }

    // End of Saturate: '<S164>/Yaw Moment Saturation'

    // Product: '<S164>/Product2'
    rtb_Product2_n = AHV_Model_B.Merge9_l * rtb_costheta_1;

    // Product: '<S111>/G*eta' incorporates:
    //   Constant: '<S111>/Spring stiffness'

    for (i = 0; i < 6; i++) {
      rtb_Integrator_f[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Integrator_f[i] += AHV_Model_ConstP.pooled55[6 * i_0 + i] *
          AHV_Model_B.Integrator1_b[i_0];
      }
    }

    // End of Product: '<S111>/G*eta'

    // Sum: '<S111>/Sum1' incorporates:
    //   Constant: '<S108>/Constant'

    tmp_8 = rtb_Integrator_f[0];
    tmp = rtb_Integrator_f[1];
    tmp_0 = rtb_Integrator_f[2];
    tmp_1 = rtb_Integrator_f[3];
    tmp_2 = rtb_Integrator_f[4];
    tmp_3 = rtb_Integrator_f[5];
    rtb_Integrator_f[0] = rtb_Row2_e - tmp_8;
    rtb_Integrator_f[1] = rtb_Product1_e3 - tmp;
    rtb_Integrator_f[2] = 0.0 - tmp_0;
    rtb_Integrator_f[3] = 0.0 - tmp_1;
    rtb_Integrator_f[4] = 0.0 - tmp_2;
    rtb_Integrator_f[5] = rtb_Product2_n - tmp_3;
    for (i = 0; i < 3; i++) {
      // Product: '<S107>/Product1' incorporates:
      //   Inport: '<Root>/tau_cable2'

      rtb_costheta_1 = rtb_TmpSignalConversionAtProd_n[i + 6] * tau_cable2[2] +
        (rtb_TmpSignalConversionAtProd_n[i + 3] * tau_cable2[1] +
         rtb_TmpSignalConversionAtProd_n[i] * tau_cable2[0]);

      // Gain: '<S106>/tau_cableGain1'
      rtb_costheta_0[i] = rtb_costheta_1;

      // Product: '<S107>/Product1'
      rtb_costheta[i] = rtb_costheta_1;
    }

    // Gain: '<S106>/tau_cableGain1' incorporates:
    //   Inport: '<Root>/ahv_fairlead2'
    //   Product: '<S158>/ae'
    //   Product: '<S158>/af'
    //   Product: '<S158>/bd'
    //   Product: '<S158>/bf'
    //   Product: '<S158>/cd'
    //   Product: '<S158>/ce'
    //   Sum: '<S158>/Sum'
    //   Sum: '<S158>/Sum1'
    //   Sum: '<S158>/Sum2'

    rtb_costheta_0[3] = ahv_fairlead2[1] * rtb_costheta[2] - ahv_fairlead2[2] *
      rtb_costheta[1];
    rtb_costheta_0[4] = ahv_fairlead2[2] * rtb_costheta[0] - ahv_fairlead2[0] *
      rtb_costheta[2];
    rtb_costheta_0[5] = ahv_fairlead2[0] * rtb_costheta[1] - ahv_fairlead2[1] *
      rtb_costheta[0];

    // Sum: '<S111>/Sum7' incorporates:
    //   Gain: '<S106>/tau_cableGain1'
    //   Sum: '<S106>/Sum'
    //   Sum: '<S111>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Integrator_f[i] = ((rtb_tau_WF_em[i] + rtb_tau_WD_n[i]) +
        rtb_Integrator_f[i]) + AHV_Model_ConstP.pooled77[i] * rtb_costheta_0[i];
    }

    // End of Sum: '<S111>/Sum7'

    // Sum: '<S111>/Sum5' incorporates:
    //   Abs: '<S117>/Abs1'
    //   Gain: '<S117>/Gain3'
    //   Product: '<S117>/Product7'
    //   Product: '<S117>/Product8'

    tmp_8 = rtb_Integrator_f[0];
    tmp = rtb_Integrator_f[1];
    tmp_0 = rtb_Integrator_f[2];
    tmp_1 = rtb_Integrator_f[3];
    tmp_2 = rtb_Integrator_f[4];
    tmp_3 = rtb_Integrator_f[5];
    rtb_Integrator_f[0] = -(std::abs(AHV_Model_B.nu_r_j[0]) *
      AHV_Model_B.nu_r_j[0] * AHV_Model_ConstB.Product5_h) + tmp_8;
    rtb_Integrator_f[1] = tmp + AHV_Model_B.dx1_a0;
    rtb_Integrator_f[2] = tmp_0;
    rtb_Integrator_f[3] = tmp_1;
    rtb_Integrator_f[4] = tmp_2;
    rtb_Integrator_f[5] = tmp_3 + AHV_Model_B.dx2_f;

    // Sum: '<S111>/Sum4' incorporates:
    //   Constant: '<S111>/damping'
    //   Product: '<S111>/D*eta '

    for (i = 0; i < 6; i++) {
      tmp_8 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        tmp_8 += AHV_Model_ConstP.pooled58[6 * i_0 + i] * AHV_Model_B.nu_r_j[i_0];
      }

      rtb_Integrator_f[i] -= tmp_8;
    }

    // End of Sum: '<S111>/Sum4'

    // Product: '<S130>/Product2' incorporates:
    //   Integrator: '<S130>/Integrator'

    rtb_Sum_j = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Sum_j += AHV_Model_ConstB.MathFunction_l[i] *
        AHV_Model_X.Integrator_CSTATE_k[i];
    }

    // Sum: '<S111>/Sum' incorporates:
    //   Constant: '<S130>/B44_inf'
    //   Constant: '<S130>/D44'
    //   Product: '<S130>/Product2'
    //   Product: '<S130>/Product3'
    //   Product: '<S130>/Product4'
    //   StateSpace: '<S118>/Dp(1,1)'
    //   StateSpace: '<S118>/Dp(2,2)'
    //   StateSpace: '<S118>/Dp(2,4)'
    //   StateSpace: '<S118>/Dp(2,6)'
    //   StateSpace: '<S118>/Dp(3,3)'
    //   StateSpace: '<S118>/Dp(3,5)'
    //   StateSpace: '<S118>/Dp(4,2)'
    //   StateSpace: '<S118>/Dp(4,6)'
    //   StateSpace: '<S118>/Dp(5,3)'
    //   StateSpace: '<S118>/Dp(5,5)'
    //   StateSpace: '<S118>/Dp(6,2)'
    //   StateSpace: '<S118>/Dp(6,4)'
    //   StateSpace: '<S118>/Dp(6,6)'
    //   Sum: '<S118>/Sum'
    //   Sum: '<S118>/Sum1'
    //   Sum: '<S118>/Sum3'
    //   Sum: '<S118>/Sum4'
    //   Sum: '<S118>/Sum5'
    //   Sum: '<S130>/Sum1'
    //   Sum: '<S130>/Sum2'

    tmp_8 = rtb_Integrator_f[0];
    tmp = rtb_Integrator_f[1];
    tmp_0 = rtb_Integrator_f[2];
    tmp_1 = rtb_Integrator_f[3];
    tmp_2 = rtb_Integrator_f[4];
    tmp_3 = rtb_Integrator_f[5];
    rtb_Integrator_f[0] = tmp_8 - (((((-192192.528024618 *
      AHV_Model_X.Dp11_CSTATE_o[0] + -4350.6445868200408 *
      AHV_Model_X.Dp11_CSTATE_o[1]) + -26196.671588656758 *
      AHV_Model_X.Dp11_CSTATE_o[2]) + -69550.629589629229 *
      AHV_Model_X.Dp11_CSTATE_o[3]) + -6667.9011063139569 *
      AHV_Model_X.Dp11_CSTATE_o[4]) + 28360.700724836512 * AHV_Model_B.nu_r_j[0]);
    rtb_Integrator_f[1] = tmp - (((((((-1.9219252802461751E+6 *
      AHV_Model_X.Dp22_CSTATE_jf[0] + -43506.445868266834 *
      AHV_Model_X.Dp22_CSTATE_jf[1]) + -261966.71588657616 *
      AHV_Model_X.Dp22_CSTATE_jf[2]) + -695506.29589629406 *
      AHV_Model_X.Dp22_CSTATE_jf[3]) + -66679.0110631567 *
      AHV_Model_X.Dp22_CSTATE_jf[4]) + 283607.00724836387 * AHV_Model_B.nu_r_j[1])
      + (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_j[0] +
             -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_j[1]) +
            782719.9224317756 * AHV_Model_X.Dp24_CSTATE_j[2]) +
           906439.02406890341 * AHV_Model_X.Dp24_CSTATE_j[3]) +
          48702.05056677026 * AHV_Model_X.Dp24_CSTATE_j[4]) + 342290.49224698928
         * AHV_Model_B.nu_r_j[3])) + (((((-377217.71116477274 *
      AHV_Model_X.Dp26_CSTATE_d[0] + 9.7759770791847613E+6 *
      AHV_Model_X.Dp26_CSTATE_d[1]) + 1.2719330237464791E+6 *
      AHV_Model_X.Dp26_CSTATE_d[2]) + 3.67278239005485E+6 *
      AHV_Model_X.Dp26_CSTATE_d[3]) + -523575.51382008579 *
      AHV_Model_X.Dp26_CSTATE_d[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_j[5]));
    rtb_Integrator_f[2] = tmp_0 - ((((((-2.8230572248621574E+6 *
      AHV_Model_X.Dp33_CSTATE_g[0] + 2971.0437052436932 *
      AHV_Model_X.Dp33_CSTATE_g[1]) + -112975.07052200216 *
      AHV_Model_X.Dp33_CSTATE_g[2]) + -459867.38057091518 *
      AHV_Model_X.Dp33_CSTATE_g[3]) + 35755.997807249405 *
      AHV_Model_X.Dp33_CSTATE_g[4]) + 425986.43085235963 * AHV_Model_B.nu_r_j[2])
      + (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE_e[0] +
             -189169.67892356269 * AHV_Model_X.Dp35_CSTATE_e[1]) +
            1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE_e[2]) +
           1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE_e[3]) +
          2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE_e[4]) +
         6.2116453321941122E+6 * AHV_Model_B.nu_r_j[4]));
    rtb_Integrator_f[3] = tmp_1 - (((((((163377.73376110662 *
      AHV_Model_X.Dp42_CSTATE_g[0] + -2.24222697524951E+6 *
      AHV_Model_X.Dp42_CSTATE_g[1]) + 782719.9224317756 *
      AHV_Model_X.Dp42_CSTATE_g[2]) + 906439.02406890341 *
      AHV_Model_X.Dp42_CSTATE_g[3]) + 48702.05056677026 *
      AHV_Model_X.Dp42_CSTATE_g[4]) + 342290.49224698928 * AHV_Model_B.nu_r_j[1])
      + ((2.3074903324953854E+7 * AHV_Model_B.nu_r_j[3] + rtb_Sum_j) +
         4.97124674352808E+7 * AHV_Model_B.nu_r_j[3])) +
      (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE_e[0] +
           2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE_e[1]) +
          -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE_e[2]) +
         1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE_e[3]) +
        -70143.531785937739 * AHV_Model_X.Dp46_CSTATE_e[4]) +
       -8.0522169282609783E+6 * AHV_Model_B.nu_r_j[5]));
    rtb_Integrator_f[4] = tmp_2 - ((((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp53_CSTATE_k[0] + -189169.67892356269 *
      AHV_Model_X.Dp53_CSTATE_k[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp53_CSTATE_k[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp53_CSTATE_k[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp53_CSTATE_k[4]) + 6.2116453321941122E+6 *
      AHV_Model_B.nu_r_j[2]) + (((((-1.3510265271416564E+7 *
      AHV_Model_X.Dp55_CSTATE_i[0] + 2.1371969925573766E+9 *
      AHV_Model_X.Dp55_CSTATE_i[1]) + 9.8710379288756877E+7 *
      AHV_Model_X.Dp55_CSTATE_i[2]) + 6.1487736935078049E+8 *
      AHV_Model_X.Dp55_CSTATE_i[3]) + 5.7220438213245332E+8 *
      AHV_Model_X.Dp55_CSTATE_i[4]) + 3.0813057877131975E+8 *
      AHV_Model_B.nu_r_j[4]));
    rtb_Integrator_f[5] = tmp_3 - (((((((-377217.71116477274 *
      AHV_Model_X.Dp62_CSTATE_h[0] + 9.7759770791847613E+6 *
      AHV_Model_X.Dp62_CSTATE_h[1]) + 1.2719330237464791E+6 *
      AHV_Model_X.Dp62_CSTATE_h[2]) + 3.67278239005485E+6 *
      AHV_Model_X.Dp62_CSTATE_h[3]) + -523575.51382008579 *
      AHV_Model_X.Dp62_CSTATE_h[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_j[1]) + (((((-1.5746048421550707E+6 *
      AHV_Model_X.Dp64_CSTATE_h[0] + 2.9884305406230576E+7 *
      AHV_Model_X.Dp64_CSTATE_h[1]) + -4.3953040127794016E+6 *
      AHV_Model_X.Dp64_CSTATE_h[2]) + 1.205872403100216E+7 *
      AHV_Model_X.Dp64_CSTATE_h[3]) + -70143.531785937739 *
      AHV_Model_X.Dp64_CSTATE_h[4]) + -8.0522169282609783E+6 *
      AHV_Model_B.nu_r_j[3])) + (((((-8.07757267677547E+8 *
      AHV_Model_X.Dp66_CSTATE_n[0] + -1.0798841515040772E+7 *
      AHV_Model_X.Dp66_CSTATE_n[1]) + -1.27285250509213E+8 *
      AHV_Model_X.Dp66_CSTATE_n[2]) + -3.0027175358715254E+8 *
      AHV_Model_X.Dp66_CSTATE_n[3]) + 1.7969070654437996E+7 *
      AHV_Model_X.Dp66_CSTATE_n[4]) + 1.1943182502908507E+8 *
      AHV_Model_B.nu_r_j[5]));

    // Product: '<S111>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_p, rtb_Integrator_f,
      AHV_Model_B.Minvtau_f);

    // Outport: '<Root>/nu2' incorporates:
    //   Integrator: '<S111>/Integrator'

    for (i = 0; i < 6; i++) {
      nu2[i] = AHV_Model_X.Integrator_CSTATE_o[i];
    }

    // End of Outport: '<Root>/nu2'

    // Gain: '<S14>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_Sum_j = 0.017453292519943295 * Current_direction;

    // Product: '<S14>/Product' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S14>/cos'

    rtb_rxpi_n = std::cos(rtb_Sum_j) * Current_speed;

    // Product: '<S14>/Product1' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Trigonometry: '<S14>/sin'

    rtb_Sum_j = std::sin(rtb_Sum_j) * Current_speed;

    // Sum: '<S15>/Sum6' incorporates:
    //   Gain: '<S10>/Current on//off '
    //   Integrator: '<S15>/Integrator'

    AHV_Model_B.nu_r_c[0] = AHV_Model_X.Integrator_CSTATE_e[0] - rtb_rxpi_n;
    AHV_Model_B.nu_r_c[1] = AHV_Model_X.Integrator_CSTATE_e[1] - rtb_Sum_j;
    AHV_Model_B.nu_r_c[2] = AHV_Model_X.Integrator_CSTATE_e[2];
    AHV_Model_B.nu_r_c[3] = AHV_Model_X.Integrator_CSTATE_e[3];
    AHV_Model_B.nu_r_c[4] = AHV_Model_X.Integrator_CSTATE_e[4];
    AHV_Model_B.nu_r_c[5] = AHV_Model_X.Integrator_CSTATE_e[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S21>/Cross-flow drag trapezoidal integration' 
      // Constant: '<S21>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx,
        AHV_Model_B.nu_r_c[1], AHV_Model_B.nu_r_c[5], &AHV_Model_B.Sum_f,
        &AHV_Model_B.Sum2_n, 82.800003);

      // End of Outputs for SubSystem: '<S21>/Cross-flow drag trapezoidal integration' 

      // Product: '<S21>/dx1' incorporates:
      //   Constant: '<S21>/2D drag coefficient '
      //   Constant: '<S21>/Transversal area//Lpp'
      //   Constant: '<S21>/rho'
      //   Gain: '<S21>/Gain1'

      AHV_Model_B.dx1_o = -0.5 * AHV_Model_B.Sum_f * 1025.0 * 5.4 *
        0.69954741847648816;

      // Product: '<S21>/dx2' incorporates:
      //   Constant: '<S21>/2D drag coefficient '
      //   Constant: '<S21>/Transversal area//Lpp'
      //   Constant: '<S21>/rho'
      //   Gain: '<S21>/Gain'

      AHV_Model_B.dx2_g = -0.5 * AHV_Model_B.Sum2_n * 1025.0 * 5.4 *
        0.69954741847648816;
    }

    // Integrator: '<S65>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init1'
    //   SignalConversion generated from: '<S65>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_l != 0) {
      AHV_Model_X.Integrator3_CSTATE_n[0] = Vessel_init1[0];
      AHV_Model_X.Integrator3_CSTATE_n[1] = Vessel_init1[1];
      AHV_Model_X.Integrator3_CSTATE_n[2] = Vessel_init1[5];
    }

    // Saturate: '<S73>/x_Saturation' incorporates:
    //   Integrator: '<S65>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_n[2] > 1.0E+10) {
      rtb_rxpi_n = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_n[2] < -1.0E+10) {
      rtb_rxpi_n = -1.0E+10;
    } else {
      rtb_rxpi_n = AHV_Model_X.Integrator3_CSTATE_n[2];
    }

    // End of Saturate: '<S73>/x_Saturation'

    // Signum: '<S73>/x_Sign'
    if (rtb_rxpi_n < 0.0) {
      rtb_Sum_j = -1.0;
    } else if (rtb_rxpi_n > 0.0) {
      rtb_Sum_j = 1.0;
    } else if (rtb_rxpi_n == 0.0) {
      rtb_Sum_j = 0.0;
    } else {
      rtb_Sum_j = (rtNaN);
    }

    // End of Signum: '<S73>/x_Sign'

    // Gain: '<S73>/pi'
    rtb_Sum_j *= 3.1415926535897931;

    // Sum: '<S73>/Sum1'
    rtb_rxpi_n += rtb_Sum_j;

    // Math: '<S73>/Math Function' incorporates:
    //   Constant: '<S73>/Constant'

    rtb_rxpi_n = rt_remd_snf(rtb_rxpi_n, 6.2831853071795862);

    // Sum: '<S73>/Sum'
    rtb_rxpi_n -= rtb_Sum_j;

    // Sum: '<S67>/Sum2' incorporates:
    //   Integrator: '<S65>/Integrator3'

    AHV_Model_B.regulation_error_i[0] = AHV_Model_B.x_ref_rel_m -
      AHV_Model_X.Integrator3_CSTATE_n[0];
    AHV_Model_B.regulation_error_i[1] = AHV_Model_B.y_ref_rel_d -
      AHV_Model_X.Integrator3_CSTATE_n[1];
    AHV_Model_B.regulation_error_i[2] = AHV_Model_B.yaw_ref_rel_c - rtb_rxpi_n;

    // Signum: '<S82>/x_Sign' incorporates:
    //   Fcn: '<S67>/yaw_angle'

    if (rtb_rxpi_n < 0.0) {
      rtb_Sum_j = -1.0;
    } else if (rtb_rxpi_n > 0.0) {
      rtb_Sum_j = 1.0;
    } else if (rtb_rxpi_n == 0.0) {
      rtb_Sum_j = 0.0;
    } else {
      rtb_Sum_j = (rtNaN);
    }

    // End of Signum: '<S82>/x_Sign'

    // Gain: '<S82>/pi'
    rtb_Product1_c = 3.1415926535897931 * rtb_Sum_j;

    // Sum: '<S82>/Sum' incorporates:
    //   Constant: '<S82>/Constant'
    //   Fcn: '<S67>/yaw_angle'
    //   Math: '<S82>/Math Function'
    //   Sum: '<S82>/Sum1'

    rtb_Sum_j = rt_remd_snf(rtb_rxpi_n + rtb_Product1_c, 6.2831853071795862) -
      rtb_Product1_c;

    // Saturate: '<S81>/x_Saturation'
    if (AHV_Model_B.regulation_error_i[2] > 1.0E+10) {
      rtb_Product1_c = 1.0E+10;
    } else if (AHV_Model_B.regulation_error_i[2] < -1.0E+10) {
      rtb_Product1_c = -1.0E+10;
    } else {
      rtb_Product1_c = AHV_Model_B.regulation_error_i[2];
    }

    // End of Saturate: '<S81>/x_Saturation'

    // Signum: '<S81>/x_Sign'
    if (rtb_Product1_c < 0.0) {
      rtb_Product2_a = -1.0;
    } else if (rtb_Product1_c > 0.0) {
      rtb_Product2_a = 1.0;
    } else if (rtb_Product1_c == 0.0) {
      rtb_Product2_a = 0.0;
    } else {
      rtb_Product2_a = (rtNaN);
    }

    // End of Signum: '<S81>/x_Sign'

    // Gain: '<S81>/pi'
    rtb_Product2_a *= 3.1415926535897931;

    // Sum: '<S81>/Sum1'
    rtb_Product1_c += rtb_Product2_a;

    // Math: '<S81>/Math Function' incorporates:
    //   Constant: '<S81>/Constant'

    rtb_Product1_c = rt_remd_snf(rtb_Product1_c, 6.2831853071795862);

    // Sum: '<S81>/Sum'
    AHV_Model_B.rxpi_e = rtb_Product1_c - rtb_Product2_a;

    // Fcn: '<S80>/Row3'
    AHV_Model_B.Row3_l = AHV_Model_B.rxpi_e;

    // Fcn: '<S80>/Row1' incorporates:
    //   Fcn: '<S80>/Row2'

    rtb_Product1_c = std::sin(rtb_Sum_j);
    rtb_Sum_j = std::cos(rtb_Sum_j);
    rtb_Product2_a = rtb_Sum_j * AHV_Model_B.regulation_error_i[0] +
      rtb_Product1_c * AHV_Model_B.regulation_error_i[1];

    // If: '<S85>/If' incorporates:
    //   Abs: '<S85>/Abs'
    //   Inport: '<S88>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_m5 = static_cast<int8_T>(!(std::abs
        (rtb_Product2_a) <= rtb_Merge3_c));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_m5) {
     case 0:
      // Outputs for IfAction SubSystem: '<S85>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S88>/Action Port'

      AHV_Model_B.Merge_k = rtb_Product2_a;

      // End of Outputs for SubSystem: '<S85>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S85>/If Action Subsystem' incorporates:
      //   ActionPort: '<S87>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge3_c, rtb_Product2_a,
        &AHV_Model_B.Merge_k);

      // End of Outputs for SubSystem: '<S85>/If Action Subsystem'
      break;
    }

    // End of If: '<S85>/If'

    // Fcn: '<S80>/Row2'
    rtb_Sum_j = -rtb_Product1_c * AHV_Model_B.regulation_error_i[0] + rtb_Sum_j *
      AHV_Model_B.regulation_error_i[1];

    // If: '<S86>/If' incorporates:
    //   Abs: '<S86>/Abs'
    //   Inport: '<S90>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_ey = static_cast<int8_T>(!(std::abs
        (rtb_Sum_j) <= rtb_Merge4_g));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_ey) {
     case 0:
      // Outputs for IfAction SubSystem: '<S86>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S90>/Action Port'

      AHV_Model_B.Merge_b = rtb_Sum_j;

      // End of Outputs for SubSystem: '<S86>/If Action Subsystem1'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S86>/If Action Subsystem' incorporates:
      //   ActionPort: '<S89>/Action Port'

      AHV_Model_IfActionSubsystem_f(rtb_Merge4_g, rtb_Sum_j,
        &AHV_Model_B.Merge_b);

      // End of Outputs for SubSystem: '<S86>/If Action Subsystem'
      break;
    }

    // End of If: '<S86>/If'

    // If: '<S78>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_h0 = static_cast<int8_T>
        (!AHV_Model_B.Merge_g);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_h0) {
     case 0:
      // Outputs for IfAction SubSystem: '<S78>/If Action Subsystem' incorporates:
      //   ActionPort: '<S83>/Action Port'

      AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_k, AHV_Model_B.Merge_b,
        AHV_Model_B.Row3_l, rtb_v, AHV_Model_ConstP.pooled2);

      // End of Outputs for SubSystem: '<S78>/If Action Subsystem'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S78>/If Action Subsystem1' incorporates:
      //   ActionPort: '<S84>/Action Port'

      AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_k, AHV_Model_B.Merge_b,
        AHV_Model_B.Row3_l, rtb_v, AHV_Model_ConstP.pooled1);

      // End of Outputs for SubSystem: '<S78>/If Action Subsystem1'
      break;
    }

    // End of If: '<S78>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S67>/Ki' incorporates:
      //   DiscreteIntegrator: '<S67>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki_m[i] = 0.0;
        AHV_Model_B.Ki_m[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_a[0];
        AHV_Model_B.Ki_m[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_a[1];
        AHV_Model_B.Ki_m[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_a[2];
      }

      // End of Gain: '<S67>/Ki'
    }

    // Sum: '<S67>/Sum3'
    rtb_sintheta[0] = rtb_v[0] + AHV_Model_B.Ki_m[0];

    // Integrator: '<S65>/Integrator4'
    rtb_v[0] = AHV_Model_X.Integrator4_CSTATE_c[0];

    // Sum: '<S67>/Sum3'
    rtb_sintheta[1] = rtb_v[1] + AHV_Model_B.Ki_m[1];

    // Integrator: '<S65>/Integrator4'
    rtb_v[1] = AHV_Model_X.Integrator4_CSTATE_c[1];

    // Sum: '<S67>/Sum3'
    rtb_sintheta[2] = rtb_v[2] + AHV_Model_B.Ki_m[2];

    // Integrator: '<S65>/Integrator4'
    rtb_v[2] = AHV_Model_X.Integrator4_CSTATE_c[2];

    // Sum: '<S67>/Sum1' incorporates:
    //   Gain: '<S67>/Kd'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = rtb_sintheta[i] - ((AHV_Model_ConstP.pooled70[i + 3] *
        rtb_v[1] + AHV_Model_ConstP.pooled70[i] * rtb_v[0]) +
        AHV_Model_ConstP.pooled70[i + 6] * rtb_v[2]);
    }

    // End of Sum: '<S67>/Sum1'

    // Saturate: '<S68>/Surge Force Saturation'
    if (rtb_costheta[0] > 2.0E+6) {
      rtb_costheta_1 = 2.0E+6;
    } else if (rtb_costheta[0] < -2.0E+6) {
      rtb_costheta_1 = -2.0E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[0];
    }

    // End of Saturate: '<S68>/Surge Force Saturation'

    // Product: '<S68>/Product'
    rtb_Sum_j = AHV_Model_B.Merge7_e * rtb_costheta_1;

    // Saturate: '<S68>/Sway Force Saturation'
    if (rtb_costheta[1] > 1.5E+6) {
      rtb_costheta_1 = 1.5E+6;
    } else if (rtb_costheta[1] < -1.5E+6) {
      rtb_costheta_1 = -1.5E+6;
    } else {
      rtb_costheta_1 = rtb_costheta[1];
    }

    // End of Saturate: '<S68>/Sway Force Saturation'

    // Product: '<S68>/Product1'
    rtb_Product1_c = AHV_Model_B.Merge8_mo * rtb_costheta_1;

    // Saturate: '<S68>/Yaw Moment Saturation'
    if (rtb_costheta[2] > 2.0E+7) {
      rtb_costheta_1 = 2.0E+7;
    } else if (rtb_costheta[2] < -2.0E+7) {
      rtb_costheta_1 = -2.0E+7;
    } else {
      rtb_costheta_1 = rtb_costheta[2];
    }

    // End of Saturate: '<S68>/Yaw Moment Saturation'

    // Product: '<S68>/Product2'
    rtb_Product2_a = AHV_Model_B.Merge9_p * rtb_costheta_1;

    // Product: '<S15>/G*eta' incorporates:
    //   Constant: '<S15>/Spring stiffness'

    for (i = 0; i < 6; i++) {
      rtb_Sum[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Sum[i] += AHV_Model_ConstP.pooled55[6 * i_0 + i] *
          AHV_Model_B.Integrator1[i_0];
      }
    }

    // End of Product: '<S15>/G*eta'

    // Sum: '<S15>/Sum1' incorporates:
    //   Constant: '<S12>/Constant'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = rtb_Sum_j - tmp_8;
    rtb_Sum[1] = rtb_Product1_c - tmp;
    rtb_Sum[2] = 0.0 - tmp_0;
    rtb_Sum[3] = 0.0 - tmp_1;
    rtb_Sum[4] = 0.0 - tmp_2;
    rtb_Sum[5] = rtb_Product2_a - tmp_3;
    for (i = 0; i < 3; i++) {
      // Product: '<S11>/Product1' incorporates:
      //   Inport: '<Root>/tau_cable1'

      rtb_costheta_1 = rtb_TmpSignalConversionAtProd_i[i + 6] * tau_cable1[2] +
        (rtb_TmpSignalConversionAtProd_i[i + 3] * tau_cable1[1] +
         rtb_TmpSignalConversionAtProd_i[i] * tau_cable1[0]);

      // Gain: '<S10>/tau_cableGain1'
      rtb_costheta_0[i] = rtb_costheta_1;

      // Product: '<S11>/Product1'
      rtb_costheta[i] = rtb_costheta_1;
    }

    // Gain: '<S10>/tau_cableGain1' incorporates:
    //   Inport: '<Root>/ahv_fairlead1'
    //   Product: '<S62>/ae'
    //   Product: '<S62>/af'
    //   Product: '<S62>/bd'
    //   Product: '<S62>/bf'
    //   Product: '<S62>/cd'
    //   Product: '<S62>/ce'
    //   Sum: '<S62>/Sum'
    //   Sum: '<S62>/Sum1'
    //   Sum: '<S62>/Sum2'

    rtb_costheta_0[3] = ahv_fairlead1[1] * rtb_costheta[2] - ahv_fairlead1[2] *
      rtb_costheta[1];
    rtb_costheta_0[4] = ahv_fairlead1[2] * rtb_costheta[0] - ahv_fairlead1[0] *
      rtb_costheta[2];
    rtb_costheta_0[5] = ahv_fairlead1[0] * rtb_costheta[1] - ahv_fairlead1[1] *
      rtb_costheta[0];

    // Sum: '<S15>/Sum7' incorporates:
    //   Gain: '<S10>/tau_cableGain1'
    //   Sum: '<S10>/Sum'
    //   Sum: '<S15>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Sum[i] = ((rtb_tau_WF_o[i] + rtb_tau_WD_p[i]) + rtb_Sum[i]) +
        AHV_Model_ConstP.pooled77[i] * rtb_costheta_0[i];
    }

    // End of Sum: '<S15>/Sum7'

    // Sum: '<S15>/Sum5' incorporates:
    //   Abs: '<S21>/Abs1'
    //   Gain: '<S21>/Gain3'
    //   Product: '<S21>/Product7'
    //   Product: '<S21>/Product8'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = -(std::abs(AHV_Model_B.nu_r_c[0]) * AHV_Model_B.nu_r_c[0] *
                   AHV_Model_ConstB.Product5) + tmp_8;
    rtb_Sum[1] = tmp + AHV_Model_B.dx1_o;
    rtb_Sum[2] = tmp_0;
    rtb_Sum[3] = tmp_1;
    rtb_Sum[4] = tmp_2;
    rtb_Sum[5] = tmp_3 + AHV_Model_B.dx2_g;

    // Sum: '<S15>/Sum4' incorporates:
    //   Constant: '<S15>/damping'
    //   Product: '<S15>/D*eta '

    for (i = 0; i < 6; i++) {
      tmp_8 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        tmp_8 += AHV_Model_ConstP.pooled58[6 * i_0 + i] * AHV_Model_B.nu_r_c[i_0];
      }

      rtb_Sum[i] -= tmp_8;
    }

    // End of Sum: '<S15>/Sum4'

    // Product: '<S34>/Product2' incorporates:
    //   Integrator: '<S34>/Integrator'

    rtb_costheta_1 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_costheta_1 += AHV_Model_ConstB.MathFunction[i] *
        AHV_Model_X.Integrator_CSTATE_p[i];
    }

    // Sum: '<S15>/Sum' incorporates:
    //   Constant: '<S34>/B44_inf'
    //   Constant: '<S34>/D44'
    //   Product: '<S34>/Product2'
    //   Product: '<S34>/Product3'
    //   Product: '<S34>/Product4'
    //   StateSpace: '<S22>/Dp(1,1)'
    //   StateSpace: '<S22>/Dp(2,2)'
    //   StateSpace: '<S22>/Dp(2,4)'
    //   StateSpace: '<S22>/Dp(2,6)'
    //   StateSpace: '<S22>/Dp(3,3)'
    //   StateSpace: '<S22>/Dp(3,5)'
    //   StateSpace: '<S22>/Dp(4,2)'
    //   StateSpace: '<S22>/Dp(4,6)'
    //   StateSpace: '<S22>/Dp(5,3)'
    //   StateSpace: '<S22>/Dp(5,5)'
    //   StateSpace: '<S22>/Dp(6,2)'
    //   StateSpace: '<S22>/Dp(6,4)'
    //   StateSpace: '<S22>/Dp(6,6)'
    //   Sum: '<S22>/Sum'
    //   Sum: '<S22>/Sum1'
    //   Sum: '<S22>/Sum3'
    //   Sum: '<S22>/Sum4'
    //   Sum: '<S22>/Sum5'
    //   Sum: '<S34>/Sum1'
    //   Sum: '<S34>/Sum2'

    tmp_8 = rtb_Sum[0];
    tmp = rtb_Sum[1];
    tmp_0 = rtb_Sum[2];
    tmp_1 = rtb_Sum[3];
    tmp_2 = rtb_Sum[4];
    tmp_3 = rtb_Sum[5];
    rtb_Sum[0] = tmp_8 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE_d[0] +
      -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE_d[1]) + -26196.671588656758 *
      AHV_Model_X.Dp11_CSTATE_d[2]) + -69550.629589629229 *
      AHV_Model_X.Dp11_CSTATE_d[3]) + -6667.9011063139569 *
      AHV_Model_X.Dp11_CSTATE_d[4]) + 28360.700724836512 * AHV_Model_B.nu_r_c[0]);
    rtb_Sum[1] = tmp - (((((((-1.9219252802461751E+6 *
      AHV_Model_X.Dp22_CSTATE_h[0] + -43506.445868266834 *
      AHV_Model_X.Dp22_CSTATE_h[1]) + -261966.71588657616 *
      AHV_Model_X.Dp22_CSTATE_h[2]) + -695506.29589629406 *
      AHV_Model_X.Dp22_CSTATE_h[3]) + -66679.0110631567 *
      AHV_Model_X.Dp22_CSTATE_h[4]) + 283607.00724836387 * AHV_Model_B.nu_r_c[1])
                         + (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_c[0]
      + -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_c[1]) + 782719.9224317756
      * AHV_Model_X.Dp24_CSTATE_c[2]) + 906439.02406890341 *
      AHV_Model_X.Dp24_CSTATE_c[3]) + 48702.05056677026 *
      AHV_Model_X.Dp24_CSTATE_c[4]) + 342290.49224698928 * AHV_Model_B.nu_r_c[3]))
                        + (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE_h[0]
      + 9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE_h[1]) +
      1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE_h[2]) +
      3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE_h[3]) + -523575.51382008579 *
      AHV_Model_X.Dp26_CSTATE_h[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_c[5]));
    rtb_Sum[2] = tmp_0 - ((((((-2.8230572248621574E+6 *
      AHV_Model_X.Dp33_CSTATE_o[0] + 2971.0437052436932 *
      AHV_Model_X.Dp33_CSTATE_o[1]) + -112975.07052200216 *
      AHV_Model_X.Dp33_CSTATE_o[2]) + -459867.38057091518 *
      AHV_Model_X.Dp33_CSTATE_o[3]) + 35755.997807249405 *
      AHV_Model_X.Dp33_CSTATE_o[4]) + 425986.43085235963 * AHV_Model_B.nu_r_c[2])
                          + (((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp35_CSTATE_h[0] + -189169.67892356269 *
      AHV_Model_X.Dp35_CSTATE_h[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp35_CSTATE_h[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp35_CSTATE_h[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp35_CSTATE_h[4]) + 6.2116453321941122E+6 *
      AHV_Model_B.nu_r_c[4]));
    rtb_Sum[3] = tmp_1 - (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE_e[0]
      + -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE_e[1]) + 782719.9224317756
      * AHV_Model_X.Dp42_CSTATE_e[2]) + 906439.02406890341 *
      AHV_Model_X.Dp42_CSTATE_e[3]) + 48702.05056677026 *
      AHV_Model_X.Dp42_CSTATE_e[4]) + 342290.49224698928 * AHV_Model_B.nu_r_c[1])
      + ((2.3074903324953854E+7 * AHV_Model_B.nu_r_c[3] + rtb_costheta_1) +
         4.97124674352808E+7 * AHV_Model_B.nu_r_c[3])) +
                          (((((-1.5746048421550707E+6 *
      AHV_Model_X.Dp46_CSTATE_h[0] + 2.9884305406230576E+7 *
      AHV_Model_X.Dp46_CSTATE_h[1]) + -4.3953040127794016E+6 *
      AHV_Model_X.Dp46_CSTATE_h[2]) + 1.205872403100216E+7 *
      AHV_Model_X.Dp46_CSTATE_h[3]) + -70143.531785937739 *
      AHV_Model_X.Dp46_CSTATE_h[4]) + -8.0522169282609783E+6 *
      AHV_Model_B.nu_r_c[5]));
    rtb_Sum[4] = tmp_2 - ((((((-4.3928932046187744E+7 *
      AHV_Model_X.Dp53_CSTATE_m[0] + -189169.67892356269 *
      AHV_Model_X.Dp53_CSTATE_m[1]) + 1.5146159421991596E+7 *
      AHV_Model_X.Dp53_CSTATE_m[2]) + 1.2040006787420809E+7 *
      AHV_Model_X.Dp53_CSTATE_m[3]) + 2.1294858953942191E+6 *
      AHV_Model_X.Dp53_CSTATE_m[4]) + 6.2116453321941122E+6 *
      AHV_Model_B.nu_r_c[2]) + (((((-1.3510265271416564E+7 *
      AHV_Model_X.Dp55_CSTATE_n[0] + 2.1371969925573766E+9 *
      AHV_Model_X.Dp55_CSTATE_n[1]) + 9.8710379288756877E+7 *
      AHV_Model_X.Dp55_CSTATE_n[2]) + 6.1487736935078049E+8 *
      AHV_Model_X.Dp55_CSTATE_n[3]) + 5.7220438213245332E+8 *
      AHV_Model_X.Dp55_CSTATE_n[4]) + 3.0813057877131975E+8 *
      AHV_Model_B.nu_r_c[4]));
    rtb_Sum[5] = tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE_a
      [0] + 9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE_a[1]) +
      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE_a[2]) +
      3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE_a[3]) + -523575.51382008579 *
      AHV_Model_X.Dp62_CSTATE_a[4]) + 1.4324880870688609E+6 *
      AHV_Model_B.nu_r_c[1]) + (((((-1.5746048421550707E+6 *
      AHV_Model_X.Dp64_CSTATE_n[0] + 2.9884305406230576E+7 *
      AHV_Model_X.Dp64_CSTATE_n[1]) + -4.3953040127794016E+6 *
      AHV_Model_X.Dp64_CSTATE_n[2]) + 1.205872403100216E+7 *
      AHV_Model_X.Dp64_CSTATE_n[3]) + -70143.531785937739 *
      AHV_Model_X.Dp64_CSTATE_n[4]) + -8.0522169282609783E+6 *
      AHV_Model_B.nu_r_c[3])) + (((((-8.07757267677547E+8 *
      AHV_Model_X.Dp66_CSTATE_fj[0] + -1.0798841515040772E+7 *
      AHV_Model_X.Dp66_CSTATE_fj[1]) + -1.27285250509213E+8 *
      AHV_Model_X.Dp66_CSTATE_fj[2]) + -3.0027175358715254E+8 *
      AHV_Model_X.Dp66_CSTATE_fj[3]) + 1.7969070654437996E+7 *
      AHV_Model_X.Dp66_CSTATE_fj[4]) + 1.1943182502908507E+8 *
      AHV_Model_B.nu_r_c[5]));

    // Product: '<S15>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2, rtb_Sum,
      AHV_Model_B.Minvtau_o);
    for (i = 0; i < 6; i++) {
      // Outport: '<Root>/nu1' incorporates:
      //   Integrator: '<S15>/Integrator'

      nu1[i] = AHV_Model_X.Integrator_CSTATE_e[i];

      // Outport: '<Root>/eta_AHV4'
      eta_AHV4[i] = AHV_Model_B.Integrator1_n[i];

      // Outport: '<Root>/eta_AHV3'
      eta_AHV3[i] = AHV_Model_B.Integrator1_p[i];

      // Outport: '<Root>/eta_AHV2'
      eta_AHV2[i] = AHV_Model_B.Integrator1_b[i];

      // Outport: '<Root>/eta_AHV1'
      eta_AHV1[i] = AHV_Model_B.Integrator1[i];
    }

    // Fcn: '<S24>/T21 ' incorporates:
    //   Fcn: '<S24>/T31 '

    rtb_costheta_1 = std::tan(AHV_Model_B.Integrator1[4]);

    // Product: '<S20>/Product1' incorporates:
    //   Constant: '<S24>/Constant'
    //   Constant: '<S24>/Constant '
    //   Fcn: '<S24>/T21 '
    //   Fcn: '<S24>/T23'
    //   Fcn: '<S24>/T31 '
    //   Fcn: '<S24>/T32'
    //   Fcn: '<S24>/T33'
    //   Reshape: '<S24>/Reshape 9x1->3x3'

    rtb_TmpSignalConversionAtProd_i[0] = 1.0;
    rtb_TmpSignalConversionAtProd_i[3] = rtb_y_dot * rtb_costheta_1;
    rtb_TmpSignalConversionAtProd_i[6] = rtb_costheta_tmp * rtb_costheta_1;
    rtb_TmpSignalConversionAtProd_i[1] = 0.0;
    rtb_TmpSignalConversionAtProd_i[4] = rtb_costheta_tmp;
    rtb_TmpSignalConversionAtProd_i[7] = -rtb_y_dot;
    rtb_TmpSignalConversionAtProd_i[2] = 0.0;
    rtb_TmpSignalConversionAtProd_i[5] = rtb_y_dot / rtb_costheta_tmp_0;
    rtb_TmpSignalConversionAtProd_i[8] = rtb_costheta_tmp / rtb_costheta_tmp_0;
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S15>/Integrator1' incorporates:
      //   Integrator: '<S15>/Integrator'
      //   Product: '<S20>/Product'
      //   Product: '<S20>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegrat_o[i] =
        rtb_TmpSignalConversionAtProduc[i + 6] *
        AHV_Model_X.Integrator_CSTATE_e[2] + (rtb_TmpSignalConversionAtProduc[i
        + 3] * AHV_Model_X.Integrator_CSTATE_e[1] +
        rtb_TmpSignalConversionAtProduc[i] * AHV_Model_X.Integrator_CSTATE_e[0]);
      AHV_Model_B.TmpSignalConversionAtIntegrat_o[i + 3] =
        rtb_TmpSignalConversionAtProd_i[i + 6] *
        AHV_Model_X.Integrator_CSTATE_e[5] + (rtb_TmpSignalConversionAtProd_i[i
        + 3] * AHV_Model_X.Integrator_CSTATE_e[4] +
        rtb_TmpSignalConversionAtProd_i[i] * AHV_Model_X.Integrator_CSTATE_e[3]);
    }

    // Fcn: '<S120>/T21 ' incorporates:
    //   Fcn: '<S120>/T31 '

    rtb_y_dot = std::tan(AHV_Model_B.Integrator1_b[4]);

    // Product: '<S116>/Product1' incorporates:
    //   Constant: '<S120>/Constant'
    //   Constant: '<S120>/Constant '
    //   Fcn: '<S120>/T21 '
    //   Fcn: '<S120>/T23'
    //   Fcn: '<S120>/T31 '
    //   Fcn: '<S120>/T32'
    //   Fcn: '<S120>/T33'
    //   Reshape: '<S120>/Reshape 9x1->3x3'

    rtb_TmpSignalConversionAtProd_i[0] = 1.0;
    rtb_TmpSignalConversionAtProd_i[3] = rtb_TransferFcn * rtb_y_dot;
    rtb_TmpSignalConversionAtProd_i[6] = rtb_costheta_tmp_1 * rtb_y_dot;
    rtb_TmpSignalConversionAtProd_i[1] = 0.0;
    rtb_TmpSignalConversionAtProd_i[4] = rtb_costheta_tmp_1;
    rtb_TmpSignalConversionAtProd_i[7] = -rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[2] = 0.0;
    rtb_TmpSignalConversionAtProd_i[5] = rtb_TransferFcn / rtb_costheta_tmp_2;
    rtb_TmpSignalConversionAtProd_i[8] = rtb_costheta_tmp_1 / rtb_costheta_tmp_2;
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S111>/Integrator1' incorporates:
      //   Integrator: '<S111>/Integrator'
      //   Product: '<S116>/Product'
      //   Product: '<S116>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegra_o3[i] =
        rtb_TmpSignalConversionAtProd_o[i + 6] *
        AHV_Model_X.Integrator_CSTATE_o[2] + (rtb_TmpSignalConversionAtProd_o[i
        + 3] * AHV_Model_X.Integrator_CSTATE_o[1] +
        rtb_TmpSignalConversionAtProd_o[i] * AHV_Model_X.Integrator_CSTATE_o[0]);
      AHV_Model_B.TmpSignalConversionAtIntegra_o3[i + 3] =
        rtb_TmpSignalConversionAtProd_i[i + 6] *
        AHV_Model_X.Integrator_CSTATE_o[5] + (rtb_TmpSignalConversionAtProd_i[i
        + 3] * AHV_Model_X.Integrator_CSTATE_o[4] +
        rtb_TmpSignalConversionAtProd_i[i] * AHV_Model_X.Integrator_CSTATE_o[3]);
    }

    // Fcn: '<S216>/T21 ' incorporates:
    //   Fcn: '<S216>/T31 '

    rtb_TransferFcn = std::tan(AHV_Model_B.Integrator1_p[4]);

    // Product: '<S212>/Product1' incorporates:
    //   Constant: '<S216>/Constant'
    //   Constant: '<S216>/Constant '
    //   Fcn: '<S216>/T21 '
    //   Fcn: '<S216>/T23'
    //   Fcn: '<S216>/T31 '
    //   Fcn: '<S216>/T32'
    //   Fcn: '<S216>/T33'
    //   Reshape: '<S216>/Reshape 9x1->3x3'

    rtb_TmpSignalConversionAtProd_i[0] = 1.0;
    rtb_TmpSignalConversionAtProd_i[3] = rtb_sintheta_tmp_0 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[6] = rtb_costheta_tmp_4 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[1] = 0.0;
    rtb_TmpSignalConversionAtProd_i[4] = rtb_costheta_tmp_4;
    rtb_TmpSignalConversionAtProd_i[7] = -rtb_sintheta_tmp_0;
    rtb_TmpSignalConversionAtProd_i[2] = 0.0;
    rtb_TmpSignalConversionAtProd_i[5] = rtb_sintheta_tmp_0 / rtb_costheta_tmp_5;
    rtb_TmpSignalConversionAtProd_i[8] = rtb_costheta_tmp_4 / rtb_costheta_tmp_5;
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S207>/Integrator1' incorporates:
      //   Integrator: '<S207>/Integrator'
      //   Product: '<S212>/Product'
      //   Product: '<S212>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i] =
        rtb_TmpSignalConversionAtProd_f[i + 6] *
        AHV_Model_X.Integrator_CSTATE_h[2] + (rtb_TmpSignalConversionAtProd_f[i
        + 3] * AHV_Model_X.Integrator_CSTATE_h[1] +
        rtb_TmpSignalConversionAtProd_f[i] * AHV_Model_X.Integrator_CSTATE_h[0]);
      AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i + 3] =
        rtb_TmpSignalConversionAtProd_i[i + 6] *
        AHV_Model_X.Integrator_CSTATE_h[5] + (rtb_TmpSignalConversionAtProd_i[i
        + 3] * AHV_Model_X.Integrator_CSTATE_h[4] +
        rtb_TmpSignalConversionAtProd_i[i] * AHV_Model_X.Integrator_CSTATE_h[3]);
    }

    // Fcn: '<S312>/T21 ' incorporates:
    //   Fcn: '<S312>/T31 '

    rtb_TransferFcn = std::tan(AHV_Model_B.Integrator1_n[4]);

    // Product: '<S308>/Product1' incorporates:
    //   Constant: '<S312>/Constant'
    //   Constant: '<S312>/Constant '
    //   Fcn: '<S312>/T21 '
    //   Fcn: '<S312>/T23'
    //   Fcn: '<S312>/T31 '
    //   Fcn: '<S312>/T32'
    //   Fcn: '<S312>/T33'
    //   Reshape: '<S312>/Reshape 9x1->3x3'

    rtb_TmpSignalConversionAtProd_i[0] = 1.0;
    rtb_TmpSignalConversionAtProd_i[3] = rtb_sintheta_tmp_2 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[6] = rtb_costheta_tmp_7 * rtb_TransferFcn;
    rtb_TmpSignalConversionAtProd_i[1] = 0.0;
    rtb_TmpSignalConversionAtProd_i[4] = rtb_costheta_tmp_7;
    rtb_TmpSignalConversionAtProd_i[7] = -rtb_sintheta_tmp_2;
    rtb_TmpSignalConversionAtProd_i[2] = 0.0;
    rtb_TmpSignalConversionAtProd_i[5] = rtb_sintheta_tmp_2 / rtb_costheta_tmp_8;
    rtb_TmpSignalConversionAtProd_i[8] = rtb_costheta_tmp_7 / rtb_costheta_tmp_8;
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S303>/Integrator1' incorporates:
      //   Integrator: '<S303>/Integrator'
      //   Product: '<S308>/Product'
      //   Product: '<S308>/Product1'

      AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i] =
        rtb_TmpSignalConversionAtProd_j[i + 6] * AHV_Model_X.Integrator_CSTATE[2]
        + (rtb_TmpSignalConversionAtProd_j[i + 3] *
           AHV_Model_X.Integrator_CSTATE[1] + rtb_TmpSignalConversionAtProd_j[i]
           * AHV_Model_X.Integrator_CSTATE[0]);
      AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i + 3] =
        rtb_TmpSignalConversionAtProd_i[i + 6] * AHV_Model_X.Integrator_CSTATE[5]
        + (rtb_TmpSignalConversionAtProd_i[i + 3] *
           AHV_Model_X.Integrator_CSTATE[4] + rtb_TmpSignalConversionAtProd_i[i]
           * AHV_Model_X.Integrator_CSTATE[3]);
    }

    for (i = 0; i < 5; i++) {
      // Sum: '<S34>/Sum'
      tmp_8 = 0.0;

      // Sum: '<S130>/Sum'
      tmp = 0.0;

      // Sum: '<S226>/Sum'
      tmp_0 = 0.0;

      // Sum: '<S322>/Sum'
      tmp_1 = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        // Sum: '<S34>/Sum' incorporates:
        //   Constant: '<S130>/A44'
        //   Constant: '<S226>/A44'
        //   Constant: '<S322>/A44'
        //   Constant: '<S34>/A44'
        //   Integrator: '<S34>/Integrator'
        //   Product: '<S130>/Product'
        //   Product: '<S226>/Product'
        //   Product: '<S322>/Product'
        //   Product: '<S34>/Product'
        //   Sum: '<S130>/Sum'
        //   Sum: '<S226>/Sum'
        //   Sum: '<S322>/Sum'

        tmp_9 = 5 * i_0 + i;
        tmp_8 += AHV_Model_ConstP.pooled59[tmp_9] *
          AHV_Model_X.Integrator_CSTATE_p[i_0];

        // Sum: '<S130>/Sum' incorporates:
        //   Constant: '<S130>/A44'
        //   Integrator: '<S130>/Integrator'
        //   Product: '<S130>/Product'

        tmp += AHV_Model_ConstP.pooled59[tmp_9] *
          AHV_Model_X.Integrator_CSTATE_k[i_0];

        // Sum: '<S226>/Sum' incorporates:
        //   Constant: '<S226>/A44'
        //   Integrator: '<S226>/Integrator'
        //   Product: '<S226>/Product'

        tmp_0 += AHV_Model_ConstP.pooled59[tmp_9] *
          AHV_Model_X.Integrator_CSTATE_c[i_0];

        // Sum: '<S322>/Sum' incorporates:
        //   Constant: '<S322>/A44'
        //   Integrator: '<S322>/Integrator'
        //   Product: '<S322>/Product'

        tmp_1 += AHV_Model_ConstP.pooled59[tmp_9] *
          AHV_Model_X.Integrator_CSTATE_d[i_0];
      }

      // Sum: '<S34>/Sum' incorporates:
      //   Constant: '<S34>/A44'
      //   Constant: '<S34>/B44'
      //   Product: '<S34>/Product'
      //   Product: '<S34>/Product1'

      AHV_Model_B.Sum[i] = AHV_Model_ConstP.pooled60[i] * AHV_Model_B.nu_r_c[3]
        + tmp_8;

      // Sum: '<S130>/Sum' incorporates:
      //   Constant: '<S130>/A44'
      //   Constant: '<S130>/B44'
      //   Product: '<S130>/Product'
      //   Product: '<S130>/Product1'

      AHV_Model_B.Sum_h[i] = AHV_Model_ConstP.pooled60[i] * AHV_Model_B.nu_r_j[3]
        + tmp;

      // Sum: '<S226>/Sum' incorporates:
      //   Constant: '<S226>/A44'
      //   Constant: '<S226>/B44'
      //   Product: '<S226>/Product'
      //   Product: '<S226>/Product1'

      AHV_Model_B.Sum_n[i] = AHV_Model_ConstP.pooled60[i] * AHV_Model_B.nu_r_m[3]
        + tmp_0;

      // Sum: '<S322>/Sum' incorporates:
      //   Constant: '<S322>/A44'
      //   Constant: '<S322>/B44'
      //   Product: '<S322>/Product'
      //   Product: '<S322>/Product1'

      AHV_Model_B.Sum_d[i] = AHV_Model_ConstP.pooled60[i] * AHV_Model_B.nu_r[3]
        + tmp_1;
    }

    // Sum: '<S65>/Sum2' incorporates:
    //   Integrator: '<S65>/Integrator1'
    //   Integrator: '<S65>/Integrator3'
    //   Sum: '<S65>/Sum4'

    rtb_costheta_tmp = AHV_Model_B.Integrator1[0] -
      (AHV_Model_X.Integrator1_CSTATE_a[0] + AHV_Model_X.Integrator3_CSTATE_n[0]);
    rtb_costheta_tmp_0 = AHV_Model_B.Integrator1[1] -
      (AHV_Model_X.Integrator1_CSTATE_a[1] + AHV_Model_X.Integrator3_CSTATE_n[1]);
    rtb_sintheta[2] = AHV_Model_B.Integrator1[5] -
      (AHV_Model_X.Integrator1_CSTATE_a[2] + rtb_rxpi_n);

    // Saturate: '<S74>/x_Saturation'
    if (rtb_sintheta[2] > 1.0E+10) {
      rtb_TransferFcn = 1.0E+10;
    } else if (rtb_sintheta[2] < -1.0E+10) {
      rtb_TransferFcn = -1.0E+10;
    } else {
      rtb_TransferFcn = rtb_sintheta[2];
    }

    // End of Saturate: '<S74>/x_Saturation'

    // Signum: '<S74>/x_Sign'
    if (rtb_TransferFcn < 0.0) {
      rtb_y_dot = -1.0;
    } else if (rtb_TransferFcn > 0.0) {
      rtb_y_dot = 1.0;
    } else if (rtb_TransferFcn == 0.0) {
      rtb_y_dot = 0.0;
    } else {
      rtb_y_dot = (rtNaN);
    }

    // End of Signum: '<S74>/x_Sign'

    // Gain: '<S74>/pi'
    rtb_y_dot *= 3.1415926535897931;

    // Sum: '<S74>/Sum1'
    rtb_TransferFcn += rtb_y_dot;

    // Math: '<S74>/Math Function' incorporates:
    //   Constant: '<S74>/Constant'

    rtb_TransferFcn = rt_remd_snf(rtb_TransferFcn, 6.2831853071795862);

    // Sum: '<S74>/Sum'
    rtb_sintheta_tmp_0 = rtb_TransferFcn - rtb_y_dot;

    // Gain: '<S65>/K4' incorporates:
    //   Sum: '<S74>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_sintheta_tmp_0 +
        (AHV_Model_ConstP.pooled114[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled114[i] * rtb_costheta_tmp);
    }

    // End of Gain: '<S65>/K4'

    // Fcn: '<S71>/Row1'
    rtb_TransferFcn = psi_mean_h * rtb_costheta[0] + rtb_sintheta_tmp *
      rtb_costheta[1];

    // Fcn: '<S71>/Row2'
    rtb_y_dot = -rtb_sintheta_tmp * rtb_costheta[0] + psi_mean_h * rtb_costheta
      [1];

    // Fcn: '<S71>/Row3'
    rtb_costheta_1 = rtb_costheta[2];

    // Gain: '<S65>/Gain6'
    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled115[i + 6] * rtb_v[2] +
        (AHV_Model_ConstP.pooled115[i + 3] * rtb_v[1] +
         AHV_Model_ConstP.pooled115[i] * rtb_v[0]);
    }

    // End of Gain: '<S65>/Gain6'

    // Sum: '<S65>/Sum8' incorporates:
    //   Fcn: '<S72>/Row1'
    //   Fcn: '<S72>/Row2'
    //   Integrator: '<S65>/Integrator6'
    //   Sum: '<S65>/Sum1'

    rtb_Gain1_n[0] = (((psi_mean_h * AHV_Model_X.Integrator6_CSTATE[0] +
                        rtb_sintheta_tmp * AHV_Model_X.Integrator6_CSTATE[1]) +
                       rtb_Sum_j) + rtb_TransferFcn) - rtb_costheta[0];
    rtb_Gain1_n[1] = (((-std::sin(AHV_Model_B.Integrator1[5]) *
                        AHV_Model_X.Integrator6_CSTATE[0] + psi_mean_h *
                        AHV_Model_X.Integrator6_CSTATE[1]) + rtb_Product1_c) +
                      rtb_y_dot) - rtb_costheta[1];
    rtb_Gain1_n[2] = ((AHV_Model_X.Integrator6_CSTATE[2] + rtb_Product2_a) +
                      rtb_costheta_1) - rtb_costheta[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S65>/Gain3'
      AHV_Model_B.M_u[i] = 0.0;
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled116[i] * rtb_Gain1_n[0];
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled116[i + 3] * rtb_Gain1_n[1];
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled116[i + 6] * rtb_Gain1_n[2];

      // Sum: '<S65>/Sum5' incorporates:
      //   Gain: '<S65>/K11'
      //   Integrator: '<S65>/Integrator1'

      AHV_Model_B.psi_WF[i] = ((AHV_Model_ConstP.pooled117[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled117[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled117[i + 6] * rtb_sintheta_tmp_0) +
        AHV_Model_X.Integrator1_CSTATE_a[i];

      // Sum: '<S65>/Sum6' incorporates:
      //   Gain: '<S65>/Gain1'
      //   Gain: '<S65>/Gain2'
      //   Gain: '<S65>/K12'
      //   Integrator: '<S65>/Integrator1'
      //   Integrator: '<S65>/Integrator2'
      //   Sum: '<S65>/Sum2'

      AHV_Model_B.Sum6[i] = ((AHV_Model_ConstP.pooled118[i + 6] *
        rtb_sintheta_tmp_0 + (AHV_Model_ConstP.pooled118[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled118[i] * rtb_costheta_tmp))
        - (AHV_Model_ConstP.pooled120[i + 6] * AHV_Model_X.Integrator2_CSTATE[2]
           + (AHV_Model_ConstP.pooled120[i + 3] *
              AHV_Model_X.Integrator2_CSTATE[1] + AHV_Model_ConstP.pooled120[i] *
              AHV_Model_X.Integrator2_CSTATE[0]))) -
        ((AHV_Model_ConstP.pooled119[i + 3] * AHV_Model_X.Integrator1_CSTATE_a[1]
          + AHV_Model_ConstP.pooled119[i] * AHV_Model_X.Integrator1_CSTATE_a[0])
         + AHV_Model_ConstP.pooled119[i + 6] * AHV_Model_X.Integrator1_CSTATE_a
         [2]);
    }

    // Fcn: '<S70>/Fcn'
    rtb_TransferFcn = psi_mean_h * rtb_v[0] - rtb_sintheta_tmp * rtb_v[1];

    // Fcn: '<S70>/Fcn1'
    rtb_y_dot = rtb_sintheta_tmp * rtb_v[0] + psi_mean_h * rtb_v[1];

    // Fcn: '<S70>/Fcn2'
    rtb_costheta_1 = rtb_v[2];

    // Fcn: '<S168>/Row1'
    tmp_4[0] = AHV_Model_B.Integrator1_b[5];
    tmp_6[0] = AHV_Model_B.Integrator1_b[5];
    for (i = 0; i < 3; i++) {
      // Gain: '<S65>/K3'
      rtb_v[i] = 0.0;
      rtb_v[i] += AHV_Model_ConstP.pooled122[i] * rtb_costheta_tmp;
      rtb_v[i] += AHV_Model_ConstP.pooled122[i + 3] * rtb_costheta_tmp_0;
      rtb_v[i] += AHV_Model_ConstP.pooled122[i + 6] * rtb_sintheta_tmp_0;

      // Sum: '<S65>/Sum7' incorporates:
      //   Gain: '<S65>/inv(T_b)'
      //   Integrator: '<S65>/Integrator6'

      AHV_Model_B.Sum7[i] = rtb_v[i] - ((AHV_Model_ConstP.pooled123[i + 3] *
        AHV_Model_X.Integrator6_CSTATE[1] + AHV_Model_ConstP.pooled123[i] *
        AHV_Model_X.Integrator6_CSTATE[0]) + AHV_Model_ConstP.pooled123[i + 6] *
        AHV_Model_X.Integrator6_CSTATE[2]);

      // Integrator: '<S161>/Integrator6'
      rtb_v[i] = AHV_Model_X.Integrator6_CSTATE_p[i];

      // Fcn: '<S168>/Row1'
      tmp_4[i + 1] = rtb_v[i];
      tmp_5[i + 1] = rtb_v[i];
      tmp_6[i + 1] = rtb_v[i];
      tmp_7[i + 1] = rtb_v[i];

      // Gain: '<S65>/K2' incorporates:
      //   Sum: '<S65>/Sum2'

      rtb_Gain1_n[i] = AHV_Model_ConstP.pooled121[i + 6] * rtb_sintheta_tmp_0 +
        (AHV_Model_ConstP.pooled121[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled121[i] * rtb_costheta_tmp);
    }

    // Sum: '<S65>/Sum3'
    AHV_Model_B.sun_k2[0] = rtb_Gain1_n[0] + rtb_TransferFcn;
    AHV_Model_B.sun_k2[1] = rtb_Gain1_n[1] + rtb_y_dot;
    AHV_Model_B.sun_k2[2] = rtb_Gain1_n[2] + rtb_costheta_1;

    // Sum: '<S161>/Sum2' incorporates:
    //   Integrator: '<S161>/Integrator1'
    //   Integrator: '<S161>/Integrator3'
    //   Sum: '<S161>/Sum4'

    rtb_costheta_tmp = AHV_Model_B.Integrator1_b[0] -
      (AHV_Model_X.Integrator1_CSTATE_h0[0] + AHV_Model_X.Integrator3_CSTATE_l[0]);
    rtb_costheta_tmp_0 = AHV_Model_B.Integrator1_b[1] -
      (AHV_Model_X.Integrator1_CSTATE_h0[1] + AHV_Model_X.Integrator3_CSTATE_l[1]);
    rtb_costheta[2] = AHV_Model_B.Integrator1_b[5] -
      (AHV_Model_X.Integrator1_CSTATE_h0[2] + rtb_rxpi_c);

    // Saturate: '<S170>/x_Saturation'
    if (rtb_costheta[2] > 1.0E+10) {
      rtb_TransferFcn = 1.0E+10;
    } else if (rtb_costheta[2] < -1.0E+10) {
      rtb_TransferFcn = -1.0E+10;
    } else {
      rtb_TransferFcn = rtb_costheta[2];
    }

    // End of Saturate: '<S170>/x_Saturation'

    // Signum: '<S170>/x_Sign'
    if (rtb_TransferFcn < 0.0) {
      rtb_y_dot = -1.0;
    } else if (rtb_TransferFcn > 0.0) {
      rtb_y_dot = 1.0;
    } else if (rtb_TransferFcn == 0.0) {
      rtb_y_dot = 0.0;
    } else {
      rtb_y_dot = (rtNaN);
    }

    // End of Signum: '<S170>/x_Sign'

    // Gain: '<S170>/pi'
    rtb_y_dot *= 3.1415926535897931;

    // Sum: '<S170>/Sum1'
    rtb_TransferFcn += rtb_y_dot;

    // Math: '<S170>/Math Function' incorporates:
    //   Constant: '<S170>/Constant'

    rtb_TransferFcn = rt_remd_snf(rtb_TransferFcn, 6.2831853071795862);

    // Sum: '<S170>/Sum'
    psi_mean_h = rtb_TransferFcn - rtb_y_dot;

    // Gain: '<S161>/K4' incorporates:
    //   Sum: '<S170>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_Gain1_n[i] = AHV_Model_ConstP.pooled114[i + 6] * psi_mean_h +
        (AHV_Model_ConstP.pooled114[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled114[i] * rtb_costheta_tmp);
    }

    // End of Gain: '<S161>/K4'

    // Fcn: '<S167>/Row1'
    rtb_TransferFcn = rtb_costheta_tmp_3 * rtb_Gain1_n[0] + rtb_sintheta_tmp_1 *
      rtb_Gain1_n[1];

    // Fcn: '<S167>/Row2'
    rtb_y_dot = -rtb_sintheta_tmp_1 * rtb_Gain1_n[0] + rtb_costheta_tmp_3 *
      rtb_Gain1_n[1];

    // Fcn: '<S167>/Row3'
    rtb_costheta_1 = rtb_Gain1_n[2];

    // Gain: '<S161>/Gain6'
    for (i = 0; i < 3; i++) {
      rtb_Gain1_n[i] = AHV_Model_ConstP.pooled115[i + 6] * rtb_k_u[2] +
        (AHV_Model_ConstP.pooled115[i + 3] * rtb_k_u[1] +
         AHV_Model_ConstP.pooled115[i] * rtb_k_u[0]);
    }

    // End of Gain: '<S161>/Gain6'

    // Sum: '<S161>/Sum8' incorporates:
    //   Fcn: '<S168>/Row1'
    //   Fcn: '<S168>/Row2'
    //   Fcn: '<S168>/Row3'
    //   Sum: '<S161>/Sum1'

    rtb_sintheta[0] = (((std::cos(tmp_4[0]) * tmp_5[1] + std::sin(tmp_6[0]) *
                         tmp_7[2]) + rtb_Row2_e) + rtb_TransferFcn) -
      rtb_Gain1_n[0];
    rtb_sintheta[1] = (((-std::sin(AHV_Model_B.Integrator1_b[5]) * rtb_v[0] +
                         rtb_costheta_tmp_3 * rtb_v[1]) + rtb_Product1_e3) +
                       rtb_y_dot) - rtb_Gain1_n[1];
    rtb_sintheta[2] = ((rtb_v[2] + rtb_Product2_n) + rtb_costheta_1) -
      rtb_Gain1_n[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S161>/Gain3'
      AHV_Model_B.M_u_o[i] = 0.0;
      AHV_Model_B.M_u_o[i] += AHV_Model_ConstP.pooled116[i] * rtb_sintheta[0];
      AHV_Model_B.M_u_o[i] += AHV_Model_ConstP.pooled116[i + 3] * rtb_sintheta[1];
      AHV_Model_B.M_u_o[i] += AHV_Model_ConstP.pooled116[i + 6] * rtb_sintheta[2];

      // Sum: '<S161>/Sum5' incorporates:
      //   Gain: '<S161>/K11'
      //   Integrator: '<S161>/Integrator1'

      AHV_Model_B.psi_WF_o[i] = ((AHV_Model_ConstP.pooled117[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled117[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled117[i + 6] * psi_mean_h) +
        AHV_Model_X.Integrator1_CSTATE_h0[i];

      // Sum: '<S161>/Sum6' incorporates:
      //   Gain: '<S161>/Gain1'
      //   Gain: '<S161>/Gain2'
      //   Gain: '<S161>/K12'
      //   Integrator: '<S161>/Integrator1'
      //   Integrator: '<S161>/Integrator2'
      //   Sum: '<S161>/Sum2'

      AHV_Model_B.Sum6_d[i] = ((AHV_Model_ConstP.pooled118[i + 6] * psi_mean_h +
        (AHV_Model_ConstP.pooled118[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled118[i] * rtb_costheta_tmp)) -
        (AHV_Model_ConstP.pooled120[i + 6] * AHV_Model_X.Integrator2_CSTATE_j[2]
         + (AHV_Model_ConstP.pooled120[i + 3] *
            AHV_Model_X.Integrator2_CSTATE_j[1] + AHV_Model_ConstP.pooled120[i] *
            AHV_Model_X.Integrator2_CSTATE_j[0]))) -
        ((AHV_Model_ConstP.pooled119[i + 3] * AHV_Model_X.Integrator1_CSTATE_h0
          [1] + AHV_Model_ConstP.pooled119[i] *
          AHV_Model_X.Integrator1_CSTATE_h0[0]) + AHV_Model_ConstP.pooled119[i +
         6] * AHV_Model_X.Integrator1_CSTATE_h0[2]);
    }

    // Gain: '<S161>/K2' incorporates:
    //   Sum: '<S170>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_sintheta[i] = AHV_Model_ConstP.pooled121[i + 6] * psi_mean_h +
        (AHV_Model_ConstP.pooled121[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled121[i] * rtb_costheta_tmp);
    }

    // End of Gain: '<S161>/K2'

    // Sum: '<S161>/Sum3' incorporates:
    //   Fcn: '<S166>/Fcn'
    //   Fcn: '<S166>/Fcn1'
    //   Fcn: '<S166>/Fcn2'

    AHV_Model_B.sun_k2_e[0] = (rtb_costheta_tmp_3 * rtb_k_u[0] -
      rtb_sintheta_tmp_1 * rtb_k_u[1]) + rtb_sintheta[0];
    AHV_Model_B.sun_k2_e[1] = (rtb_sintheta_tmp_1 * rtb_k_u[0] +
      rtb_costheta_tmp_3 * rtb_k_u[1]) + rtb_sintheta[1];
    AHV_Model_B.sun_k2_e[2] = rtb_sintheta[2] + rtb_k_u[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S161>/inv(T_b)'
      rtb_k_u[i] = 0.0;
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i] * rtb_v[0];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 3] * rtb_v[1];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 6] * rtb_v[2];

      // Sum: '<S161>/Sum7' incorporates:
      //   Gain: '<S161>/K3'
      //   Sum: '<S170>/Sum'

      AHV_Model_B.Sum7_o[i] = ((AHV_Model_ConstP.pooled122[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled122[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled122[i + 6] * psi_mean_h) - rtb_k_u[i];
    }

    // Integrator: '<S257>/Integrator6'
    rtb_v[0] = AHV_Model_X.Integrator6_CSTATE_i[0];
    rtb_v[1] = AHV_Model_X.Integrator6_CSTATE_i[1];
    rtb_v[2] = AHV_Model_X.Integrator6_CSTATE_i[2];

    // Integrator: '<S257>/Integrator1'
    rtb_k_u[0] = AHV_Model_X.Integrator1_CSTATE_hl[0];
    rtb_k_u[1] = AHV_Model_X.Integrator1_CSTATE_hl[1];
    rtb_k_u[2] = AHV_Model_X.Integrator1_CSTATE_hl[2];

    // Sum: '<S257>/Sum2' incorporates:
    //   Sum: '<S257>/Sum4'

    rtb_Gain1_n[2] = AHV_Model_B.Integrator1_p[5] - (rtb_k_u[2] + rtb_Switch2);

    // Saturate: '<S266>/x_Saturation'
    if (rtb_Gain1_n[2] > 1.0E+10) {
      rtb_TransferFcn = 1.0E+10;
    } else if (rtb_Gain1_n[2] < -1.0E+10) {
      rtb_TransferFcn = -1.0E+10;
    } else {
      rtb_TransferFcn = rtb_Gain1_n[2];
    }

    // End of Saturate: '<S266>/x_Saturation'

    // Signum: '<S266>/x_Sign'
    if (rtb_TransferFcn < 0.0) {
      rtb_y_dot = -1.0;
    } else if (rtb_TransferFcn > 0.0) {
      rtb_y_dot = 1.0;
    } else if (rtb_TransferFcn == 0.0) {
      rtb_y_dot = 0.0;
    } else {
      rtb_y_dot = (rtNaN);
    }

    // End of Signum: '<S266>/x_Sign'

    // Gain: '<S266>/pi'
    rtb_y_dot *= 3.1415926535897931;

    // Sum: '<S266>/Sum1'
    rtb_TransferFcn += rtb_y_dot;

    // Math: '<S266>/Math Function' incorporates:
    //   Constant: '<S266>/Constant'

    rtb_TransferFcn = rt_remd_snf(rtb_TransferFcn, 6.2831853071795862);

    // Sum: '<S257>/Sum2' incorporates:
    //   Integrator: '<S257>/Integrator3'
    //   Sum: '<S257>/Sum4'

    psi_mean_h = AHV_Model_B.Integrator1_p[0] - (rtb_k_u[0] +
      AHV_Model_X.Integrator3_CSTATE_c[0]);

    // SignalConversion generated from: '<S257>/K11' incorporates:
    //   Integrator: '<S257>/Integrator3'
    //   Sum: '<S257>/Sum2'
    //   Sum: '<S257>/Sum4'
    //   Sum: '<S266>/Sum'

    rtb_sintheta_tmp = AHV_Model_B.Integrator1_p[1] - (rtb_k_u[1] +
      AHV_Model_X.Integrator3_CSTATE_c[1]);
    rtb_sintheta_tmp_1 = rtb_TransferFcn - rtb_y_dot;

    // Gain: '<S257>/K4' incorporates:
    //   SignalConversion generated from: '<S257>/K11'
    //   Sum: '<S257>/Sum2'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_sintheta_tmp_1 +
        (AHV_Model_ConstP.pooled114[i + 3] * rtb_sintheta_tmp +
         AHV_Model_ConstP.pooled114[i] * psi_mean_h);
    }

    // End of Gain: '<S257>/K4'

    // Fcn: '<S263>/Row1'
    rtb_TransferFcn = rtb_costheta_tmp_6 * rtb_costheta[0] +
      rtb_TmpSignalConversionAtProd_9 * rtb_costheta[1];

    // Fcn: '<S263>/Row2'
    rtb_y_dot = -rtb_TmpSignalConversionAtProd_9 * rtb_costheta[0] +
      rtb_costheta_tmp_6 * rtb_costheta[1];

    // Fcn: '<S263>/Row3'
    rtb_costheta_1 = rtb_costheta[2];

    // Gain: '<S257>/Gain6'
    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled115[i + 6] * rtb_Gain2_g[2] +
        (AHV_Model_ConstP.pooled115[i + 3] * rtb_Gain2_g[1] +
         AHV_Model_ConstP.pooled115[i] * rtb_Gain2_g[0]);
    }

    // End of Gain: '<S257>/Gain6'

    // Sum: '<S257>/Sum8' incorporates:
    //   Fcn: '<S264>/Row1'
    //   Fcn: '<S264>/Row2'
    //   Fcn: '<S264>/Row3'
    //   Sum: '<S257>/Sum1'

    rtb_sintheta[0] = (((rtb_costheta_tmp_6 * rtb_v[0] +
                         rtb_TmpSignalConversionAtProd_9 * rtb_v[1]) +
                        rtb_Switch3) + rtb_TransferFcn) - rtb_costheta[0];
    rtb_sintheta[1] = (((-std::sin(AHV_Model_B.Integrator1_p[5]) * rtb_v[0] +
                         rtb_costheta_tmp_6 * rtb_v[1]) + rtb_Product1_lf) +
                       rtb_y_dot) - rtb_costheta[1];
    rtb_sintheta[2] = ((rtb_v[2] + SwayFailure) + rtb_costheta_1) -
      rtb_costheta[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S257>/Gain3'
      AHV_Model_B.M_u_i[i] = 0.0;
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled116[i] * rtb_sintheta[0];
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled116[i + 3] * rtb_sintheta[1];
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled116[i + 6] * rtb_sintheta[2];

      // Sum: '<S257>/Sum5' incorporates:
      //   Gain: '<S257>/K11'
      //   SignalConversion generated from: '<S257>/K11'

      AHV_Model_B.psi_WF_m[i] = ((AHV_Model_ConstP.pooled117[i + 3] *
        rtb_sintheta_tmp + AHV_Model_ConstP.pooled117[i] * psi_mean_h) +
        AHV_Model_ConstP.pooled117[i + 6] * rtb_sintheta_tmp_1) + rtb_k_u[i];

      // Gain: '<S257>/Gain1'
      rtb_costheta[i] = AHV_Model_ConstP.pooled119[i + 6] * rtb_k_u[2] +
        (AHV_Model_ConstP.pooled119[i + 3] * rtb_k_u[1] +
         AHV_Model_ConstP.pooled119[i] * rtb_k_u[0]);
    }

    for (i = 0; i < 3; i++) {
      // Integrator: '<S257>/Integrator2'
      rtb_k_u[i] = AHV_Model_X.Integrator2_CSTATE_f[i];

      // Gain: '<S257>/K12' incorporates:
      //   SignalConversion generated from: '<S257>/K11'
      //   Sum: '<S257>/Sum2'

      rtb_sintheta[i] = AHV_Model_ConstP.pooled118[i + 6] * rtb_sintheta_tmp_1 +
        (AHV_Model_ConstP.pooled118[i + 3] * rtb_sintheta_tmp +
         AHV_Model_ConstP.pooled118[i] * psi_mean_h);
    }

    for (i = 0; i < 3; i++) {
      // Sum: '<S257>/Sum6' incorporates:
      //   Gain: '<S257>/Gain2'

      AHV_Model_B.Sum6_j[i] = (rtb_sintheta[i] - (AHV_Model_ConstP.pooled120[i +
        6] * rtb_k_u[2] + (AHV_Model_ConstP.pooled120[i + 3] * rtb_k_u[1] +
                           AHV_Model_ConstP.pooled120[i] * rtb_k_u[0]))) -
        rtb_costheta[i];
    }

    for (i = 0; i < 3; i++) {
      // Gain: '<S257>/inv(T_b)'
      rtb_k_u[i] = 0.0;
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i] * rtb_v[0];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 3] * rtb_v[1];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 6] * rtb_v[2];

      // Sum: '<S257>/Sum7' incorporates:
      //   Gain: '<S257>/K3'
      //   SignalConversion generated from: '<S257>/K11'

      AHV_Model_B.Sum7_n[i] = ((AHV_Model_ConstP.pooled122[i + 3] *
        rtb_sintheta_tmp + AHV_Model_ConstP.pooled122[i] * psi_mean_h) +
        AHV_Model_ConstP.pooled122[i + 6] * rtb_sintheta_tmp_1) - rtb_k_u[i];

      // Gain: '<S257>/K2' incorporates:
      //   SignalConversion generated from: '<S257>/K11'
      //   Sum: '<S257>/Sum2'

      rtb_sintheta[i] = AHV_Model_ConstP.pooled121[i + 6] * rtb_sintheta_tmp_1 +
        (AHV_Model_ConstP.pooled121[i + 3] * rtb_sintheta_tmp +
         AHV_Model_ConstP.pooled121[i] * psi_mean_h);
    }

    // Sum: '<S257>/Sum3' incorporates:
    //   Fcn: '<S262>/Fcn'
    //   Fcn: '<S262>/Fcn1'
    //   Fcn: '<S262>/Fcn2'

    AHV_Model_B.sun_k2_e3[0] = (rtb_costheta_tmp_6 * rtb_Gain2_g[0] -
      rtb_TmpSignalConversionAtProd_9 * rtb_Gain2_g[1]) + rtb_sintheta[0];
    AHV_Model_B.sun_k2_e3[1] = (rtb_TmpSignalConversionAtProd_9 * rtb_Gain2_g[0]
      + rtb_costheta_tmp_6 * rtb_Gain2_g[1]) + rtb_sintheta[1];
    AHV_Model_B.sun_k2_e3[2] = rtb_sintheta[2] + rtb_Gain2_g[2];

    // Integrator: '<S353>/Integrator6'
    rtb_v[0] = AHV_Model_X.Integrator6_CSTATE_ij[0];
    rtb_v[1] = AHV_Model_X.Integrator6_CSTATE_ij[1];
    rtb_v[2] = AHV_Model_X.Integrator6_CSTATE_ij[2];

    // Integrator: '<S353>/Integrator1'
    rtb_k_u[0] = AHV_Model_X.Integrator1_CSTATE_nd[0];
    rtb_k_u[1] = AHV_Model_X.Integrator1_CSTATE_nd[1];
    rtb_k_u[2] = AHV_Model_X.Integrator1_CSTATE_nd[2];

    // Sum: '<S353>/Sum2' incorporates:
    //   Integrator: '<S353>/Integrator3'
    //   Sum: '<S353>/Sum4'

    rtb_costheta_tmp = AHV_Model_B.Integrator1_n[0] - (rtb_k_u[0] +
      AHV_Model_X.Integrator3_CSTATE[0]);
    rtb_costheta_tmp_0 = AHV_Model_B.Integrator1_n[1] - (rtb_k_u[1] +
      AHV_Model_X.Integrator3_CSTATE[1]);
    rtb_Gain1_n[2] = AHV_Model_B.Integrator1_n[5] - (rtb_k_u[2] + rtb_rxpi_h);

    // Saturate: '<S362>/x_Saturation'
    if (rtb_Gain1_n[2] > 1.0E+10) {
      rtb_TransferFcn = 1.0E+10;
    } else if (rtb_Gain1_n[2] < -1.0E+10) {
      rtb_TransferFcn = -1.0E+10;
    } else {
      rtb_TransferFcn = rtb_Gain1_n[2];
    }

    // End of Saturate: '<S362>/x_Saturation'

    // Signum: '<S362>/x_Sign'
    if (rtb_TransferFcn < 0.0) {
      rtb_y_dot = -1.0;
    } else if (rtb_TransferFcn > 0.0) {
      rtb_y_dot = 1.0;
    } else if (rtb_TransferFcn == 0.0) {
      rtb_y_dot = 0.0;
    } else {
      rtb_y_dot = (rtNaN);
    }

    // End of Signum: '<S362>/x_Sign'

    // Gain: '<S362>/pi'
    rtb_y_dot *= 3.1415926535897931;

    // Sum: '<S362>/Sum1'
    rtb_TransferFcn += rtb_y_dot;

    // Math: '<S362>/Math Function' incorporates:
    //   Constant: '<S362>/Constant'

    rtb_TransferFcn = rt_remd_snf(rtb_TransferFcn, 6.2831853071795862);

    // Sum: '<S362>/Sum'
    psi_mean_h = rtb_TransferFcn - rtb_y_dot;

    // Gain: '<S353>/K4' incorporates:
    //   Sum: '<S362>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled114[i + 6] * psi_mean_h +
        (AHV_Model_ConstP.pooled114[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled114[i] * rtb_costheta_tmp);
    }

    // End of Gain: '<S353>/K4'

    // Fcn: '<S359>/Row1'
    rtb_TransferFcn = rtb_costheta_tmp_9 * rtb_costheta[0] +
      rtb_TmpSignalConversionAtProd_0 * rtb_costheta[1];

    // Fcn: '<S359>/Row2'
    rtb_y_dot = -rtb_TmpSignalConversionAtProd_0 * rtb_costheta[0] +
      rtb_costheta_tmp_9 * rtb_costheta[1];

    // Fcn: '<S359>/Row3'
    rtb_costheta_1 = rtb_costheta[2];

    // Gain: '<S353>/Gain6'
    for (i = 0; i < 3; i++) {
      rtb_costheta[i] = AHV_Model_ConstP.pooled115[i + 6] * rtb_nu[2] +
        (AHV_Model_ConstP.pooled115[i + 3] * rtb_nu[1] +
         AHV_Model_ConstP.pooled115[i] * rtb_nu[0]);
    }

    // End of Gain: '<S353>/Gain6'

    // Sum: '<S353>/Sum8' incorporates:
    //   Fcn: '<S360>/Row1'
    //   Fcn: '<S360>/Row2'
    //   Fcn: '<S360>/Row3'
    //   Sum: '<S353>/Sum1'

    rtb_sintheta[0] = (((rtb_costheta_tmp_9 * rtb_v[0] +
                         rtb_TmpSignalConversionAtProd_0 * rtb_v[1]) +
                        rtb_Product_g) + rtb_TransferFcn) - rtb_costheta[0];
    rtb_sintheta[1] = (((-std::sin(AHV_Model_B.Integrator1_n[5]) * rtb_v[0] +
                         rtb_costheta_tmp_9 * rtb_v[1]) + rtb_Product1_lq) +
                       rtb_y_dot) - rtb_costheta[1];
    rtb_sintheta[2] = ((rtb_v[2] + rtb_Product2_b) + rtb_costheta_1) -
      rtb_costheta[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S353>/Gain3'
      AHV_Model_B.M_u_n[i] = 0.0;
      AHV_Model_B.M_u_n[i] += AHV_Model_ConstP.pooled116[i] * rtb_sintheta[0];
      AHV_Model_B.M_u_n[i] += AHV_Model_ConstP.pooled116[i + 3] * rtb_sintheta[1];
      AHV_Model_B.M_u_n[i] += AHV_Model_ConstP.pooled116[i + 6] * rtb_sintheta[2];

      // Sum: '<S353>/Sum5' incorporates:
      //   Gain: '<S353>/K11'
      //   Sum: '<S362>/Sum'

      AHV_Model_B.psi_WF_i[i] = ((AHV_Model_ConstP.pooled117[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled117[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled117[i + 6] * psi_mean_h) + rtb_k_u[i];

      // Gain: '<S353>/Gain1'
      rtb_costheta[i] = AHV_Model_ConstP.pooled119[i + 6] * rtb_k_u[2] +
        (AHV_Model_ConstP.pooled119[i + 3] * rtb_k_u[1] +
         AHV_Model_ConstP.pooled119[i] * rtb_k_u[0]);
    }

    // Integrator: '<S353>/Integrator2'
    rtb_k_u[0] = AHV_Model_X.Integrator2_CSTATE_p[0];
    rtb_k_u[1] = AHV_Model_X.Integrator2_CSTATE_p[1];
    rtb_k_u[2] = AHV_Model_X.Integrator2_CSTATE_p[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S353>/Gain2'
      rtb_Gain2_g[i] = 0.0;
      rtb_Gain2_g[i] += AHV_Model_ConstP.pooled120[i] * rtb_k_u[0];
      rtb_Gain2_g[i] += AHV_Model_ConstP.pooled120[i + 3] * rtb_k_u[1];
      rtb_Gain2_g[i] += AHV_Model_ConstP.pooled120[i + 6] * rtb_k_u[2];

      // Sum: '<S353>/Sum6' incorporates:
      //   Gain: '<S353>/K12'

      AHV_Model_B.Sum6_i[i] = (((AHV_Model_ConstP.pooled118[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled118[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled118[i + 6] * psi_mean_h) - rtb_Gain2_g[i]) -
        rtb_costheta[i];

      // Gain: '<S353>/K2' incorporates:
      //   Sum: '<S353>/Sum2'

      rtb_sintheta[i] = AHV_Model_ConstP.pooled121[i + 6] * psi_mean_h +
        (AHV_Model_ConstP.pooled121[i + 3] * rtb_costheta_tmp_0 +
         AHV_Model_ConstP.pooled121[i] * rtb_costheta_tmp);
    }

    // Sum: '<S353>/Sum3' incorporates:
    //   Fcn: '<S358>/Fcn'
    //   Fcn: '<S358>/Fcn1'
    //   Fcn: '<S358>/Fcn2'

    AHV_Model_B.sun_k2_a[0] = (rtb_costheta_tmp_9 * rtb_nu[0] -
      rtb_TmpSignalConversionAtProd_0 * rtb_nu[1]) + rtb_sintheta[0];
    AHV_Model_B.sun_k2_a[1] = (rtb_TmpSignalConversionAtProd_0 * rtb_nu[0] +
      rtb_costheta_tmp_9 * rtb_nu[1]) + rtb_sintheta[1];
    AHV_Model_B.sun_k2_a[2] = rtb_sintheta[2] + rtb_nu[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S353>/inv(T_b)'
      rtb_k_u[i] = 0.0;
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i] * rtb_v[0];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 3] * rtb_v[1];
      rtb_k_u[i] += AHV_Model_ConstP.pooled123[i + 6] * rtb_v[2];

      // Sum: '<S353>/Sum7' incorporates:
      //   Gain: '<S353>/K3'
      //   Sum: '<S362>/Sum'

      AHV_Model_B.Sum7_n4[i] = ((AHV_Model_ConstP.pooled122[i + 3] *
        rtb_costheta_tmp_0 + AHV_Model_ConstP.pooled122[i] * rtb_costheta_tmp) +
        AHV_Model_ConstP.pooled122[i + 6] * psi_mean_h) - rtb_k_u[i];
    }

    // Fcn: '<S18>/Fcn' incorporates:
    //   Integrator: '<S15>/Integrator'

    tmp_8 = AHV_Model_X.Integrator_CSTATE_e[0] *
      AHV_Model_X.Integrator_CSTATE_e[0] + AHV_Model_X.Integrator_CSTATE_e[1] *
      AHV_Model_X.Integrator_CSTATE_e[1];
    if (tmp_8 < 0.0) {
      AHV_Model_B.Fcn = -std::sqrt(-tmp_8);
    } else {
      AHV_Model_B.Fcn = std::sqrt(tmp_8);
    }

    // End of Fcn: '<S18>/Fcn'

    // Fcn: '<S114>/Fcn' incorporates:
    //   Integrator: '<S111>/Integrator'

    tmp_8 = AHV_Model_X.Integrator_CSTATE_o[0] *
      AHV_Model_X.Integrator_CSTATE_o[0] + AHV_Model_X.Integrator_CSTATE_o[1] *
      AHV_Model_X.Integrator_CSTATE_o[1];
    if (tmp_8 < 0.0) {
      AHV_Model_B.Fcn_j = -std::sqrt(-tmp_8);
    } else {
      AHV_Model_B.Fcn_j = std::sqrt(tmp_8);
    }

    // End of Fcn: '<S114>/Fcn'

    // Fcn: '<S210>/Fcn' incorporates:
    //   Integrator: '<S207>/Integrator'

    tmp_8 = AHV_Model_X.Integrator_CSTATE_h[0] *
      AHV_Model_X.Integrator_CSTATE_h[0] + AHV_Model_X.Integrator_CSTATE_h[1] *
      AHV_Model_X.Integrator_CSTATE_h[1];
    if (tmp_8 < 0.0) {
      AHV_Model_B.Fcn_b = -std::sqrt(-tmp_8);
    } else {
      AHV_Model_B.Fcn_b = std::sqrt(tmp_8);
    }

    // End of Fcn: '<S210>/Fcn'

    // Fcn: '<S306>/Fcn' incorporates:
    //   Integrator: '<S303>/Integrator'

    tmp_8 = AHV_Model_X.Integrator_CSTATE[0] * AHV_Model_X.Integrator_CSTATE[0]
      + AHV_Model_X.Integrator_CSTATE[1] * AHV_Model_X.Integrator_CSTATE[1];
    if (tmp_8 < 0.0) {
      AHV_Model_B.Fcn_e = -std::sqrt(-tmp_8);
    } else {
      AHV_Model_B.Fcn_e = std::sqrt(tmp_8);
    }

    // End of Fcn: '<S306>/Fcn'

    // Outport: '<Root>/AHV_speed1' incorporates:
    //   Gain: '<S18>/Gain'
    //   Gain: '<S18>/Gain1'
    //   Rounding: '<S18>/Rounding Function'
    //   TransferFcn: '<S18>/Transfer Fcn'

    AHV_speed1 = std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE * 10.0) * 0.1;

    // Outport: '<Root>/AHV_speed2' incorporates:
    //   Gain: '<S114>/Gain'
    //   Gain: '<S114>/Gain1'
    //   Rounding: '<S114>/Rounding Function'
    //   TransferFcn: '<S114>/Transfer Fcn'

    AHV_speed2 = std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_c * 10.0) * 0.1;

    // Outport: '<Root>/AHV_speed3' incorporates:
    //   Gain: '<S210>/Gain'
    //   Gain: '<S210>/Gain1'
    //   Rounding: '<S210>/Rounding Function'
    //   TransferFcn: '<S210>/Transfer Fcn'

    AHV_speed3 = std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_d * 10.0) * 0.1;

    // Outport: '<Root>/AHV_speed4' incorporates:
    //   Gain: '<S306>/Gain'
    //   Gain: '<S306>/Gain1'
    //   Rounding: '<S306>/Rounding Function'
    //   TransferFcn: '<S306>/Transfer Fcn'

    AHV_speed4 = std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_p * 10.0) * 0.1;
  }

  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    real_T tmp;
    real_T tmp_0;
    real_T tmp_1;

    // Update for Integrator: '<S15>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK = 0;

    // Update for Integrator: '<S111>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_o = 0;

    // Update for Integrator: '<S207>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_b = 0;

    // Update for Integrator: '<S303>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_c = 0;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for UnitDelay: '<S9>/Unit Delay'
      AHV_Model_DW.UnitDelay_DSTATE[0] = AHV_Model_B.regulation_error_i[0];
      AHV_Model_DW.UnitDelay_DSTATE[1] = AHV_Model_B.regulation_error_i[1];
      AHV_Model_DW.UnitDelay_DSTATE[2] = AHV_Model_B.rxpi_e;
    }

    // Update for If: '<S8>/If'
    if (AHV_Model_DW.If_ActiveSubsystem == 0) {
      // Update for IfAction SubSystem: '<S8>/If Action Subsystem' incorporates:
      //   ActionPort: '<S91>/Action Port'

      // Update for RateLimiter: '<S94>/Rate Limiter1' incorporates:
      //   RateLimiter: '<S94>/Rate Limiter2'

      AHV_Model_DW.PrevY_a = AHV_Model_B.RateLimiter1;
      AHV_Model_DW.LastMajorTime_e = (&AHV_Model_M)->Timing.t[0];

      // Update for RateLimiter: '<S94>/Rate Limiter2'
      AHV_Model_DW.PrevY_ic = AHV_Model_B.RateLimiter2;
      AHV_Model_DW.LastMajorTime_h = AHV_Model_DW.LastMajorTime_e;

      // End of Update for SubSystem: '<S8>/If Action Subsystem'
    }

    // End of Update for If: '<S8>/If'

    // Update for RateLimiter: '<S66>/Rate Limiter' incorporates:
    //   RateLimiter: '<S162>/Rate Limiter'
    //   RateLimiter: '<S258>/Rate Limiter'
    //   RateLimiter: '<S354>/Rate Limiter'

    AHV_Model_DW.PrevY = AHV_Model_B.RateLimiter;
    AHV_Model_DW.LastMajorTime = (&AHV_Model_M)->Timing.t[0];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for UnitDelay: '<S105>/Unit Delay'
      AHV_Model_DW.UnitDelay_DSTATE_h[0] = AHV_Model_B.regulation_error_k[0];
      AHV_Model_DW.UnitDelay_DSTATE_h[1] = AHV_Model_B.regulation_error_k[1];
      AHV_Model_DW.UnitDelay_DSTATE_h[2] = AHV_Model_B.rxpi_dyz;
    }

    // Update for If: '<S104>/If'
    if (AHV_Model_DW.If_ActiveSubsystem_l == 0) {
      // Update for IfAction SubSystem: '<S104>/If Action Subsystem' incorporates:
      //   ActionPort: '<S187>/Action Port'

      AHV_Mo_IfActionSubsystem_Update(&AHV_Model_B.IfActionSubsystem_c,
        &AHV_Model_DW.IfActionSubsystem_c);

      // End of Update for SubSystem: '<S104>/If Action Subsystem'
    }

    // End of Update for If: '<S104>/If'

    // Update for RateLimiter: '<S162>/Rate Limiter'
    AHV_Model_DW.PrevY_p = AHV_Model_B.RateLimiter_b;
    AHV_Model_DW.LastMajorTime_g = AHV_Model_DW.LastMajorTime;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for UnitDelay: '<S201>/Unit Delay'
      AHV_Model_DW.UnitDelay_DSTATE_c[0] = AHV_Model_B.regulation_error_d[0];
      AHV_Model_DW.UnitDelay_DSTATE_c[1] = AHV_Model_B.regulation_error_d[1];
      AHV_Model_DW.UnitDelay_DSTATE_c[2] = AHV_Model_B.rxpi_k;
    }

    // Update for If: '<S200>/If'
    if (AHV_Model_DW.If_ActiveSubsystem_a == 0) {
      // Update for IfAction SubSystem: '<S200>/If Action Subsystem' incorporates:
      //   ActionPort: '<S283>/Action Port'

      AHV_Mo_IfActionSubsystem_Update(&AHV_Model_B.IfActionSubsystem_e,
        &AHV_Model_DW.IfActionSubsystem_e);

      // End of Update for SubSystem: '<S200>/If Action Subsystem'
    }

    // End of Update for If: '<S200>/If'

    // Update for RateLimiter: '<S258>/Rate Limiter'
    AHV_Model_DW.PrevY_i = AHV_Model_B.RateLimiter_c;
    AHV_Model_DW.LastMajorTime_a = AHV_Model_DW.LastMajorTime;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for UnitDelay: '<S297>/Unit Delay'
      AHV_Model_DW.UnitDelay_DSTATE_a[0] = AHV_Model_B.regulation_error[0];
      AHV_Model_DW.UnitDelay_DSTATE_a[1] = AHV_Model_B.regulation_error[1];
      AHV_Model_DW.UnitDelay_DSTATE_a[2] = AHV_Model_B.rxpi_j5;
    }

    // Update for If: '<S296>/If'
    if (AHV_Model_DW.If_ActiveSubsystem_k == 0) {
      // Update for IfAction SubSystem: '<S296>/If Action Subsystem' incorporates:
      //   ActionPort: '<S379>/Action Port'

      AHV_Mo_IfActionSubsystem_Update(&AHV_Model_B.IfActionSubsystem_eo,
        &AHV_Model_DW.IfActionSubsystem_eo);

      // End of Update for SubSystem: '<S296>/If Action Subsystem'
    }

    // End of Update for If: '<S296>/If'

    // Update for Integrator: '<S353>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK = 0;

    // Update for RateLimiter: '<S354>/Rate Limiter'
    AHV_Model_DW.PrevY_o = AHV_Model_B.RateLimiter_m;
    AHV_Model_DW.LastMajorTime_c = AHV_Model_DW.LastMajorTime;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for DiscreteIntegrator: '<S355>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE[2];
      AHV_Model_DW.Integrator1_DSTATE[0] = 0.05 * AHV_Model_B.Merge + tmp;
      AHV_Model_DW.Integrator1_DSTATE[1] = 0.05 * AHV_Model_B.Merge_e + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE[2] = 0.05 * AHV_Model_B.Row3 + tmp_1;

      // Update for DiscreteIntegrator: '<S259>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE_d[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE_d[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE_d[2];
      AHV_Model_DW.Integrator1_DSTATE_d[0] = 0.05 * AHV_Model_B.Merge_o + tmp;
      AHV_Model_DW.Integrator1_DSTATE_d[1] = 0.05 * AHV_Model_B.Merge_p + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE_d[2] = 0.05 * AHV_Model_B.Row3_j + tmp_1;
    }

    // Update for Integrator: '<S257>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_i = 0;

    // Update for Integrator: '<S161>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_c = 0;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for DiscreteIntegrator: '<S163>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE_p[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE_p[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE_p[2];
      AHV_Model_DW.Integrator1_DSTATE_p[0] = 0.05 * AHV_Model_B.Merge_j + tmp;
      AHV_Model_DW.Integrator1_DSTATE_p[1] = 0.05 * AHV_Model_B.Merge_i + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE_p[2] = 0.05 * AHV_Model_B.Row3_h + tmp_1;

      // Update for DiscreteIntegrator: '<S67>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE_a[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE_a[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE_a[2];
      AHV_Model_DW.Integrator1_DSTATE_a[0] = 0.05 * AHV_Model_B.Merge_k + tmp;
      AHV_Model_DW.Integrator1_DSTATE_a[1] = 0.05 * AHV_Model_B.Merge_b + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE_a[2] = 0.05 * AHV_Model_B.Row3_l + tmp_1;
    }

    // Update for Integrator: '<S65>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_l = 0;
  }                                    // end MajorTimeStep

  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    rt_ertODEUpdateContinuousStates(&(&AHV_Model_M)->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++(&AHV_Model_M)->Timing.clockTick0;
    (&AHV_Model_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&AHV_Model_M)
      ->solverInfo);

    {
      // Update absolute timer for sample time: [0.05s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.05, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      (&AHV_Model_M)->Timing.clockTick1++;
    }
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void AH_Model_v1ModelClass::AHV_Model_derivatives()
{
  int32_T i;
  XDot_AHV_Model_T *_rtXdot;
  _rtXdot = ((XDot_AHV_Model_T *) (&AHV_Model_M)->derivs);
  for (i = 0; i < 6; i++) {
    // Derivatives for Integrator: '<S15>/Integrator1'
    _rtXdot->Integrator1_CSTATE[i] =
      AHV_Model_B.TmpSignalConversionAtIntegrat_o[i];

    // Derivatives for Integrator: '<S111>/Integrator1'
    _rtXdot->Integrator1_CSTATE_n[i] =
      AHV_Model_B.TmpSignalConversionAtIntegra_o3[i];

    // Derivatives for Integrator: '<S207>/Integrator1'
    _rtXdot->Integrator1_CSTATE_h[i] =
      AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i];

    // Derivatives for Integrator: '<S303>/Integrator1'
    _rtXdot->Integrator1_CSTATE_hn[i] =
      AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i];

    // Derivatives for Integrator: '<S303>/Integrator'
    _rtXdot->Integrator_CSTATE[i] = AHV_Model_B.Minvtau[i];
  }

  // Derivatives for Integrator: '<S353>/Integrator3'
  _rtXdot->Integrator3_CSTATE[0] = AHV_Model_B.sun_k2_a[0];

  // Derivatives for Integrator: '<S353>/Integrator4'
  _rtXdot->Integrator4_CSTATE[0] = AHV_Model_B.M_u_n[0];

  // Derivatives for Integrator: '<S353>/Integrator3'
  _rtXdot->Integrator3_CSTATE[1] = AHV_Model_B.sun_k2_a[1];

  // Derivatives for Integrator: '<S353>/Integrator4'
  _rtXdot->Integrator4_CSTATE[1] = AHV_Model_B.M_u_n[1];

  // Derivatives for Integrator: '<S353>/Integrator3'
  _rtXdot->Integrator3_CSTATE[2] = AHV_Model_B.sun_k2_a[2];

  // Derivatives for Integrator: '<S353>/Integrator4'
  _rtXdot->Integrator4_CSTATE[2] = AHV_Model_B.M_u_n[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S310>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE[i] = 0.0;

    // Derivatives for Integrator: '<S322>/Integrator'
    _rtXdot->Integrator_CSTATE_d[i] = AHV_Model_B.Sum_d[i];

    // Derivatives for StateSpace: '<S310>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S310>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S310>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[0] += -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[0] += -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[0] += -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[0] += -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[1] += 1.2054566685310679 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[1] += -0.00075311499999405624 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[1] += -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[1] += -0.023103995048734959 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[1] += -0.0023400735275004216 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[2] += -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[2] += 0.010562919761939651 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[2] += -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[2] += -0.089746811088073059 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[3] += 1.0500590385158888 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[3] += -0.023103995048766656 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[3] += 1.4618515443990421 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[3] += -4.560672872596915 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[3] += -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[4] += -0.095423424399719625 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[4] += 0.0023400735275035814 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[4] += -0.089746811088075445 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[4] += -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[0] += -3.3881001186476949 * AHV_Model_B.nu_r[0];
  _rtXdot->Dp11_CSTATE[1] += 0.07669611088588818 * AHV_Model_B.nu_r[0];
  _rtXdot->Dp11_CSTATE[2] += -0.46181267830731043 * AHV_Model_B.nu_r[0];
  _rtXdot->Dp11_CSTATE[3] += 1.226085627712131 * AHV_Model_B.nu_r[0];
  _rtXdot->Dp11_CSTATE[4] += -0.11754627904442222 * AHV_Model_B.nu_r[0];

  // Derivatives for StateSpace: '<S310>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE[0] += -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[0] += -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[0] += -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[0] += -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[0] += -0.095423424399736376 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[1] += 1.2054566685310135 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[1] += -0.00075311499999501316 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[1] += -0.010562919761942575 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[1] += -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[1] += -0.002340073527504722 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[2] += -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[2] += 0.010562919761928094 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[2] += -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[2] += -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[2] += -0.089746811088091308 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[3] += 1.0500590385158797 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[3] += -0.023103995048724405 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[3] += 1.4618515443989983 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[3] += -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[3] += -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[4] += -0.095423424399735959 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[4] += 0.00234007352750089 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[4] += -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[4] += -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[0] += -3.3881001186476967 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp22_CSTATE[1] += 0.076696110885761712 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp22_CSTATE[2] += -0.46181267830733025 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp22_CSTATE[3] += 1.2260856277121315 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp22_CSTATE[4] += -0.11754627904444465 * AHV_Model_B.nu_r[1];

  // Derivatives for StateSpace: '<S310>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE[0] += -0.0069962390232986621 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[0] += -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[0] += -0.074692039015985229 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[0] += -0.0042999629253483171 * AHV_Model_X.Dp24_CSTATE[4];
  _rtXdot->Dp24_CSTATE[1] += -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[1] += -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[1] += 0.95129605479096324 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[1] += 1.3710797263241543 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[1] += 0.067632454414034648 * AHV_Model_X.Dp24_CSTATE[4];
  _rtXdot->Dp24_CSTATE[2] += 0.078929006579498431 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[2] += 0.9512960547909608 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[2] += -1.065705011247746 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[2] += -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[2] += -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE[4];
  _rtXdot->Dp24_CSTATE[3] += -0.074692039015982259 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[3] += -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[3] += 3.3410455288652465 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[3] += -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[3] += -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE[4];
  _rtXdot->Dp24_CSTATE[4] += 0.0042999629253487083 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[4] += 0.067632454414049734 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[4] += -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[4] += -0.020742386481004484 * AHV_Model_X.Dp24_CSTATE[4];
  _rtXdot->Dp24_CSTATE[0] += -0.23872786278308805 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp24_CSTATE[1] += -3.2763464234276056 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp24_CSTATE[2] += 1.1437118751635387 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp24_CSTATE[3] += -1.3244904674438165 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp24_CSTATE[4] += 0.071163531145327891 * AHV_Model_B.nu_r[3];

  // Derivatives for StateSpace: '<S310>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE[0] += -0.0017008072822443973 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[0] += 1.2361995937266759 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[0] += 0.012808448456780379 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[0] += 0.031758809083141569 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[0] += -0.0047828760333047419 * AHV_Model_X.Dp26_CSTATE[4];
  _rtXdot->Dp26_CSTATE[1] += -1.236199593726705 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[1] += -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[1] += -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[1] += -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[1] += 0.12995794495601007 * AHV_Model_X.Dp26_CSTATE[4];
  _rtXdot->Dp26_CSTATE[2] += -0.012808448456765344 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[2] += -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[2] += -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[2] += -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[2] += 0.1356516649888091 * AHV_Model_X.Dp26_CSTATE[4];
  _rtXdot->Dp26_CSTATE[3] += 0.031758809083164766 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[3] += 0.96895750705376427 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[3] += 1.8113135858263727 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[3] += -3.762337151306161 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[3] += 1.5320613622870083 * AHV_Model_X.Dp26_CSTATE[4];
  _rtXdot->Dp26_CSTATE[4] += 0.0047828760333030167 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[4] += 0.12995794495600341 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[4] += 0.13565166498879488 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[4] += -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[4] += -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE[4];
  _rtXdot->Dp26_CSTATE[0] += 0.13167647316600378 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp26_CSTATE[1] += 3.4125284827245985 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp26_CSTATE[2] += 0.44399732492152644 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp26_CSTATE[3] += -1.2820687298457427 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp26_CSTATE[4] += -0.18276601298221137 * AHV_Model_B.nu_r[5];

  // Derivatives for StateSpace: '<S310>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE[0] += -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[0] += -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[0] += -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[0] += 0.033960047754095488 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[1] += -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[1] += -1.538297528844471E-6 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[1] += 0.00010972576839633958 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[1] += 0.00050341083172827721 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[1] += -3.8193652830025795E-5 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[2] += 0.11835430912143848 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[2] += 0.0001097257683747581 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[2] += -0.033603468126584782 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[2] += -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[2] += 0.039518459830251443 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[3] += -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[3] += -0.00050341083159401254 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[3] += 1.4905346203698526 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[3] += -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[3] += 0.067736905936386635 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[4] += 0.033960047754138974 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[4] += 3.8193652822306167E-5 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[4] += -0.039518459830266792 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[4] += 0.067736905936572681 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[4] += -0.0072896299361453719 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[0] += -3.3182263091979736 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[1] += -0.0034921698713633515 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[2] += 0.13279109186599056 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[3] += -0.54052890870000092 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[4] += 0.042027661214546957 * AHV_Model_B.nu_r[2];

  // Derivatives for StateSpace: '<S310>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE[0] += -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[0] += -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[0] += 2.1365278323686749 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[0] += 1.6362305234027079 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[1] += 0.9998076276438177 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[1] += -5.3039287993722015E-5 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[1] += 0.0077107435119831529 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[1] += 0.0063309425482545112 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[1] += 0.0012479014747342309 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[2] += -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[2] += 0.0077107435118432491 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[2] += -3.350088004888828 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[2] += -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[2] += -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[3] += -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[3] += 0.0063309425481496455 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[3] += -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[3] += -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[3] += -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[4] += 0.2597101217072279 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[4] += -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[4] += 1.6370488166464041 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[4] += 3.2557107662584994 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[4] += -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[0] += -3.5374929885435322 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp35_CSTATE[1] += 0.015233386783081164 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp35_CSTATE[2] += -1.2196843916515756 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp35_CSTATE[3] += -0.96955326725759594 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp35_CSTATE[4] += 0.17148246208757431 * AHV_Model_B.nu_r[4];

  // Derivatives for StateSpace: '<S310>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE[0] += -0.0069962390232986621 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[0] += -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[0] += -0.074692039015985229 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[0] += -0.0042999629253483171 * AHV_Model_X.Dp42_CSTATE[4];
  _rtXdot->Dp42_CSTATE[1] += -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[1] += -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[1] += 0.95129605479096324 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[1] += 1.3710797263241543 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[1] += 0.067632454414034648 * AHV_Model_X.Dp42_CSTATE[4];
  _rtXdot->Dp42_CSTATE[2] += 0.078929006579498431 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[2] += 0.9512960547909608 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[2] += -1.065705011247746 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[2] += -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[2] += -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE[4];
  _rtXdot->Dp42_CSTATE[3] += -0.074692039015982259 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[3] += -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[3] += 3.3410455288652465 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[3] += -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[3] += -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE[4];
  _rtXdot->Dp42_CSTATE[4] += 0.0042999629253487083 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[4] += 0.067632454414049734 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[4] += -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[4] += -0.020742386481004484 * AHV_Model_X.Dp42_CSTATE[4];
  _rtXdot->Dp42_CSTATE[0] += -0.23872786278308805 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp42_CSTATE[1] += -3.2763464234276056 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp42_CSTATE[2] += 1.1437118751635387 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp42_CSTATE[3] += -1.3244904674438165 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp42_CSTATE[4] += 0.071163531145327891 * AHV_Model_B.nu_r[1];

  // Derivatives for StateSpace: '<S310>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE[0] += -0.004396104715914141 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[0] += 1.3927195985501903 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[0] += -0.027511960605097124 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[0] += 0.065100301994703638 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[0] += -0.0003953747088775357 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[1] += -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[1] += -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[1] += 0.47138330474887147 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[1] += -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[1] += 0.008356128979928646 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[2] += 0.02751196060511972 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[2] += 0.47138330474889245 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[2] += -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[2] += 2.5521376035497743 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[2] += -0.0093184009624265578 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[3] += 0.06510030199471134 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[3] += 1.5107242348458856 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[3] += -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[3] += -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[3] += 0.12040468647842413 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[4] += 0.00039537470888004991 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[4] += 0.0083561289799246388 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[4] += -0.0093184009624202816 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[4] += -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[4] += -0.00092959984608804514 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[0] += -0.33524288945878539 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[1] += -6.3625492730866693 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[2] += 0.93578679415036736 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[3] += 2.5673752417819116 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[4] += 0.014933981938330865 * AHV_Model_B.nu_r[5];

  // Derivatives for StateSpace: '<S310>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE[0] += -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[0] += -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[0] += 2.1365278323686749 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[0] += 1.6362305234027079 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[1] += 0.9998076276438177 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[1] += -5.3039287993722015E-5 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[1] += 0.0077107435119831529 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[1] += 0.0063309425482545112 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[1] += 0.0012479014747342309 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[2] += -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[2] += 0.0077107435118432491 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[2] += -3.350088004888828 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[2] += -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[2] += -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[3] += -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[3] += 0.0063309425481496455 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[3] += -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[3] += -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[3] += -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[4] += 0.2597101217072279 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[4] += -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[4] += 1.6370488166464041 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[4] += 3.2557107662584994 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[4] += -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[0] += -3.5374929885435322 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp53_CSTATE[1] += 0.015233386783081164 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp53_CSTATE[2] += -1.2196843916515756 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp53_CSTATE[3] += -0.96955326725759594 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp53_CSTATE[4] += 0.17148246208757431 * AHV_Model_B.nu_r[2];

  // Derivatives for StateSpace: '<S310>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE[0] += -8.7466634083809418E-5 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[0] += 0.70495353490101365 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[0] += 0.0013720743421809071 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[0] += 0.0075482812368147384 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[0] += 0.00709599540037472 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[1] += -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[1] += -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[1] += -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[1] += -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[1] += -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[2] += -0.0013720743421854747 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[2] += -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[2] += -0.068186331750790974 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[2] += -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[3] += 0.007548281236770463 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[3] += 1.3901163749547758 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[3] += 4.2372658827383063 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[3] += -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[3] += -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[4] += 0.0070959954003254081 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[4] += 1.2786273904776095 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[4] += 2.2216067374575723 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[4] += -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[4] += -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[0] += 0.021904981880023978 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp55_CSTATE[1] += 3.4651622640804733 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp55_CSTATE[2] += 0.16004490113712252 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp55_CSTATE[3] += -0.9969365784863482 * AHV_Model_B.nu_r[4];
  _rtXdot->Dp55_CSTATE[4] += -0.92774837285087253 * AHV_Model_B.nu_r[4];

  // Derivatives for StateSpace: '<S310>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE[0] += -0.0017008072822443973 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[0] += 1.2361995937266759 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[0] += 0.012808448456780379 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[0] += 0.031758809083141569 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[0] += -0.0047828760333047419 * AHV_Model_X.Dp62_CSTATE[4];
  _rtXdot->Dp62_CSTATE[1] += -1.236199593726705 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[1] += -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[1] += -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[1] += -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[1] += 0.12995794495601007 * AHV_Model_X.Dp62_CSTATE[4];
  _rtXdot->Dp62_CSTATE[2] += -0.012808448456765344 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[2] += -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[2] += -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[2] += -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[2] += 0.1356516649888091 * AHV_Model_X.Dp62_CSTATE[4];
  _rtXdot->Dp62_CSTATE[3] += 0.031758809083164766 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[3] += 0.96895750705376427 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[3] += 1.8113135858263727 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[3] += -3.762337151306161 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[3] += 1.5320613622870083 * AHV_Model_X.Dp62_CSTATE[4];
  _rtXdot->Dp62_CSTATE[4] += 0.0047828760333030167 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[4] += 0.12995794495600341 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[4] += 0.13565166498879488 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[4] += -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[4] += -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE[4];
  _rtXdot->Dp62_CSTATE[0] += 0.13167647316600378 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp62_CSTATE[1] += 3.4125284827245985 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp62_CSTATE[2] += 0.44399732492152644 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp62_CSTATE[3] += -1.2820687298457427 * AHV_Model_B.nu_r[1];
  _rtXdot->Dp62_CSTATE[4] += -0.18276601298221137 * AHV_Model_B.nu_r[1];

  // Derivatives for StateSpace: '<S310>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE[0] += -0.004396104715914141 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[0] += 1.3927195985501903 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[0] += -0.027511960605097124 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[0] += 0.065100301994703638 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[0] += -0.0003953747088775357 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[1] += -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[1] += -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[1] += 0.47138330474887147 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[1] += -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[1] += 0.008356128979928646 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[2] += 0.02751196060511972 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[2] += 0.47138330474889245 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[2] += -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[2] += 2.5521376035497743 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[2] += -0.0093184009624265578 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[3] += 0.06510030199471134 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[3] += 1.5107242348458856 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[3] += -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[3] += -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[3] += 0.12040468647842413 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[4] += 0.00039537470888004991 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[4] += 0.0083561289799246388 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[4] += -0.0093184009624202816 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[4] += -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[4] += -0.00092959984608804514 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[0] += -0.33524288945878539 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[1] += -6.3625492730866693 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[2] += 0.93578679415036736 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[3] += 2.5673752417819116 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[4] += 0.014933981938330865 * AHV_Model_B.nu_r[3];

  // Derivatives for StateSpace: '<S310>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE[0] += -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[0] += -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[0] += -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[0] += -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[0] += 0.069477290559771143 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[1] += 1.3217886535783951 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[1] += -0.00029199832472599611 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[1] += -0.0080078570280257711 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[1] += -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[1] += 0.00098472174210678556 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[2] += -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[2] += 0.00800785702802341 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[2] += -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[2] += -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[2] += 0.074587161707531782 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[3] += 1.2174625158903034 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[3] += -0.015684008067190565 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[3] += 1.8221092762782234 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[3] += -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[4] += 0.069477290559780386 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[4] += -0.000984721742106249 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[4] += 0.074587161707547089 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[4] += -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[4] += -0.061420217921035136 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[0] += -3.3844821061206916 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp66_CSTATE[1] += 0.045246867266885989 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp66_CSTATE[2] += -0.53332191483729374 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp66_CSTATE[3] += 1.2581308985447515 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp66_CSTATE[4] += 0.075289942321940723 * AHV_Model_B.nu_r[5];

  // Derivatives for Integrator: '<S207>/Integrator'
  for (i = 0; i < 6; i++) {
    _rtXdot->Integrator_CSTATE_h[i] = AHV_Model_B.Minvtau_n[i];
  }

  // End of Derivatives for Integrator: '<S207>/Integrator'

  // Derivatives for Integrator: '<S257>/Integrator3'
  _rtXdot->Integrator3_CSTATE_c[0] = AHV_Model_B.sun_k2_e3[0];

  // Derivatives for Integrator: '<S257>/Integrator4'
  _rtXdot->Integrator4_CSTATE_f[0] = AHV_Model_B.M_u_i[0];

  // Derivatives for Integrator: '<S257>/Integrator3'
  _rtXdot->Integrator3_CSTATE_c[1] = AHV_Model_B.sun_k2_e3[1];

  // Derivatives for Integrator: '<S257>/Integrator4'
  _rtXdot->Integrator4_CSTATE_f[1] = AHV_Model_B.M_u_i[1];

  // Derivatives for Integrator: '<S257>/Integrator3'
  _rtXdot->Integrator3_CSTATE_c[2] = AHV_Model_B.sun_k2_e3[2];

  // Derivatives for Integrator: '<S257>/Integrator4'
  _rtXdot->Integrator4_CSTATE_f[2] = AHV_Model_B.M_u_i[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S214>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_c[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_j[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_b[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_n[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_m[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_b[i] = 0.0;

    // Derivatives for Integrator: '<S226>/Integrator'
    _rtXdot->Integrator_CSTATE_c[i] = AHV_Model_B.Sum_n[i];

    // Derivatives for StateSpace: '<S214>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_o[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_g[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S214>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_f[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S214>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_l[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_l[0];
  _rtXdot->Dp11_CSTATE_l[0] += -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_l[1];
  _rtXdot->Dp11_CSTATE_l[0] += -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_l
    [2];
  _rtXdot->Dp11_CSTATE_l[0] += -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_l[3];
  _rtXdot->Dp11_CSTATE_l[0] += -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_l[4];
  _rtXdot->Dp11_CSTATE_l[1] += 1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_l[0];
  _rtXdot->Dp11_CSTATE_l[1] += -0.00075311499999405624 *
    AHV_Model_X.Dp11_CSTATE_l[1];
  _rtXdot->Dp11_CSTATE_l[1] += -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_l
    [2];
  _rtXdot->Dp11_CSTATE_l[1] += -0.023103995048734959 *
    AHV_Model_X.Dp11_CSTATE_l[3];
  _rtXdot->Dp11_CSTATE_l[1] += -0.0023400735275004216 *
    AHV_Model_X.Dp11_CSTATE_l[4];
  _rtXdot->Dp11_CSTATE_l[2] += -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_l
    [0];
  _rtXdot->Dp11_CSTATE_l[2] += 0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_l
    [1];
  _rtXdot->Dp11_CSTATE_l[2] += -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_l
    [2];
  _rtXdot->Dp11_CSTATE_l[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_l[3];
  _rtXdot->Dp11_CSTATE_l[2] += -0.089746811088073059 *
    AHV_Model_X.Dp11_CSTATE_l[4];
  _rtXdot->Dp11_CSTATE_l[3] += 1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_l[0];
  _rtXdot->Dp11_CSTATE_l[3] += -0.023103995048766656 *
    AHV_Model_X.Dp11_CSTATE_l[1];
  _rtXdot->Dp11_CSTATE_l[3] += 1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_l[2];
  _rtXdot->Dp11_CSTATE_l[3] += -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_l[3];
  _rtXdot->Dp11_CSTATE_l[3] += -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_l[4];
  _rtXdot->Dp11_CSTATE_l[4] += -0.095423424399719625 *
    AHV_Model_X.Dp11_CSTATE_l[0];
  _rtXdot->Dp11_CSTATE_l[4] += 0.0023400735275035814 *
    AHV_Model_X.Dp11_CSTATE_l[1];
  _rtXdot->Dp11_CSTATE_l[4] += -0.089746811088075445 *
    AHV_Model_X.Dp11_CSTATE_l[2];
  _rtXdot->Dp11_CSTATE_l[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_l[3];
  _rtXdot->Dp11_CSTATE_l[4] += -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_l
    [4];
  _rtXdot->Dp11_CSTATE_l[0] += -3.3881001186476949 * AHV_Model_B.nu_r_m[0];
  _rtXdot->Dp11_CSTATE_l[1] += 0.07669611088588818 * AHV_Model_B.nu_r_m[0];
  _rtXdot->Dp11_CSTATE_l[2] += -0.46181267830731043 * AHV_Model_B.nu_r_m[0];
  _rtXdot->Dp11_CSTATE_l[3] += 1.226085627712131 * AHV_Model_B.nu_r_m[0];
  _rtXdot->Dp11_CSTATE_l[4] += -0.11754627904442222 * AHV_Model_B.nu_r_m[0];

  // Derivatives for StateSpace: '<S214>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_j[0] += -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE_j[0];
  _rtXdot->Dp22_CSTATE_j[0] += -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE_j[1];
  _rtXdot->Dp22_CSTATE_j[0] += -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE_j
    [2];
  _rtXdot->Dp22_CSTATE_j[0] += -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE_j[3];
  _rtXdot->Dp22_CSTATE_j[0] += -0.095423424399736376 *
    AHV_Model_X.Dp22_CSTATE_j[4];
  _rtXdot->Dp22_CSTATE_j[1] += 1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_j[0];
  _rtXdot->Dp22_CSTATE_j[1] += -0.00075311499999501316 *
    AHV_Model_X.Dp22_CSTATE_j[1];
  _rtXdot->Dp22_CSTATE_j[1] += -0.010562919761942575 *
    AHV_Model_X.Dp22_CSTATE_j[2];
  _rtXdot->Dp22_CSTATE_j[1] += -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE_j[3];
  _rtXdot->Dp22_CSTATE_j[1] += -0.002340073527504722 *
    AHV_Model_X.Dp22_CSTATE_j[4];
  _rtXdot->Dp22_CSTATE_j[2] += -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE_j
    [0];
  _rtXdot->Dp22_CSTATE_j[2] += 0.010562919761928094 * AHV_Model_X.Dp22_CSTATE_j
    [1];
  _rtXdot->Dp22_CSTATE_j[2] += -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE_j
    [2];
  _rtXdot->Dp22_CSTATE_j[2] += -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE_j[3];
  _rtXdot->Dp22_CSTATE_j[2] += -0.089746811088091308 *
    AHV_Model_X.Dp22_CSTATE_j[4];
  _rtXdot->Dp22_CSTATE_j[3] += 1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_j[0];
  _rtXdot->Dp22_CSTATE_j[3] += -0.023103995048724405 *
    AHV_Model_X.Dp22_CSTATE_j[1];
  _rtXdot->Dp22_CSTATE_j[3] += 1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_j[2];
  _rtXdot->Dp22_CSTATE_j[3] += -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE_j[3];
  _rtXdot->Dp22_CSTATE_j[3] += -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE_j[4];
  _rtXdot->Dp22_CSTATE_j[4] += -0.095423424399735959 *
    AHV_Model_X.Dp22_CSTATE_j[0];
  _rtXdot->Dp22_CSTATE_j[4] += 0.00234007352750089 * AHV_Model_X.Dp22_CSTATE_j[1];
  _rtXdot->Dp22_CSTATE_j[4] += -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE_j[2];
  _rtXdot->Dp22_CSTATE_j[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_j[3];
  _rtXdot->Dp22_CSTATE_j[4] += -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE_j
    [4];
  _rtXdot->Dp22_CSTATE_j[0] += -3.3881001186476967 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp22_CSTATE_j[1] += 0.076696110885761712 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp22_CSTATE_j[2] += -0.46181267830733025 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp22_CSTATE_j[3] += 1.2260856277121315 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp22_CSTATE_j[4] += -0.11754627904444465 * AHV_Model_B.nu_r_m[1];

  // Derivatives for StateSpace: '<S214>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_b[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp24_CSTATE_b[0];
  _rtXdot->Dp24_CSTATE_b[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_b[1];
  _rtXdot->Dp24_CSTATE_b[0] += -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_b
    [2];
  _rtXdot->Dp24_CSTATE_b[0] += -0.074692039015985229 *
    AHV_Model_X.Dp24_CSTATE_b[3];
  _rtXdot->Dp24_CSTATE_b[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp24_CSTATE_b[4];
  _rtXdot->Dp24_CSTATE_b[1] += -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_b[0];
  _rtXdot->Dp24_CSTATE_b[1] += -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_b[1];
  _rtXdot->Dp24_CSTATE_b[1] += 0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_b[2];
  _rtXdot->Dp24_CSTATE_b[1] += 1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_b[3];
  _rtXdot->Dp24_CSTATE_b[1] += 0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_b
    [4];
  _rtXdot->Dp24_CSTATE_b[2] += 0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_b
    [0];
  _rtXdot->Dp24_CSTATE_b[2] += 0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_b[1];
  _rtXdot->Dp24_CSTATE_b[2] += -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_b[2];
  _rtXdot->Dp24_CSTATE_b[2] += -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_b[3];
  _rtXdot->Dp24_CSTATE_b[2] += -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_b
    [4];
  _rtXdot->Dp24_CSTATE_b[3] += -0.074692039015982259 *
    AHV_Model_X.Dp24_CSTATE_b[0];
  _rtXdot->Dp24_CSTATE_b[3] += -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_b[1];
  _rtXdot->Dp24_CSTATE_b[3] += 3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_b[2];
  _rtXdot->Dp24_CSTATE_b[3] += -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_b[3];
  _rtXdot->Dp24_CSTATE_b[3] += -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_b[4];
  _rtXdot->Dp24_CSTATE_b[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp24_CSTATE_b[0];
  _rtXdot->Dp24_CSTATE_b[4] += 0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_b
    [1];
  _rtXdot->Dp24_CSTATE_b[4] += -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_b
    [2];
  _rtXdot->Dp24_CSTATE_b[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_b[3];
  _rtXdot->Dp24_CSTATE_b[4] += -0.020742386481004484 *
    AHV_Model_X.Dp24_CSTATE_b[4];
  _rtXdot->Dp24_CSTATE_b[0] += -0.23872786278308805 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp24_CSTATE_b[1] += -3.2763464234276056 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp24_CSTATE_b[2] += 1.1437118751635387 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp24_CSTATE_b[3] += -1.3244904674438165 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp24_CSTATE_b[4] += 0.071163531145327891 * AHV_Model_B.nu_r_m[3];

  // Derivatives for StateSpace: '<S214>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_n[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp26_CSTATE_n[0];
  _rtXdot->Dp26_CSTATE_n[0] += 1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_n[1];
  _rtXdot->Dp26_CSTATE_n[0] += 0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_n
    [2];
  _rtXdot->Dp26_CSTATE_n[0] += 0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_n
    [3];
  _rtXdot->Dp26_CSTATE_n[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp26_CSTATE_n[4];
  _rtXdot->Dp26_CSTATE_n[1] += -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_n[0];
  _rtXdot->Dp26_CSTATE_n[1] += -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_n[1];
  _rtXdot->Dp26_CSTATE_n[1] += -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_n
    [2];
  _rtXdot->Dp26_CSTATE_n[1] += -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_n
    [3];
  _rtXdot->Dp26_CSTATE_n[1] += 0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_n[4];
  _rtXdot->Dp26_CSTATE_n[2] += -0.012808448456765344 *
    AHV_Model_X.Dp26_CSTATE_n[0];
  _rtXdot->Dp26_CSTATE_n[2] += -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_n
    [1];
  _rtXdot->Dp26_CSTATE_n[2] += -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_n
    [2];
  _rtXdot->Dp26_CSTATE_n[2] += -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_n[3];
  _rtXdot->Dp26_CSTATE_n[2] += 0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_n[4];
  _rtXdot->Dp26_CSTATE_n[3] += 0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_n
    [0];
  _rtXdot->Dp26_CSTATE_n[3] += 0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_n[1];
  _rtXdot->Dp26_CSTATE_n[3] += 1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_n[2];
  _rtXdot->Dp26_CSTATE_n[3] += -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_n[3];
  _rtXdot->Dp26_CSTATE_n[3] += 1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_n[4];
  _rtXdot->Dp26_CSTATE_n[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp26_CSTATE_n[0];
  _rtXdot->Dp26_CSTATE_n[4] += 0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_n[1];
  _rtXdot->Dp26_CSTATE_n[4] += 0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_n[2];
  _rtXdot->Dp26_CSTATE_n[4] += -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_n[3];
  _rtXdot->Dp26_CSTATE_n[4] += -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_n
    [4];
  _rtXdot->Dp26_CSTATE_n[0] += 0.13167647316600378 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp26_CSTATE_n[1] += 3.4125284827245985 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp26_CSTATE_n[2] += 0.44399732492152644 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp26_CSTATE_n[3] += -1.2820687298457427 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp26_CSTATE_n[4] += -0.18276601298221137 * AHV_Model_B.nu_r_m[5];

  // Derivatives for StateSpace: '<S214>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_m[0] += -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[0] += -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_m
    [2];
  _rtXdot->Dp33_CSTATE_m[0] += -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_m
    [3];
  _rtXdot->Dp33_CSTATE_m[0] += 0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_m
    [4];
  _rtXdot->Dp33_CSTATE_m[1] += -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_m
    [0];
  _rtXdot->Dp33_CSTATE_m[1] += -1.538297528844471E-6 *
    AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[1] += 0.00010972576839633958 *
    AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[1] += 0.00050341083172827721 *
    AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[1] += -3.8193652830025795E-5 *
    AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[2] += 0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[2] += 0.0001097257683747581 *
    AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[2] += -0.033603468126584782 *
    AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[2] += -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[2] += 0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_m
    [4];
  _rtXdot->Dp33_CSTATE_m[3] += -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_m
    [0];
  _rtXdot->Dp33_CSTATE_m[3] += -0.00050341083159401254 *
    AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[3] += 1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[3] += -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_m
    [3];
  _rtXdot->Dp33_CSTATE_m[3] += 0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_m
    [4];
  _rtXdot->Dp33_CSTATE_m[4] += 0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_m
    [0];
  _rtXdot->Dp33_CSTATE_m[4] += 3.8193652822306167E-5 *
    AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[4] += -0.039518459830266792 *
    AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[4] += 0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_m
    [3];
  _rtXdot->Dp33_CSTATE_m[4] += -0.0072896299361453719 *
    AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[0] += -3.3182263091979736 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp33_CSTATE_m[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp33_CSTATE_m[2] += 0.13279109186599056 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp33_CSTATE_m[3] += -0.54052890870000092 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp33_CSTATE_m[4] += 0.042027661214546957 * AHV_Model_B.nu_r_m[2];

  // Derivatives for StateSpace: '<S214>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_f[0] += -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_f[0];
  _rtXdot->Dp35_CSTATE_f[0] += -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_f
    [1];
  _rtXdot->Dp35_CSTATE_f[0] += 2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_f[2];
  _rtXdot->Dp35_CSTATE_f[0] += 1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_f[3];
  _rtXdot->Dp35_CSTATE_f[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_f[4];
  _rtXdot->Dp35_CSTATE_f[1] += 0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_f[0];
  _rtXdot->Dp35_CSTATE_f[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp35_CSTATE_f[1];
  _rtXdot->Dp35_CSTATE_f[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp35_CSTATE_f[2];
  _rtXdot->Dp35_CSTATE_f[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp35_CSTATE_f[3];
  _rtXdot->Dp35_CSTATE_f[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp35_CSTATE_f[4];
  _rtXdot->Dp35_CSTATE_f[2] += -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_f[0];
  _rtXdot->Dp35_CSTATE_f[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp35_CSTATE_f[1];
  _rtXdot->Dp35_CSTATE_f[2] += -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_f[2];
  _rtXdot->Dp35_CSTATE_f[2] += -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_f[3];
  _rtXdot->Dp35_CSTATE_f[2] += -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_f[4];
  _rtXdot->Dp35_CSTATE_f[3] += -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_f[0];
  _rtXdot->Dp35_CSTATE_f[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp35_CSTATE_f[1];
  _rtXdot->Dp35_CSTATE_f[3] += -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_f[2];
  _rtXdot->Dp35_CSTATE_f[3] += -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_f[3];
  _rtXdot->Dp35_CSTATE_f[3] += -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_f[4];
  _rtXdot->Dp35_CSTATE_f[4] += 0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_f[0];
  _rtXdot->Dp35_CSTATE_f[4] += -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_f
    [1];
  _rtXdot->Dp35_CSTATE_f[4] += 1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_f[2];
  _rtXdot->Dp35_CSTATE_f[4] += 3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_f[3];
  _rtXdot->Dp35_CSTATE_f[4] += -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_f
    [4];
  _rtXdot->Dp35_CSTATE_f[0] += -3.5374929885435322 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp35_CSTATE_f[1] += 0.015233386783081164 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp35_CSTATE_f[2] += -1.2196843916515756 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp35_CSTATE_f[3] += -0.96955326725759594 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp35_CSTATE_f[4] += 0.17148246208757431 * AHV_Model_B.nu_r_m[4];

  // Derivatives for StateSpace: '<S214>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_b[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp42_CSTATE_b[0];
  _rtXdot->Dp42_CSTATE_b[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_b[1];
  _rtXdot->Dp42_CSTATE_b[0] += -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_b
    [2];
  _rtXdot->Dp42_CSTATE_b[0] += -0.074692039015985229 *
    AHV_Model_X.Dp42_CSTATE_b[3];
  _rtXdot->Dp42_CSTATE_b[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp42_CSTATE_b[4];
  _rtXdot->Dp42_CSTATE_b[1] += -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_b[0];
  _rtXdot->Dp42_CSTATE_b[1] += -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_b[1];
  _rtXdot->Dp42_CSTATE_b[1] += 0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_b[2];
  _rtXdot->Dp42_CSTATE_b[1] += 1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_b[3];
  _rtXdot->Dp42_CSTATE_b[1] += 0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_b
    [4];
  _rtXdot->Dp42_CSTATE_b[2] += 0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_b
    [0];
  _rtXdot->Dp42_CSTATE_b[2] += 0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_b[1];
  _rtXdot->Dp42_CSTATE_b[2] += -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_b[2];
  _rtXdot->Dp42_CSTATE_b[2] += -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_b[3];
  _rtXdot->Dp42_CSTATE_b[2] += -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_b
    [4];
  _rtXdot->Dp42_CSTATE_b[3] += -0.074692039015982259 *
    AHV_Model_X.Dp42_CSTATE_b[0];
  _rtXdot->Dp42_CSTATE_b[3] += -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_b[1];
  _rtXdot->Dp42_CSTATE_b[3] += 3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_b[2];
  _rtXdot->Dp42_CSTATE_b[3] += -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_b[3];
  _rtXdot->Dp42_CSTATE_b[3] += -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_b[4];
  _rtXdot->Dp42_CSTATE_b[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp42_CSTATE_b[0];
  _rtXdot->Dp42_CSTATE_b[4] += 0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_b
    [1];
  _rtXdot->Dp42_CSTATE_b[4] += -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_b
    [2];
  _rtXdot->Dp42_CSTATE_b[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_b[3];
  _rtXdot->Dp42_CSTATE_b[4] += -0.020742386481004484 *
    AHV_Model_X.Dp42_CSTATE_b[4];
  _rtXdot->Dp42_CSTATE_b[0] += -0.23872786278308805 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp42_CSTATE_b[1] += -3.2763464234276056 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp42_CSTATE_b[2] += 1.1437118751635387 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp42_CSTATE_b[3] += -1.3244904674438165 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp42_CSTATE_b[4] += 0.071163531145327891 * AHV_Model_B.nu_r_m[1];

  // Derivatives for StateSpace: '<S214>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_o[0] += -0.004396104715914141 *
    AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[0] += 1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[0] += -0.027511960605097124 *
    AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[0] += 0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_o
    [3];
  _rtXdot->Dp46_CSTATE_o[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[1] += -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[1] += -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[1] += 0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[1] += -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[1] += 0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_o
    [4];
  _rtXdot->Dp46_CSTATE_o[2] += 0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[2] += 0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[2] += -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_o
    [2];
  _rtXdot->Dp46_CSTATE_o[2] += 2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[3] += 0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[3] += 1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[3] += -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[3] += -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[3] += 0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[4] += -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_o
    [3];
  _rtXdot->Dp46_CSTATE_o[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[0] += -0.33524288945878539 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp46_CSTATE_o[1] += -6.3625492730866693 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp46_CSTATE_o[2] += 0.93578679415036736 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp46_CSTATE_o[3] += 2.5673752417819116 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp46_CSTATE_o[4] += 0.014933981938330865 * AHV_Model_B.nu_r_m[5];

  // Derivatives for StateSpace: '<S214>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_i[0] += -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_i[0];
  _rtXdot->Dp53_CSTATE_i[0] += -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_i
    [1];
  _rtXdot->Dp53_CSTATE_i[0] += 2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_i[2];
  _rtXdot->Dp53_CSTATE_i[0] += 1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_i[3];
  _rtXdot->Dp53_CSTATE_i[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE_i[4];
  _rtXdot->Dp53_CSTATE_i[1] += 0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_i[0];
  _rtXdot->Dp53_CSTATE_i[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp53_CSTATE_i[1];
  _rtXdot->Dp53_CSTATE_i[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp53_CSTATE_i[2];
  _rtXdot->Dp53_CSTATE_i[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp53_CSTATE_i[3];
  _rtXdot->Dp53_CSTATE_i[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp53_CSTATE_i[4];
  _rtXdot->Dp53_CSTATE_i[2] += -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_i[0];
  _rtXdot->Dp53_CSTATE_i[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp53_CSTATE_i[1];
  _rtXdot->Dp53_CSTATE_i[2] += -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_i[2];
  _rtXdot->Dp53_CSTATE_i[2] += -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_i[3];
  _rtXdot->Dp53_CSTATE_i[2] += -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_i[4];
  _rtXdot->Dp53_CSTATE_i[3] += -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_i[0];
  _rtXdot->Dp53_CSTATE_i[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp53_CSTATE_i[1];
  _rtXdot->Dp53_CSTATE_i[3] += -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_i[2];
  _rtXdot->Dp53_CSTATE_i[3] += -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_i[3];
  _rtXdot->Dp53_CSTATE_i[3] += -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_i[4];
  _rtXdot->Dp53_CSTATE_i[4] += 0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_i[0];
  _rtXdot->Dp53_CSTATE_i[4] += -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_i
    [1];
  _rtXdot->Dp53_CSTATE_i[4] += 1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_i[2];
  _rtXdot->Dp53_CSTATE_i[4] += 3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_i[3];
  _rtXdot->Dp53_CSTATE_i[4] += -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_i
    [4];
  _rtXdot->Dp53_CSTATE_i[0] += -3.5374929885435322 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp53_CSTATE_i[1] += 0.015233386783081164 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp53_CSTATE_i[2] += -1.2196843916515756 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp53_CSTATE_i[3] += -0.96955326725759594 * AHV_Model_B.nu_r_m[2];
  _rtXdot->Dp53_CSTATE_i[4] += 0.17148246208757431 * AHV_Model_B.nu_r_m[2];

  // Derivatives for StateSpace: '<S214>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_k[0] += -8.7466634083809418E-5 *
    AHV_Model_X.Dp55_CSTATE_k[0];
  _rtXdot->Dp55_CSTATE_k[0] += 0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_k[1];
  _rtXdot->Dp55_CSTATE_k[0] += 0.0013720743421809071 *
    AHV_Model_X.Dp55_CSTATE_k[2];
  _rtXdot->Dp55_CSTATE_k[0] += 0.0075482812368147384 *
    AHV_Model_X.Dp55_CSTATE_k[3];
  _rtXdot->Dp55_CSTATE_k[0] += 0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_k[4];
  _rtXdot->Dp55_CSTATE_k[1] += -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_k
    [0];
  _rtXdot->Dp55_CSTATE_k[1] += -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_k[1];
  _rtXdot->Dp55_CSTATE_k[1] += -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_k
    [2];
  _rtXdot->Dp55_CSTATE_k[1] += -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_k[3];
  _rtXdot->Dp55_CSTATE_k[1] += -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_k[4];
  _rtXdot->Dp55_CSTATE_k[2] += -0.0013720743421854747 *
    AHV_Model_X.Dp55_CSTATE_k[0];
  _rtXdot->Dp55_CSTATE_k[2] += -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_k
    [1];
  _rtXdot->Dp55_CSTATE_k[2] += -0.068186331750790974 *
    AHV_Model_X.Dp55_CSTATE_k[2];
  _rtXdot->Dp55_CSTATE_k[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_k[3];
  _rtXdot->Dp55_CSTATE_k[2] += -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_k[4];
  _rtXdot->Dp55_CSTATE_k[3] += 0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_k
    [0];
  _rtXdot->Dp55_CSTATE_k[3] += 1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_k[1];
  _rtXdot->Dp55_CSTATE_k[3] += 4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_k[2];
  _rtXdot->Dp55_CSTATE_k[3] += -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_k[3];
  _rtXdot->Dp55_CSTATE_k[3] += -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_k[4];
  _rtXdot->Dp55_CSTATE_k[4] += 0.0070959954003254081 *
    AHV_Model_X.Dp55_CSTATE_k[0];
  _rtXdot->Dp55_CSTATE_k[4] += 1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_k[1];
  _rtXdot->Dp55_CSTATE_k[4] += 2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_k[2];
  _rtXdot->Dp55_CSTATE_k[4] += -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_k[3];
  _rtXdot->Dp55_CSTATE_k[4] += -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_k[4];
  _rtXdot->Dp55_CSTATE_k[0] += 0.021904981880023978 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp55_CSTATE_k[1] += 3.4651622640804733 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp55_CSTATE_k[2] += 0.16004490113712252 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp55_CSTATE_k[3] += -0.9969365784863482 * AHV_Model_B.nu_r_m[4];
  _rtXdot->Dp55_CSTATE_k[4] += -0.92774837285087253 * AHV_Model_B.nu_r_m[4];

  // Derivatives for StateSpace: '<S214>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_g[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp62_CSTATE_g[0];
  _rtXdot->Dp62_CSTATE_g[0] += 1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_g[1];
  _rtXdot->Dp62_CSTATE_g[0] += 0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_g
    [2];
  _rtXdot->Dp62_CSTATE_g[0] += 0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_g
    [3];
  _rtXdot->Dp62_CSTATE_g[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp62_CSTATE_g[4];
  _rtXdot->Dp62_CSTATE_g[1] += -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_g[0];
  _rtXdot->Dp62_CSTATE_g[1] += -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_g[1];
  _rtXdot->Dp62_CSTATE_g[1] += -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_g
    [2];
  _rtXdot->Dp62_CSTATE_g[1] += -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_g
    [3];
  _rtXdot->Dp62_CSTATE_g[1] += 0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_g[4];
  _rtXdot->Dp62_CSTATE_g[2] += -0.012808448456765344 *
    AHV_Model_X.Dp62_CSTATE_g[0];
  _rtXdot->Dp62_CSTATE_g[2] += -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_g
    [1];
  _rtXdot->Dp62_CSTATE_g[2] += -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_g
    [2];
  _rtXdot->Dp62_CSTATE_g[2] += -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_g[3];
  _rtXdot->Dp62_CSTATE_g[2] += 0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_g[4];
  _rtXdot->Dp62_CSTATE_g[3] += 0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_g
    [0];
  _rtXdot->Dp62_CSTATE_g[3] += 0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_g[1];
  _rtXdot->Dp62_CSTATE_g[3] += 1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_g[2];
  _rtXdot->Dp62_CSTATE_g[3] += -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_g[3];
  _rtXdot->Dp62_CSTATE_g[3] += 1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_g[4];
  _rtXdot->Dp62_CSTATE_g[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp62_CSTATE_g[0];
  _rtXdot->Dp62_CSTATE_g[4] += 0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_g[1];
  _rtXdot->Dp62_CSTATE_g[4] += 0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_g[2];
  _rtXdot->Dp62_CSTATE_g[4] += -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_g[3];
  _rtXdot->Dp62_CSTATE_g[4] += -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_g
    [4];
  _rtXdot->Dp62_CSTATE_g[0] += 0.13167647316600378 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp62_CSTATE_g[1] += 3.4125284827245985 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp62_CSTATE_g[2] += 0.44399732492152644 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp62_CSTATE_g[3] += -1.2820687298457427 * AHV_Model_B.nu_r_m[1];
  _rtXdot->Dp62_CSTATE_g[4] += -0.18276601298221137 * AHV_Model_B.nu_r_m[1];

  // Derivatives for StateSpace: '<S214>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_p[0] += -0.004396104715914141 *
    AHV_Model_X.Dp64_CSTATE_p[0];
  _rtXdot->Dp64_CSTATE_p[0] += 1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_p[1];
  _rtXdot->Dp64_CSTATE_p[0] += -0.027511960605097124 *
    AHV_Model_X.Dp64_CSTATE_p[2];
  _rtXdot->Dp64_CSTATE_p[0] += 0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_p
    [3];
  _rtXdot->Dp64_CSTATE_p[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp64_CSTATE_p[4];
  _rtXdot->Dp64_CSTATE_p[1] += -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_p[0];
  _rtXdot->Dp64_CSTATE_p[1] += -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_p[1];
  _rtXdot->Dp64_CSTATE_p[1] += 0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_p[2];
  _rtXdot->Dp64_CSTATE_p[1] += -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_p[3];
  _rtXdot->Dp64_CSTATE_p[1] += 0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_p
    [4];
  _rtXdot->Dp64_CSTATE_p[2] += 0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_p[0];
  _rtXdot->Dp64_CSTATE_p[2] += 0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_p[1];
  _rtXdot->Dp64_CSTATE_p[2] += -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_p
    [2];
  _rtXdot->Dp64_CSTATE_p[2] += 2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_p[3];
  _rtXdot->Dp64_CSTATE_p[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp64_CSTATE_p[4];
  _rtXdot->Dp64_CSTATE_p[3] += 0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_p[0];
  _rtXdot->Dp64_CSTATE_p[3] += 1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_p[1];
  _rtXdot->Dp64_CSTATE_p[3] += -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_p[2];
  _rtXdot->Dp64_CSTATE_p[3] += -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_p[3];
  _rtXdot->Dp64_CSTATE_p[3] += 0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_p[4];
  _rtXdot->Dp64_CSTATE_p[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp64_CSTATE_p[0];
  _rtXdot->Dp64_CSTATE_p[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp64_CSTATE_p[1];
  _rtXdot->Dp64_CSTATE_p[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp64_CSTATE_p[2];
  _rtXdot->Dp64_CSTATE_p[4] += -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_p
    [3];
  _rtXdot->Dp64_CSTATE_p[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp64_CSTATE_p[4];
  _rtXdot->Dp64_CSTATE_p[0] += -0.33524288945878539 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp64_CSTATE_p[1] += -6.3625492730866693 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp64_CSTATE_p[2] += 0.93578679415036736 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp64_CSTATE_p[3] += 2.5673752417819116 * AHV_Model_B.nu_r_m[3];
  _rtXdot->Dp64_CSTATE_p[4] += 0.014933981938330865 * AHV_Model_B.nu_r_m[3];

  // Derivatives for StateSpace: '<S214>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_f[0] += -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE_f[0];
  _rtXdot->Dp66_CSTATE_f[0] += -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE_f[1];
  _rtXdot->Dp66_CSTATE_f[0] += -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE_f
    [2];
  _rtXdot->Dp66_CSTATE_f[0] += -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE_f[3];
  _rtXdot->Dp66_CSTATE_f[0] += 0.069477290559771143 * AHV_Model_X.Dp66_CSTATE_f
    [4];
  _rtXdot->Dp66_CSTATE_f[1] += 1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_f[0];
  _rtXdot->Dp66_CSTATE_f[1] += -0.00029199832472599611 *
    AHV_Model_X.Dp66_CSTATE_f[1];
  _rtXdot->Dp66_CSTATE_f[1] += -0.0080078570280257711 *
    AHV_Model_X.Dp66_CSTATE_f[2];
  _rtXdot->Dp66_CSTATE_f[1] += -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE_f
    [3];
  _rtXdot->Dp66_CSTATE_f[1] += 0.00098472174210678556 *
    AHV_Model_X.Dp66_CSTATE_f[4];
  _rtXdot->Dp66_CSTATE_f[2] += -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE_f
    [0];
  _rtXdot->Dp66_CSTATE_f[2] += 0.00800785702802341 * AHV_Model_X.Dp66_CSTATE_f[1];
  _rtXdot->Dp66_CSTATE_f[2] += -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE_f
    [2];
  _rtXdot->Dp66_CSTATE_f[2] += -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE_f[3];
  _rtXdot->Dp66_CSTATE_f[2] += 0.074587161707531782 * AHV_Model_X.Dp66_CSTATE_f
    [4];
  _rtXdot->Dp66_CSTATE_f[3] += 1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_f[0];
  _rtXdot->Dp66_CSTATE_f[3] += -0.015684008067190565 *
    AHV_Model_X.Dp66_CSTATE_f[1];
  _rtXdot->Dp66_CSTATE_f[3] += 1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_f[2];
  _rtXdot->Dp66_CSTATE_f[3] += -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE_f[3];
  _rtXdot->Dp66_CSTATE_f[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE_f[4];
  _rtXdot->Dp66_CSTATE_f[4] += 0.069477290559780386 * AHV_Model_X.Dp66_CSTATE_f
    [0];
  _rtXdot->Dp66_CSTATE_f[4] += -0.000984721742106249 *
    AHV_Model_X.Dp66_CSTATE_f[1];
  _rtXdot->Dp66_CSTATE_f[4] += 0.074587161707547089 * AHV_Model_X.Dp66_CSTATE_f
    [2];
  _rtXdot->Dp66_CSTATE_f[4] += -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE_f[3];
  _rtXdot->Dp66_CSTATE_f[4] += -0.061420217921035136 *
    AHV_Model_X.Dp66_CSTATE_f[4];
  _rtXdot->Dp66_CSTATE_f[0] += -3.3844821061206916 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp66_CSTATE_f[1] += 0.045246867266885989 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp66_CSTATE_f[2] += -0.53332191483729374 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp66_CSTATE_f[3] += 1.2581308985447515 * AHV_Model_B.nu_r_m[5];
  _rtXdot->Dp66_CSTATE_f[4] += 0.075289942321940723 * AHV_Model_B.nu_r_m[5];

  // Derivatives for Integrator: '<S111>/Integrator'
  for (i = 0; i < 6; i++) {
    _rtXdot->Integrator_CSTATE_o[i] = AHV_Model_B.Minvtau_f[i];
  }

  // End of Derivatives for Integrator: '<S111>/Integrator'

  // Derivatives for Integrator: '<S161>/Integrator3'
  _rtXdot->Integrator3_CSTATE_l[0] = AHV_Model_B.sun_k2_e[0];

  // Derivatives for Integrator: '<S161>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[0] = AHV_Model_B.M_u_o[0];

  // Derivatives for Integrator: '<S161>/Integrator3'
  _rtXdot->Integrator3_CSTATE_l[1] = AHV_Model_B.sun_k2_e[1];

  // Derivatives for Integrator: '<S161>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[1] = AHV_Model_B.M_u_o[1];

  // Derivatives for Integrator: '<S161>/Integrator3'
  _rtXdot->Integrator3_CSTATE_l[2] = AHV_Model_B.sun_k2_e[2];

  // Derivatives for Integrator: '<S161>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[2] = AHV_Model_B.M_u_o[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S118>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_o[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_cn[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_g[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_jf[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_j[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_d[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_j[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_g[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_g[i] = 0.0;

    // Derivatives for Integrator: '<S130>/Integrator'
    _rtXdot->Integrator_CSTATE_k[i] = AHV_Model_B.Sum_h[i];

    // Derivatives for StateSpace: '<S118>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S118>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_n[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S118>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_o[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_o[0];
  _rtXdot->Dp11_CSTATE_o[0] += -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_o[1];
  _rtXdot->Dp11_CSTATE_o[0] += -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_o
    [2];
  _rtXdot->Dp11_CSTATE_o[0] += -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_o[3];
  _rtXdot->Dp11_CSTATE_o[0] += -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_o[4];
  _rtXdot->Dp11_CSTATE_o[1] += 1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_o[0];
  _rtXdot->Dp11_CSTATE_o[1] += -0.00075311499999405624 *
    AHV_Model_X.Dp11_CSTATE_o[1];
  _rtXdot->Dp11_CSTATE_o[1] += -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_o
    [2];
  _rtXdot->Dp11_CSTATE_o[1] += -0.023103995048734959 *
    AHV_Model_X.Dp11_CSTATE_o[3];
  _rtXdot->Dp11_CSTATE_o[1] += -0.0023400735275004216 *
    AHV_Model_X.Dp11_CSTATE_o[4];
  _rtXdot->Dp11_CSTATE_o[2] += -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_o
    [0];
  _rtXdot->Dp11_CSTATE_o[2] += 0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_o
    [1];
  _rtXdot->Dp11_CSTATE_o[2] += -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_o
    [2];
  _rtXdot->Dp11_CSTATE_o[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_o[3];
  _rtXdot->Dp11_CSTATE_o[2] += -0.089746811088073059 *
    AHV_Model_X.Dp11_CSTATE_o[4];
  _rtXdot->Dp11_CSTATE_o[3] += 1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_o[0];
  _rtXdot->Dp11_CSTATE_o[3] += -0.023103995048766656 *
    AHV_Model_X.Dp11_CSTATE_o[1];
  _rtXdot->Dp11_CSTATE_o[3] += 1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_o[2];
  _rtXdot->Dp11_CSTATE_o[3] += -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_o[3];
  _rtXdot->Dp11_CSTATE_o[3] += -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_o[4];
  _rtXdot->Dp11_CSTATE_o[4] += -0.095423424399719625 *
    AHV_Model_X.Dp11_CSTATE_o[0];
  _rtXdot->Dp11_CSTATE_o[4] += 0.0023400735275035814 *
    AHV_Model_X.Dp11_CSTATE_o[1];
  _rtXdot->Dp11_CSTATE_o[4] += -0.089746811088075445 *
    AHV_Model_X.Dp11_CSTATE_o[2];
  _rtXdot->Dp11_CSTATE_o[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_o[3];
  _rtXdot->Dp11_CSTATE_o[4] += -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_o
    [4];
  _rtXdot->Dp11_CSTATE_o[0] += -3.3881001186476949 * AHV_Model_B.nu_r_j[0];
  _rtXdot->Dp11_CSTATE_o[1] += 0.07669611088588818 * AHV_Model_B.nu_r_j[0];
  _rtXdot->Dp11_CSTATE_o[2] += -0.46181267830731043 * AHV_Model_B.nu_r_j[0];
  _rtXdot->Dp11_CSTATE_o[3] += 1.226085627712131 * AHV_Model_B.nu_r_j[0];
  _rtXdot->Dp11_CSTATE_o[4] += -0.11754627904442222 * AHV_Model_B.nu_r_j[0];

  // Derivatives for StateSpace: '<S118>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_jf[0] += -1.3928141301132206 *
    AHV_Model_X.Dp22_CSTATE_jf[0];
  _rtXdot->Dp22_CSTATE_jf[0] += -1.2054566685310828 *
    AHV_Model_X.Dp22_CSTATE_jf[1];
  _rtXdot->Dp22_CSTATE_jf[0] += -0.33483039198973219 *
    AHV_Model_X.Dp22_CSTATE_jf[2];
  _rtXdot->Dp22_CSTATE_jf[0] += -1.0500590385158863 *
    AHV_Model_X.Dp22_CSTATE_jf[3];
  _rtXdot->Dp22_CSTATE_jf[0] += -0.095423424399736376 *
    AHV_Model_X.Dp22_CSTATE_jf[4];
  _rtXdot->Dp22_CSTATE_jf[1] += 1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_jf
    [0];
  _rtXdot->Dp22_CSTATE_jf[1] += -0.00075311499999501316 *
    AHV_Model_X.Dp22_CSTATE_jf[1];
  _rtXdot->Dp22_CSTATE_jf[1] += -0.010562919761942575 *
    AHV_Model_X.Dp22_CSTATE_jf[2];
  _rtXdot->Dp22_CSTATE_jf[1] += -0.0231039950487718 *
    AHV_Model_X.Dp22_CSTATE_jf[3];
  _rtXdot->Dp22_CSTATE_jf[1] += -0.002340073527504722 *
    AHV_Model_X.Dp22_CSTATE_jf[4];
  _rtXdot->Dp22_CSTATE_jf[2] += -0.33483039198972253 *
    AHV_Model_X.Dp22_CSTATE_jf[0];
  _rtXdot->Dp22_CSTATE_jf[2] += 0.010562919761928094 *
    AHV_Model_X.Dp22_CSTATE_jf[1];
  _rtXdot->Dp22_CSTATE_jf[2] += -0.19313064085990478 *
    AHV_Model_X.Dp22_CSTATE_jf[2];
  _rtXdot->Dp22_CSTATE_jf[2] += -1.4618515443989919 *
    AHV_Model_X.Dp22_CSTATE_jf[3];
  _rtXdot->Dp22_CSTATE_jf[2] += -0.089746811088091308 *
    AHV_Model_X.Dp22_CSTATE_jf[4];
  _rtXdot->Dp22_CSTATE_jf[3] += 1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_jf
    [0];
  _rtXdot->Dp22_CSTATE_jf[3] += -0.023103995048724405 *
    AHV_Model_X.Dp22_CSTATE_jf[1];
  _rtXdot->Dp22_CSTATE_jf[3] += 1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_jf
    [2];
  _rtXdot->Dp22_CSTATE_jf[3] += -4.5606728725966441 *
    AHV_Model_X.Dp22_CSTATE_jf[3];
  _rtXdot->Dp22_CSTATE_jf[3] += -1.2857667586959629 *
    AHV_Model_X.Dp22_CSTATE_jf[4];
  _rtXdot->Dp22_CSTATE_jf[4] += -0.095423424399735959 *
    AHV_Model_X.Dp22_CSTATE_jf[0];
  _rtXdot->Dp22_CSTATE_jf[4] += 0.00234007352750089 *
    AHV_Model_X.Dp22_CSTATE_jf[1];
  _rtXdot->Dp22_CSTATE_jf[4] += -0.0897468110880912 *
    AHV_Model_X.Dp22_CSTATE_jf[2];
  _rtXdot->Dp22_CSTATE_jf[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_jf[3];
  _rtXdot->Dp22_CSTATE_jf[4] += -0.13104378720667384 *
    AHV_Model_X.Dp22_CSTATE_jf[4];
  _rtXdot->Dp22_CSTATE_jf[0] += -3.3881001186476967 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp22_CSTATE_jf[1] += 0.076696110885761712 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp22_CSTATE_jf[2] += -0.46181267830733025 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp22_CSTATE_jf[3] += 1.2260856277121315 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp22_CSTATE_jf[4] += -0.11754627904444465 * AHV_Model_B.nu_r_j[1];

  // Derivatives for StateSpace: '<S118>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_j[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp24_CSTATE_j[0];
  _rtXdot->Dp24_CSTATE_j[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_j[1];
  _rtXdot->Dp24_CSTATE_j[0] += -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_j
    [2];
  _rtXdot->Dp24_CSTATE_j[0] += -0.074692039015985229 *
    AHV_Model_X.Dp24_CSTATE_j[3];
  _rtXdot->Dp24_CSTATE_j[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp24_CSTATE_j[4];
  _rtXdot->Dp24_CSTATE_j[1] += -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_j[0];
  _rtXdot->Dp24_CSTATE_j[1] += -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_j[1];
  _rtXdot->Dp24_CSTATE_j[1] += 0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_j[2];
  _rtXdot->Dp24_CSTATE_j[1] += 1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_j[3];
  _rtXdot->Dp24_CSTATE_j[1] += 0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_j
    [4];
  _rtXdot->Dp24_CSTATE_j[2] += 0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_j
    [0];
  _rtXdot->Dp24_CSTATE_j[2] += 0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_j[1];
  _rtXdot->Dp24_CSTATE_j[2] += -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_j[2];
  _rtXdot->Dp24_CSTATE_j[2] += -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_j[3];
  _rtXdot->Dp24_CSTATE_j[2] += -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_j
    [4];
  _rtXdot->Dp24_CSTATE_j[3] += -0.074692039015982259 *
    AHV_Model_X.Dp24_CSTATE_j[0];
  _rtXdot->Dp24_CSTATE_j[3] += -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_j[1];
  _rtXdot->Dp24_CSTATE_j[3] += 3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_j[2];
  _rtXdot->Dp24_CSTATE_j[3] += -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_j[3];
  _rtXdot->Dp24_CSTATE_j[3] += -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_j[4];
  _rtXdot->Dp24_CSTATE_j[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp24_CSTATE_j[0];
  _rtXdot->Dp24_CSTATE_j[4] += 0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_j
    [1];
  _rtXdot->Dp24_CSTATE_j[4] += -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_j
    [2];
  _rtXdot->Dp24_CSTATE_j[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_j[3];
  _rtXdot->Dp24_CSTATE_j[4] += -0.020742386481004484 *
    AHV_Model_X.Dp24_CSTATE_j[4];
  _rtXdot->Dp24_CSTATE_j[0] += -0.23872786278308805 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp24_CSTATE_j[1] += -3.2763464234276056 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp24_CSTATE_j[2] += 1.1437118751635387 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp24_CSTATE_j[3] += -1.3244904674438165 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp24_CSTATE_j[4] += 0.071163531145327891 * AHV_Model_B.nu_r_j[3];

  // Derivatives for StateSpace: '<S118>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_d[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp26_CSTATE_d[0];
  _rtXdot->Dp26_CSTATE_d[0] += 1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_d[1];
  _rtXdot->Dp26_CSTATE_d[0] += 0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_d
    [2];
  _rtXdot->Dp26_CSTATE_d[0] += 0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_d
    [3];
  _rtXdot->Dp26_CSTATE_d[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp26_CSTATE_d[4];
  _rtXdot->Dp26_CSTATE_d[1] += -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_d[0];
  _rtXdot->Dp26_CSTATE_d[1] += -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_d[1];
  _rtXdot->Dp26_CSTATE_d[1] += -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_d
    [2];
  _rtXdot->Dp26_CSTATE_d[1] += -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_d
    [3];
  _rtXdot->Dp26_CSTATE_d[1] += 0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_d[4];
  _rtXdot->Dp26_CSTATE_d[2] += -0.012808448456765344 *
    AHV_Model_X.Dp26_CSTATE_d[0];
  _rtXdot->Dp26_CSTATE_d[2] += -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_d
    [1];
  _rtXdot->Dp26_CSTATE_d[2] += -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_d
    [2];
  _rtXdot->Dp26_CSTATE_d[2] += -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_d[3];
  _rtXdot->Dp26_CSTATE_d[2] += 0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_d[4];
  _rtXdot->Dp26_CSTATE_d[3] += 0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_d
    [0];
  _rtXdot->Dp26_CSTATE_d[3] += 0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_d[1];
  _rtXdot->Dp26_CSTATE_d[3] += 1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_d[2];
  _rtXdot->Dp26_CSTATE_d[3] += -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_d[3];
  _rtXdot->Dp26_CSTATE_d[3] += 1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_d[4];
  _rtXdot->Dp26_CSTATE_d[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp26_CSTATE_d[0];
  _rtXdot->Dp26_CSTATE_d[4] += 0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_d[1];
  _rtXdot->Dp26_CSTATE_d[4] += 0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_d[2];
  _rtXdot->Dp26_CSTATE_d[4] += -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_d[3];
  _rtXdot->Dp26_CSTATE_d[4] += -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_d
    [4];
  _rtXdot->Dp26_CSTATE_d[0] += 0.13167647316600378 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp26_CSTATE_d[1] += 3.4125284827245985 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp26_CSTATE_d[2] += 0.44399732492152644 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp26_CSTATE_d[3] += -1.2820687298457427 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp26_CSTATE_d[4] += -0.18276601298221137 * AHV_Model_B.nu_r_j[5];

  // Derivatives for StateSpace: '<S118>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_g[0] += -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_g[0];
  _rtXdot->Dp33_CSTATE_g[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_g[1];
  _rtXdot->Dp33_CSTATE_g[0] += -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_g
    [2];
  _rtXdot->Dp33_CSTATE_g[0] += -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_g
    [3];
  _rtXdot->Dp33_CSTATE_g[0] += 0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_g
    [4];
  _rtXdot->Dp33_CSTATE_g[1] += -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_g
    [0];
  _rtXdot->Dp33_CSTATE_g[1] += -1.538297528844471E-6 *
    AHV_Model_X.Dp33_CSTATE_g[1];
  _rtXdot->Dp33_CSTATE_g[1] += 0.00010972576839633958 *
    AHV_Model_X.Dp33_CSTATE_g[2];
  _rtXdot->Dp33_CSTATE_g[1] += 0.00050341083172827721 *
    AHV_Model_X.Dp33_CSTATE_g[3];
  _rtXdot->Dp33_CSTATE_g[1] += -3.8193652830025795E-5 *
    AHV_Model_X.Dp33_CSTATE_g[4];
  _rtXdot->Dp33_CSTATE_g[2] += 0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_g[0];
  _rtXdot->Dp33_CSTATE_g[2] += 0.0001097257683747581 *
    AHV_Model_X.Dp33_CSTATE_g[1];
  _rtXdot->Dp33_CSTATE_g[2] += -0.033603468126584782 *
    AHV_Model_X.Dp33_CSTATE_g[2];
  _rtXdot->Dp33_CSTATE_g[2] += -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_g[3];
  _rtXdot->Dp33_CSTATE_g[2] += 0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_g
    [4];
  _rtXdot->Dp33_CSTATE_g[3] += -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_g
    [0];
  _rtXdot->Dp33_CSTATE_g[3] += -0.00050341083159401254 *
    AHV_Model_X.Dp33_CSTATE_g[1];
  _rtXdot->Dp33_CSTATE_g[3] += 1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_g[2];
  _rtXdot->Dp33_CSTATE_g[3] += -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_g
    [3];
  _rtXdot->Dp33_CSTATE_g[3] += 0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_g
    [4];
  _rtXdot->Dp33_CSTATE_g[4] += 0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_g
    [0];
  _rtXdot->Dp33_CSTATE_g[4] += 3.8193652822306167E-5 *
    AHV_Model_X.Dp33_CSTATE_g[1];
  _rtXdot->Dp33_CSTATE_g[4] += -0.039518459830266792 *
    AHV_Model_X.Dp33_CSTATE_g[2];
  _rtXdot->Dp33_CSTATE_g[4] += 0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_g
    [3];
  _rtXdot->Dp33_CSTATE_g[4] += -0.0072896299361453719 *
    AHV_Model_X.Dp33_CSTATE_g[4];
  _rtXdot->Dp33_CSTATE_g[0] += -3.3182263091979736 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp33_CSTATE_g[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp33_CSTATE_g[2] += 0.13279109186599056 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp33_CSTATE_g[3] += -0.54052890870000092 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp33_CSTATE_g[4] += 0.042027661214546957 * AHV_Model_B.nu_r_j[2];

  // Derivatives for StateSpace: '<S118>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_e[0] += -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_e[0];
  _rtXdot->Dp35_CSTATE_e[0] += -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_e
    [1];
  _rtXdot->Dp35_CSTATE_e[0] += 2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_e[2];
  _rtXdot->Dp35_CSTATE_e[0] += 1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_e[3];
  _rtXdot->Dp35_CSTATE_e[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_e[4];
  _rtXdot->Dp35_CSTATE_e[1] += 0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_e[0];
  _rtXdot->Dp35_CSTATE_e[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp35_CSTATE_e[1];
  _rtXdot->Dp35_CSTATE_e[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp35_CSTATE_e[2];
  _rtXdot->Dp35_CSTATE_e[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp35_CSTATE_e[3];
  _rtXdot->Dp35_CSTATE_e[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp35_CSTATE_e[4];
  _rtXdot->Dp35_CSTATE_e[2] += -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_e[0];
  _rtXdot->Dp35_CSTATE_e[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp35_CSTATE_e[1];
  _rtXdot->Dp35_CSTATE_e[2] += -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_e[2];
  _rtXdot->Dp35_CSTATE_e[2] += -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_e[3];
  _rtXdot->Dp35_CSTATE_e[2] += -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_e[4];
  _rtXdot->Dp35_CSTATE_e[3] += -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_e[0];
  _rtXdot->Dp35_CSTATE_e[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp35_CSTATE_e[1];
  _rtXdot->Dp35_CSTATE_e[3] += -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_e[2];
  _rtXdot->Dp35_CSTATE_e[3] += -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_e[3];
  _rtXdot->Dp35_CSTATE_e[3] += -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_e[4];
  _rtXdot->Dp35_CSTATE_e[4] += 0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_e[0];
  _rtXdot->Dp35_CSTATE_e[4] += -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_e
    [1];
  _rtXdot->Dp35_CSTATE_e[4] += 1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_e[2];
  _rtXdot->Dp35_CSTATE_e[4] += 3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_e[3];
  _rtXdot->Dp35_CSTATE_e[4] += -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_e
    [4];
  _rtXdot->Dp35_CSTATE_e[0] += -3.5374929885435322 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp35_CSTATE_e[1] += 0.015233386783081164 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp35_CSTATE_e[2] += -1.2196843916515756 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp35_CSTATE_e[3] += -0.96955326725759594 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp35_CSTATE_e[4] += 0.17148246208757431 * AHV_Model_B.nu_r_j[4];

  // Derivatives for StateSpace: '<S118>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_g[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp42_CSTATE_g[0];
  _rtXdot->Dp42_CSTATE_g[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_g[1];
  _rtXdot->Dp42_CSTATE_g[0] += -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_g
    [2];
  _rtXdot->Dp42_CSTATE_g[0] += -0.074692039015985229 *
    AHV_Model_X.Dp42_CSTATE_g[3];
  _rtXdot->Dp42_CSTATE_g[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp42_CSTATE_g[4];
  _rtXdot->Dp42_CSTATE_g[1] += -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_g[0];
  _rtXdot->Dp42_CSTATE_g[1] += -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_g[1];
  _rtXdot->Dp42_CSTATE_g[1] += 0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_g[2];
  _rtXdot->Dp42_CSTATE_g[1] += 1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_g[3];
  _rtXdot->Dp42_CSTATE_g[1] += 0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_g
    [4];
  _rtXdot->Dp42_CSTATE_g[2] += 0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_g
    [0];
  _rtXdot->Dp42_CSTATE_g[2] += 0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_g[1];
  _rtXdot->Dp42_CSTATE_g[2] += -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_g[2];
  _rtXdot->Dp42_CSTATE_g[2] += -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_g[3];
  _rtXdot->Dp42_CSTATE_g[2] += -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_g
    [4];
  _rtXdot->Dp42_CSTATE_g[3] += -0.074692039015982259 *
    AHV_Model_X.Dp42_CSTATE_g[0];
  _rtXdot->Dp42_CSTATE_g[3] += -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_g[1];
  _rtXdot->Dp42_CSTATE_g[3] += 3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_g[2];
  _rtXdot->Dp42_CSTATE_g[3] += -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_g[3];
  _rtXdot->Dp42_CSTATE_g[3] += -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_g[4];
  _rtXdot->Dp42_CSTATE_g[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp42_CSTATE_g[0];
  _rtXdot->Dp42_CSTATE_g[4] += 0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_g
    [1];
  _rtXdot->Dp42_CSTATE_g[4] += -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_g
    [2];
  _rtXdot->Dp42_CSTATE_g[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_g[3];
  _rtXdot->Dp42_CSTATE_g[4] += -0.020742386481004484 *
    AHV_Model_X.Dp42_CSTATE_g[4];
  _rtXdot->Dp42_CSTATE_g[0] += -0.23872786278308805 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp42_CSTATE_g[1] += -3.2763464234276056 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp42_CSTATE_g[2] += 1.1437118751635387 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp42_CSTATE_g[3] += -1.3244904674438165 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp42_CSTATE_g[4] += 0.071163531145327891 * AHV_Model_B.nu_r_j[1];

  // Derivatives for StateSpace: '<S118>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_e[0] += -0.004396104715914141 *
    AHV_Model_X.Dp46_CSTATE_e[0];
  _rtXdot->Dp46_CSTATE_e[0] += 1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_e[1];
  _rtXdot->Dp46_CSTATE_e[0] += -0.027511960605097124 *
    AHV_Model_X.Dp46_CSTATE_e[2];
  _rtXdot->Dp46_CSTATE_e[0] += 0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_e
    [3];
  _rtXdot->Dp46_CSTATE_e[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp46_CSTATE_e[4];
  _rtXdot->Dp46_CSTATE_e[1] += -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_e[0];
  _rtXdot->Dp46_CSTATE_e[1] += -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_e[1];
  _rtXdot->Dp46_CSTATE_e[1] += 0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_e[2];
  _rtXdot->Dp46_CSTATE_e[1] += -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_e[3];
  _rtXdot->Dp46_CSTATE_e[1] += 0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_e
    [4];
  _rtXdot->Dp46_CSTATE_e[2] += 0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_e[0];
  _rtXdot->Dp46_CSTATE_e[2] += 0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_e[1];
  _rtXdot->Dp46_CSTATE_e[2] += -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_e
    [2];
  _rtXdot->Dp46_CSTATE_e[2] += 2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_e[3];
  _rtXdot->Dp46_CSTATE_e[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp46_CSTATE_e[4];
  _rtXdot->Dp46_CSTATE_e[3] += 0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_e[0];
  _rtXdot->Dp46_CSTATE_e[3] += 1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_e[1];
  _rtXdot->Dp46_CSTATE_e[3] += -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_e[2];
  _rtXdot->Dp46_CSTATE_e[3] += -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_e[3];
  _rtXdot->Dp46_CSTATE_e[3] += 0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_e[4];
  _rtXdot->Dp46_CSTATE_e[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp46_CSTATE_e[0];
  _rtXdot->Dp46_CSTATE_e[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp46_CSTATE_e[1];
  _rtXdot->Dp46_CSTATE_e[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp46_CSTATE_e[2];
  _rtXdot->Dp46_CSTATE_e[4] += -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_e
    [3];
  _rtXdot->Dp46_CSTATE_e[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp46_CSTATE_e[4];
  _rtXdot->Dp46_CSTATE_e[0] += -0.33524288945878539 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp46_CSTATE_e[1] += -6.3625492730866693 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp46_CSTATE_e[2] += 0.93578679415036736 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp46_CSTATE_e[3] += 2.5673752417819116 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp46_CSTATE_e[4] += 0.014933981938330865 * AHV_Model_B.nu_r_j[5];

  // Derivatives for StateSpace: '<S118>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_k[0] += -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_k[0];
  _rtXdot->Dp53_CSTATE_k[0] += -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_k
    [1];
  _rtXdot->Dp53_CSTATE_k[0] += 2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_k[2];
  _rtXdot->Dp53_CSTATE_k[0] += 1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_k[3];
  _rtXdot->Dp53_CSTATE_k[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE_k[4];
  _rtXdot->Dp53_CSTATE_k[1] += 0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_k[0];
  _rtXdot->Dp53_CSTATE_k[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp53_CSTATE_k[1];
  _rtXdot->Dp53_CSTATE_k[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp53_CSTATE_k[2];
  _rtXdot->Dp53_CSTATE_k[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp53_CSTATE_k[3];
  _rtXdot->Dp53_CSTATE_k[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp53_CSTATE_k[4];
  _rtXdot->Dp53_CSTATE_k[2] += -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_k[0];
  _rtXdot->Dp53_CSTATE_k[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp53_CSTATE_k[1];
  _rtXdot->Dp53_CSTATE_k[2] += -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_k[2];
  _rtXdot->Dp53_CSTATE_k[2] += -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_k[3];
  _rtXdot->Dp53_CSTATE_k[2] += -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_k[4];
  _rtXdot->Dp53_CSTATE_k[3] += -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_k[0];
  _rtXdot->Dp53_CSTATE_k[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp53_CSTATE_k[1];
  _rtXdot->Dp53_CSTATE_k[3] += -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_k[2];
  _rtXdot->Dp53_CSTATE_k[3] += -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_k[3];
  _rtXdot->Dp53_CSTATE_k[3] += -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_k[4];
  _rtXdot->Dp53_CSTATE_k[4] += 0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_k[0];
  _rtXdot->Dp53_CSTATE_k[4] += -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_k
    [1];
  _rtXdot->Dp53_CSTATE_k[4] += 1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_k[2];
  _rtXdot->Dp53_CSTATE_k[4] += 3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_k[3];
  _rtXdot->Dp53_CSTATE_k[4] += -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_k
    [4];
  _rtXdot->Dp53_CSTATE_k[0] += -3.5374929885435322 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp53_CSTATE_k[1] += 0.015233386783081164 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp53_CSTATE_k[2] += -1.2196843916515756 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp53_CSTATE_k[3] += -0.96955326725759594 * AHV_Model_B.nu_r_j[2];
  _rtXdot->Dp53_CSTATE_k[4] += 0.17148246208757431 * AHV_Model_B.nu_r_j[2];

  // Derivatives for StateSpace: '<S118>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_i[0] += -8.7466634083809418E-5 *
    AHV_Model_X.Dp55_CSTATE_i[0];
  _rtXdot->Dp55_CSTATE_i[0] += 0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_i[1];
  _rtXdot->Dp55_CSTATE_i[0] += 0.0013720743421809071 *
    AHV_Model_X.Dp55_CSTATE_i[2];
  _rtXdot->Dp55_CSTATE_i[0] += 0.0075482812368147384 *
    AHV_Model_X.Dp55_CSTATE_i[3];
  _rtXdot->Dp55_CSTATE_i[0] += 0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_i[4];
  _rtXdot->Dp55_CSTATE_i[1] += -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_i
    [0];
  _rtXdot->Dp55_CSTATE_i[1] += -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_i[1];
  _rtXdot->Dp55_CSTATE_i[1] += -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_i
    [2];
  _rtXdot->Dp55_CSTATE_i[1] += -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_i[3];
  _rtXdot->Dp55_CSTATE_i[1] += -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_i[4];
  _rtXdot->Dp55_CSTATE_i[2] += -0.0013720743421854747 *
    AHV_Model_X.Dp55_CSTATE_i[0];
  _rtXdot->Dp55_CSTATE_i[2] += -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_i
    [1];
  _rtXdot->Dp55_CSTATE_i[2] += -0.068186331750790974 *
    AHV_Model_X.Dp55_CSTATE_i[2];
  _rtXdot->Dp55_CSTATE_i[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_i[3];
  _rtXdot->Dp55_CSTATE_i[2] += -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_i[4];
  _rtXdot->Dp55_CSTATE_i[3] += 0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_i
    [0];
  _rtXdot->Dp55_CSTATE_i[3] += 1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_i[1];
  _rtXdot->Dp55_CSTATE_i[3] += 4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_i[2];
  _rtXdot->Dp55_CSTATE_i[3] += -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_i[3];
  _rtXdot->Dp55_CSTATE_i[3] += -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_i[4];
  _rtXdot->Dp55_CSTATE_i[4] += 0.0070959954003254081 *
    AHV_Model_X.Dp55_CSTATE_i[0];
  _rtXdot->Dp55_CSTATE_i[4] += 1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_i[1];
  _rtXdot->Dp55_CSTATE_i[4] += 2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_i[2];
  _rtXdot->Dp55_CSTATE_i[4] += -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_i[3];
  _rtXdot->Dp55_CSTATE_i[4] += -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_i[4];
  _rtXdot->Dp55_CSTATE_i[0] += 0.021904981880023978 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp55_CSTATE_i[1] += 3.4651622640804733 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp55_CSTATE_i[2] += 0.16004490113712252 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp55_CSTATE_i[3] += -0.9969365784863482 * AHV_Model_B.nu_r_j[4];
  _rtXdot->Dp55_CSTATE_i[4] += -0.92774837285087253 * AHV_Model_B.nu_r_j[4];

  // Derivatives for StateSpace: '<S118>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_h[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp62_CSTATE_h[0];
  _rtXdot->Dp62_CSTATE_h[0] += 1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_h[1];
  _rtXdot->Dp62_CSTATE_h[0] += 0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_h
    [2];
  _rtXdot->Dp62_CSTATE_h[0] += 0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_h
    [3];
  _rtXdot->Dp62_CSTATE_h[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp62_CSTATE_h[4];
  _rtXdot->Dp62_CSTATE_h[1] += -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_h[0];
  _rtXdot->Dp62_CSTATE_h[1] += -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_h[1];
  _rtXdot->Dp62_CSTATE_h[1] += -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_h
    [2];
  _rtXdot->Dp62_CSTATE_h[1] += -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_h
    [3];
  _rtXdot->Dp62_CSTATE_h[1] += 0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_h[4];
  _rtXdot->Dp62_CSTATE_h[2] += -0.012808448456765344 *
    AHV_Model_X.Dp62_CSTATE_h[0];
  _rtXdot->Dp62_CSTATE_h[2] += -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_h
    [1];
  _rtXdot->Dp62_CSTATE_h[2] += -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_h
    [2];
  _rtXdot->Dp62_CSTATE_h[2] += -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_h[3];
  _rtXdot->Dp62_CSTATE_h[2] += 0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_h[4];
  _rtXdot->Dp62_CSTATE_h[3] += 0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_h
    [0];
  _rtXdot->Dp62_CSTATE_h[3] += 0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_h[1];
  _rtXdot->Dp62_CSTATE_h[3] += 1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_h[2];
  _rtXdot->Dp62_CSTATE_h[3] += -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_h[3];
  _rtXdot->Dp62_CSTATE_h[3] += 1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_h[4];
  _rtXdot->Dp62_CSTATE_h[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp62_CSTATE_h[0];
  _rtXdot->Dp62_CSTATE_h[4] += 0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_h[1];
  _rtXdot->Dp62_CSTATE_h[4] += 0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_h[2];
  _rtXdot->Dp62_CSTATE_h[4] += -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_h[3];
  _rtXdot->Dp62_CSTATE_h[4] += -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_h
    [4];
  _rtXdot->Dp62_CSTATE_h[0] += 0.13167647316600378 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp62_CSTATE_h[1] += 3.4125284827245985 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp62_CSTATE_h[2] += 0.44399732492152644 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp62_CSTATE_h[3] += -1.2820687298457427 * AHV_Model_B.nu_r_j[1];
  _rtXdot->Dp62_CSTATE_h[4] += -0.18276601298221137 * AHV_Model_B.nu_r_j[1];

  // Derivatives for StateSpace: '<S118>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_h[0] += -0.004396104715914141 *
    AHV_Model_X.Dp64_CSTATE_h[0];
  _rtXdot->Dp64_CSTATE_h[0] += 1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_h[1];
  _rtXdot->Dp64_CSTATE_h[0] += -0.027511960605097124 *
    AHV_Model_X.Dp64_CSTATE_h[2];
  _rtXdot->Dp64_CSTATE_h[0] += 0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_h
    [3];
  _rtXdot->Dp64_CSTATE_h[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp64_CSTATE_h[4];
  _rtXdot->Dp64_CSTATE_h[1] += -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_h[0];
  _rtXdot->Dp64_CSTATE_h[1] += -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_h[1];
  _rtXdot->Dp64_CSTATE_h[1] += 0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_h[2];
  _rtXdot->Dp64_CSTATE_h[1] += -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_h[3];
  _rtXdot->Dp64_CSTATE_h[1] += 0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_h
    [4];
  _rtXdot->Dp64_CSTATE_h[2] += 0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_h[0];
  _rtXdot->Dp64_CSTATE_h[2] += 0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_h[1];
  _rtXdot->Dp64_CSTATE_h[2] += -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_h
    [2];
  _rtXdot->Dp64_CSTATE_h[2] += 2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_h[3];
  _rtXdot->Dp64_CSTATE_h[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp64_CSTATE_h[4];
  _rtXdot->Dp64_CSTATE_h[3] += 0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_h[0];
  _rtXdot->Dp64_CSTATE_h[3] += 1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_h[1];
  _rtXdot->Dp64_CSTATE_h[3] += -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_h[2];
  _rtXdot->Dp64_CSTATE_h[3] += -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_h[3];
  _rtXdot->Dp64_CSTATE_h[3] += 0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_h[4];
  _rtXdot->Dp64_CSTATE_h[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp64_CSTATE_h[0];
  _rtXdot->Dp64_CSTATE_h[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp64_CSTATE_h[1];
  _rtXdot->Dp64_CSTATE_h[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp64_CSTATE_h[2];
  _rtXdot->Dp64_CSTATE_h[4] += -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_h
    [3];
  _rtXdot->Dp64_CSTATE_h[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp64_CSTATE_h[4];
  _rtXdot->Dp64_CSTATE_h[0] += -0.33524288945878539 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp64_CSTATE_h[1] += -6.3625492730866693 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp64_CSTATE_h[2] += 0.93578679415036736 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp64_CSTATE_h[3] += 2.5673752417819116 * AHV_Model_B.nu_r_j[3];
  _rtXdot->Dp64_CSTATE_h[4] += 0.014933981938330865 * AHV_Model_B.nu_r_j[3];

  // Derivatives for StateSpace: '<S118>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_n[0] += -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[0] += -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[0] += -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE_n
    [2];
  _rtXdot->Dp66_CSTATE_n[0] += -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[0] += 0.069477290559771143 * AHV_Model_X.Dp66_CSTATE_n
    [4];
  _rtXdot->Dp66_CSTATE_n[1] += 1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[1] += -0.00029199832472599611 *
    AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[1] += -0.0080078570280257711 *
    AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[1] += -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE_n
    [3];
  _rtXdot->Dp66_CSTATE_n[1] += 0.00098472174210678556 *
    AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[2] += -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE_n
    [0];
  _rtXdot->Dp66_CSTATE_n[2] += 0.00800785702802341 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[2] += -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE_n
    [2];
  _rtXdot->Dp66_CSTATE_n[2] += -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[2] += 0.074587161707531782 * AHV_Model_X.Dp66_CSTATE_n
    [4];
  _rtXdot->Dp66_CSTATE_n[3] += 1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[3] += -0.015684008067190565 *
    AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[3] += 1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[3] += -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[4] += 0.069477290559780386 * AHV_Model_X.Dp66_CSTATE_n
    [0];
  _rtXdot->Dp66_CSTATE_n[4] += -0.000984721742106249 *
    AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[4] += 0.074587161707547089 * AHV_Model_X.Dp66_CSTATE_n
    [2];
  _rtXdot->Dp66_CSTATE_n[4] += -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[4] += -0.061420217921035136 *
    AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[0] += -3.3844821061206916 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp66_CSTATE_n[1] += 0.045246867266885989 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp66_CSTATE_n[2] += -0.53332191483729374 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp66_CSTATE_n[3] += 1.2581308985447515 * AHV_Model_B.nu_r_j[5];
  _rtXdot->Dp66_CSTATE_n[4] += 0.075289942321940723 * AHV_Model_B.nu_r_j[5];

  // Derivatives for Integrator: '<S15>/Integrator'
  for (i = 0; i < 6; i++) {
    _rtXdot->Integrator_CSTATE_e[i] = AHV_Model_B.Minvtau_o[i];
  }

  // End of Derivatives for Integrator: '<S15>/Integrator'

  // Derivatives for Integrator: '<S65>/Integrator3'
  _rtXdot->Integrator3_CSTATE_n[0] = AHV_Model_B.sun_k2[0];

  // Derivatives for Integrator: '<S65>/Integrator4'
  _rtXdot->Integrator4_CSTATE_c[0] = AHV_Model_B.M_u[0];

  // Derivatives for Integrator: '<S65>/Integrator3'
  _rtXdot->Integrator3_CSTATE_n[1] = AHV_Model_B.sun_k2[1];

  // Derivatives for Integrator: '<S65>/Integrator4'
  _rtXdot->Integrator4_CSTATE_c[1] = AHV_Model_B.M_u[1];

  // Derivatives for Integrator: '<S65>/Integrator3'
  _rtXdot->Integrator3_CSTATE_n[2] = AHV_Model_B.sun_k2[2];

  // Derivatives for Integrator: '<S65>/Integrator4'
  _rtXdot->Integrator4_CSTATE_c[2] = AHV_Model_B.M_u[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S22>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_d[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_d[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_m[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_c[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_b[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_o[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_e[i] = 0.0;

    // Derivatives for Integrator: '<S34>/Integrator'
    _rtXdot->Integrator_CSTATE_p[i] = AHV_Model_B.Sum[i];

    // Derivatives for StateSpace: '<S22>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_g[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_m[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_n[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_n[i] = 0.0;

    // Derivatives for StateSpace: '<S22>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_fj[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S22>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_d[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_d[0];
  _rtXdot->Dp11_CSTATE_d[0] += -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_d[1];
  _rtXdot->Dp11_CSTATE_d[0] += -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_d
    [2];
  _rtXdot->Dp11_CSTATE_d[0] += -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_d[3];
  _rtXdot->Dp11_CSTATE_d[0] += -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_d[4];
  _rtXdot->Dp11_CSTATE_d[1] += 1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_d[0];
  _rtXdot->Dp11_CSTATE_d[1] += -0.00075311499999405624 *
    AHV_Model_X.Dp11_CSTATE_d[1];
  _rtXdot->Dp11_CSTATE_d[1] += -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_d
    [2];
  _rtXdot->Dp11_CSTATE_d[1] += -0.023103995048734959 *
    AHV_Model_X.Dp11_CSTATE_d[3];
  _rtXdot->Dp11_CSTATE_d[1] += -0.0023400735275004216 *
    AHV_Model_X.Dp11_CSTATE_d[4];
  _rtXdot->Dp11_CSTATE_d[2] += -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_d
    [0];
  _rtXdot->Dp11_CSTATE_d[2] += 0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_d
    [1];
  _rtXdot->Dp11_CSTATE_d[2] += -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_d
    [2];
  _rtXdot->Dp11_CSTATE_d[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_d[3];
  _rtXdot->Dp11_CSTATE_d[2] += -0.089746811088073059 *
    AHV_Model_X.Dp11_CSTATE_d[4];
  _rtXdot->Dp11_CSTATE_d[3] += 1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_d[0];
  _rtXdot->Dp11_CSTATE_d[3] += -0.023103995048766656 *
    AHV_Model_X.Dp11_CSTATE_d[1];
  _rtXdot->Dp11_CSTATE_d[3] += 1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_d[2];
  _rtXdot->Dp11_CSTATE_d[3] += -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_d[3];
  _rtXdot->Dp11_CSTATE_d[3] += -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_d[4];
  _rtXdot->Dp11_CSTATE_d[4] += -0.095423424399719625 *
    AHV_Model_X.Dp11_CSTATE_d[0];
  _rtXdot->Dp11_CSTATE_d[4] += 0.0023400735275035814 *
    AHV_Model_X.Dp11_CSTATE_d[1];
  _rtXdot->Dp11_CSTATE_d[4] += -0.089746811088075445 *
    AHV_Model_X.Dp11_CSTATE_d[2];
  _rtXdot->Dp11_CSTATE_d[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_d[3];
  _rtXdot->Dp11_CSTATE_d[4] += -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_d
    [4];
  _rtXdot->Dp11_CSTATE_d[0] += -3.3881001186476949 * AHV_Model_B.nu_r_c[0];
  _rtXdot->Dp11_CSTATE_d[1] += 0.07669611088588818 * AHV_Model_B.nu_r_c[0];
  _rtXdot->Dp11_CSTATE_d[2] += -0.46181267830731043 * AHV_Model_B.nu_r_c[0];
  _rtXdot->Dp11_CSTATE_d[3] += 1.226085627712131 * AHV_Model_B.nu_r_c[0];
  _rtXdot->Dp11_CSTATE_d[4] += -0.11754627904442222 * AHV_Model_B.nu_r_c[0];

  // Derivatives for StateSpace: '<S22>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_h[0] += -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE_h[0];
  _rtXdot->Dp22_CSTATE_h[0] += -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE_h[1];
  _rtXdot->Dp22_CSTATE_h[0] += -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE_h
    [2];
  _rtXdot->Dp22_CSTATE_h[0] += -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE_h[3];
  _rtXdot->Dp22_CSTATE_h[0] += -0.095423424399736376 *
    AHV_Model_X.Dp22_CSTATE_h[4];
  _rtXdot->Dp22_CSTATE_h[1] += 1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_h[0];
  _rtXdot->Dp22_CSTATE_h[1] += -0.00075311499999501316 *
    AHV_Model_X.Dp22_CSTATE_h[1];
  _rtXdot->Dp22_CSTATE_h[1] += -0.010562919761942575 *
    AHV_Model_X.Dp22_CSTATE_h[2];
  _rtXdot->Dp22_CSTATE_h[1] += -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE_h[3];
  _rtXdot->Dp22_CSTATE_h[1] += -0.002340073527504722 *
    AHV_Model_X.Dp22_CSTATE_h[4];
  _rtXdot->Dp22_CSTATE_h[2] += -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE_h
    [0];
  _rtXdot->Dp22_CSTATE_h[2] += 0.010562919761928094 * AHV_Model_X.Dp22_CSTATE_h
    [1];
  _rtXdot->Dp22_CSTATE_h[2] += -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE_h
    [2];
  _rtXdot->Dp22_CSTATE_h[2] += -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE_h[3];
  _rtXdot->Dp22_CSTATE_h[2] += -0.089746811088091308 *
    AHV_Model_X.Dp22_CSTATE_h[4];
  _rtXdot->Dp22_CSTATE_h[3] += 1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_h[0];
  _rtXdot->Dp22_CSTATE_h[3] += -0.023103995048724405 *
    AHV_Model_X.Dp22_CSTATE_h[1];
  _rtXdot->Dp22_CSTATE_h[3] += 1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_h[2];
  _rtXdot->Dp22_CSTATE_h[3] += -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE_h[3];
  _rtXdot->Dp22_CSTATE_h[3] += -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE_h[4];
  _rtXdot->Dp22_CSTATE_h[4] += -0.095423424399735959 *
    AHV_Model_X.Dp22_CSTATE_h[0];
  _rtXdot->Dp22_CSTATE_h[4] += 0.00234007352750089 * AHV_Model_X.Dp22_CSTATE_h[1];
  _rtXdot->Dp22_CSTATE_h[4] += -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE_h[2];
  _rtXdot->Dp22_CSTATE_h[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_h[3];
  _rtXdot->Dp22_CSTATE_h[4] += -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE_h
    [4];
  _rtXdot->Dp22_CSTATE_h[0] += -3.3881001186476967 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp22_CSTATE_h[1] += 0.076696110885761712 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp22_CSTATE_h[2] += -0.46181267830733025 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp22_CSTATE_h[3] += 1.2260856277121315 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp22_CSTATE_h[4] += -0.11754627904444465 * AHV_Model_B.nu_r_c[1];

  // Derivatives for StateSpace: '<S22>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_c[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp24_CSTATE_c[0];
  _rtXdot->Dp24_CSTATE_c[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_c[1];
  _rtXdot->Dp24_CSTATE_c[0] += -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_c
    [2];
  _rtXdot->Dp24_CSTATE_c[0] += -0.074692039015985229 *
    AHV_Model_X.Dp24_CSTATE_c[3];
  _rtXdot->Dp24_CSTATE_c[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp24_CSTATE_c[4];
  _rtXdot->Dp24_CSTATE_c[1] += -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_c[0];
  _rtXdot->Dp24_CSTATE_c[1] += -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_c[1];
  _rtXdot->Dp24_CSTATE_c[1] += 0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_c[2];
  _rtXdot->Dp24_CSTATE_c[1] += 1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_c[3];
  _rtXdot->Dp24_CSTATE_c[1] += 0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_c
    [4];
  _rtXdot->Dp24_CSTATE_c[2] += 0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_c
    [0];
  _rtXdot->Dp24_CSTATE_c[2] += 0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_c[1];
  _rtXdot->Dp24_CSTATE_c[2] += -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_c[2];
  _rtXdot->Dp24_CSTATE_c[2] += -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_c[3];
  _rtXdot->Dp24_CSTATE_c[2] += -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_c
    [4];
  _rtXdot->Dp24_CSTATE_c[3] += -0.074692039015982259 *
    AHV_Model_X.Dp24_CSTATE_c[0];
  _rtXdot->Dp24_CSTATE_c[3] += -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_c[1];
  _rtXdot->Dp24_CSTATE_c[3] += 3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_c[2];
  _rtXdot->Dp24_CSTATE_c[3] += -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_c[3];
  _rtXdot->Dp24_CSTATE_c[3] += -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_c[4];
  _rtXdot->Dp24_CSTATE_c[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp24_CSTATE_c[0];
  _rtXdot->Dp24_CSTATE_c[4] += 0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_c
    [1];
  _rtXdot->Dp24_CSTATE_c[4] += -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_c
    [2];
  _rtXdot->Dp24_CSTATE_c[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_c[3];
  _rtXdot->Dp24_CSTATE_c[4] += -0.020742386481004484 *
    AHV_Model_X.Dp24_CSTATE_c[4];
  _rtXdot->Dp24_CSTATE_c[0] += -0.23872786278308805 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp24_CSTATE_c[1] += -3.2763464234276056 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp24_CSTATE_c[2] += 1.1437118751635387 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp24_CSTATE_c[3] += -1.3244904674438165 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp24_CSTATE_c[4] += 0.071163531145327891 * AHV_Model_B.nu_r_c[3];

  // Derivatives for StateSpace: '<S22>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_h[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp26_CSTATE_h[0];
  _rtXdot->Dp26_CSTATE_h[0] += 1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_h[1];
  _rtXdot->Dp26_CSTATE_h[0] += 0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_h
    [2];
  _rtXdot->Dp26_CSTATE_h[0] += 0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_h
    [3];
  _rtXdot->Dp26_CSTATE_h[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp26_CSTATE_h[4];
  _rtXdot->Dp26_CSTATE_h[1] += -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_h[0];
  _rtXdot->Dp26_CSTATE_h[1] += -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_h[1];
  _rtXdot->Dp26_CSTATE_h[1] += -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_h
    [2];
  _rtXdot->Dp26_CSTATE_h[1] += -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_h
    [3];
  _rtXdot->Dp26_CSTATE_h[1] += 0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_h[4];
  _rtXdot->Dp26_CSTATE_h[2] += -0.012808448456765344 *
    AHV_Model_X.Dp26_CSTATE_h[0];
  _rtXdot->Dp26_CSTATE_h[2] += -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_h
    [1];
  _rtXdot->Dp26_CSTATE_h[2] += -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_h
    [2];
  _rtXdot->Dp26_CSTATE_h[2] += -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_h[3];
  _rtXdot->Dp26_CSTATE_h[2] += 0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_h[4];
  _rtXdot->Dp26_CSTATE_h[3] += 0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_h
    [0];
  _rtXdot->Dp26_CSTATE_h[3] += 0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_h[1];
  _rtXdot->Dp26_CSTATE_h[3] += 1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_h[2];
  _rtXdot->Dp26_CSTATE_h[3] += -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_h[3];
  _rtXdot->Dp26_CSTATE_h[3] += 1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_h[4];
  _rtXdot->Dp26_CSTATE_h[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp26_CSTATE_h[0];
  _rtXdot->Dp26_CSTATE_h[4] += 0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_h[1];
  _rtXdot->Dp26_CSTATE_h[4] += 0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_h[2];
  _rtXdot->Dp26_CSTATE_h[4] += -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_h[3];
  _rtXdot->Dp26_CSTATE_h[4] += -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_h
    [4];
  _rtXdot->Dp26_CSTATE_h[0] += 0.13167647316600378 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp26_CSTATE_h[1] += 3.4125284827245985 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp26_CSTATE_h[2] += 0.44399732492152644 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp26_CSTATE_h[3] += -1.2820687298457427 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp26_CSTATE_h[4] += -0.18276601298221137 * AHV_Model_B.nu_r_c[5];

  // Derivatives for StateSpace: '<S22>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_o[0] += -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[0] += -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_o
    [2];
  _rtXdot->Dp33_CSTATE_o[0] += -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_o
    [3];
  _rtXdot->Dp33_CSTATE_o[0] += 0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_o
    [4];
  _rtXdot->Dp33_CSTATE_o[1] += -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_o
    [0];
  _rtXdot->Dp33_CSTATE_o[1] += -1.538297528844471E-6 *
    AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[1] += 0.00010972576839633958 *
    AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[1] += 0.00050341083172827721 *
    AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[1] += -3.8193652830025795E-5 *
    AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[2] += 0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[2] += 0.0001097257683747581 *
    AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[2] += -0.033603468126584782 *
    AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[2] += -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[2] += 0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_o
    [4];
  _rtXdot->Dp33_CSTATE_o[3] += -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_o
    [0];
  _rtXdot->Dp33_CSTATE_o[3] += -0.00050341083159401254 *
    AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[3] += 1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[3] += -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_o
    [3];
  _rtXdot->Dp33_CSTATE_o[3] += 0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_o
    [4];
  _rtXdot->Dp33_CSTATE_o[4] += 0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_o
    [0];
  _rtXdot->Dp33_CSTATE_o[4] += 3.8193652822306167E-5 *
    AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[4] += -0.039518459830266792 *
    AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[4] += 0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_o
    [3];
  _rtXdot->Dp33_CSTATE_o[4] += -0.0072896299361453719 *
    AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[0] += -3.3182263091979736 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp33_CSTATE_o[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp33_CSTATE_o[2] += 0.13279109186599056 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp33_CSTATE_o[3] += -0.54052890870000092 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp33_CSTATE_o[4] += 0.042027661214546957 * AHV_Model_B.nu_r_c[2];

  // Derivatives for StateSpace: '<S22>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_h[0] += -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_h[0];
  _rtXdot->Dp35_CSTATE_h[0] += -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_h
    [1];
  _rtXdot->Dp35_CSTATE_h[0] += 2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_h[2];
  _rtXdot->Dp35_CSTATE_h[0] += 1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_h[3];
  _rtXdot->Dp35_CSTATE_h[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_h[4];
  _rtXdot->Dp35_CSTATE_h[1] += 0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_h[0];
  _rtXdot->Dp35_CSTATE_h[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp35_CSTATE_h[1];
  _rtXdot->Dp35_CSTATE_h[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp35_CSTATE_h[2];
  _rtXdot->Dp35_CSTATE_h[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp35_CSTATE_h[3];
  _rtXdot->Dp35_CSTATE_h[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp35_CSTATE_h[4];
  _rtXdot->Dp35_CSTATE_h[2] += -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_h[0];
  _rtXdot->Dp35_CSTATE_h[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp35_CSTATE_h[1];
  _rtXdot->Dp35_CSTATE_h[2] += -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_h[2];
  _rtXdot->Dp35_CSTATE_h[2] += -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_h[3];
  _rtXdot->Dp35_CSTATE_h[2] += -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_h[4];
  _rtXdot->Dp35_CSTATE_h[3] += -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_h[0];
  _rtXdot->Dp35_CSTATE_h[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp35_CSTATE_h[1];
  _rtXdot->Dp35_CSTATE_h[3] += -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_h[2];
  _rtXdot->Dp35_CSTATE_h[3] += -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_h[3];
  _rtXdot->Dp35_CSTATE_h[3] += -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_h[4];
  _rtXdot->Dp35_CSTATE_h[4] += 0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_h[0];
  _rtXdot->Dp35_CSTATE_h[4] += -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_h
    [1];
  _rtXdot->Dp35_CSTATE_h[4] += 1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_h[2];
  _rtXdot->Dp35_CSTATE_h[4] += 3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_h[3];
  _rtXdot->Dp35_CSTATE_h[4] += -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_h
    [4];
  _rtXdot->Dp35_CSTATE_h[0] += -3.5374929885435322 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp35_CSTATE_h[1] += 0.015233386783081164 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp35_CSTATE_h[2] += -1.2196843916515756 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp35_CSTATE_h[3] += -0.96955326725759594 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp35_CSTATE_h[4] += 0.17148246208757431 * AHV_Model_B.nu_r_c[4];

  // Derivatives for StateSpace: '<S22>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_e[0] += -0.0069962390232986621 *
    AHV_Model_X.Dp42_CSTATE_e[0];
  _rtXdot->Dp42_CSTATE_e[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_e[1];
  _rtXdot->Dp42_CSTATE_e[0] += -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_e
    [2];
  _rtXdot->Dp42_CSTATE_e[0] += -0.074692039015985229 *
    AHV_Model_X.Dp42_CSTATE_e[3];
  _rtXdot->Dp42_CSTATE_e[0] += -0.0042999629253483171 *
    AHV_Model_X.Dp42_CSTATE_e[4];
  _rtXdot->Dp42_CSTATE_e[1] += -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_e[0];
  _rtXdot->Dp42_CSTATE_e[1] += -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_e[1];
  _rtXdot->Dp42_CSTATE_e[1] += 0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_e[2];
  _rtXdot->Dp42_CSTATE_e[1] += 1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_e[3];
  _rtXdot->Dp42_CSTATE_e[1] += 0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_e
    [4];
  _rtXdot->Dp42_CSTATE_e[2] += 0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_e
    [0];
  _rtXdot->Dp42_CSTATE_e[2] += 0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_e[1];
  _rtXdot->Dp42_CSTATE_e[2] += -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_e[2];
  _rtXdot->Dp42_CSTATE_e[2] += -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_e[3];
  _rtXdot->Dp42_CSTATE_e[2] += -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_e
    [4];
  _rtXdot->Dp42_CSTATE_e[3] += -0.074692039015982259 *
    AHV_Model_X.Dp42_CSTATE_e[0];
  _rtXdot->Dp42_CSTATE_e[3] += -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_e[1];
  _rtXdot->Dp42_CSTATE_e[3] += 3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_e[2];
  _rtXdot->Dp42_CSTATE_e[3] += -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_e[3];
  _rtXdot->Dp42_CSTATE_e[3] += -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_e[4];
  _rtXdot->Dp42_CSTATE_e[4] += 0.0042999629253487083 *
    AHV_Model_X.Dp42_CSTATE_e[0];
  _rtXdot->Dp42_CSTATE_e[4] += 0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_e
    [1];
  _rtXdot->Dp42_CSTATE_e[4] += -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_e
    [2];
  _rtXdot->Dp42_CSTATE_e[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_e[3];
  _rtXdot->Dp42_CSTATE_e[4] += -0.020742386481004484 *
    AHV_Model_X.Dp42_CSTATE_e[4];
  _rtXdot->Dp42_CSTATE_e[0] += -0.23872786278308805 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp42_CSTATE_e[1] += -3.2763464234276056 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp42_CSTATE_e[2] += 1.1437118751635387 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp42_CSTATE_e[3] += -1.3244904674438165 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp42_CSTATE_e[4] += 0.071163531145327891 * AHV_Model_B.nu_r_c[1];

  // Derivatives for StateSpace: '<S22>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_h[0] += -0.004396104715914141 *
    AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[0] += 1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[0] += -0.027511960605097124 *
    AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[0] += 0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_h
    [3];
  _rtXdot->Dp46_CSTATE_h[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[1] += -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[1] += -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[1] += 0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[1] += -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[1] += 0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_h
    [4];
  _rtXdot->Dp46_CSTATE_h[2] += 0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[2] += 0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[2] += -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_h
    [2];
  _rtXdot->Dp46_CSTATE_h[2] += 2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[3] += 0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[3] += 1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[3] += -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[3] += -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[3] += 0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[4] += -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_h
    [3];
  _rtXdot->Dp46_CSTATE_h[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[0] += -0.33524288945878539 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp46_CSTATE_h[1] += -6.3625492730866693 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp46_CSTATE_h[2] += 0.93578679415036736 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp46_CSTATE_h[3] += 2.5673752417819116 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp46_CSTATE_h[4] += 0.014933981938330865 * AHV_Model_B.nu_r_c[5];

  // Derivatives for StateSpace: '<S22>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_m[0] += -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_m[0];
  _rtXdot->Dp53_CSTATE_m[0] += -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_m
    [1];
  _rtXdot->Dp53_CSTATE_m[0] += 2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_m[2];
  _rtXdot->Dp53_CSTATE_m[0] += 1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_m[3];
  _rtXdot->Dp53_CSTATE_m[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE_m[4];
  _rtXdot->Dp53_CSTATE_m[1] += 0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_m[0];
  _rtXdot->Dp53_CSTATE_m[1] += -5.3039287993722015E-5 *
    AHV_Model_X.Dp53_CSTATE_m[1];
  _rtXdot->Dp53_CSTATE_m[1] += 0.0077107435119831529 *
    AHV_Model_X.Dp53_CSTATE_m[2];
  _rtXdot->Dp53_CSTATE_m[1] += 0.0063309425482545112 *
    AHV_Model_X.Dp53_CSTATE_m[3];
  _rtXdot->Dp53_CSTATE_m[1] += 0.0012479014747342309 *
    AHV_Model_X.Dp53_CSTATE_m[4];
  _rtXdot->Dp53_CSTATE_m[2] += -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_m[0];
  _rtXdot->Dp53_CSTATE_m[2] += 0.0077107435118432491 *
    AHV_Model_X.Dp53_CSTATE_m[1];
  _rtXdot->Dp53_CSTATE_m[2] += -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_m[2];
  _rtXdot->Dp53_CSTATE_m[2] += -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_m[3];
  _rtXdot->Dp53_CSTATE_m[2] += -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_m[4];
  _rtXdot->Dp53_CSTATE_m[3] += -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_m[0];
  _rtXdot->Dp53_CSTATE_m[3] += 0.0063309425481496455 *
    AHV_Model_X.Dp53_CSTATE_m[1];
  _rtXdot->Dp53_CSTATE_m[3] += -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_m[2];
  _rtXdot->Dp53_CSTATE_m[3] += -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_m[3];
  _rtXdot->Dp53_CSTATE_m[3] += -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_m[4];
  _rtXdot->Dp53_CSTATE_m[4] += 0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_m[0];
  _rtXdot->Dp53_CSTATE_m[4] += -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_m
    [1];
  _rtXdot->Dp53_CSTATE_m[4] += 1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_m[2];
  _rtXdot->Dp53_CSTATE_m[4] += 3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_m[3];
  _rtXdot->Dp53_CSTATE_m[4] += -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_m
    [4];
  _rtXdot->Dp53_CSTATE_m[0] += -3.5374929885435322 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp53_CSTATE_m[1] += 0.015233386783081164 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp53_CSTATE_m[2] += -1.2196843916515756 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp53_CSTATE_m[3] += -0.96955326725759594 * AHV_Model_B.nu_r_c[2];
  _rtXdot->Dp53_CSTATE_m[4] += 0.17148246208757431 * AHV_Model_B.nu_r_c[2];

  // Derivatives for StateSpace: '<S22>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_n[0] += -8.7466634083809418E-5 *
    AHV_Model_X.Dp55_CSTATE_n[0];
  _rtXdot->Dp55_CSTATE_n[0] += 0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_n[1];
  _rtXdot->Dp55_CSTATE_n[0] += 0.0013720743421809071 *
    AHV_Model_X.Dp55_CSTATE_n[2];
  _rtXdot->Dp55_CSTATE_n[0] += 0.0075482812368147384 *
    AHV_Model_X.Dp55_CSTATE_n[3];
  _rtXdot->Dp55_CSTATE_n[0] += 0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_n[4];
  _rtXdot->Dp55_CSTATE_n[1] += -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_n
    [0];
  _rtXdot->Dp55_CSTATE_n[1] += -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_n[1];
  _rtXdot->Dp55_CSTATE_n[1] += -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_n
    [2];
  _rtXdot->Dp55_CSTATE_n[1] += -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_n[3];
  _rtXdot->Dp55_CSTATE_n[1] += -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_n[4];
  _rtXdot->Dp55_CSTATE_n[2] += -0.0013720743421854747 *
    AHV_Model_X.Dp55_CSTATE_n[0];
  _rtXdot->Dp55_CSTATE_n[2] += -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_n
    [1];
  _rtXdot->Dp55_CSTATE_n[2] += -0.068186331750790974 *
    AHV_Model_X.Dp55_CSTATE_n[2];
  _rtXdot->Dp55_CSTATE_n[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_n[3];
  _rtXdot->Dp55_CSTATE_n[2] += -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_n[4];
  _rtXdot->Dp55_CSTATE_n[3] += 0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_n
    [0];
  _rtXdot->Dp55_CSTATE_n[3] += 1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_n[1];
  _rtXdot->Dp55_CSTATE_n[3] += 4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_n[2];
  _rtXdot->Dp55_CSTATE_n[3] += -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_n[3];
  _rtXdot->Dp55_CSTATE_n[3] += -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_n[4];
  _rtXdot->Dp55_CSTATE_n[4] += 0.0070959954003254081 *
    AHV_Model_X.Dp55_CSTATE_n[0];
  _rtXdot->Dp55_CSTATE_n[4] += 1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_n[1];
  _rtXdot->Dp55_CSTATE_n[4] += 2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_n[2];
  _rtXdot->Dp55_CSTATE_n[4] += -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_n[3];
  _rtXdot->Dp55_CSTATE_n[4] += -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_n[4];
  _rtXdot->Dp55_CSTATE_n[0] += 0.021904981880023978 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp55_CSTATE_n[1] += 3.4651622640804733 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp55_CSTATE_n[2] += 0.16004490113712252 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp55_CSTATE_n[3] += -0.9969365784863482 * AHV_Model_B.nu_r_c[4];
  _rtXdot->Dp55_CSTATE_n[4] += -0.92774837285087253 * AHV_Model_B.nu_r_c[4];

  // Derivatives for StateSpace: '<S22>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_a[0] += -0.0017008072822443973 *
    AHV_Model_X.Dp62_CSTATE_a[0];
  _rtXdot->Dp62_CSTATE_a[0] += 1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_a[1];
  _rtXdot->Dp62_CSTATE_a[0] += 0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_a
    [2];
  _rtXdot->Dp62_CSTATE_a[0] += 0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_a
    [3];
  _rtXdot->Dp62_CSTATE_a[0] += -0.0047828760333047419 *
    AHV_Model_X.Dp62_CSTATE_a[4];
  _rtXdot->Dp62_CSTATE_a[1] += -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_a[0];
  _rtXdot->Dp62_CSTATE_a[1] += -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_a[1];
  _rtXdot->Dp62_CSTATE_a[1] += -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_a
    [2];
  _rtXdot->Dp62_CSTATE_a[1] += -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_a
    [3];
  _rtXdot->Dp62_CSTATE_a[1] += 0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_a[4];
  _rtXdot->Dp62_CSTATE_a[2] += -0.012808448456765344 *
    AHV_Model_X.Dp62_CSTATE_a[0];
  _rtXdot->Dp62_CSTATE_a[2] += -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_a
    [1];
  _rtXdot->Dp62_CSTATE_a[2] += -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_a
    [2];
  _rtXdot->Dp62_CSTATE_a[2] += -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_a[3];
  _rtXdot->Dp62_CSTATE_a[2] += 0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_a[4];
  _rtXdot->Dp62_CSTATE_a[3] += 0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_a
    [0];
  _rtXdot->Dp62_CSTATE_a[3] += 0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_a[1];
  _rtXdot->Dp62_CSTATE_a[3] += 1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_a[2];
  _rtXdot->Dp62_CSTATE_a[3] += -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_a[3];
  _rtXdot->Dp62_CSTATE_a[3] += 1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_a[4];
  _rtXdot->Dp62_CSTATE_a[4] += 0.0047828760333030167 *
    AHV_Model_X.Dp62_CSTATE_a[0];
  _rtXdot->Dp62_CSTATE_a[4] += 0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_a[1];
  _rtXdot->Dp62_CSTATE_a[4] += 0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_a[2];
  _rtXdot->Dp62_CSTATE_a[4] += -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_a[3];
  _rtXdot->Dp62_CSTATE_a[4] += -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_a
    [4];
  _rtXdot->Dp62_CSTATE_a[0] += 0.13167647316600378 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp62_CSTATE_a[1] += 3.4125284827245985 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp62_CSTATE_a[2] += 0.44399732492152644 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp62_CSTATE_a[3] += -1.2820687298457427 * AHV_Model_B.nu_r_c[1];
  _rtXdot->Dp62_CSTATE_a[4] += -0.18276601298221137 * AHV_Model_B.nu_r_c[1];

  // Derivatives for StateSpace: '<S22>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_n[0] += -0.004396104715914141 *
    AHV_Model_X.Dp64_CSTATE_n[0];
  _rtXdot->Dp64_CSTATE_n[0] += 1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_n[1];
  _rtXdot->Dp64_CSTATE_n[0] += -0.027511960605097124 *
    AHV_Model_X.Dp64_CSTATE_n[2];
  _rtXdot->Dp64_CSTATE_n[0] += 0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_n
    [3];
  _rtXdot->Dp64_CSTATE_n[0] += -0.0003953747088775357 *
    AHV_Model_X.Dp64_CSTATE_n[4];
  _rtXdot->Dp64_CSTATE_n[1] += -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_n[0];
  _rtXdot->Dp64_CSTATE_n[1] += -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_n[1];
  _rtXdot->Dp64_CSTATE_n[1] += 0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_n[2];
  _rtXdot->Dp64_CSTATE_n[1] += -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_n[3];
  _rtXdot->Dp64_CSTATE_n[1] += 0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_n
    [4];
  _rtXdot->Dp64_CSTATE_n[2] += 0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_n[0];
  _rtXdot->Dp64_CSTATE_n[2] += 0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_n[1];
  _rtXdot->Dp64_CSTATE_n[2] += -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_n
    [2];
  _rtXdot->Dp64_CSTATE_n[2] += 2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_n[3];
  _rtXdot->Dp64_CSTATE_n[2] += -0.0093184009624265578 *
    AHV_Model_X.Dp64_CSTATE_n[4];
  _rtXdot->Dp64_CSTATE_n[3] += 0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_n[0];
  _rtXdot->Dp64_CSTATE_n[3] += 1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_n[1];
  _rtXdot->Dp64_CSTATE_n[3] += -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_n[2];
  _rtXdot->Dp64_CSTATE_n[3] += -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_n[3];
  _rtXdot->Dp64_CSTATE_n[3] += 0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_n[4];
  _rtXdot->Dp64_CSTATE_n[4] += 0.00039537470888004991 *
    AHV_Model_X.Dp64_CSTATE_n[0];
  _rtXdot->Dp64_CSTATE_n[4] += 0.0083561289799246388 *
    AHV_Model_X.Dp64_CSTATE_n[1];
  _rtXdot->Dp64_CSTATE_n[4] += -0.0093184009624202816 *
    AHV_Model_X.Dp64_CSTATE_n[2];
  _rtXdot->Dp64_CSTATE_n[4] += -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_n
    [3];
  _rtXdot->Dp64_CSTATE_n[4] += -0.00092959984608804514 *
    AHV_Model_X.Dp64_CSTATE_n[4];
  _rtXdot->Dp64_CSTATE_n[0] += -0.33524288945878539 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp64_CSTATE_n[1] += -6.3625492730866693 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp64_CSTATE_n[2] += 0.93578679415036736 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp64_CSTATE_n[3] += 2.5673752417819116 * AHV_Model_B.nu_r_c[3];
  _rtXdot->Dp64_CSTATE_n[4] += 0.014933981938330865 * AHV_Model_B.nu_r_c[3];

  // Derivatives for StateSpace: '<S22>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_fj[0] += -1.5814921936539046 *
    AHV_Model_X.Dp66_CSTATE_fj[0];
  _rtXdot->Dp66_CSTATE_fj[0] += -1.3217886535784036 *
    AHV_Model_X.Dp66_CSTATE_fj[1];
  _rtXdot->Dp66_CSTATE_fj[0] += -0.43878308112001396 *
    AHV_Model_X.Dp66_CSTATE_fj[2];
  _rtXdot->Dp66_CSTATE_fj[0] += -1.2174625158903036 *
    AHV_Model_X.Dp66_CSTATE_fj[3];
  _rtXdot->Dp66_CSTATE_fj[0] += 0.069477290559771143 *
    AHV_Model_X.Dp66_CSTATE_fj[4];
  _rtXdot->Dp66_CSTATE_fj[1] += 1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_fj
    [0];
  _rtXdot->Dp66_CSTATE_fj[1] += -0.00029199832472599611 *
    AHV_Model_X.Dp66_CSTATE_fj[1];
  _rtXdot->Dp66_CSTATE_fj[1] += -0.0080078570280257711 *
    AHV_Model_X.Dp66_CSTATE_fj[2];
  _rtXdot->Dp66_CSTATE_fj[1] += -0.01568400806720225 *
    AHV_Model_X.Dp66_CSTATE_fj[3];
  _rtXdot->Dp66_CSTATE_fj[1] += 0.00098472174210678556 *
    AHV_Model_X.Dp66_CSTATE_fj[4];
  _rtXdot->Dp66_CSTATE_fj[2] += -0.43878308111999881 *
    AHV_Model_X.Dp66_CSTATE_fj[0];
  _rtXdot->Dp66_CSTATE_fj[2] += 0.00800785702802341 *
    AHV_Model_X.Dp66_CSTATE_fj[1];
  _rtXdot->Dp66_CSTATE_fj[2] += -0.28893903358444806 *
    AHV_Model_X.Dp66_CSTATE_fj[2];
  _rtXdot->Dp66_CSTATE_fj[2] += -1.8221092762782123 *
    AHV_Model_X.Dp66_CSTATE_fj[3];
  _rtXdot->Dp66_CSTATE_fj[2] += 0.074587161707531782 *
    AHV_Model_X.Dp66_CSTATE_fj[4];
  _rtXdot->Dp66_CSTATE_fj[3] += 1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_fj
    [0];
  _rtXdot->Dp66_CSTATE_fj[3] += -0.015684008067190565 *
    AHV_Model_X.Dp66_CSTATE_fj[1];
  _rtXdot->Dp66_CSTATE_fj[3] += 1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_fj
    [2];
  _rtXdot->Dp66_CSTATE_fj[3] += -6.3850648427943675 *
    AHV_Model_X.Dp66_CSTATE_fj[3];
  _rtXdot->Dp66_CSTATE_fj[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE_fj[4];
  _rtXdot->Dp66_CSTATE_fj[4] += 0.069477290559780386 *
    AHV_Model_X.Dp66_CSTATE_fj[0];
  _rtXdot->Dp66_CSTATE_fj[4] += -0.000984721742106249 *
    AHV_Model_X.Dp66_CSTATE_fj[1];
  _rtXdot->Dp66_CSTATE_fj[4] += 0.074587161707547089 *
    AHV_Model_X.Dp66_CSTATE_fj[2];
  _rtXdot->Dp66_CSTATE_fj[4] += -1.2174309773027607 *
    AHV_Model_X.Dp66_CSTATE_fj[3];
  _rtXdot->Dp66_CSTATE_fj[4] += -0.061420217921035136 *
    AHV_Model_X.Dp66_CSTATE_fj[4];
  _rtXdot->Dp66_CSTATE_fj[0] += -3.3844821061206916 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp66_CSTATE_fj[1] += 0.045246867266885989 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp66_CSTATE_fj[2] += -0.53332191483729374 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp66_CSTATE_fj[3] += 1.2581308985447515 * AHV_Model_B.nu_r_c[5];
  _rtXdot->Dp66_CSTATE_fj[4] += 0.075289942321940723 * AHV_Model_B.nu_r_c[5];

  // Derivatives for Integrator: '<S65>/Integrator6'
  _rtXdot->Integrator6_CSTATE[0] = AHV_Model_B.Sum7[0];

  // Derivatives for Integrator: '<S65>/Integrator1'
  _rtXdot->Integrator1_CSTATE_a[0] = AHV_Model_B.Sum6[0];

  // Derivatives for Integrator: '<S65>/Integrator2'
  _rtXdot->Integrator2_CSTATE[0] = AHV_Model_B.psi_WF[0];

  // Derivatives for Integrator: '<S161>/Integrator6'
  _rtXdot->Integrator6_CSTATE_p[0] = AHV_Model_B.Sum7_o[0];

  // Derivatives for Integrator: '<S161>/Integrator1'
  _rtXdot->Integrator1_CSTATE_h0[0] = AHV_Model_B.Sum6_d[0];

  // Derivatives for Integrator: '<S161>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[0] = AHV_Model_B.psi_WF_o[0];

  // Derivatives for Integrator: '<S257>/Integrator6'
  _rtXdot->Integrator6_CSTATE_i[0] = AHV_Model_B.Sum7_n[0];

  // Derivatives for Integrator: '<S257>/Integrator1'
  _rtXdot->Integrator1_CSTATE_hl[0] = AHV_Model_B.Sum6_j[0];

  // Derivatives for Integrator: '<S257>/Integrator2'
  _rtXdot->Integrator2_CSTATE_f[0] = AHV_Model_B.psi_WF_m[0];

  // Derivatives for Integrator: '<S353>/Integrator6'
  _rtXdot->Integrator6_CSTATE_ij[0] = AHV_Model_B.Sum7_n4[0];

  // Derivatives for Integrator: '<S353>/Integrator1'
  _rtXdot->Integrator1_CSTATE_nd[0] = AHV_Model_B.Sum6_i[0];

  // Derivatives for Integrator: '<S353>/Integrator2'
  _rtXdot->Integrator2_CSTATE_p[0] = AHV_Model_B.psi_WF_i[0];

  // Derivatives for Integrator: '<S65>/Integrator6'
  _rtXdot->Integrator6_CSTATE[1] = AHV_Model_B.Sum7[1];

  // Derivatives for Integrator: '<S65>/Integrator1'
  _rtXdot->Integrator1_CSTATE_a[1] = AHV_Model_B.Sum6[1];

  // Derivatives for Integrator: '<S65>/Integrator2'
  _rtXdot->Integrator2_CSTATE[1] = AHV_Model_B.psi_WF[1];

  // Derivatives for Integrator: '<S161>/Integrator6'
  _rtXdot->Integrator6_CSTATE_p[1] = AHV_Model_B.Sum7_o[1];

  // Derivatives for Integrator: '<S161>/Integrator1'
  _rtXdot->Integrator1_CSTATE_h0[1] = AHV_Model_B.Sum6_d[1];

  // Derivatives for Integrator: '<S161>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[1] = AHV_Model_B.psi_WF_o[1];

  // Derivatives for Integrator: '<S257>/Integrator6'
  _rtXdot->Integrator6_CSTATE_i[1] = AHV_Model_B.Sum7_n[1];

  // Derivatives for Integrator: '<S257>/Integrator1'
  _rtXdot->Integrator1_CSTATE_hl[1] = AHV_Model_B.Sum6_j[1];

  // Derivatives for Integrator: '<S257>/Integrator2'
  _rtXdot->Integrator2_CSTATE_f[1] = AHV_Model_B.psi_WF_m[1];

  // Derivatives for Integrator: '<S353>/Integrator6'
  _rtXdot->Integrator6_CSTATE_ij[1] = AHV_Model_B.Sum7_n4[1];

  // Derivatives for Integrator: '<S353>/Integrator1'
  _rtXdot->Integrator1_CSTATE_nd[1] = AHV_Model_B.Sum6_i[1];

  // Derivatives for Integrator: '<S353>/Integrator2'
  _rtXdot->Integrator2_CSTATE_p[1] = AHV_Model_B.psi_WF_i[1];

  // Derivatives for Integrator: '<S65>/Integrator6'
  _rtXdot->Integrator6_CSTATE[2] = AHV_Model_B.Sum7[2];

  // Derivatives for Integrator: '<S65>/Integrator1'
  _rtXdot->Integrator1_CSTATE_a[2] = AHV_Model_B.Sum6[2];

  // Derivatives for Integrator: '<S65>/Integrator2'
  _rtXdot->Integrator2_CSTATE[2] = AHV_Model_B.psi_WF[2];

  // Derivatives for Integrator: '<S161>/Integrator6'
  _rtXdot->Integrator6_CSTATE_p[2] = AHV_Model_B.Sum7_o[2];

  // Derivatives for Integrator: '<S161>/Integrator1'
  _rtXdot->Integrator1_CSTATE_h0[2] = AHV_Model_B.Sum6_d[2];

  // Derivatives for Integrator: '<S161>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[2] = AHV_Model_B.psi_WF_o[2];

  // Derivatives for Integrator: '<S257>/Integrator6'
  _rtXdot->Integrator6_CSTATE_i[2] = AHV_Model_B.Sum7_n[2];

  // Derivatives for Integrator: '<S257>/Integrator1'
  _rtXdot->Integrator1_CSTATE_hl[2] = AHV_Model_B.Sum6_j[2];

  // Derivatives for Integrator: '<S257>/Integrator2'
  _rtXdot->Integrator2_CSTATE_f[2] = AHV_Model_B.psi_WF_m[2];

  // Derivatives for Integrator: '<S353>/Integrator6'
  _rtXdot->Integrator6_CSTATE_ij[2] = AHV_Model_B.Sum7_n4[2];

  // Derivatives for Integrator: '<S353>/Integrator1'
  _rtXdot->Integrator1_CSTATE_nd[2] = AHV_Model_B.Sum6_i[2];

  // Derivatives for Integrator: '<S353>/Integrator2'
  _rtXdot->Integrator2_CSTATE_p[2] = AHV_Model_B.psi_WF_i[2];

  // Derivatives for TransferFcn: '<S18>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += -0.2 * AHV_Model_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += AHV_Model_B.Fcn;

  // Derivatives for TransferFcn: '<S114>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_c = 0.0;
  _rtXdot->TransferFcn_CSTATE_c += -0.2 * AHV_Model_X.TransferFcn_CSTATE_c;
  _rtXdot->TransferFcn_CSTATE_c += AHV_Model_B.Fcn_j;

  // Derivatives for TransferFcn: '<S210>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_d = 0.0;
  _rtXdot->TransferFcn_CSTATE_d += -0.2 * AHV_Model_X.TransferFcn_CSTATE_d;
  _rtXdot->TransferFcn_CSTATE_d += AHV_Model_B.Fcn_b;

  // Derivatives for TransferFcn: '<S306>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_p = 0.0;
  _rtXdot->TransferFcn_CSTATE_p += -0.2 * AHV_Model_X.TransferFcn_CSTATE_p;
  _rtXdot->TransferFcn_CSTATE_p += AHV_Model_B.Fcn_e;
}

// Model initialize function
void AH_Model_v1ModelClass::initialize(double dtime)
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&AHV_Model_M)->solverInfo, &rtmGetTPtr((&AHV_Model_M)));
    rtsiSetStepSizePtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)
                       ->Timing.stepSize0);
    rtsiSetdXPtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)->derivs);
    rtsiSetContStatesPtr(&(&AHV_Model_M)->solverInfo, (real_T **) &(&AHV_Model_M)
                         ->contStates);
    rtsiSetNumContStatesPtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M
      )->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&AHV_Model_M)->solverInfo,
      &(&AHV_Model_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&AHV_Model_M)->solverInfo,
      &(&AHV_Model_M)->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&(&AHV_Model_M)->solverInfo, (&rtmGetErrorStatus
      ((&AHV_Model_M))));
    rtsiSetRTModelPtr(&(&AHV_Model_M)->solverInfo, (&AHV_Model_M));
  }

  rtsiSetSimTimeStep(&(&AHV_Model_M)->solverInfo, MAJOR_TIME_STEP);
  (&AHV_Model_M)->intgData.y = (&AHV_Model_M)->odeY;
  (&AHV_Model_M)->intgData.f[0] = (&AHV_Model_M)->odeF[0];
  (&AHV_Model_M)->intgData.f[1] = (&AHV_Model_M)->odeF[1];
  (&AHV_Model_M)->intgData.f[2] = (&AHV_Model_M)->odeF[2];
  (&AHV_Model_M)->intgData.f[3] = (&AHV_Model_M)->odeF[3];
  (&AHV_Model_M)->contStates = ((X_AHV_Model_T *) &AHV_Model_X);
  rtsiSetSolverData(&(&AHV_Model_M)->solverInfo, static_cast<void *>
                    (&(&AHV_Model_M)->intgData));
  rtsiSetSolverName(&(&AHV_Model_M)->solverInfo,"dtime");
  rtmSetTPtr((&AHV_Model_M), &(&AHV_Model_M)->Timing.tArray[0]);
  (&AHV_Model_M)->Timing.stepSize0 = 0.05;
  rtmSetFirstInitCond((&AHV_Model_M), 1);

  {
    real_T heading_deg;
    int32_T i;

    // Start for If: '<S285>/If1'
    AHV_Model_DW.If1_ActiveSubsystem = -1;

    // Start for If: '<S189>/If1'
    AHV_Model_DW.If1_ActiveSubsystem_p = -1;

    // Start for If: '<S93>/If1'
    AHV_Model_DW.If1_ActiveSubsystem_g = -1;

    // Start for If: '<S8>/If'
    AHV_Model_DW.If_ActiveSubsystem = -1;

    // Start for If: '<S104>/If'
    AHV_Model_DW.If_ActiveSubsystem_l = -1;

    // Start for If: '<S200>/If'
    AHV_Model_DW.If_ActiveSubsystem_a = -1;

    // Start for If: '<S381>/If1'
    AHV_Model_DW.If1_ActiveSubsystem_l = -1;

    // Start for If: '<S296>/If'
    AHV_Model_DW.If_ActiveSubsystem_k = -1;

    // Start for If: '<S373>/If'
    AHV_Model_DW.If_ActiveSubsystem_i = -1;

    // Start for If: '<S374>/If'
    AHV_Model_DW.If_ActiveSubsystem_j = -1;

    // Start for If: '<S366>/If'
    AHV_Model_DW.If_ActiveSubsystem_d = -1;

    // Start for If: '<S277>/If'
    AHV_Model_DW.If_ActiveSubsystem_h = -1;

    // Start for If: '<S278>/If'
    AHV_Model_DW.If_ActiveSubsystem_m = -1;

    // Start for If: '<S270>/If'
    AHV_Model_DW.If_ActiveSubsystem_c = -1;

    // Start for If: '<S181>/If'
    AHV_Model_DW.If_ActiveSubsystem_e = -1;

    // Start for If: '<S182>/If'
    AHV_Model_DW.If_ActiveSubsystem_g = -1;

    // Start for If: '<S174>/If'
    AHV_Model_DW.If_ActiveSubsystem_dx = -1;

    // Start for If: '<S85>/If'
    AHV_Model_DW.If_ActiveSubsystem_m5 = -1;

    // Start for If: '<S86>/If'
    AHV_Model_DW.If_ActiveSubsystem_ey = -1;

    // Start for If: '<S78>/If'
    AHV_Model_DW.If_ActiveSubsystem_h0 = -1;

    // InitializeConditions for Integrator: '<S15>/Integrator1' incorporates:
    //   Integrator: '<S111>/Integrator1'

    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[5] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_n[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK = 1;

    // End of InitializeConditions for Integrator: '<S15>/Integrator1'

    // InitializeConditions for Integrator: '<S111>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_o = 1;

    // InitializeConditions for Integrator: '<S207>/Integrator1' incorporates:
    //   Integrator: '<S303>/Integrator1'

    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE_h[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_h[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_h[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_h[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_h[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_h[5] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_hn[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK_b = 1;

    // End of InitializeConditions for Integrator: '<S207>/Integrator1'

    // InitializeConditions for Integrator: '<S303>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_c = 1;

    // InitializeConditions for Integrator: '<S303>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S303>/Integrator'

    // InitializeConditions for RateLimiter: '<S66>/Rate Limiter'
    AHV_Model_DW.LastMajorTime = (rtInf);

    // InitializeConditions for RateLimiter: '<S162>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_g = (rtInf);

    // InitializeConditions for RateLimiter: '<S258>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_a = (rtInf);

    // InitializeConditions for UnitDelay: '<S9>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE[0] = 1.0;

    // InitializeConditions for UnitDelay: '<S105>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_h[0] = 1.0;

    // InitializeConditions for UnitDelay: '<S201>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_c[0] = 1.0;

    // InitializeConditions for UnitDelay: '<S297>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_a[0] = 1.0;

    // InitializeConditions for UnitDelay: '<S9>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE[1] = 1.0;

    // InitializeConditions for UnitDelay: '<S105>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_h[1] = 1.0;

    // InitializeConditions for UnitDelay: '<S201>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_c[1] = 1.0;

    // InitializeConditions for UnitDelay: '<S297>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_a[1] = 1.0;

    // InitializeConditions for UnitDelay: '<S9>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE[2] = 1.0;

    // InitializeConditions for UnitDelay: '<S105>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_h[2] = 1.0;

    // InitializeConditions for UnitDelay: '<S201>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_c[2] = 1.0;

    // InitializeConditions for UnitDelay: '<S297>/Unit Delay'
    AHV_Model_DW.UnitDelay_DSTATE_a[2] = 1.0;

    // InitializeConditions for Integrator: '<S353>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK = 1;

    // End of InitializeConditions for Integrator: '<S353>/Integrator3'

    // InitializeConditions for RateLimiter: '<S354>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_c = (rtInf);

    // InitializeConditions for Integrator: '<S353>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S310>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE[i] = 0.0;

      // InitializeConditions for Integrator: '<S322>/Integrator'
      AHV_Model_X.Integrator_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S310>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S207>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_h[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S207>/Integrator'

    // InitializeConditions for Integrator: '<S257>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_c[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_c[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_c[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_i = 1;

    // End of InitializeConditions for Integrator: '<S257>/Integrator3'

    // InitializeConditions for Integrator: '<S257>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_f[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_f[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_f[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S214>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_j[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_b[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_n[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_m[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_b[i] = 0.0;

      // InitializeConditions for Integrator: '<S226>/Integrator'
      AHV_Model_X.Integrator_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_o[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_g[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S214>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_f[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S111>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_o[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S111>/Integrator'

    // InitializeConditions for Integrator: '<S161>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_l[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_l[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_l[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_c = 1;

    // End of InitializeConditions for Integrator: '<S161>/Integrator3'

    // InitializeConditions for Integrator: '<S161>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_p[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_p[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_p[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S118>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_o[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_cn[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_g[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_jf[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_j[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_j[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_g[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_g[i] = 0.0;

      // InitializeConditions for Integrator: '<S130>/Integrator'
      AHV_Model_X.Integrator_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S118>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_n[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S15>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_e[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S15>/Integrator'

    // InitializeConditions for Integrator: '<S65>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_n[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_n[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_n[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_l = 1;

    // End of InitializeConditions for Integrator: '<S65>/Integrator3'

    // InitializeConditions for Integrator: '<S65>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_c[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_c[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_c[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S22>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_m[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_b[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_o[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_e[i] = 0.0;

      // InitializeConditions for Integrator: '<S34>/Integrator'
      AHV_Model_X.Integrator_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_g[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_m[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_n[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_n[i] = 0.0;

      // InitializeConditions for StateSpace: '<S22>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_fj[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S65>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[0] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_a[0] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[0] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_p[0] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_h0[0] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[0] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_i[0] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_hl[0] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_f[0] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_ij[0] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_nd[0] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_p[0] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[1] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_a[1] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[1] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_p[1] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_h0[1] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[1] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_i[1] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_hl[1] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_f[1] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_ij[1] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_nd[1] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_p[1] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[2] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_a[2] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[2] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_p[2] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_h0[2] = 0.0;

    // InitializeConditions for Integrator: '<S161>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[2] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_i[2] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_hl[2] = 0.0;

    // InitializeConditions for Integrator: '<S257>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_f[2] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_ij[2] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_nd[2] = 0.0;

    // InitializeConditions for Integrator: '<S353>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_p[2] = 0.0;

    // InitializeConditions for TransferFcn: '<S18>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S114>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_c = 0.0;

    // InitializeConditions for TransferFcn: '<S210>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_d = 0.0;

    // InitializeConditions for TransferFcn: '<S306>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_p = 0.0;

    // SystemInitialize for Atomic SubSystem: '<S10>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3);

    // End of SystemInitialize for SubSystem: '<S10>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S106>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_d);

    // End of SystemInitialize for SubSystem: '<S106>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S202>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_h);

    // End of SystemInitialize for SubSystem: '<S202>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S298>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_o);

    // End of SystemInitialize for SubSystem: '<S298>/Subsystem3'

    // SystemInitialize for IfAction SubSystem: '<S8>/If Action Subsystem'
    // InitializeConditions for RateLimiter: '<S94>/Rate Limiter1'
    AHV_Model_DW.LastMajorTime_e = (rtInf);

    // InitializeConditions for RateLimiter: '<S94>/Rate Limiter2'
    AHV_Model_DW.LastMajorTime_h = (rtInf);

    // SystemInitialize for Chart: '<S94>/Chart2'
    AHV_Model_Chart1_Init(AHV_Model_B.RateLimiter1, AHV_Model_B.Gain4,
                          &AHV_Model_B.heading_deg_g, &AHV_Model_DW.sf_Chart2);

    // SystemInitialize for Chart: '<S94>/Chart3'
    AHV_Model_Chart3_Init(AHV_Model_B.Integrator1[0], AHV_Model_B.Product,
                          &heading_deg, &AHV_Model_DW.sf_Chart3);

    // SystemInitialize for Chart: '<S94>/Chart1'
    AHV_Model_Chart1_Init(AHV_Model_B.Product1, AHV_Model_B.Integrator1[1],
                          &heading_deg, &AHV_Model_DW.sf_Chart1);

    // End of SystemInitialize for SubSystem: '<S8>/If Action Subsystem'

    // SystemInitialize for IfAction SubSystem: '<S104>/If Action Subsystem'
    AHV_Mo_IfActionSubsystem_j_Init(AHV_Model_B.Integrator1_b,
      &AHV_Model_B.IfActionSubsystem_c, &AHV_Model_DW.IfActionSubsystem_c);

    // End of SystemInitialize for SubSystem: '<S104>/If Action Subsystem'

    // SystemInitialize for IfAction SubSystem: '<S200>/If Action Subsystem'
    AHV_Mo_IfActionSubsystem_j_Init(AHV_Model_B.Integrator1_p,
      &AHV_Model_B.IfActionSubsystem_e, &AHV_Model_DW.IfActionSubsystem_e);

    // End of SystemInitialize for SubSystem: '<S200>/If Action Subsystem'

    // SystemInitialize for IfAction SubSystem: '<S296>/If Action Subsystem'
    AHV_Mo_IfActionSubsystem_j_Init(AHV_Model_B.Integrator1_n,
      &AHV_Model_B.IfActionSubsystem_eo, &AHV_Model_DW.IfActionSubsystem_eo);

    // End of SystemInitialize for SubSystem: '<S296>/If Action Subsystem'
  }

  // set "at time zero" to false
  if (rtmIsFirstInitCond((&AHV_Model_M))) {
    rtmSetFirstInitCond((&AHV_Model_M), 0);
  }
}

// Model terminate function
void AH_Model_v1ModelClass::terminate()
{
  // (no terminate code required)
}

// Constructor
AH_Model_v1ModelClass::AH_Model_v1ModelClass():
  AHV_Model_B()
  ,AHV_Model_DW()
  ,AHV_Model_X()
  ,AHV_Model_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
AH_Model_v1ModelClass::~AH_Model_v1ModelClass()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
AH_Model_v1ModelClass::RT_MODEL_AHV_Model_T * AH_Model_v1ModelClass::getRTM()
{
  return (&AHV_Model_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
