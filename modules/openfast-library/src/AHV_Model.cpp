//
// File: AHV_Model.cpp
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
#include "AHV_Model.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "AHV_Model_private.h"
using namespace std;
std::vector<std::string> entries;
std::ofstream fout;
int write_local = 0;

int fileload = 1;

// Named constants for Chart: '<S64>/Chart'
const uint8_T AHV_Model_IN_hold = 1U;
const uint8_T AHV_Model_IN_unhold = 2U;
//--------��������---------
real_T Featured = 1;
real_T SidePush = 10;
real_T MainThrust = 1;
real_T FuThrust = 1;
real_T AzimuthThrust = 10;

// Exported block signals
real_T spectrum_type;       // '<Root>/spectrum_type'
real_T hs;                  // '<Root>/hs'
real_T omega_peak;          // '<Root>/omega_peak'
real_T psi_mean;            // '<Root>/psi_mean'
real_T gamma_value;         // '<Root>/gamma_value'
real_T spread;              // '<Root>/spread'
real_T depth;               // '<Root>/depth'
real_T nfreq;               // '<Root>/nfreq'
real_T ndir;                // '<Root>/ndir'
real_T energylim;           // '<Root>/energylim'
real_T freq_cutoff;         // '<Root>/freq_cutoff'
real_T dir_cutoff;          // '<Root>/dir_cutoff'
real_T rand_freq;           // '<Root>/rand_freq'
real_T rand_dir;            // '<Root>/rand_dir'
real_T Current_direction;   // '<Root>/Current_direction'
real_T Current_speed;       // '<Root>/Current_speed'
real_T Vessel_init1[6];     // '<Root>/Vessel_init1'
boolean_T Hold_Position1;   // '<Root>/Hold_Position1'
real_T heading_mode1;       // '<Root>/heading_mode1'
real_T heading_angle_ref1;  // '<Root>/heading_angle_ref1'
real_T Vessel_X_Ref1;       // '<Root>/Vessel_X_Ref1'
real_T Vessel_Y_Ref1;       // '<Root>/Vessel_Y_Ref1'
real_T ahv_fairlead1[3];    // '<Root>/ahv_fairlead1'
real_T tau_cable1[3];       // '<Root>/tau_cable1'
real_T Vessel_init2[6];     // '<Root>/Vessel_init2'
boolean_T Hold_Position2;   // '<Root>/Hold_Position2'
real_T heading_mode2;       // '<Root>/heading_mode2'
real_T heading_angle_ref2;  // '<Root>/heading_angle_ref2'
real_T Vessel_X_Ref2;       // '<Root>/Vessel_X_Ref2'
real_T Vessel_Y_Ref2;       // '<Root>/Vessel_Y_Ref2'
real_T ahv_fairlead2[3];    // '<Root>/ahv_fairlead2'
real_T tau_cable2[3];       // '<Root>/tau_cable2'
real_T Vessel_init3[6];     // '<Root>/Vessel_init3'
boolean_T Hold_Position3;   // '<Root>/Hold_Position3'
real_T heading_mode3;       // '<Root>/heading_mode3'
real_T heading_angle_ref3;  // '<Root>/heading_angle_ref3'
real_T Vessel_X_Ref3;       // '<Root>/Vessel_X_Ref3'
real_T Vessel_Y_Ref3;       // '<Root>/Vessel_Y_Ref3'
real_T ahv_fairlead3[3];    // '<Root>/ahv_fairlead3'
real_T tau_cable3[3];       // '<Root>/tau_cable3'
real_T Vessel_init4[6];     // '<Root>/Vessel_init4'
boolean_T Hold_Position4;   // '<Root>/Hold_Position4'
real_T heading_mode4;       // '<Root>/heading_mode4'
real_T heading_angle_ref4;  // '<Root>/heading_angle_ref4'
real_T Vessel_X_Ref4;       // '<Root>/Vessel_X_Ref4'
real_T Vessel_Y_Ref4;       // '<Root>/Vessel_Y_Ref4'
real_T ahv_fairlead4[3];    // '<Root>/ahv_fairlead4'
real_T tau_cable4[3];       // '<Root>/tau_cable4'
real_T eta_AHV1[6];         // '<Root>/eta_AHV1'
real_T nu1[6];              // '<Root>/nu1'
real_T AHV_speed1;          // '<Root>/AHV_speed1'
real_T eta_AHV2[6];         // '<Root>/eta_AHV2'
real_T nu2[6];              // '<Root>/nu2'
real_T AHV_speed2;          // '<Root>/AHV_speed2'
real_T eta_AHV3[6];         // '<Root>/eta_AHV3'
real_T nu3[6];              // '<Root>/nu3'
real_T AHV_speed3;          // '<Root>/AHV_speed3'
real_T eta_AHV4[6];         // '<Root>/eta_AHV4'
real_T nu4[6];              // '<Root>/nu4'
real_T AHV_speed4;          // '<Root>/AHV_speed4'

string ReadLine(string filename, int line) {
  int i = 0;
  string temp;
  fstream file;
  file.open(filename, ios::in);

  if (line <= 0) {
    fileload = 0;
    // cout << "Error 1: �������󣬲���Ϊ0��������" <<
    // endl;
  }

  if (file.fail()) {
    fileload = 0;
    // cout << "Error 2: �ļ������ڡ�" << endl;
  }

  while (getline(file, temp) && i < line - 1) {
    i++;
  }

  file.close();
  return temp;
}

vector<string> split(const string& s,
                     char delim)  // TODO: remove delim arg, it's unused!
{
  vector<string> elems;  // the vector of words to return

  char str[100];                 // this gives some memory for the char array
  strncpy(str, s.c_str(), 100);  // copy input string to str
  char* pch;
  pch = strtok(str, " \t");  // give strtok the c string of s
  while (pch != NULL) {
    elems.push_back(pch);
    pch = strtok(NULL, " \t");  // split by spaces or tabs
  }
  return elems;
}

uint32_T plook_u32d_evencka(real_T u, real_T bp0, real_T bpSpace,
                            uint32_T maxIndex) {
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

uint32_T plook_lincp(real_T u, const real_T bp[], uint32_T maxIndex,
                     real_T* fraction, uint32_T* prevIndex) {
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

real_T intrp3d_l_pw(const uint32_T bpIndex[], const real_T frac[],
                    const real_T table[], const uint32_T stride[]) {
  real_T yL_2d;
  uint32_T offset_2d;
  real_T yL_1d;
  uint32_T offset_0d;

  // Column-major Interpolation 3-D
  // Interpolation method: 'Linear point-slope'
  // Use last breakpoint for index at or above upper limit: 'off'
  // Overflow mode: 'portable wrapping'

  offset_2d =
      (bpIndex[2U] * stride[2U] + bpIndex[1U] * stride[1U]) + bpIndex[0U];
  yL_1d =
      (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] + table[offset_2d];
  offset_0d = offset_2d + stride[1U];
  yL_2d = (((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
            table[offset_0d]) -
           yL_1d) *
              frac[1U] +
          yL_1d;
  offset_2d += stride[2U];
  yL_1d =
      (table[offset_2d + 1U] - table[offset_2d]) * frac[0U] + table[offset_2d];
  offset_0d = offset_2d + stride[1U];
  return (((((table[offset_0d + 1U] - table[offset_0d]) * frac[0U] +
             table[offset_0d]) -
            yL_1d) *
               frac[1U] +
           yL_1d) -
          yL_2d) *
             frac[2U] +
         yL_2d;
}

uint32_T linsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex) {
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
                     real_T* fraction) {
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
void AH_Model_v1ModelClass::rt_ertODEUpdateContinuousStates(RTWSolverInfo* si) {
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T* x = rtsiGetContStates(si);
  ODE4_IntgData* id = static_cast<ODE4_IntgData*>(rtsiGetSolverData(si));
  real_T* y = id->y;
  real_T* f0 = id->f[0];
  real_T* f1 = id->f[1];
  real_T* f2 = id->f[2];
  real_T* f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 475;
  rtsiSetSimTimeStep(si, MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void)std::memcpy(y, x, static_cast<uint_T>(nXc) * sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  AHV_Model_derivatives();

  // f1 = f(t + (h/2), y + (h/2)*f0)
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp * f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  this->step();
  AHV_Model_derivatives();

  // f2 = f(t + (h/2), y + (h/2)*f1)
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp * f1[i]);
  }

  rtsiSetdX(si, f2);
  this->step();
  AHV_Model_derivatives();

  // f3 = f(t + h, y + h*f2)
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h * f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  this->step();
  AHV_Model_derivatives();

  // tnew = t + h
  // ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3)
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si, MAJOR_TIME_STEP);
}

//
// System reset for iterator system:
//    '<S19>/Cross-flow drag trapezoidal integration'
//    '<S103>/Cross-flow drag trapezoidal integration'
//    '<S187>/Cross-flow drag trapezoidal integration'
//    '<S271>/Cross-flow drag trapezoidal integration'
//
void AH_Model_v1ModelClass::Crossflowdragtrapezoidali_Reset(
    real_T* memory2_PreviousInput, real_T* memory1_PreviousInput) {
  // InitializeConditions for Memory: '<S23>/memory2'
  *memory2_PreviousInput = 0.0;

  // InitializeConditions for Memory: '<S23>/memory1'
  *memory1_PreviousInput = 0.0;
}

//
// Output and update for iterator system:
//    '<S19>/Cross-flow drag trapezoidal integration'
//    '<S103>/Cross-flow drag trapezoidal integration'
//    '<S187>/Cross-flow drag trapezoidal integration'
//    '<S271>/Cross-flow drag trapezoidal integration'
//
void AH_Model_v1ModelClass::Crossflowdragtrapezoidalintegra(
    real_T rtu_N, real_T rtu_dx, real_T rtu_v_r, real_T rtu_r, real_T* rty_sum1,
    real_T* rty_sum2, real_T rtp_Lpp) {
  int32_T tmp;
  int32_T i;
  real_T rtb_x;
  real_T rtb_memory1;
  real_T memory2_PreviousInput;
  real_T memory1_PreviousInput;

  // Outputs for Iterator SubSystem: '<S19>/Cross-flow drag trapezoidal
  // integration' incorporates:
  //   ForIterator: '<S23>/For Iterator'

  Crossflowdragtrapezoidali_Reset(&memory2_PreviousInput,
                                  &memory1_PreviousInput);
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
    // Sum: '<S23>/Sum3' incorporates:
    //   Constant: '<S23>/Lpp//2'
    //   Constant: '<S23>/constant'
    //   Product: '<S23>/Product1'
    //   Sum: '<S23>/Sum4'

    rtb_x = (static_cast<real_T>(i) - 1.0) * rtu_dx - rtp_Lpp / 2.0;

    // Sum: '<S23>/Sum1' incorporates:
    //   Product: '<S23>/x * r'

    rtb_memory1 = rtu_r * rtb_x + rtu_v_r;

    // If: '<S23>/If i=1 or i=N' incorporates:
    //   Abs: '<S23>/Abs'
    //   Inport: '<S25>/In1'
    //   Product: '<S23>/Product'
    //   Product: '<S23>/Product3'

    if ((i == 1) || (i == rtu_N)) {
      // Outputs for IfAction SubSystem: '<S23>/multiply with 0.5' incorporates:
      //   ActionPort: '<S24>/Action Port'

      // Gain: '<S24>/Gain' incorporates:
      //   Abs: '<S23>/Abs'
      //   Product: '<S23>/Product'
      //   Product: '<S23>/Product3'

      rtb_memory1 = 0.5 * (std::abs(rtb_memory1) * rtb_memory1 * rtu_dx);

      // End of Outputs for SubSystem: '<S23>/multiply with 0.5'
    } else {
      // Outputs for IfAction SubSystem: '<S23>/multiply with 1' incorporates:
      //   ActionPort: '<S25>/Action Port'

      rtb_memory1 = std::abs(rtb_memory1) * rtb_memory1 * rtu_dx;

      // End of Outputs for SubSystem: '<S23>/multiply with 1'
    }

    // End of If: '<S23>/If i=1 or i=N'

    // Sum: '<S23>/Sum2' incorporates:
    //   Memory: '<S23>/memory2'
    //   Product: '<S23>/Product4'

    *rty_sum2 = rtb_memory1 * rtb_x + memory2_PreviousInput;

    // Sum: '<S23>/Sum' incorporates:
    //   Memory: '<S23>/memory1'

    *rty_sum1 = memory1_PreviousInput + rtb_memory1;

    // Update for Memory: '<S23>/memory2'
    memory2_PreviousInput = *rty_sum2;

    // Update for Memory: '<S23>/memory1'
    memory1_PreviousInput = *rty_sum1;
  }

  // End of Outputs for SubSystem: '<S19>/Cross-flow drag trapezoidal
  // integration'
}

//
// Output and update for atomic system:
//    '<S64>/Chart'
//    '<S148>/Chart'
//    '<S232>/Chart'
//    '<S316>/Chart'
//
void AH_Model_v1ModelClass::AHV_Model_Chart(
    boolean_T rtu_hold, real_T rtu_x_ref, real_T rtu_y_ref, real_T rtu_x_hold,
    real_T rtu_y_hold, real_T* rty_x_ref_rel, real_T* rty_y_ref_rel,
    DW_Chart_AHV_Model_T* localDW) {
  // Chart: '<S64>/Chart'
  if (localDW->is_active_c129_AHV_Model == 0U) {
    localDW->is_active_c129_AHV_Model = 1U;
    if (rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_hold;
      localDW->x_local = rtu_x_hold;
      localDW->y_local = rtu_y_hold;
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    } else {
      localDW->is_c129_AHV_Model = AHV_Model_IN_unhold;
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    }
  } else if (localDW->is_c129_AHV_Model == AHV_Model_IN_hold) {
    if (!rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_unhold;
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    } else {
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    }
  } else {
    // case IN_unhold:
    if (rtu_hold) {
      localDW->is_c129_AHV_Model = AHV_Model_IN_hold;
      localDW->x_local = rtu_x_hold;
      localDW->y_local = rtu_y_hold;
      *rty_x_ref_rel = localDW->x_local;
      *rty_y_ref_rel = localDW->y_local;
    } else {
      *rty_x_ref_rel = rtu_x_ref;
      *rty_y_ref_rel = rtu_y_ref;
    }
  }

  // End of Chart: '<S64>/Chart'
}

//
// Output and update for action system:
//    '<S76>/If Action Subsystem'
//    '<S160>/If Action Subsystem'
//    '<S244>/If Action Subsystem'
//    '<S328>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem(
    real_T rtu_In1, real_T rtu_In1_a, real_T rtu_In1_d, real_T rty_Kp[3],
    const real_T rtp_Kp_A[9]) {
  int32_T i;

  // Gain: '<S81>/Kp_u' incorporates:
  //   SignalConversion generated from: '<S81>/Kp_u'

  for (i = 0; i < 3; i++) {
    rty_Kp[i] = 0.0;
    rty_Kp[i] += rtp_Kp_A[i] * rtu_In1;
    rty_Kp[i] += rtp_Kp_A[i + 3] * rtu_In1_a;
    rty_Kp[i] += rtp_Kp_A[i + 6] * rtu_In1_d;
  }

  // End of Gain: '<S81>/Kp_u'
}

//
// Output and update for action system:
//    '<S76>/If Action Subsystem1'
//    '<S160>/If Action Subsystem1'
//    '<S244>/If Action Subsystem1'
//    '<S328>/If Action Subsystem1'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem1(
    real_T rtu_In1, real_T rtu_In1_g, real_T rtu_In1_l, real_T rty_Kp[3],
    const real_T rtp_Kp[9]) {
  int32_T i;

  // Gain: '<S82>/Kp_u' incorporates:
  //   SignalConversion generated from: '<S82>/Kp_u'

  for (i = 0; i < 3; i++) {
    rty_Kp[i] = 0.0;
    rty_Kp[i] += rtp_Kp[i] * rtu_In1;
    rty_Kp[i] += rtp_Kp[i + 3] * rtu_In1_g;
    rty_Kp[i] += rtp_Kp[i + 6] * rtu_In1_l;
  }

  // End of Gain: '<S82>/Kp_u'
}

//
// Output and update for action system:
//    '<S83>/If Action Subsystem'
//    '<S84>/If Action Subsystem'
//    '<S167>/If Action Subsystem'
//    '<S168>/If Action Subsystem'
//    '<S251>/If Action Subsystem'
//    '<S252>/If Action Subsystem'
//    '<S335>/If Action Subsystem'
//    '<S336>/If Action Subsystem'
//
void AH_Model_v1ModelClass::AHV_Model_IfActionSubsystem_l(
    real_T rtu_Increment, real_T rtu_main_T, real_T* rty_Out1,
    const ConstB_IfActionSubsystem_AH_p_T* localC) {
  // Switch: '<S85>/Switch'
  if (rtu_main_T > 0.0) {
    *rty_Out1 = rtu_Increment;
  } else {
    *rty_Out1 = localC->Gain;
  }

  // End of Switch: '<S85>/Switch'
}

real_T rt_atan2d_snf(real_T u0, real_T u1) {
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
                                     real_T y[6]) {
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
void AH_Model_v1ModelClass::step() {
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // set solver stop time
    rtsiSetSolverStopTime(&(&AHV_Model_M)->solverInfo,
                          (((&AHV_Model_M)->Timing.clockTick0 + 1) *
                           (&AHV_Model_M)->Timing.stepSize0));
  }  // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if (rtmIsMinorTimeStep((&AHV_Model_M))) {
    (&AHV_Model_M)->Timing.t[0] = rtsiGetT(&(&AHV_Model_M)->solverInfo);
  }

  {
    real_T rtb_T33;
    real_T rtb_K3_u[3];
    real_T rtb_x_dot;
    real_T rtb_Sum1_di[6];
    real_T rtb_psi_dot;
    real_T rtb_Switch_a;
    real_T rtb_K2_u[3];
    real_T rtb_Row3;
    real_T rtb_Switch1;
    real_T rtb_Gain2[3];
    int8_T rtAction;
    real_T rtb_TmpSignalConversionAtProd_k[9];
    real_T rtb_TmpSignalConversionAtProduc[9];
    real_T rtb_tau_WD[6];
    real_T rtb_Integrator[6];
    real_T rtb_Row1;
    real_T rtb_tau_WF[6];
    real_T rtb_Row1_e;
    real_T psi_mean_m;
    int32_T i;
    real_T rtb_K3_u_0[9];
    real_T tmp;
    real_T tmp_0;
    real_T tmp_1;
    real_T tmp_2;
    real_T tmp_3;
    real_T rtb_Sum5_0;
    int32_T i_0;
    real_T rtb_K3_u_tmp;
    real_T rtb_K2_u_tmp;
    if (write_local == 0) {
      write_local = 1;
      fout.open("Thrust.txt");

      fout << "x "
           << "\t";
      fout << "y "
           << "\t";
      fout << "heave "
           << "\t";
      fout << "roll "
           << "\t";
      fout << "pitch "
           << "\t";
      fout << "yaw "
           << "\t";

      fout << "velocities1 "
           << "\t";
      fout << "velocities2 "
           << "\t";
      fout << "velocities3 "
           << "\t";
      fout << "velocities4 "
           << "\t";
      fout << "velocities5 "
           << "\t";
      fout << "velocities6 "
           << "\t";
    }

    //��
    //./5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth/5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth.fst
    // string path = "../openfast
    // /home/zzzzzsh/Documents/source/openfast2/reg_tests/r-test/glue-codes/openfast/5MW_OC4Semi_Linear/ATH_SET.txt";
    string path = "ATH_SET.txt";
    string exePath[8];
    for (int i = 1; i < 8; i++) {
      exePath[i] = ReadLine(path, i);
    }
    if (fileload == 1)  //���سɹ��ٸ�ֵ
    {
      entries = split(exePath[2], ' ');
      Featured = atof(entries[1].c_str());

      entries = split(exePath[3], ' ');
      SidePush = atof(entries[1].c_str());

      entries = split(exePath[5], ' ');
      MainThrust = atof(entries[1].c_str());

      entries = split(exePath[6], ' ');
      FuThrust = atof(entries[1].c_str());

      entries = split(exePath[7], ' ');
      AzimuthThrust = atof(entries[1].c_str());
    } else {
      Featured = 1;
      SidePush = 1;
      MainThrust = 1;
      FuThrust = 1;
      AzimuthThrust = 1;
    }

    for (i = 0; i < 6; i++) {
      // Integrator: '<S13>/Integrator1' incorporates:
      //   Inport: '<Root>/Vessel_init1'

      if (AHV_Model_DW.Integrator1_IWORK != 0) {
        AHV_Model_X.Integrator1_CSTATE[i] = Vessel_init1[i];
      }

      // Outport: '<Root>/eta_AHV1' incorporates:
      //   Integrator: '<S13>/Integrator1'

      eta_AHV1[i] = AHV_Model_X.Integrator1_CSTATE[i];

      // Outport: '<Root>/nu1' incorporates:
      //   Integrator: '<S13>/Integrator'

      nu1[i] = AHV_Model_X.Integrator_CSTATE[i];
    }
    fout << eta_AHV1[0] << "\t";
    fout << eta_AHV1[1] << "\t";
    fout << eta_AHV1[2] << "\t";
    fout << eta_AHV1[3] << "\t";
    fout << eta_AHV1[4] << "\t";
    fout << eta_AHV1[5] << "\t";

    fout << nu1[0] << "\t";
    fout << nu1[1] << "\t";
    fout << nu1[2] << "\t";
    fout << nu1[3] << "\t";
    fout << nu1[4] << "\t";
    fout << nu1[5] << "\n";
    // Outport: '<Root>/AHV_speed1' incorporates:
    //   Gain: '<S16>/Gain'
    //   Gain: '<S16>/Gain1'
    //   Rounding: '<S16>/Rounding Function'
    //   TransferFcn: '<S16>/Transfer Fcn'

    AHV_speed1 = std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE * 10.0) * 0.1;
    // cout << AHV_speed1 << endl;
    //  Gain: '<S12>/Gain' incorporates:
    //    Inport: '<Root>/Current_direction'

    rtb_T33 = 0.017453292519943295 * Current_direction;

    // Trigonometry: '<S61>/sin(theta)' incorporates:
    //   Fcn: '<S22>/T21 '
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/sin(theta)'

    rtb_x_dot = std::sin(AHV_Model_X.Integrator1_CSTATE[3]);

    // Trigonometry: '<S61>/cos(theta)' incorporates:
    //   Fcn: '<S22>/T31 '
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/cos(theta)'

    psi_mean_m = std::cos(AHV_Model_X.Integrator1_CSTATE[3]);

    // Trigonometry: '<S61>/sin(theta)' incorporates:
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/sin(theta)'

    rtb_psi_dot = std::sin(AHV_Model_X.Integrator1_CSTATE[4]);

    // Trigonometry: '<S61>/cos(theta)' incorporates:
    //   Fcn: '<S22>/T23'
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/cos(theta)'

    rtb_Row1_e = std::cos(AHV_Model_X.Integrator1_CSTATE[4]);

    // Trigonometry: '<S61>/sin(theta)' incorporates:
    //   Fcn: '<S68>/Fcn'
    //   Fcn: '<S68>/Fcn1'
    //   Fcn: '<S69>/Row1'
    //   Fcn: '<S70>/Row1'
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/sin(theta)'

    rtb_K3_u_tmp = std::sin(AHV_Model_X.Integrator1_CSTATE[5]);

    // Trigonometry: '<S61>/cos(theta)' incorporates:
    //   Fcn: '<S68>/Fcn'
    //   Fcn: '<S68>/Fcn1'
    //   Fcn: '<S69>/Row1'
    //   Fcn: '<S69>/Row2'
    //   Fcn: '<S70>/Row1'
    //   Fcn: '<S70>/Row2'
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S21>/cos(theta)'

    rtb_K2_u_tmp = std::cos(AHV_Model_X.Integrator1_CSTATE[5]);

    // Product: '<S9>/Product1' incorporates:
    //   Fcn: '<S61>/R11'
    //   Fcn: '<S61>/R12'
    //   Fcn: '<S61>/R13'
    //   Fcn: '<S61>/R21 '
    //   Fcn: '<S61>/R22'
    //   Fcn: '<S61>/R23'
    //   Fcn: '<S61>/R31 '
    //   Fcn: '<S61>/R32'
    //   Fcn: '<S61>/R33'
    //   Inport: '<Root>/ahv_fairlead1'
    //   Trigonometry: '<S61>/cos(theta)'
    //   Trigonometry: '<S61>/sin(theta)'

    rtb_K3_u_0[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[3] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[6] = -rtb_psi_dot;
    rtb_K3_u_0[1] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_x_dot + -rtb_K3_u_tmp * psi_mean_m;
    rtb_K3_u_0[4] =
        rtb_K3_u_tmp * rtb_psi_dot * rtb_x_dot + rtb_K2_u_tmp * psi_mean_m;
    rtb_K3_u_0[7] = rtb_Row1_e * rtb_x_dot;
    rtb_K3_u_0[2] =
        rtb_K2_u_tmp * rtb_psi_dot * psi_mean_m + rtb_K3_u_tmp * rtb_x_dot;
    rtb_K3_u_0[5] =
        rtb_K3_u_tmp * rtb_psi_dot * psi_mean_m + -rtb_K2_u_tmp * rtb_x_dot;
    rtb_K3_u_0[8] = rtb_Row1_e * psi_mean_m;
    for (i = 0; i < 3; i++) {
      rtb_Switch_a = rtb_K3_u_0[i + 6] * ahv_fairlead1[2] +
                     (rtb_K3_u_0[i + 3] * ahv_fairlead1[1] +
                      rtb_K3_u_0[i] * ahv_fairlead1[0]);

      // Gain: '<S8>/Gain1' incorporates:
      //   Inport: '<Root>/ahv_fairlead1'

      rtb_tau_WD[i] = AHV_Model_ConstP.pooled60[i] * rtb_Switch_a;
      rtb_K3_u[i] = rtb_Switch_a;
    }

    // End of Product: '<S9>/Product1'

    // Gain: '<S8>/Gain1' incorporates:
    //   Inport: '<Root>/tau_cable1'
    //   Product: '<S60>/ae'
    //   Product: '<S60>/af'
    //   Product: '<S60>/bd'
    //   Product: '<S60>/bf'
    //   Product: '<S60>/cd'
    //   Product: '<S60>/ce'
    //   Sum: '<S60>/Sum'
    //   Sum: '<S60>/Sum1'
    //   Sum: '<S60>/Sum2'

    rtb_tau_WD[3] = tau_cable1[1] * rtb_K3_u[2] - tau_cable1[2] * rtb_K3_u[1];
    rtb_tau_WD[4] =
        (tau_cable1[2] * rtb_K3_u[0] - tau_cable1[0] * rtb_K3_u[2]) * 0.5;
    rtb_tau_WD[5] = tau_cable1[0] * rtb_K3_u[1] - tau_cable1[1] * rtb_K3_u[0];

    // SignalConversion generated from: '<S18>/Product' incorporates:
    //   Fcn: '<S21>/R11'
    //   Fcn: '<S21>/R12'
    //   Fcn: '<S21>/R21 '
    //   Fcn: '<S21>/R31 '

    rtb_TmpSignalConversionAtProduc[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProduc[1] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProduc[2] = -rtb_psi_dot;
    rtb_TmpSignalConversionAtProduc[3] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_x_dot - rtb_K3_u_tmp * psi_mean_m;

    // Fcn: '<S21>/R22' incorporates:
    //   Fcn: '<S21>/R23'

    rtb_Switch_a = rtb_K3_u_tmp * rtb_psi_dot;

    // SignalConversion generated from: '<S18>/Product' incorporates:
    //   Fcn: '<S21>/R13'
    //   Fcn: '<S21>/R22'
    //   Fcn: '<S21>/R23'
    //   Fcn: '<S21>/R32'
    //   Fcn: '<S21>/R33'

    rtb_TmpSignalConversionAtProduc[4] =
        rtb_Switch_a * rtb_x_dot + rtb_K2_u_tmp * psi_mean_m;
    rtb_TmpSignalConversionAtProduc[5] = rtb_Row1_e * rtb_x_dot;
    rtb_TmpSignalConversionAtProduc[6] =
        psi_mean_m * rtb_psi_dot * rtb_K2_u_tmp + rtb_K3_u_tmp * rtb_x_dot;
    rtb_TmpSignalConversionAtProduc[7] =
        rtb_Switch_a * psi_mean_m - rtb_K2_u_tmp * rtb_x_dot;
    rtb_TmpSignalConversionAtProduc[8] = rtb_Row1_e * psi_mean_m;

    // Sum: '<S13>/Sum6' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Integrator: '<S13>/Integrator'
    //   Product: '<S12>/Product'
    //   Product: '<S12>/Product1'
    //   Trigonometry: '<S12>/cos'
    //   Trigonometry: '<S12>/sin'

    AHV_Model_B.nu_r[0] =
        AHV_Model_X.Integrator_CSTATE[0] - std::cos(rtb_T33) * Current_speed;
    AHV_Model_B.nu_r[1] =
        AHV_Model_X.Integrator_CSTATE[1] - std::sin(rtb_T33) * Current_speed;
    AHV_Model_B.nu_r[2] = AHV_Model_X.Integrator_CSTATE[2];
    AHV_Model_B.nu_r[3] = AHV_Model_X.Integrator_CSTATE[3];
    AHV_Model_B.nu_r[4] = AHV_Model_X.Integrator_CSTATE[4];
    AHV_Model_B.nu_r[5] = AHV_Model_X.Integrator_CSTATE[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S19>/Cross-flow drag trapezoidal
      // integration' Constant: '<S19>/N'
      Crossflowdragtrapezoidalintegra(
          20.0, AHV_Model_ConstB.dx, AHV_Model_B.nu_r[1], AHV_Model_B.nu_r[5],
          &AHV_Model_B.Sum_n, &AHV_Model_B.Sum2_pj, 82.800003);

      // End of Outputs for SubSystem: '<S19>/Cross-flow drag trapezoidal
      // integration'

      // Product: '<S19>/dx1' incorporates:
      //   Constant: '<S19>/2D drag coefficient '
      //   Constant: '<S19>/Transversal area//Lpp'
      //   Constant: '<S19>/rho'
      //   Gain: '<S19>/Gain1'

      AHV_Model_B.dx1 =
          -0.5 * AHV_Model_B.Sum_n * 1025.0 * 5.4 * 0.69954741847648816;

      // Product: '<S19>/dx2' incorporates:
      //   Constant: '<S19>/2D drag coefficient '
      //   Constant: '<S19>/Transversal area//Lpp'
      //   Constant: '<S19>/rho'
      //   Gain: '<S19>/Gain'

      AHV_Model_B.dx2 =
          -0.5 * AHV_Model_B.Sum2_pj * 1025.0 * 5.4 * 0.69954741847648816;
    }

    // Product: '<S13>/G*eta' incorporates:
    //   Constant: '<S13>/Spring stiffness'
    //   Integrator: '<S13>/Integrator1'

    for (i = 0; i < 6; i++) {
      rtb_Sum1_di[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Sum1_di[i] += AHV_Model_ConstP.pooled48[6 * i_0 + i] *
                          AHV_Model_X.Integrator1_CSTATE[i_0];
      }
    }

    // End of Product: '<S13>/G*eta'

    // Switch: '<S6>/Switch' incorporates:
    //   Constant: '<S6>/Constant'
    //   Constant: '<S6>/Constant1'
    //   Inport: '<Root>/Hold_Position1'

    if (Hold_Position1) {
      rtb_T33 = 2.0;
    } else {
      rtb_T33 = 1.0;
    }

    // End of Switch: '<S6>/Switch'

    // Integrator: '<S63>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init1'
    //   SignalConversion generated from: '<S63>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK != 0) {
      AHV_Model_X.Integrator3_CSTATE[0] = Vessel_init1[0];
      AHV_Model_X.Integrator3_CSTATE[1] = Vessel_init1[1];
      AHV_Model_X.Integrator3_CSTATE[2] = Vessel_init1[5];
    }

    // Saturate: '<S71>/x_Saturation' incorporates:
    //   Integrator: '<S63>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE[2] > 1.0E+10) {
      rtb_psi_dot = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE[2] < -1.0E+10) {
      rtb_psi_dot = -1.0E+10;
    } else {
      rtb_psi_dot = AHV_Model_X.Integrator3_CSTATE[2];
    }

    // End of Saturate: '<S71>/x_Saturation'

    // Signum: '<S71>/x_Sign'
    if (rtb_psi_dot < 0.0) {
      rtb_Switch_a = -1.0;
    } else if (rtb_psi_dot > 0.0) {
      rtb_Switch_a = 1.0;
    } else if (rtb_psi_dot == 0.0) {
      rtb_Switch_a = 0.0;
    } else {
      rtb_Switch_a = (rtNaN);
    }

    // End of Signum: '<S71>/x_Sign'

    // Gain: '<S71>/pi'
    rtb_Switch_a *= 3.1415926535897931;

    // Sum: '<S71>/Sum1'
    rtb_psi_dot += rtb_Switch_a;

    // Math: '<S71>/Math Function' incorporates:
    //   Constant: '<S71>/Constant'

    rtb_psi_dot = rt_remd_snf(rtb_psi_dot, 6.2831853071795862);

    // Sum: '<S71>/Sum'
    rtb_psi_dot -= rtb_Switch_a;

    // Signum: '<S80>/x_Sign' incorporates:
    //   Fcn: '<S65>/yaw_angle'

    if (rtb_psi_dot < 0.0) {
      rtb_Switch_a = -1.0;
    } else if (rtb_psi_dot > 0.0) {
      rtb_Switch_a = 1.0;
    } else if (rtb_psi_dot == 0.0) {
      rtb_Switch_a = 0.0;
    } else {
      rtb_Switch_a = (rtNaN);
    }

    // End of Signum: '<S80>/x_Sign'

    // Gain: '<S80>/pi'
    rtb_Switch1 = 3.1415926535897931 * rtb_Switch_a;

    // Sum: '<S80>/Sum' incorporates:
    //   Constant: '<S80>/Constant'
    //   Fcn: '<S65>/yaw_angle'
    //   Math: '<S80>/Math Function'
    //   Sum: '<S80>/Sum1'

    rtb_Switch_a = rt_remd_snf(rtb_psi_dot + rtb_Switch1, 6.2831853071795862) -
                   rtb_Switch1;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S64>/Chart' incorporates:
      //   Inport: '<Root>/Hold_Position1'
      //   Inport: '<Root>/Vessel_X_Ref1'
      //   Inport: '<Root>/Vessel_Y_Ref1'
      //   Integrator: '<S13>/Integrator1'

      AHV_Model_Chart(Hold_Position1, Vessel_X_Ref1, Vessel_Y_Ref1,
                      AHV_Model_X.Integrator1_CSTATE[0],
                      AHV_Model_X.Integrator1_CSTATE[1],
                      &AHV_Model_B.x_ref_rel_eg, &AHV_Model_B.y_ref_rel_no,
                      &AHV_Model_DW.sf_Chart);
    }

    // Switch: '<S64>/Switch' incorporates:
    //   Gain: '<S67>/Gain'
    //   Inport: '<Root>/heading_angle_ref1'
    //   Inport: '<Root>/heading_mode1'
    //   Integrator: '<S13>/Integrator1'
    //   Trigonometry: '<S62>/Trigonometric Function'

    if (heading_mode1 > 0.5) {
      AHV_Model_B.RateLimiter =
          57.295779513082323 * rt_atan2d_snf(AHV_Model_X.Integrator1_CSTATE[1],
                                             AHV_Model_X.Integrator1_CSTATE[0]);
    } else {
      AHV_Model_B.RateLimiter = heading_angle_ref1;
    }

    // End of Switch: '<S64>/Switch'

    // RateLimiter: '<S64>/Rate Limiter'
    if (!(AHV_Model_DW.LastMajorTime == (rtInf))) {
      rtb_Row3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime;
      rtb_Row1 = rtb_Row3 * 10.0;
      rtb_Switch1 = AHV_Model_B.RateLimiter - AHV_Model_DW.PrevY;
      if (rtb_Switch1 > rtb_Row1) {
        AHV_Model_B.RateLimiter = AHV_Model_DW.PrevY + rtb_Row1;
      } else {
        rtb_Row3 *= -10.0;
        if (rtb_Switch1 < rtb_Row3) {
          AHV_Model_B.RateLimiter = AHV_Model_DW.PrevY + rtb_Row3;
        }
      }
    }

    // End of RateLimiter: '<S64>/Rate Limiter'

    // Gain: '<S74>/Gain1'
    rtb_Switch1 = 0.017453292519943295 * AHV_Model_B.RateLimiter;

    // Saturate: '<S75>/x_Saturation'
    if (rtb_Switch1 > 1.0E+10) {
      rtb_Switch1 = 1.0E+10;
    } else {
      if (rtb_Switch1 < -1.0E+10) {
        rtb_Switch1 = -1.0E+10;
      }
    }

    // End of Saturate: '<S75>/x_Saturation'

    // Signum: '<S75>/x_Sign'
    if (rtb_Switch1 < 0.0) {
      rtb_Row3 = -1.0;
    } else if (rtb_Switch1 > 0.0) {
      rtb_Row3 = 1.0;
    } else if (rtb_Switch1 == 0.0) {
      rtb_Row3 = 0.0;
    } else {
      rtb_Row3 = (rtNaN);
    }

    // End of Signum: '<S75>/x_Sign'

    // Gain: '<S75>/pi'
    rtb_Row3 *= 3.1415926535897931;

    // Sum: '<S75>/Sum1'
    rtb_Switch1 += rtb_Row3;

    // Math: '<S75>/Math Function' incorporates:
    //   Constant: '<S75>/Constant'

    rtb_Switch1 = rt_remd_snf(rtb_Switch1, 6.2831853071795862);

    // Sum: '<S75>/Sum'
    rtb_Switch1 -= rtb_Row3;

    // Sum: '<S65>/Sum2' incorporates:
    //   Integrator: '<S63>/Integrator3'

    rtb_K2_u[0] = AHV_Model_B.x_ref_rel_eg - AHV_Model_X.Integrator3_CSTATE[0];
    rtb_K2_u[1] = AHV_Model_B.y_ref_rel_no - AHV_Model_X.Integrator3_CSTATE[1];

    // Saturate: '<S79>/x_Saturation' incorporates:
    //   Sum: '<S65>/Sum2'

    rtb_Switch1 -= rtb_psi_dot;

    // Signum: '<S79>/x_Sign' incorporates:
    //   Saturate: '<S79>/x_Saturation'

    if (rtb_Switch1 < 0.0) {
      rtb_Row3 = -1.0;
    } else if (rtb_Switch1 > 0.0) {
      rtb_Row3 = 1.0;
    } else if (rtb_Switch1 == 0.0) {
      rtb_Row3 = 0.0;
    } else {
      rtb_Row3 = (rtNaN);
    }

    // End of Signum: '<S79>/x_Sign'

    // Gain: '<S79>/pi'
    rtb_Row3 *= 3.1415926535897931;

    // Sum: '<S79>/Sum' incorporates:
    //   Constant: '<S79>/Constant'
    //   Math: '<S79>/Math Function'
    //   Saturate: '<S79>/x_Saturation'
    //   Sum: '<S79>/Sum1'

    rtb_Switch1 =
        rt_remd_snf(rtb_Switch1 + rtb_Row3, 6.2831853071795862) - rtb_Row3;

    // SignalConversion generated from: '<S65>/Integrator' incorporates:
    //   Fcn: '<S78>/Row3'

    AHV_Model_B.TmpSignalConversionAtIntegrat_p[2] = rtb_Switch1;

    // Fcn: '<S78>/Row1' incorporates:
    //   Fcn: '<S78>/Row2'

    rtb_Row3 = std::sin(rtb_Switch_a);
    rtb_Switch_a = std::cos(rtb_Switch_a);
    rtb_Row1 = rtb_Switch_a * rtb_K2_u[0] + rtb_Row3 * rtb_K2_u[1];

    // If: '<S83>/If' incorporates:
    //   Abs: '<S83>/Abs'
    //   Constant: '<Root>/Constant41'
    //   Inport: '<S86>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem =
          static_cast<int8_T>(!(std::abs(rtb_Row1) <= Featured));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem) {
      case 0:
        // Outputs for IfAction SubSystem: '<S83>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S86>/Action Port'

        AHV_Model_B.Merge = rtb_Row1;

        // End of Outputs for SubSystem: '<S83>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S83>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S85>/Action Port'

        AHV_Model_IfActionSubsystem_l(1.0, rtb_Row1, &AHV_Model_B.Merge,
                                      &AHV_Model_ConstB.IfActionSubsystem_l);

        // End of Outputs for SubSystem: '<S83>/If Action Subsystem'
        break;
    }

    // End of If: '<S83>/If'

    // Fcn: '<S78>/Row2'
    rtb_Switch_a = -rtb_Row3 * rtb_K2_u[0] + rtb_Switch_a * rtb_K2_u[1];

    // If: '<S84>/If' incorporates:
    //   Abs: '<S84>/Abs'
    //   Constant: '<Root>/Constant42'
    //   Inport: '<S88>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_a =
          static_cast<int8_T>(!(std::abs(rtb_Switch_a) <= SidePush));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_a) {
      case 0:
        // Outputs for IfAction SubSystem: '<S84>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S88>/Action Port'

        AHV_Model_B.Merge_m = rtb_Switch_a;

        // End of Outputs for SubSystem: '<S84>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S84>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S87>/Action Port'

        AHV_Model_IfActionSubsystem_l(10.0, rtb_Switch_a, &AHV_Model_B.Merge_m,
                                      &AHV_Model_ConstB.IfActionSubsystem_m);

        // End of Outputs for SubSystem: '<S84>/If Action Subsystem'
        break;
    }

    // End of If: '<S84>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // If: '<S76>/If' incorporates:
      //   Fcn: '<S78>/Row3'

      if (rtmIsMajorTimeStep((&AHV_Model_M))) {
        rtAction = 1;
        AHV_Model_DW.If_ActiveSubsystem_i = 1;
      } else {
        rtAction = AHV_Model_DW.If_ActiveSubsystem_i;
      }

      switch (rtAction) {
        case 0:
          // Outputs for IfAction SubSystem: '<S76>/If Action Subsystem'
          // incorporates:
          //   ActionPort: '<S81>/Action Port'

          AHV_Model_IfActionSubsystem(AHV_Model_B.Merge, AHV_Model_B.Merge_m,
                                      rtb_Switch1, AHV_Model_B.Kp,
                                      AHV_Model_ConstP.pooled2);

          // End of Outputs for SubSystem: '<S76>/If Action Subsystem'
          break;

        case 1:
          // Outputs for IfAction SubSystem: '<S76>/If Action Subsystem1'
          // incorporates:
          //   ActionPort: '<S82>/Action Port'

          AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge, AHV_Model_B.Merge_m,
                                       rtb_Switch1, AHV_Model_B.Kp,
                                       AHV_Model_ConstP.pooled1);

          // End of Outputs for SubSystem: '<S76>/If Action Subsystem1'
          break;
      }

      // End of If: '<S76>/If'
    }

    for (i = 0; i < 3; i++) {
      // Sum: '<S65>/Sum1' incorporates:
      //   Gain: '<S65>/Kd'
      //   Gain: '<S65>/Ki'
      //   Integrator: '<S63>/Integrator4'
      //   Integrator: '<S65>/Integrator'
      //   Sum: '<S65>/Sum3'

      rtb_K2_u[i] =
          (((0.0 * AHV_Model_X.Integrator_CSTATE_b[0] +
             0.0 * AHV_Model_X.Integrator_CSTATE_b[1]) +
            0.0 * AHV_Model_X.Integrator_CSTATE_b[2]) +
           AHV_Model_B.Kp[i]) -
          (AHV_Model_ConstP.pooled64[i + 6] *
               AHV_Model_X.Integrator4_CSTATE[2] +
           (AHV_Model_ConstP.pooled64[i + 3] *
                AHV_Model_X.Integrator4_CSTATE[1] +
            AHV_Model_ConstP.pooled64[i] * AHV_Model_X.Integrator4_CSTATE[0]));
    }

    // Saturate: '<S66>/Surge Force Saturation'
    if (rtb_K2_u[0] > 2.0E+6) {
      rtb_Switch_a = 2.0E+6;
    } else if (rtb_K2_u[0] < -2.0E+6) {
      rtb_Switch_a = -2.0E+6;
    } else {
      rtb_Switch_a = rtb_K2_u[0];
    }

    // End of Saturate: '<S66>/Surge Force Saturation'

    // Saturate: '<S66>/Sway Force Saturation'
    if (rtb_K2_u[1] > 1.5E+6) {
      rtb_Switch1 = 1.5E+6;
    } else if (rtb_K2_u[1] < -1.5E+6) {
      rtb_Switch1 = -1.5E+6;
    } else {
      rtb_Switch1 = rtb_K2_u[1];
    }

    // End of Saturate: '<S66>/Sway Force Saturation'

    // Saturate: '<S66>/Yaw Moment Saturation'
    if (rtb_K2_u[2] > 2.0E+7) {
      rtb_Row3 = 2.0E+7;
    } else if (rtb_K2_u[2] < -2.0E+7) {
      rtb_Row3 = -2.0E+7;
    } else {
      rtb_Row3 = rtb_K2_u[2];
    }

    // End of Saturate: '<S66>/Yaw Moment Saturation'

    // Sum: '<S13>/Sum1' incorporates:
    //   Constant: '<S10>/Constant'
    //   Product: '<S6>/Product'

    rtb_Integrator[0] = rtb_T33 * rtb_Switch_a - rtb_Sum1_di[0];
    rtb_Integrator[1] = rtb_T33 * rtb_Switch1 - rtb_Sum1_di[1];
    rtb_Integrator[2] = 0.0 - rtb_Sum1_di[2];
    rtb_Integrator[3] = 0.0 - rtb_Sum1_di[3];
    rtb_Integrator[4] = 0.0 - rtb_Sum1_di[4];
    rtb_Integrator[5] = rtb_T33 * rtb_Row3 - rtb_Sum1_di[5];
    for (i = 0; i < 6; i++) {
      rtb_Sum1_di[i] = rtb_Integrator[i];
    }

    // End of Sum: '<S13>/Sum1'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S8>/Subsystem3'
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

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread,
                depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff,
                rand_freq, rand_dir, AHV_Model_B.Zeta_a_k, AHV_Model_B.Omega_b,
                AHV_Model_B.Phase_kf, AHV_Model_B.Wavenum_d, AHV_Model_B.Psi_l,
                &rtb_T33, &AHV_Model_B.Subsystem3, &AHV_Model_DW.Subsystem3,
                20.0, 10.0);

      // End of Outputs for SubSystem: '<S8>/Subsystem3'
    }

    // Outputs for Atomic SubSystem: '<S8>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_X.Integrator1_CSTATE, AHV_Model_B.Psi_l,
                   AHV_Model_B.Wavenum_d, AHV_Model_B.Omega_b,
                   AHV_Model_B.Phase_kf, AHV_Model_B.Zeta_a_k, rtb_tau_WF,
                   rtb_Integrator, &AHV_Model_B.WaveloadsU0,
                   &AHV_Model_ConstB.WaveloadsU0, &AHV_Model_DW.WaveloadsU0);

    // End of Outputs for SubSystem: '<S8>/Wave loads (U=0)'

    // Sum: '<S13>/Sum7' incorporates:
    //   Sum: '<S13>/Sum3'
    //   Sum: '<S8>/Sum'

    for (i = 0; i < 6; i++) {
      rtb_tau_WD[i] += (rtb_tau_WF[i] + rtb_Integrator[i]) + rtb_Sum1_di[i];
    }

    // End of Sum: '<S13>/Sum7'

    // Sum: '<S13>/Sum5' incorporates:
    //   Abs: '<S19>/Abs1'
    //   Gain: '<S19>/Gain3'
    //   Product: '<S19>/Product7'
    //   Product: '<S19>/Product8'

    rtb_Row1 = rtb_tau_WD[0];
    tmp = rtb_tau_WD[1];
    tmp_0 = rtb_tau_WD[2];
    tmp_1 = rtb_tau_WD[3];
    tmp_2 = rtb_tau_WD[4];
    tmp_3 = rtb_tau_WD[5];
    rtb_tau_WD[0] = -(std::abs(AHV_Model_B.nu_r[0]) * AHV_Model_B.nu_r[0] *
                      AHV_Model_ConstB.Product5) +
                    rtb_Row1;
    rtb_tau_WD[1] = tmp + AHV_Model_B.dx1;
    rtb_tau_WD[2] = tmp_0;
    rtb_tau_WD[3] = tmp_1;
    rtb_tau_WD[4] = tmp_2;
    rtb_tau_WD[5] = tmp_3 + AHV_Model_B.dx2;

    // Sum: '<S13>/Sum4' incorporates:
    //   Constant: '<S13>/damping'
    //   Product: '<S13>/D*eta '

    for (i = 0; i < 6; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Row1 +=
            AHV_Model_ConstP.pooled51[6 * i_0 + i] * AHV_Model_B.nu_r[i_0];
      }

      rtb_tau_WD[i] -= rtb_Row1;
    }

    // End of Sum: '<S13>/Sum4'

    // Product: '<S32>/Product2' incorporates:
    //   Integrator: '<S32>/Integrator'

    rtb_Sum5_0 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Sum5_0 +=
          AHV_Model_ConstB.MathFunction[i] * AHV_Model_X.Integrator_CSTATE_l[i];
    }

    // Sum: '<S13>/Sum' incorporates:
    //   Constant: '<S32>/B44_inf'
    //   Constant: '<S32>/D44'
    //   Product: '<S32>/Product2'
    //   Product: '<S32>/Product3'
    //   Product: '<S32>/Product4'
    //   StateSpace: '<S20>/Dp(1,1)'
    //   StateSpace: '<S20>/Dp(2,2)'
    //   StateSpace: '<S20>/Dp(2,4)'
    //   StateSpace: '<S20>/Dp(2,6)'
    //   StateSpace: '<S20>/Dp(3,3)'
    //   StateSpace: '<S20>/Dp(3,5)'
    //   StateSpace: '<S20>/Dp(4,2)'
    //   StateSpace: '<S20>/Dp(4,6)'
    //   StateSpace: '<S20>/Dp(5,3)'
    //   StateSpace: '<S20>/Dp(5,5)'
    //   StateSpace: '<S20>/Dp(6,2)'
    //   StateSpace: '<S20>/Dp(6,4)'
    //   StateSpace: '<S20>/Dp(6,6)'
    //   Sum: '<S20>/Sum'
    //   Sum: '<S20>/Sum1'
    //   Sum: '<S20>/Sum3'
    //   Sum: '<S20>/Sum4'
    //   Sum: '<S20>/Sum5'
    //   Sum: '<S32>/Sum1'
    //   Sum: '<S32>/Sum2'

    rtb_Row1 = rtb_tau_WD[0];
    tmp = rtb_tau_WD[1];
    tmp_0 = rtb_tau_WD[2];
    tmp_1 = rtb_tau_WD[3];
    tmp_2 = rtb_tau_WD[4];
    tmp_3 = rtb_tau_WD[5];
    rtb_tau_WD[0] =
        rtb_Row1 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE[0] +
                        -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE[1]) +
                       -26196.671588656758 * AHV_Model_X.Dp11_CSTATE[2]) +
                      -69550.629589629229 * AHV_Model_X.Dp11_CSTATE[3]) +
                     -6667.9011063139569 * AHV_Model_X.Dp11_CSTATE[4]) +
                    28360.700724836512 * AHV_Model_B.nu_r[0]);
    rtb_tau_WD[1] =
        tmp - (((((((-1.9219252802461751E+6 * AHV_Model_X.Dp22_CSTATE[0] +
                     -43506.445868266834 * AHV_Model_X.Dp22_CSTATE[1]) +
                    -261966.71588657616 * AHV_Model_X.Dp22_CSTATE[2]) +
                   -695506.29589629406 * AHV_Model_X.Dp22_CSTATE[3]) +
                  -66679.0110631567 * AHV_Model_X.Dp22_CSTATE[4]) +
                 283607.00724836387 * AHV_Model_B.nu_r[1]) +
                (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE[0] +
                     -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE[1]) +
                    782719.9224317756 * AHV_Model_X.Dp24_CSTATE[2]) +
                   906439.02406890341 * AHV_Model_X.Dp24_CSTATE[3]) +
                  48702.05056677026 * AHV_Model_X.Dp24_CSTATE[4]) +
                 342290.49224698928 * AHV_Model_B.nu_r[3])) +
               (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE[0] +
                    9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE[1]) +
                   1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE[2]) +
                  3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE[3]) +
                 -523575.51382008579 * AHV_Model_X.Dp26_CSTATE[4]) +
                1.4324880870688609E+6 * AHV_Model_B.nu_r[5]));
    rtb_tau_WD[2] =
        tmp_0 - ((((((-2.8230572248621574E+6 * AHV_Model_X.Dp33_CSTATE[0] +
                      2971.0437052436932 * AHV_Model_X.Dp33_CSTATE[1]) +
                     -112975.07052200216 * AHV_Model_X.Dp33_CSTATE[2]) +
                    -459867.38057091518 * AHV_Model_X.Dp33_CSTATE[3]) +
                   35755.997807249405 * AHV_Model_X.Dp33_CSTATE[4]) +
                  425986.43085235963 * AHV_Model_B.nu_r[2]) +
                 (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE[0] +
                      -189169.67892356269 * AHV_Model_X.Dp35_CSTATE[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r[4]));
    rtb_tau_WD[3] =
        tmp_1 - (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE[0] +
                       -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE[1]) +
                      782719.9224317756 * AHV_Model_X.Dp42_CSTATE[2]) +
                     906439.02406890341 * AHV_Model_X.Dp42_CSTATE[3]) +
                    48702.05056677026 * AHV_Model_X.Dp42_CSTATE[4]) +
                   342290.49224698928 * AHV_Model_B.nu_r[1]) +
                  ((2.3074903324953854E+7 * AHV_Model_B.nu_r[3] + rtb_Sum5_0) +
                   4.97124674352808E+7 * AHV_Model_B.nu_r[3])) +
                 (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE[0] +
                      2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE[1]) +
                     -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE[2]) +
                    1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE[3]) +
                   -70143.531785937739 * AHV_Model_X.Dp46_CSTATE[4]) +
                  -8.0522169282609783E+6 * AHV_Model_B.nu_r[5]));
    rtb_tau_WD[4] =
        tmp_2 - ((((((-4.3928932046187744E+7 * AHV_Model_X.Dp53_CSTATE[0] +
                      -189169.67892356269 * AHV_Model_X.Dp53_CSTATE[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp53_CSTATE[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp53_CSTATE[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp53_CSTATE[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r[2]) +
                 (((((-1.3510265271416564E+7 * AHV_Model_X.Dp55_CSTATE[0] +
                      2.1371969925573766E+9 * AHV_Model_X.Dp55_CSTATE[1]) +
                     9.8710379288756877E+7 * AHV_Model_X.Dp55_CSTATE[2]) +
                    6.1487736935078049E+8 * AHV_Model_X.Dp55_CSTATE[3]) +
                   5.7220438213245332E+8 * AHV_Model_X.Dp55_CSTATE[4]) +
                  3.0813057877131975E+8 * AHV_Model_B.nu_r[4]));
    rtb_tau_WD[5] =
        tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE[0] +
                       9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE[1]) +
                      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE[2]) +
                     3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE[3]) +
                    -523575.51382008579 * AHV_Model_X.Dp62_CSTATE[4]) +
                   1.4324880870688609E+6 * AHV_Model_B.nu_r[1]) +
                  (((((-1.5746048421550707E+6 * AHV_Model_X.Dp64_CSTATE[0] +
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

    // Product: '<S13>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2, rtb_tau_WD,
                                    AHV_Model_B.Minvtau);

    // Fcn: '<S22>/T21 ' incorporates:
    //   Fcn: '<S22>/T31 '
    //   Integrator: '<S13>/Integrator1'

    rtb_T33 = std::tan(AHV_Model_X.Integrator1_CSTATE[4]);

    // Product: '<S18>/Product1' incorporates:
    //   Constant: '<S22>/Constant'
    //   Constant: '<S22>/Constant '
    //   Fcn: '<S22>/T21 '
    //   Fcn: '<S22>/T23'
    //   Fcn: '<S22>/T31 '
    //   Fcn: '<S22>/T32'
    //   Fcn: '<S22>/T33'
    //   Reshape: '<S22>/Reshape 9x1->3x3'

    rtb_K3_u_0[0] = 1.0;
    rtb_K3_u_0[3] = rtb_x_dot * rtb_T33;
    rtb_K3_u_0[6] = psi_mean_m * rtb_T33;
    rtb_K3_u_0[1] = 0.0;
    rtb_K3_u_0[4] = psi_mean_m;
    rtb_K3_u_0[7] = -rtb_x_dot;
    rtb_K3_u_0[2] = 0.0;
    rtb_K3_u_0[5] = rtb_x_dot / rtb_Row1_e;
    rtb_K3_u_0[8] = psi_mean_m / rtb_Row1_e;
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S13>/Integrator1' incorporates:
      //   Integrator: '<S13>/Integrator'
      //   Product: '<S18>/Product'
      //   Product: '<S18>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegrat_o[i] =
          rtb_TmpSignalConversionAtProduc[i + 6] *
              AHV_Model_X.Integrator_CSTATE[2] +
          (rtb_TmpSignalConversionAtProduc[i + 3] *
               AHV_Model_X.Integrator_CSTATE[1] +
           rtb_TmpSignalConversionAtProduc[i] *
               AHV_Model_X.Integrator_CSTATE[0]);
      AHV_Model_B.TmpSignalConversionAtIntegrat_o[i + 3] =
          rtb_K3_u_0[i + 6] * AHV_Model_X.Integrator_CSTATE[5] +
          (rtb_K3_u_0[i + 3] * AHV_Model_X.Integrator_CSTATE[4] +
           rtb_K3_u_0[i] * AHV_Model_X.Integrator_CSTATE[3]);
    }

    // Sum: '<S32>/Sum' incorporates:
    //   Constant: '<S32>/A44'
    //   Constant: '<S32>/B44'
    //   Integrator: '<S32>/Integrator'
    //   Product: '<S32>/Product'
    //   Product: '<S32>/Product1'

    for (i = 0; i < 5; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        rtb_Row1 += AHV_Model_ConstP.pooled52[5 * i_0 + i] *
                    AHV_Model_X.Integrator_CSTATE_l[i_0];
      }

      AHV_Model_B.Sum[i] =
          AHV_Model_ConstP.pooled53[i] * AHV_Model_B.nu_r[3] + rtb_Row1;
    }

    // End of Sum: '<S32>/Sum'

    // Fcn: '<S16>/Fcn' incorporates:
    //   Integrator: '<S13>/Integrator'

    rtb_Row1 =
        AHV_Model_X.Integrator_CSTATE[0] * AHV_Model_X.Integrator_CSTATE[0] +
        AHV_Model_X.Integrator_CSTATE[1] * AHV_Model_X.Integrator_CSTATE[1];
    if (rtb_Row1 < 0.0) {
      AHV_Model_B.Fcn = -std::sqrt(-rtb_Row1);
    } else {
      AHV_Model_B.Fcn = std::sqrt(rtb_Row1);
    }

    // End of Fcn: '<S16>/Fcn'

    // Sum: '<S63>/Sum2' incorporates:
    //   Integrator: '<S13>/Integrator1'
    //   Integrator: '<S63>/Integrator1'
    //   Integrator: '<S63>/Integrator3'
    //   Sum: '<S63>/Sum4'

    psi_mean_m = AHV_Model_X.Integrator1_CSTATE[0] -
                 (AHV_Model_X.Integrator1_CSTATE_c[0] +
                  AHV_Model_X.Integrator3_CSTATE[0]);
    rtb_Row1_e = AHV_Model_X.Integrator1_CSTATE[1] -
                 (AHV_Model_X.Integrator1_CSTATE_c[1] +
                  AHV_Model_X.Integrator3_CSTATE[1]);
    rtb_Gain2[2] = AHV_Model_X.Integrator1_CSTATE[5] -
                   (AHV_Model_X.Integrator1_CSTATE_c[2] + rtb_psi_dot);

    // Saturate: '<S72>/x_Saturation'
    if (rtb_Gain2[2] > 1.0E+10) {
      rtb_x_dot = 1.0E+10;
    } else if (rtb_Gain2[2] < -1.0E+10) {
      rtb_x_dot = -1.0E+10;
    } else {
      rtb_x_dot = rtb_Gain2[2];
    }

    // End of Saturate: '<S72>/x_Saturation'

    // Signum: '<S72>/x_Sign'
    if (rtb_x_dot < 0.0) {
      rtb_T33 = -1.0;
    } else if (rtb_x_dot > 0.0) {
      rtb_T33 = 1.0;
    } else if (rtb_x_dot == 0.0) {
      rtb_T33 = 0.0;
    } else {
      rtb_T33 = (rtNaN);
    }

    // End of Signum: '<S72>/x_Sign'

    // Gain: '<S72>/pi'
    rtb_T33 *= 3.1415926535897931;

    // Sum: '<S72>/Sum1'
    rtb_x_dot += rtb_T33;

    // Math: '<S72>/Math Function' incorporates:
    //   Constant: '<S72>/Constant'

    rtb_x_dot = rt_remd_snf(rtb_x_dot, 6.2831853071795862);

    // Sum: '<S72>/Sum'
    rtb_Row1 = rtb_x_dot - rtb_T33;

    // Gain: '<S63>/K4' incorporates:
    //   Sum: '<S72>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = AHV_Model_ConstP.pooled109[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled109[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled109[i] * psi_mean_m);
    }

    // End of Gain: '<S63>/K4'

    // Fcn: '<S69>/Row1'
    rtb_x_dot = rtb_K2_u_tmp * rtb_K2_u[0] + rtb_K3_u_tmp * rtb_K2_u[1];

    // Fcn: '<S69>/Row2' incorporates:
    //   Integrator: '<S13>/Integrator1'

    rtb_T33 = -std::sin(AHV_Model_X.Integrator1_CSTATE[5]) * rtb_K2_u[0] +
              rtb_K2_u_tmp * rtb_K2_u[1];

    // Fcn: '<S69>/Row3'
    rtb_psi_dot = rtb_K2_u[2];

    // Gain: '<S63>/Gain6' incorporates:
    //   Integrator: '<S63>/Integrator4'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] =
          AHV_Model_ConstP.pooled110[i + 6] *
              AHV_Model_X.Integrator4_CSTATE[2] +
          (AHV_Model_ConstP.pooled110[i + 3] *
               AHV_Model_X.Integrator4_CSTATE[1] +
           AHV_Model_ConstP.pooled110[i] * AHV_Model_X.Integrator4_CSTATE[0]);
    }

    // End of Gain: '<S63>/Gain6'

    // Sum: '<S63>/Sum8' incorporates:
    //   Fcn: '<S70>/Row1'
    //   Fcn: '<S70>/Row2'
    //   Integrator: '<S13>/Integrator1'
    //   Integrator: '<S63>/Integrator6'
    //   Sum: '<S63>/Sum1'

    rtb_Gain2[0] = (((rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE[0] +
                      rtb_K3_u_tmp * AHV_Model_X.Integrator6_CSTATE[1]) +
                     rtb_Switch_a) +
                    rtb_x_dot) -
                   rtb_K2_u[0];
    rtb_Gain2[1] = (((-std::sin(AHV_Model_X.Integrator1_CSTATE[5]) *
                          AHV_Model_X.Integrator6_CSTATE[0] +
                      rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE[1]) +
                     rtb_Switch1) +
                    rtb_T33) -
                   rtb_K2_u[1];
    rtb_Gain2[2] =
        ((AHV_Model_X.Integrator6_CSTATE[2] + rtb_Row3) + rtb_psi_dot) -
        rtb_K2_u[2];
    for (i = 0; i < 3; i++) {
      // Gain: '<S63>/Gain3'
      AHV_Model_B.M_u[i] = 0.0;
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled111[i] * rtb_Gain2[0];
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled111[i + 3] * rtb_Gain2[1];
      AHV_Model_B.M_u[i] += AHV_Model_ConstP.pooled111[i + 6] * rtb_Gain2[2];

      // Sum: '<S63>/Sum5' incorporates:
      //   Gain: '<S63>/K11'
      //   Integrator: '<S63>/Integrator1'

      AHV_Model_B.psi_WF[i] = ((AHV_Model_ConstP.pooled112[i + 3] * rtb_Row1_e +
                                AHV_Model_ConstP.pooled112[i] * psi_mean_m) +
                               AHV_Model_ConstP.pooled112[i + 6] * rtb_Row1) +
                              AHV_Model_X.Integrator1_CSTATE_c[i];

      // Sum: '<S63>/Sum6' incorporates:
      //   Gain: '<S63>/Gain1'
      //   Gain: '<S63>/Gain2'
      //   Gain: '<S63>/K12'
      //   Integrator: '<S63>/Integrator1'
      //   Integrator: '<S63>/Integrator2'
      //   Sum: '<S72>/Sum'

      AHV_Model_B.Sum6[i] = ((AHV_Model_ConstP.pooled113[i + 6] * rtb_Row1 +
                              (AHV_Model_ConstP.pooled113[i + 3] * rtb_Row1_e +
                               AHV_Model_ConstP.pooled113[i] * psi_mean_m)) -
                             (AHV_Model_ConstP.pooled108[i + 6] *
                                  AHV_Model_X.Integrator2_CSTATE[2] +
                              (AHV_Model_ConstP.pooled108[i + 3] *
                                   AHV_Model_X.Integrator2_CSTATE[1] +
                               AHV_Model_ConstP.pooled108[i] *
                                   AHV_Model_X.Integrator2_CSTATE[0]))) -
                            ((AHV_Model_ConstP.pooled107[i + 3] *
                                  AHV_Model_X.Integrator1_CSTATE_c[1] +
                              AHV_Model_ConstP.pooled107[i] *
                                  AHV_Model_X.Integrator1_CSTATE_c[0]) +
                             AHV_Model_ConstP.pooled107[i + 6] *
                                 AHV_Model_X.Integrator1_CSTATE_c[2]);

      // Sum: '<S63>/Sum7' incorporates:
      //   Gain: '<S63>/K3'
      //   Gain: '<S63>/inv(T_b)'
      //   Integrator: '<S63>/Integrator6'
      //   Sum: '<S72>/Sum'

      AHV_Model_B.Sum7[i] =
          (AHV_Model_ConstP.pooled115[i + 6] * rtb_Row1 +
           (AHV_Model_ConstP.pooled115[i + 3] * rtb_Row1_e +
            AHV_Model_ConstP.pooled115[i] * psi_mean_m)) -
          (AHV_Model_ConstP.pooled116[i + 6] *
               AHV_Model_X.Integrator6_CSTATE[2] +
           (AHV_Model_ConstP.pooled116[i + 3] *
                AHV_Model_X.Integrator6_CSTATE[1] +
            AHV_Model_ConstP.pooled116[i] * AHV_Model_X.Integrator6_CSTATE[0]));

      // Gain: '<S63>/K2' incorporates:
      //   Sum: '<S63>/Sum2'

      rtb_K3_u[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled114[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled114[i] * psi_mean_m);
    }

    // Sum: '<S63>/Sum3' incorporates:
    //   Fcn: '<S68>/Fcn'
    //   Fcn: '<S68>/Fcn1'
    //   Fcn: '<S68>/Fcn2'
    //   Integrator: '<S63>/Integrator4'

    AHV_Model_B.sun_k2[0] = (rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE[0] -
                             rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE[1]) +
                            rtb_K3_u[0];
    AHV_Model_B.sun_k2[1] = (rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE[0] +
                             rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE[1]) +
                            rtb_K3_u[1];
    AHV_Model_B.sun_k2[2] = rtb_K3_u[2] + AHV_Model_X.Integrator4_CSTATE[2];

    // SignalConversion generated from: '<S65>/Integrator'
    AHV_Model_B.TmpSignalConversionAtIntegrat_p[0] = AHV_Model_B.Merge;
    AHV_Model_B.TmpSignalConversionAtIntegrat_p[1] = AHV_Model_B.Merge_m;
    for (i = 0; i < 6; i++) {
      // Integrator: '<S97>/Integrator1' incorporates:
      //   Inport: '<Root>/Vessel_init2'

      if (AHV_Model_DW.Integrator1_IWORK_a != 0) {
        AHV_Model_X.Integrator1_CSTATE_l[i] = Vessel_init2[i];
      }

      // Outport: '<Root>/eta_AHV2' incorporates:
      //   Integrator: '<S97>/Integrator1'

      eta_AHV2[i] = AHV_Model_X.Integrator1_CSTATE_l[i];

      // Outport: '<Root>/nu2' incorporates:
      //   Integrator: '<S97>/Integrator'

      nu2[i] = AHV_Model_X.Integrator_CSTATE_n[i];
    }

    // Outport: '<Root>/AHV_speed2' incorporates:
    //   Gain: '<S100>/Gain'
    //   Gain: '<S100>/Gain1'
    //   Rounding: '<S100>/Rounding Function'
    //   TransferFcn: '<S100>/Transfer Fcn'

    AHV_speed2 =
        std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_e * 10.0) * 0.1;

    // Gain: '<S96>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_x_dot = 0.017453292519943295 * Current_direction;

    // Trigonometry: '<S145>/sin(theta)' incorporates:
    //   Fcn: '<S106>/T23'
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/sin(theta)'

    rtb_T33 = std::sin(AHV_Model_X.Integrator1_CSTATE_l[3]);

    // Trigonometry: '<S145>/cos(theta)' incorporates:
    //   Fcn: '<S106>/T31 '
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/cos(theta)'

    psi_mean_m = std::cos(AHV_Model_X.Integrator1_CSTATE_l[3]);

    // Trigonometry: '<S145>/sin(theta)' incorporates:
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/sin(theta)'

    rtb_psi_dot = std::sin(AHV_Model_X.Integrator1_CSTATE_l[4]);

    // Trigonometry: '<S145>/cos(theta)' incorporates:
    //   Fcn: '<S106>/T23'
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/cos(theta)'

    rtb_Row1_e = std::cos(AHV_Model_X.Integrator1_CSTATE_l[4]);

    // Trigonometry: '<S145>/sin(theta)' incorporates:
    //   Fcn: '<S152>/Fcn'
    //   Fcn: '<S152>/Fcn1'
    //   Fcn: '<S153>/Row1'
    //   Fcn: '<S154>/Row1'
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/sin(theta)'

    rtb_K3_u_tmp = std::sin(AHV_Model_X.Integrator1_CSTATE_l[5]);

    // Trigonometry: '<S145>/cos(theta)' incorporates:
    //   Fcn: '<S152>/Fcn'
    //   Fcn: '<S152>/Fcn1'
    //   Fcn: '<S153>/Row1'
    //   Fcn: '<S153>/Row2'
    //   Fcn: '<S154>/Row1'
    //   Fcn: '<S154>/Row2'
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S105>/cos(theta)'

    rtb_K2_u_tmp = std::cos(AHV_Model_X.Integrator1_CSTATE_l[5]);

    // Product: '<S93>/Product1' incorporates:
    //   Fcn: '<S145>/R11'
    //   Fcn: '<S145>/R12'
    //   Fcn: '<S145>/R13'
    //   Fcn: '<S145>/R21 '
    //   Fcn: '<S145>/R22'
    //   Fcn: '<S145>/R23'
    //   Fcn: '<S145>/R31 '
    //   Fcn: '<S145>/R32'
    //   Fcn: '<S145>/R33'
    //   Inport: '<Root>/ahv_fairlead2'
    //   Trigonometry: '<S145>/cos(theta)'
    //   Trigonometry: '<S145>/sin(theta)'

    rtb_K3_u_0[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[3] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[6] = -rtb_psi_dot;
    rtb_K3_u_0[1] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 + -rtb_K3_u_tmp * psi_mean_m;
    rtb_K3_u_0[4] =
        rtb_K3_u_tmp * rtb_psi_dot * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_K3_u_0[7] = rtb_Row1_e * rtb_T33;
    rtb_K3_u_0[2] =
        rtb_K2_u_tmp * rtb_psi_dot * psi_mean_m + rtb_K3_u_tmp * rtb_T33;
    rtb_K3_u_0[5] =
        rtb_K3_u_tmp * rtb_psi_dot * psi_mean_m + -rtb_K2_u_tmp * rtb_T33;
    rtb_K3_u_0[8] = rtb_Row1_e * psi_mean_m;
    for (i = 0; i < 3; i++) {
      rtb_Switch_a = rtb_K3_u_0[i + 6] * ahv_fairlead2[2] +
                     (rtb_K3_u_0[i + 3] * ahv_fairlead2[1] +
                      rtb_K3_u_0[i] * ahv_fairlead2[0]);

      // Gain: '<S92>/Gain1' incorporates:
      //   Inport: '<Root>/ahv_fairlead2'

      rtb_Integrator[i] = AHV_Model_ConstP.pooled60[i] * rtb_Switch_a;
      rtb_Gain2[i] = rtb_Switch_a;
    }

    // End of Product: '<S93>/Product1'

    // Gain: '<S92>/Gain1' incorporates:
    //   Inport: '<Root>/tau_cable2'
    //   Product: '<S144>/ae'
    //   Product: '<S144>/af'
    //   Product: '<S144>/bd'
    //   Product: '<S144>/bf'
    //   Product: '<S144>/cd'
    //   Product: '<S144>/ce'
    //   Sum: '<S144>/Sum'
    //   Sum: '<S144>/Sum1'
    //   Sum: '<S144>/Sum2'

    rtb_Integrator[3] =
        tau_cable2[1] * rtb_Gain2[2] - tau_cable2[2] * rtb_Gain2[1];
    rtb_Integrator[4] =
        (tau_cable2[2] * rtb_Gain2[0] - tau_cable2[0] * rtb_Gain2[2]) * 0.5;
    rtb_Integrator[5] =
        tau_cable2[0] * rtb_Gain2[1] - tau_cable2[1] * rtb_Gain2[0];

    // Fcn: '<S106>/T21 ' incorporates:
    //   Fcn: '<S106>/T31 '
    //   Integrator: '<S97>/Integrator1'

    rtb_Switch_a = std::tan(AHV_Model_X.Integrator1_CSTATE_l[4]);

    // Reshape: '<S106>/Reshape 9x1->3x3' incorporates:
    //   Constant: '<S106>/Constant'
    //   Constant: '<S106>/Constant '
    //   Fcn: '<S106>/T21 '
    //   Fcn: '<S106>/T23'
    //   Fcn: '<S106>/T31 '
    //   Fcn: '<S106>/T32'
    //   Fcn: '<S106>/T33'

    rtb_TmpSignalConversionAtProduc[0] = 1.0;
    rtb_TmpSignalConversionAtProduc[1] = 0.0;
    rtb_TmpSignalConversionAtProduc[2] = 0.0;
    rtb_TmpSignalConversionAtProduc[3] = rtb_T33 * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[4] = psi_mean_m;
    rtb_TmpSignalConversionAtProduc[5] = rtb_T33 / rtb_Row1_e;
    rtb_TmpSignalConversionAtProduc[6] = psi_mean_m * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[7] = -rtb_T33;
    rtb_TmpSignalConversionAtProduc[8] = psi_mean_m / rtb_Row1_e;

    // SignalConversion generated from: '<S102>/Product' incorporates:
    //   Fcn: '<S105>/R11'
    //   Fcn: '<S105>/R12'
    //   Fcn: '<S105>/R21 '
    //   Fcn: '<S105>/R31 '

    rtb_TmpSignalConversionAtProd_k[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[1] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[2] = -rtb_psi_dot;
    rtb_TmpSignalConversionAtProd_k[3] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 - rtb_K3_u_tmp * psi_mean_m;

    // Fcn: '<S105>/R22' incorporates:
    //   Fcn: '<S105>/R23'

    rtb_Switch_a = rtb_K3_u_tmp * rtb_psi_dot;

    // SignalConversion generated from: '<S102>/Product' incorporates:
    //   Fcn: '<S105>/R13'
    //   Fcn: '<S105>/R22'
    //   Fcn: '<S105>/R23'
    //   Fcn: '<S105>/R32'
    //   Fcn: '<S105>/R33'

    rtb_TmpSignalConversionAtProd_k[4] =
        rtb_Switch_a * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_TmpSignalConversionAtProd_k[5] = rtb_Row1_e * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[6] =
        psi_mean_m * rtb_psi_dot * rtb_K2_u_tmp + rtb_K3_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[7] =
        rtb_Switch_a * psi_mean_m - rtb_K2_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[8] = rtb_Row1_e * psi_mean_m;

    // Sum: '<S97>/Sum6' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Integrator: '<S97>/Integrator'
    //   Product: '<S96>/Product'
    //   Product: '<S96>/Product1'
    //   Trigonometry: '<S96>/cos'
    //   Trigonometry: '<S96>/sin'

    AHV_Model_B.nu_r_b[0] = AHV_Model_X.Integrator_CSTATE_n[0] -
                            std::cos(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_b[1] = AHV_Model_X.Integrator_CSTATE_n[1] -
                            std::sin(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_b[2] = AHV_Model_X.Integrator_CSTATE_n[2];
    AHV_Model_B.nu_r_b[3] = AHV_Model_X.Integrator_CSTATE_n[3];
    AHV_Model_B.nu_r_b[4] = AHV_Model_X.Integrator_CSTATE_n[4];
    AHV_Model_B.nu_r_b[5] = AHV_Model_X.Integrator_CSTATE_n[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S103>/Cross-flow drag trapezoidal
      // integration' Constant: '<S103>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_d,
                                      AHV_Model_B.nu_r_b[1],
                                      AHV_Model_B.nu_r_b[5], &AHV_Model_B.Sum_k,
                                      &AHV_Model_B.Sum2_p, 82.800003);

      // End of Outputs for SubSystem: '<S103>/Cross-flow drag trapezoidal
      // integration'

      // Product: '<S103>/dx1' incorporates:
      //   Constant: '<S103>/2D drag coefficient '
      //   Constant: '<S103>/Transversal area//Lpp'
      //   Constant: '<S103>/rho'
      //   Gain: '<S103>/Gain1'

      AHV_Model_B.dx1_g =
          -0.5 * AHV_Model_B.Sum_k * 1025.0 * 5.4 * 0.69954741847648816;

      // Product: '<S103>/dx2' incorporates:
      //   Constant: '<S103>/2D drag coefficient '
      //   Constant: '<S103>/Transversal area//Lpp'
      //   Constant: '<S103>/rho'
      //   Gain: '<S103>/Gain'

      AHV_Model_B.dx2_o =
          -0.5 * AHV_Model_B.Sum2_p * 1025.0 * 5.4 * 0.69954741847648816;
    }

    // Product: '<S97>/G*eta' incorporates:
    //   Constant: '<S97>/Spring stiffness'
    //   Integrator: '<S97>/Integrator1'

    for (i = 0; i < 6; i++) {
      rtb_tau_WF[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_tau_WF[i] += AHV_Model_ConstP.pooled48[6 * i_0 + i] *
                         AHV_Model_X.Integrator1_CSTATE_l[i_0];
      }
    }

    // End of Product: '<S97>/G*eta'

    // Integrator: '<S147>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init2'
    //   SignalConversion generated from: '<S147>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_p != 0) {
      AHV_Model_X.Integrator3_CSTATE_i[0] = Vessel_init2[0];
      AHV_Model_X.Integrator3_CSTATE_i[1] = Vessel_init2[1];
      AHV_Model_X.Integrator3_CSTATE_i[2] = Vessel_init2[5];
    }

    // Saturate: '<S155>/x_Saturation' incorporates:
    //   Integrator: '<S147>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_i[2] > 1.0E+10) {
      rtb_T33 = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_i[2] < -1.0E+10) {
      rtb_T33 = -1.0E+10;
    } else {
      rtb_T33 = AHV_Model_X.Integrator3_CSTATE_i[2];
    }

    // End of Saturate: '<S155>/x_Saturation'

    // Signum: '<S155>/x_Sign'
    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S155>/x_Sign'

    // Gain: '<S155>/pi'
    rtb_psi_dot *= 3.1415926535897931;

    // Sum: '<S155>/Sum1'
    rtb_T33 += rtb_psi_dot;

    // Math: '<S155>/Math Function' incorporates:
    //   Constant: '<S155>/Constant'

    rtb_T33 = rt_remd_snf(rtb_T33, 6.2831853071795862);

    // Sum: '<S155>/Sum'
    rtb_T33 -= rtb_psi_dot;

    // Signum: '<S164>/x_Sign' incorporates:
    //   Fcn: '<S149>/yaw_angle'

    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S164>/x_Sign'

    // Gain: '<S164>/pi'
    rtb_Switch_a = 3.1415926535897931 * rtb_psi_dot;

    // Sum: '<S164>/Sum' incorporates:
    //   Constant: '<S164>/Constant'
    //   Fcn: '<S149>/yaw_angle'
    //   Math: '<S164>/Math Function'
    //   Sum: '<S164>/Sum1'

    rtb_psi_dot =
        rt_remd_snf(rtb_T33 + rtb_Switch_a, 6.2831853071795862) - rtb_Switch_a;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S148>/Chart' incorporates:
      //   Inport: '<Root>/Hold_Position2'
      //   Inport: '<Root>/Vessel_X_Ref2'
      //   Inport: '<Root>/Vessel_Y_Ref2'
      //   Integrator: '<S97>/Integrator1'

      AHV_Model_Chart(Hold_Position2, Vessel_X_Ref2, Vessel_Y_Ref2,
                      AHV_Model_X.Integrator1_CSTATE_l[0],
                      AHV_Model_X.Integrator1_CSTATE_l[1],
                      &AHV_Model_B.x_ref_rel_e, &AHV_Model_B.y_ref_rel_n,
                      &AHV_Model_DW.sf_Chart_o);
    }

    // Switch: '<S148>/Switch' incorporates:
    //   Gain: '<S151>/Gain'
    //   Inport: '<Root>/heading_angle_ref2'
    //   Inport: '<Root>/heading_mode2'
    //   Integrator: '<S97>/Integrator1'
    //   Trigonometry: '<S146>/Trigonometric Function'

    if (heading_mode2 > 0.5) {
      AHV_Model_B.RateLimiter_j =
          57.295779513082323 *
          rt_atan2d_snf(AHV_Model_X.Integrator1_CSTATE_l[1],
                        AHV_Model_X.Integrator1_CSTATE_l[0]);
    } else {
      AHV_Model_B.RateLimiter_j = heading_angle_ref2;
    }

    // End of Switch: '<S148>/Switch'

    // RateLimiter: '<S148>/Rate Limiter'
    if (!(AHV_Model_DW.LastMajorTime_n == (rtInf))) {
      rtb_Row3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_n;
      rtb_Row1 = rtb_Row3 * 10.0;
      rtb_Switch1 = AHV_Model_B.RateLimiter_j - AHV_Model_DW.PrevY_p;
      if (rtb_Switch1 > rtb_Row1) {
        AHV_Model_B.RateLimiter_j = AHV_Model_DW.PrevY_p + rtb_Row1;
      } else {
        rtb_Row3 *= -10.0;
        if (rtb_Switch1 < rtb_Row3) {
          AHV_Model_B.RateLimiter_j = AHV_Model_DW.PrevY_p + rtb_Row3;
        }
      }
    }

    // End of RateLimiter: '<S148>/Rate Limiter'

    // Gain: '<S158>/Gain1'
    rtb_Switch_a = 0.017453292519943295 * AHV_Model_B.RateLimiter_j;

    // Saturate: '<S159>/x_Saturation'
    if (rtb_Switch_a > 1.0E+10) {
      rtb_Switch_a = 1.0E+10;
    } else {
      if (rtb_Switch_a < -1.0E+10) {
        rtb_Switch_a = -1.0E+10;
      }
    }

    // End of Saturate: '<S159>/x_Saturation'

    // Signum: '<S159>/x_Sign'
    if (rtb_Switch_a < 0.0) {
      psi_mean_m = -1.0;
    } else if (rtb_Switch_a > 0.0) {
      psi_mean_m = 1.0;
    } else if (rtb_Switch_a == 0.0) {
      psi_mean_m = 0.0;
    } else {
      psi_mean_m = (rtNaN);
    }

    // End of Signum: '<S159>/x_Sign'

    // Gain: '<S159>/pi'
    rtb_Switch1 = 3.1415926535897931 * psi_mean_m;

    // Sum: '<S159>/Sum1'
    rtb_Switch_a += rtb_Switch1;

    // Math: '<S159>/Math Function' incorporates:
    //   Constant: '<S159>/Constant'

    rtb_Switch_a = rt_remd_snf(rtb_Switch_a, 6.2831853071795862);

    // Sum: '<S159>/Sum'
    rtb_Switch_a -= rtb_Switch1;

    // Sum: '<S149>/Sum2' incorporates:
    //   Integrator: '<S147>/Integrator3'

    rtb_Gain2[0] =
        AHV_Model_B.x_ref_rel_e - AHV_Model_X.Integrator3_CSTATE_i[0];
    rtb_Gain2[1] =
        AHV_Model_B.y_ref_rel_n - AHV_Model_X.Integrator3_CSTATE_i[1];
    rtb_Gain2[2] = rtb_Switch_a - rtb_T33;

    // Signum: '<S163>/x_Sign' incorporates:
    //   Saturate: '<S163>/x_Saturation'

    if (rtb_Gain2[2] < 0.0) {
      rtb_Switch_a = -1.0;
    } else if (rtb_Gain2[2] > 0.0) {
      rtb_Switch_a = 1.0;
    } else if (rtb_Gain2[2] == 0.0) {
      rtb_Switch_a = 0.0;
    } else {
      rtb_Switch_a = (rtNaN);
    }

    // End of Signum: '<S163>/x_Sign'

    // Gain: '<S163>/pi'
    rtb_Switch1 = 3.1415926535897931 * rtb_Switch_a;

    // Fcn: '<S162>/Row3' incorporates:
    //   Constant: '<S163>/Constant'
    //   Math: '<S163>/Math Function'
    //   Saturate: '<S163>/x_Saturation'
    //   Sum: '<S163>/Sum'
    //   Sum: '<S163>/Sum1'

    AHV_Model_B.Row3_j =
        rt_remd_snf(rtb_Gain2[2] + rtb_Switch1, 6.2831853071795862) -
        rtb_Switch1;

    // Fcn: '<S162>/Row1' incorporates:
    //   Fcn: '<S162>/Row2'

    psi_mean_m = std::sin(rtb_psi_dot);
    rtb_psi_dot = std::cos(rtb_psi_dot);
    rtb_Row1_e = rtb_psi_dot * rtb_Gain2[0] + psi_mean_m * rtb_Gain2[1];

    // If: '<S167>/If' incorporates:
    //   Abs: '<S167>/Abs'
    //   Constant: '<Root>/Constant41'
    //   Inport: '<S170>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_h =
          static_cast<int8_T>(!(std::abs(rtb_Row1_e) <= Featured));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_h) {
      case 0:
        // Outputs for IfAction SubSystem: '<S167>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S170>/Action Port'

        AHV_Model_B.Merge_a = rtb_Row1_e;

        // End of Outputs for SubSystem: '<S167>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S167>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S169>/Action Port'

        AHV_Model_IfActionSubsystem_l(1.0, rtb_Row1_e, &AHV_Model_B.Merge_a,
                                      &AHV_Model_ConstB.IfActionSubsystem_d);

        // End of Outputs for SubSystem: '<S167>/If Action Subsystem'
        break;
    }

    // End of If: '<S167>/If'

    // Fcn: '<S162>/Row2'
    rtb_psi_dot = -psi_mean_m * rtb_Gain2[0] + rtb_psi_dot * rtb_Gain2[1];

    // If: '<S168>/If' incorporates:
    //   Abs: '<S168>/Abs'
    //   Constant: '<Root>/Constant42'
    //   Inport: '<S172>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_ig =
          static_cast<int8_T>(!(std::abs(rtb_psi_dot) <= SidePush));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_ig) {
      case 0:
        // Outputs for IfAction SubSystem: '<S168>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S172>/Action Port'

        AHV_Model_B.Merge_g = rtb_psi_dot;

        // End of Outputs for SubSystem: '<S168>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S168>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S171>/Action Port'

        AHV_Model_IfActionSubsystem_l(10.0, rtb_psi_dot, &AHV_Model_B.Merge_g,
                                      &AHV_Model_ConstB.IfActionSubsystem_f);

        // End of Outputs for SubSystem: '<S168>/If Action Subsystem'
        break;
    }

    // End of If: '<S168>/If'

    // If: '<S160>/If' incorporates:
    //   Inport: '<Root>/Hold_Position2'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_ab = static_cast<int8_T>(!Hold_Position2);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_ab) {
      case 0:
        // Outputs for IfAction SubSystem: '<S160>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S165>/Action Port'

        AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_a, AHV_Model_B.Merge_g,
                                    AHV_Model_B.Row3_j, rtb_Gain2,
                                    AHV_Model_ConstP.pooled2);

        // End of Outputs for SubSystem: '<S160>/If Action Subsystem'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S160>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S166>/Action Port'

        AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_a, AHV_Model_B.Merge_g,
                                     AHV_Model_B.Row3_j, rtb_Gain2,
                                     AHV_Model_ConstP.pooled1);

        // End of Outputs for SubSystem: '<S160>/If Action Subsystem1'
        break;
    }

    // End of If: '<S160>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S149>/Ki' incorporates:
      //   DiscreteIntegrator: '<S149>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki[i] = 0.0;
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[0];
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[1];
        AHV_Model_B.Ki[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE[2];
      }

      // End of Gain: '<S149>/Ki'
    }

    // Sum: '<S149>/Sum1' incorporates:
    //   Gain: '<S149>/Kd'
    //   Integrator: '<S147>/Integrator4'
    //   Sum: '<S149>/Sum3'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = (rtb_Gain2[i] + AHV_Model_B.Ki[i]) -
                    ((AHV_Model_ConstP.pooled64[i + 3] *
                          AHV_Model_X.Integrator4_CSTATE_p[1] +
                      AHV_Model_ConstP.pooled64[i] *
                          AHV_Model_X.Integrator4_CSTATE_p[0]) +
                     AHV_Model_ConstP.pooled64[i + 6] *
                         AHV_Model_X.Integrator4_CSTATE_p[2]);
    }

    // End of Sum: '<S149>/Sum1'

    // Saturate: '<S150>/Surge Force Saturation'
    if (rtb_K2_u[0] > 2.0E+6) {
      rtb_psi_dot = 2.0E+6;
    } else if (rtb_K2_u[0] < -2.0E+6) {
      rtb_psi_dot = -2.0E+6;
    } else {
      rtb_psi_dot = rtb_K2_u[0] * MainThrust;
    }

    // End of Saturate: '<S150>/Surge Force Saturation'

    // Saturate: '<S150>/Sway Force Saturation'
    if (rtb_K2_u[1] > 1.5E+6) {
      rtb_Switch_a = 1.5E+6;
    } else if (rtb_K2_u[1] < -1.5E+6) {
      rtb_Switch_a = -1.5E+6;
    } else {
      rtb_Switch_a = rtb_K2_u[1] * FuThrust;
    }

    // End of Saturate: '<S150>/Sway Force Saturation'

    // Saturate: '<S150>/Yaw Moment Saturation'
    if (rtb_K2_u[2] > 2.0E+7) {
      rtb_Switch1 = 2.0E+7;
    } else if (rtb_K2_u[2] < -2.0E+7) {
      rtb_Switch1 = -2.0E+7;
    } else {
      rtb_Switch1 = rtb_K2_u[2] * AzimuthThrust;
    }

    // End of Saturate: '<S150>/Yaw Moment Saturation'

    // Sum: '<S97>/Sum1' incorporates:
    //   Constant: '<S94>/Constant'

    rtb_Row1 = rtb_tau_WF[0];
    tmp = rtb_tau_WF[1];
    tmp_0 = rtb_tau_WF[2];
    tmp_1 = rtb_tau_WF[3];
    tmp_2 = rtb_tau_WF[4];
    tmp_3 = rtb_tau_WF[5];
    rtb_tau_WF[0] = rtb_psi_dot - rtb_Row1;
    rtb_tau_WF[1] = rtb_Switch_a - tmp;
    rtb_tau_WF[2] = 0.0 - tmp_0;
    rtb_tau_WF[3] = 0.0 - tmp_1;
    rtb_tau_WF[4] = 0.0 - tmp_2;
    rtb_tau_WF[5] = rtb_Switch1 - tmp_3;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S92>/Subsystem3'
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

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread,
                depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff,
                rand_freq, rand_dir, AHV_Model_B.Zeta_a_d, AHV_Model_B.Omega_e,
                AHV_Model_B.Phase_j, AHV_Model_B.Wavenum_m, AHV_Model_B.Psi_a,
                &psi_mean_m, &AHV_Model_B.Subsystem3_b,
                &AHV_Model_DW.Subsystem3_b, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S92>/Subsystem3'
    }

    // Outputs for Atomic SubSystem: '<S92>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_X.Integrator1_CSTATE_l, AHV_Model_B.Psi_a,
                   AHV_Model_B.Wavenum_m, AHV_Model_B.Omega_e,
                   AHV_Model_B.Phase_j, AHV_Model_B.Zeta_a_d, rtb_Sum1_di,
                   rtb_tau_WD, &AHV_Model_B.WaveloadsU0_b,
                   &AHV_Model_ConstB.WaveloadsU0_b,
                   &AHV_Model_DW.WaveloadsU0_b);

    // End of Outputs for SubSystem: '<S92>/Wave loads (U=0)'

    // Sum: '<S97>/Sum7' incorporates:
    //   Sum: '<S92>/Sum'
    //   Sum: '<S97>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Integrator[i] += (rtb_Sum1_di[i] + rtb_tau_WD[i]) + rtb_tau_WF[i];
    }

    // End of Sum: '<S97>/Sum7'

    // Sum: '<S97>/Sum5' incorporates:
    //   Abs: '<S103>/Abs1'
    //   Gain: '<S103>/Gain3'
    //   Product: '<S103>/Product7'
    //   Product: '<S103>/Product8'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] = -(std::abs(AHV_Model_B.nu_r_b[0]) *
                          AHV_Model_B.nu_r_b[0] * AHV_Model_ConstB.Product5_o) +
                        rtb_Row1;
    rtb_Integrator[1] = tmp + AHV_Model_B.dx1_g;
    rtb_Integrator[2] = tmp_0;
    rtb_Integrator[3] = tmp_1;
    rtb_Integrator[4] = tmp_2;
    rtb_Integrator[5] = tmp_3 + AHV_Model_B.dx2_o;

    // Sum: '<S97>/Sum4' incorporates:
    //   Constant: '<S97>/damping'
    //   Product: '<S97>/D*eta '

    for (i = 0; i < 6; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Row1 +=
            AHV_Model_ConstP.pooled51[6 * i_0 + i] * AHV_Model_B.nu_r_b[i_0];
      }

      rtb_Integrator[i] -= rtb_Row1;
    }

    // End of Sum: '<S97>/Sum4'

    // Product: '<S116>/Product2' incorporates:
    //   Integrator: '<S116>/Integrator'

    rtb_Sum5_0 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Sum5_0 += AHV_Model_ConstB.MathFunction_i[i] *
                    AHV_Model_X.Integrator_CSTATE_m[i];
    }

    // Sum: '<S97>/Sum' incorporates:
    //   Constant: '<S116>/B44_inf'
    //   Constant: '<S116>/D44'
    //   Product: '<S116>/Product2'
    //   Product: '<S116>/Product3'
    //   Product: '<S116>/Product4'
    //   StateSpace: '<S104>/Dp(1,1)'
    //   StateSpace: '<S104>/Dp(2,2)'
    //   StateSpace: '<S104>/Dp(2,4)'
    //   StateSpace: '<S104>/Dp(2,6)'
    //   StateSpace: '<S104>/Dp(3,3)'
    //   StateSpace: '<S104>/Dp(3,5)'
    //   StateSpace: '<S104>/Dp(4,2)'
    //   StateSpace: '<S104>/Dp(4,6)'
    //   StateSpace: '<S104>/Dp(5,3)'
    //   StateSpace: '<S104>/Dp(5,5)'
    //   StateSpace: '<S104>/Dp(6,2)'
    //   StateSpace: '<S104>/Dp(6,4)'
    //   StateSpace: '<S104>/Dp(6,6)'
    //   Sum: '<S104>/Sum'
    //   Sum: '<S104>/Sum1'
    //   Sum: '<S104>/Sum3'
    //   Sum: '<S104>/Sum4'
    //   Sum: '<S104>/Sum5'
    //   Sum: '<S116>/Sum1'
    //   Sum: '<S116>/Sum2'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] =
        rtb_Row1 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE_k[0] +
                        -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE_k[1]) +
                       -26196.671588656758 * AHV_Model_X.Dp11_CSTATE_k[2]) +
                      -69550.629589629229 * AHV_Model_X.Dp11_CSTATE_k[3]) +
                     -6667.9011063139569 * AHV_Model_X.Dp11_CSTATE_k[4]) +
                    28360.700724836512 * AHV_Model_B.nu_r_b[0]);
    rtb_Integrator[1] =
        tmp - (((((((-1.9219252802461751E+6 * AHV_Model_X.Dp22_CSTATE_l[0] +
                     -43506.445868266834 * AHV_Model_X.Dp22_CSTATE_l[1]) +
                    -261966.71588657616 * AHV_Model_X.Dp22_CSTATE_l[2]) +
                   -695506.29589629406 * AHV_Model_X.Dp22_CSTATE_l[3]) +
                  -66679.0110631567 * AHV_Model_X.Dp22_CSTATE_l[4]) +
                 283607.00724836387 * AHV_Model_B.nu_r_b[1]) +
                (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_i[0] +
                     -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_i[1]) +
                    782719.9224317756 * AHV_Model_X.Dp24_CSTATE_i[2]) +
                   906439.02406890341 * AHV_Model_X.Dp24_CSTATE_i[3]) +
                  48702.05056677026 * AHV_Model_X.Dp24_CSTATE_i[4]) +
                 342290.49224698928 * AHV_Model_B.nu_r_b[3])) +
               (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE_p[0] +
                    9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE_p[1]) +
                   1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE_p[2]) +
                  3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE_p[3]) +
                 -523575.51382008579 * AHV_Model_X.Dp26_CSTATE_p[4]) +
                1.4324880870688609E+6 * AHV_Model_B.nu_r_b[5]));
    rtb_Integrator[2] =
        tmp_0 - ((((((-2.8230572248621574E+6 * AHV_Model_X.Dp33_CSTATE_m[0] +
                      2971.0437052436932 * AHV_Model_X.Dp33_CSTATE_m[1]) +
                     -112975.07052200216 * AHV_Model_X.Dp33_CSTATE_m[2]) +
                    -459867.38057091518 * AHV_Model_X.Dp33_CSTATE_m[3]) +
                   35755.997807249405 * AHV_Model_X.Dp33_CSTATE_m[4]) +
                  425986.43085235963 * AHV_Model_B.nu_r_b[2]) +
                 (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE_n[0] +
                      -189169.67892356269 * AHV_Model_X.Dp35_CSTATE_n[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE_n[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE_n[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE_n[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_b[4]));
    rtb_Integrator[3] =
        tmp_1 -
        (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE_c[0] +
               -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE_c[1]) +
              782719.9224317756 * AHV_Model_X.Dp42_CSTATE_c[2]) +
             906439.02406890341 * AHV_Model_X.Dp42_CSTATE_c[3]) +
            48702.05056677026 * AHV_Model_X.Dp42_CSTATE_c[4]) +
           342290.49224698928 * AHV_Model_B.nu_r_b[1]) +
          ((2.3074903324953854E+7 * AHV_Model_B.nu_r_b[3] + rtb_Sum5_0) +
           4.97124674352808E+7 * AHV_Model_B.nu_r_b[3])) +
         (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE_o[0] +
              2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE_o[1]) +
             -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE_o[2]) +
            1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE_o[3]) +
           -70143.531785937739 * AHV_Model_X.Dp46_CSTATE_o[4]) +
          -8.0522169282609783E+6 * AHV_Model_B.nu_r_b[5]));
    rtb_Integrator[4] =
        tmp_2 - ((((((-4.3928932046187744E+7 * AHV_Model_X.Dp53_CSTATE_p[0] +
                      -189169.67892356269 * AHV_Model_X.Dp53_CSTATE_p[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp53_CSTATE_p[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp53_CSTATE_p[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp53_CSTATE_p[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_b[2]) +
                 (((((-1.3510265271416564E+7 * AHV_Model_X.Dp55_CSTATE_e[0] +
                      2.1371969925573766E+9 * AHV_Model_X.Dp55_CSTATE_e[1]) +
                     9.8710379288756877E+7 * AHV_Model_X.Dp55_CSTATE_e[2]) +
                    6.1487736935078049E+8 * AHV_Model_X.Dp55_CSTATE_e[3]) +
                   5.7220438213245332E+8 * AHV_Model_X.Dp55_CSTATE_e[4]) +
                  3.0813057877131975E+8 * AHV_Model_B.nu_r_b[4]));
    rtb_Integrator[5] =
        tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE_k[0] +
                       9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE_k[1]) +
                      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE_k[2]) +
                     3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE_k[3]) +
                    -523575.51382008579 * AHV_Model_X.Dp62_CSTATE_k[4]) +
                   1.4324880870688609E+6 * AHV_Model_B.nu_r_b[1]) +
                  (((((-1.5746048421550707E+6 * AHV_Model_X.Dp64_CSTATE_l[0] +
                       2.9884305406230576E+7 * AHV_Model_X.Dp64_CSTATE_l[1]) +
                      -4.3953040127794016E+6 * AHV_Model_X.Dp64_CSTATE_l[2]) +
                     1.205872403100216E+7 * AHV_Model_X.Dp64_CSTATE_l[3]) +
                    -70143.531785937739 * AHV_Model_X.Dp64_CSTATE_l[4]) +
                   -8.0522169282609783E+6 * AHV_Model_B.nu_r_b[3])) +
                 (((((-8.07757267677547E+8 * AHV_Model_X.Dp66_CSTATE_g[0] +
                      -1.0798841515040772E+7 * AHV_Model_X.Dp66_CSTATE_g[1]) +
                     -1.27285250509213E+8 * AHV_Model_X.Dp66_CSTATE_g[2]) +
                    -3.0027175358715254E+8 * AHV_Model_X.Dp66_CSTATE_g[3]) +
                   1.7969070654437996E+7 * AHV_Model_X.Dp66_CSTATE_g[4]) +
                  1.1943182502908507E+8 * AHV_Model_B.nu_r_b[5]));

    // Product: '<S97>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_o, rtb_Integrator,
                                    AHV_Model_B.Minvtau_n);
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S97>/Integrator1' incorporates:
      //   Integrator: '<S97>/Integrator'
      //   Product: '<S102>/Product'
      //   Product: '<S102>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegra_o3[i] =
          rtb_TmpSignalConversionAtProd_k[i + 6] *
              AHV_Model_X.Integrator_CSTATE_n[2] +
          (rtb_TmpSignalConversionAtProd_k[i + 3] *
               AHV_Model_X.Integrator_CSTATE_n[1] +
           rtb_TmpSignalConversionAtProd_k[i] *
               AHV_Model_X.Integrator_CSTATE_n[0]);
      AHV_Model_B.TmpSignalConversionAtIntegra_o3[i + 3] =
          rtb_TmpSignalConversionAtProduc[i + 6] *
              AHV_Model_X.Integrator_CSTATE_n[5] +
          (rtb_TmpSignalConversionAtProduc[i + 3] *
               AHV_Model_X.Integrator_CSTATE_n[4] +
           rtb_TmpSignalConversionAtProduc[i] *
               AHV_Model_X.Integrator_CSTATE_n[3]);
    }

    // Sum: '<S116>/Sum' incorporates:
    //   Constant: '<S116>/A44'
    //   Constant: '<S116>/B44'
    //   Integrator: '<S116>/Integrator'
    //   Product: '<S116>/Product'
    //   Product: '<S116>/Product1'

    for (i = 0; i < 5; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        rtb_Row1 += AHV_Model_ConstP.pooled52[5 * i_0 + i] *
                    AHV_Model_X.Integrator_CSTATE_m[i_0];
      }

      AHV_Model_B.Sum_j[i] =
          AHV_Model_ConstP.pooled53[i] * AHV_Model_B.nu_r_b[3] + rtb_Row1;
    }

    // End of Sum: '<S116>/Sum'

    // Fcn: '<S100>/Fcn' incorporates:
    //   Integrator: '<S97>/Integrator'

    rtb_Row1 =
        AHV_Model_X.Integrator_CSTATE_n[0] *
            AHV_Model_X.Integrator_CSTATE_n[0] +
        AHV_Model_X.Integrator_CSTATE_n[1] * AHV_Model_X.Integrator_CSTATE_n[1];
    if (rtb_Row1 < 0.0) {
      AHV_Model_B.Fcn_j = -std::sqrt(-rtb_Row1);
    } else {
      AHV_Model_B.Fcn_j = std::sqrt(rtb_Row1);
    }

    // End of Fcn: '<S100>/Fcn'

    // Sum: '<S147>/Sum2' incorporates:
    //   Integrator: '<S147>/Integrator1'
    //   Integrator: '<S147>/Integrator3'
    //   Integrator: '<S97>/Integrator1'
    //   Sum: '<S147>/Sum4'

    psi_mean_m = AHV_Model_X.Integrator1_CSTATE_l[0] -
                 (AHV_Model_X.Integrator1_CSTATE_d[0] +
                  AHV_Model_X.Integrator3_CSTATE_i[0]);
    rtb_Row1_e = AHV_Model_X.Integrator1_CSTATE_l[1] -
                 (AHV_Model_X.Integrator1_CSTATE_d[1] +
                  AHV_Model_X.Integrator3_CSTATE_i[1]);
    rtb_Gain2[2] = AHV_Model_X.Integrator1_CSTATE_l[5] -
                   (AHV_Model_X.Integrator1_CSTATE_d[2] + rtb_T33);

    // Saturate: '<S156>/x_Saturation'
    if (rtb_Gain2[2] > 1.0E+10) {
      rtb_x_dot = 1.0E+10;
    } else if (rtb_Gain2[2] < -1.0E+10) {
      rtb_x_dot = -1.0E+10;
    } else {
      rtb_x_dot = rtb_Gain2[2];
    }

    // End of Saturate: '<S156>/x_Saturation'

    // Signum: '<S156>/x_Sign'
    if (rtb_x_dot < 0.0) {
      rtb_T33 = -1.0;
    } else if (rtb_x_dot > 0.0) {
      rtb_T33 = 1.0;
    } else if (rtb_x_dot == 0.0) {
      rtb_T33 = 0.0;
    } else {
      rtb_T33 = (rtNaN);
    }

    // End of Signum: '<S156>/x_Sign'

    // Gain: '<S156>/pi'
    rtb_T33 *= 3.1415926535897931;

    // Sum: '<S156>/Sum1'
    rtb_x_dot += rtb_T33;

    // Math: '<S156>/Math Function' incorporates:
    //   Constant: '<S156>/Constant'

    rtb_x_dot = rt_remd_snf(rtb_x_dot, 6.2831853071795862);

    // Sum: '<S156>/Sum'
    rtb_Row1 = rtb_x_dot - rtb_T33;

    // Gain: '<S147>/K4' incorporates:
    //   Sum: '<S156>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = AHV_Model_ConstP.pooled109[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled109[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled109[i] * psi_mean_m);
    }

    // End of Gain: '<S147>/K4'

    // Fcn: '<S153>/Row1'
    rtb_x_dot = rtb_K2_u_tmp * rtb_K2_u[0] + rtb_K3_u_tmp * rtb_K2_u[1];

    // Fcn: '<S153>/Row2' incorporates:
    //   Integrator: '<S97>/Integrator1'

    rtb_T33 = -std::sin(AHV_Model_X.Integrator1_CSTATE_l[5]) * rtb_K2_u[0] +
              rtb_K2_u_tmp * rtb_K2_u[1];

    // Fcn: '<S153>/Row3'
    rtb_Row3 = rtb_K2_u[2];

    // Gain: '<S147>/Gain6' incorporates:
    //   Integrator: '<S147>/Integrator4'

    for (i = 0; i < 3; i++) {
      rtb_Gain2[i] =
          AHV_Model_ConstP.pooled110[i + 6] *
              AHV_Model_X.Integrator4_CSTATE_p[2] +
          (AHV_Model_ConstP.pooled110[i + 3] *
               AHV_Model_X.Integrator4_CSTATE_p[1] +
           AHV_Model_ConstP.pooled110[i] * AHV_Model_X.Integrator4_CSTATE_p[0]);
    }

    // End of Gain: '<S147>/Gain6'

    // Sum: '<S147>/Sum8' incorporates:
    //   Fcn: '<S154>/Row1'
    //   Fcn: '<S154>/Row2'
    //   Integrator: '<S147>/Integrator6'
    //   Integrator: '<S97>/Integrator1'
    //   Sum: '<S147>/Sum1'

    rtb_K2_u[0] = (((rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_b[0] +
                     rtb_K3_u_tmp * AHV_Model_X.Integrator6_CSTATE_b[1]) +
                    rtb_psi_dot) +
                   rtb_x_dot) -
                  rtb_Gain2[0];
    rtb_K2_u[1] = (((-std::sin(AHV_Model_X.Integrator1_CSTATE_l[5]) *
                         AHV_Model_X.Integrator6_CSTATE_b[0] +
                     rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_b[1]) +
                    rtb_Switch_a) +
                   rtb_T33) -
                  rtb_Gain2[1];
    rtb_K2_u[2] =
        ((AHV_Model_X.Integrator6_CSTATE_b[2] + rtb_Switch1) + rtb_Row3) -
        rtb_Gain2[2];

    // Gain: '<S147>/Gain3'
    for (i = 0; i < 3; i++) {
      AHV_Model_B.M_u_f[i] = 0.0;
      AHV_Model_B.M_u_f[i] += AHV_Model_ConstP.pooled111[i] * rtb_K2_u[0];
      AHV_Model_B.M_u_f[i] += AHV_Model_ConstP.pooled111[i + 3] * rtb_K2_u[1];
      AHV_Model_B.M_u_f[i] += AHV_Model_ConstP.pooled111[i + 6] * rtb_K2_u[2];
    }

    // End of Gain: '<S147>/Gain3'
    for (i = 0; i < 3; i++) {
      // Sum: '<S147>/Sum5' incorporates:
      //   Gain: '<S147>/K11'
      //   Integrator: '<S147>/Integrator1'

      AHV_Model_B.psi_WF_g[i] =
          ((AHV_Model_ConstP.pooled112[i + 3] * rtb_Row1_e +
            AHV_Model_ConstP.pooled112[i] * psi_mean_m) +
           AHV_Model_ConstP.pooled112[i + 6] * rtb_Row1) +
          AHV_Model_X.Integrator1_CSTATE_d[i];

      // Sum: '<S147>/Sum6' incorporates:
      //   Gain: '<S147>/Gain1'
      //   Gain: '<S147>/Gain2'
      //   Gain: '<S147>/K12'
      //   Integrator: '<S147>/Integrator1'
      //   Integrator: '<S147>/Integrator2'
      //   Sum: '<S156>/Sum'

      AHV_Model_B.Sum6_i[i] =
          ((AHV_Model_ConstP.pooled113[i + 6] * rtb_Row1 +
            (AHV_Model_ConstP.pooled113[i + 3] * rtb_Row1_e +
             AHV_Model_ConstP.pooled113[i] * psi_mean_m)) -
           (AHV_Model_ConstP.pooled108[i + 6] *
                AHV_Model_X.Integrator2_CSTATE_i[2] +
            (AHV_Model_ConstP.pooled108[i + 3] *
                 AHV_Model_X.Integrator2_CSTATE_i[1] +
             AHV_Model_ConstP.pooled108[i] *
                 AHV_Model_X.Integrator2_CSTATE_i[0]))) -
          ((AHV_Model_ConstP.pooled107[i + 3] *
                AHV_Model_X.Integrator1_CSTATE_d[1] +
            AHV_Model_ConstP.pooled107[i] *
                AHV_Model_X.Integrator1_CSTATE_d[0]) +
           AHV_Model_ConstP.pooled107[i + 6] *
               AHV_Model_X.Integrator1_CSTATE_d[2]);

      // Sum: '<S147>/Sum7' incorporates:
      //   Gain: '<S147>/K3'
      //   Gain: '<S147>/inv(T_b)'
      //   Integrator: '<S147>/Integrator6'
      //   Sum: '<S156>/Sum'

      AHV_Model_B.Sum7_k[i] = (AHV_Model_ConstP.pooled115[i + 6] * rtb_Row1 +
                               (AHV_Model_ConstP.pooled115[i + 3] * rtb_Row1_e +
                                AHV_Model_ConstP.pooled115[i] * psi_mean_m)) -
                              (AHV_Model_ConstP.pooled116[i + 6] *
                                   AHV_Model_X.Integrator6_CSTATE_b[2] +
                               (AHV_Model_ConstP.pooled116[i + 3] *
                                    AHV_Model_X.Integrator6_CSTATE_b[1] +
                                AHV_Model_ConstP.pooled116[i] *
                                    AHV_Model_X.Integrator6_CSTATE_b[0]));

      // Gain: '<S147>/K2' incorporates:
      //   Sum: '<S147>/Sum2'

      rtb_K2_u[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled114[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled114[i] * psi_mean_m);
    }

    // Sum: '<S147>/Sum3' incorporates:
    //   Fcn: '<S152>/Fcn'
    //   Fcn: '<S152>/Fcn1'
    //   Fcn: '<S152>/Fcn2'
    //   Integrator: '<S147>/Integrator4'

    AHV_Model_B.sun_k2_d[0] =
        (rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_p[0] -
         rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_p[1]) +
        rtb_K2_u[0];
    AHV_Model_B.sun_k2_d[1] =
        (rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_p[0] +
         rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_p[1]) +
        rtb_K2_u[1];
    AHV_Model_B.sun_k2_d[2] = rtb_K2_u[2] + AHV_Model_X.Integrator4_CSTATE_p[2];
    for (i = 0; i < 6; i++) {
      // Integrator: '<S181>/Integrator1' incorporates:
      //   Inport: '<Root>/Vessel_init3'

      if (AHV_Model_DW.Integrator1_IWORK_p != 0) {
        AHV_Model_X.Integrator1_CSTATE_lj[i] = Vessel_init3[i];
      }

      // Outport: '<Root>/eta_AHV3' incorporates:
      //   Integrator: '<S181>/Integrator1'

      eta_AHV3[i] = AHV_Model_X.Integrator1_CSTATE_lj[i];

      // Outport: '<Root>/nu3' incorporates:
      //   Integrator: '<S181>/Integrator'

      nu3[i] = AHV_Model_X.Integrator_CSTATE_o[i];
    }

    // Outport: '<Root>/AHV_speed3' incorporates:
    //   Gain: '<S184>/Gain'
    //   Gain: '<S184>/Gain1'
    //   Rounding: '<S184>/Rounding Function'
    //   TransferFcn: '<S184>/Transfer Fcn'

    AHV_speed3 =
        std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_l * 10.0) * 0.1;

    // Gain: '<S180>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_x_dot = 0.017453292519943295 * Current_direction;

    // Trigonometry: '<S229>/sin(theta)' incorporates:
    //   Fcn: '<S190>/T23'
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/sin(theta)'

    rtb_T33 = std::sin(AHV_Model_X.Integrator1_CSTATE_lj[3]);

    // Trigonometry: '<S229>/cos(theta)' incorporates:
    //   Fcn: '<S190>/T31 '
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/cos(theta)'

    psi_mean_m = std::cos(AHV_Model_X.Integrator1_CSTATE_lj[3]);

    // Trigonometry: '<S229>/sin(theta)' incorporates:
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/sin(theta)'

    rtb_psi_dot = std::sin(AHV_Model_X.Integrator1_CSTATE_lj[4]);

    // Trigonometry: '<S229>/cos(theta)' incorporates:
    //   Fcn: '<S190>/T23'
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/cos(theta)'

    rtb_Row1_e = std::cos(AHV_Model_X.Integrator1_CSTATE_lj[4]);

    // Trigonometry: '<S229>/sin(theta)' incorporates:
    //   Fcn: '<S236>/Fcn'
    //   Fcn: '<S236>/Fcn1'
    //   Fcn: '<S237>/Row1'
    //   Fcn: '<S238>/Row1'
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/sin(theta)'

    rtb_K3_u_tmp = std::sin(AHV_Model_X.Integrator1_CSTATE_lj[5]);

    // Trigonometry: '<S229>/cos(theta)' incorporates:
    //   Fcn: '<S236>/Fcn'
    //   Fcn: '<S236>/Fcn1'
    //   Fcn: '<S237>/Row1'
    //   Fcn: '<S237>/Row2'
    //   Fcn: '<S238>/Row1'
    //   Fcn: '<S238>/Row2'
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S189>/cos(theta)'

    rtb_K2_u_tmp = std::cos(AHV_Model_X.Integrator1_CSTATE_lj[5]);

    // Product: '<S177>/Product1' incorporates:
    //   Fcn: '<S229>/R11'
    //   Fcn: '<S229>/R12'
    //   Fcn: '<S229>/R13'
    //   Fcn: '<S229>/R21 '
    //   Fcn: '<S229>/R22'
    //   Fcn: '<S229>/R23'
    //   Fcn: '<S229>/R31 '
    //   Fcn: '<S229>/R32'
    //   Fcn: '<S229>/R33'
    //   Inport: '<Root>/ahv_fairlead3'
    //   Trigonometry: '<S229>/cos(theta)'
    //   Trigonometry: '<S229>/sin(theta)'

    rtb_K3_u_0[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[3] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[6] = -rtb_psi_dot;
    rtb_K3_u_0[1] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 + -rtb_K3_u_tmp * psi_mean_m;
    rtb_K3_u_0[4] =
        rtb_K3_u_tmp * rtb_psi_dot * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_K3_u_0[7] = rtb_Row1_e * rtb_T33;
    rtb_K3_u_0[2] =
        rtb_K2_u_tmp * rtb_psi_dot * psi_mean_m + rtb_K3_u_tmp * rtb_T33;
    rtb_K3_u_0[5] =
        rtb_K3_u_tmp * rtb_psi_dot * psi_mean_m + -rtb_K2_u_tmp * rtb_T33;
    rtb_K3_u_0[8] = rtb_Row1_e * psi_mean_m;
    for (i = 0; i < 3; i++) {
      rtb_Switch_a = rtb_K3_u_0[i + 6] * ahv_fairlead3[2] +
                     (rtb_K3_u_0[i + 3] * ahv_fairlead3[1] +
                      rtb_K3_u_0[i] * ahv_fairlead3[0]);

      // Gain: '<S176>/Gain1' incorporates:
      //   Inport: '<Root>/ahv_fairlead3'

      rtb_Integrator[i] = AHV_Model_ConstP.pooled60[i] * rtb_Switch_a;
      rtb_Gain2[i] = rtb_Switch_a;
    }

    // End of Product: '<S177>/Product1'

    // Gain: '<S176>/Gain1' incorporates:
    //   Inport: '<Root>/tau_cable3'
    //   Product: '<S228>/ae'
    //   Product: '<S228>/af'
    //   Product: '<S228>/bd'
    //   Product: '<S228>/bf'
    //   Product: '<S228>/cd'
    //   Product: '<S228>/ce'
    //   Sum: '<S228>/Sum'
    //   Sum: '<S228>/Sum1'
    //   Sum: '<S228>/Sum2'

    rtb_Integrator[3] =
        tau_cable3[1] * rtb_Gain2[2] - tau_cable3[2] * rtb_Gain2[1];
    rtb_Integrator[4] =
        (tau_cable3[2] * rtb_Gain2[0] - tau_cable3[0] * rtb_Gain2[2]) * 0.5;
    rtb_Integrator[5] =
        tau_cable3[0] * rtb_Gain2[1] - tau_cable3[1] * rtb_Gain2[0];

    // Fcn: '<S190>/T21 ' incorporates:
    //   Fcn: '<S190>/T31 '
    //   Integrator: '<S181>/Integrator1'

    rtb_Switch_a = std::tan(AHV_Model_X.Integrator1_CSTATE_lj[4]);

    // Reshape: '<S190>/Reshape 9x1->3x3' incorporates:
    //   Constant: '<S190>/Constant'
    //   Constant: '<S190>/Constant '
    //   Fcn: '<S190>/T21 '
    //   Fcn: '<S190>/T23'
    //   Fcn: '<S190>/T31 '
    //   Fcn: '<S190>/T32'
    //   Fcn: '<S190>/T33'

    rtb_TmpSignalConversionAtProduc[0] = 1.0;
    rtb_TmpSignalConversionAtProduc[1] = 0.0;
    rtb_TmpSignalConversionAtProduc[2] = 0.0;
    rtb_TmpSignalConversionAtProduc[3] = rtb_T33 * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[4] = psi_mean_m;
    rtb_TmpSignalConversionAtProduc[5] = rtb_T33 / rtb_Row1_e;
    rtb_TmpSignalConversionAtProduc[6] = psi_mean_m * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[7] = -rtb_T33;
    rtb_TmpSignalConversionAtProduc[8] = psi_mean_m / rtb_Row1_e;

    // SignalConversion generated from: '<S186>/Product' incorporates:
    //   Fcn: '<S189>/R11'
    //   Fcn: '<S189>/R12'
    //   Fcn: '<S189>/R21 '
    //   Fcn: '<S189>/R31 '

    rtb_TmpSignalConversionAtProd_k[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[1] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[2] = -rtb_psi_dot;
    rtb_TmpSignalConversionAtProd_k[3] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 - rtb_K3_u_tmp * psi_mean_m;

    // Fcn: '<S189>/R22' incorporates:
    //   Fcn: '<S189>/R23'

    rtb_Switch_a = rtb_K3_u_tmp * rtb_psi_dot;

    // SignalConversion generated from: '<S186>/Product' incorporates:
    //   Fcn: '<S189>/R13'
    //   Fcn: '<S189>/R22'
    //   Fcn: '<S189>/R23'
    //   Fcn: '<S189>/R32'
    //   Fcn: '<S189>/R33'

    rtb_TmpSignalConversionAtProd_k[4] =
        rtb_Switch_a * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_TmpSignalConversionAtProd_k[5] = rtb_Row1_e * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[6] =
        psi_mean_m * rtb_psi_dot * rtb_K2_u_tmp + rtb_K3_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[7] =
        rtb_Switch_a * psi_mean_m - rtb_K2_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[8] = rtb_Row1_e * psi_mean_m;

    // Sum: '<S181>/Sum6' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Integrator: '<S181>/Integrator'
    //   Product: '<S180>/Product'
    //   Product: '<S180>/Product1'
    //   Trigonometry: '<S180>/cos'
    //   Trigonometry: '<S180>/sin'

    AHV_Model_B.nu_r_d[0] = AHV_Model_X.Integrator_CSTATE_o[0] -
                            std::cos(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_d[1] = AHV_Model_X.Integrator_CSTATE_o[1] -
                            std::sin(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_d[2] = AHV_Model_X.Integrator_CSTATE_o[2];
    AHV_Model_B.nu_r_d[3] = AHV_Model_X.Integrator_CSTATE_o[3];
    AHV_Model_B.nu_r_d[4] = AHV_Model_X.Integrator_CSTATE_o[4];
    AHV_Model_B.nu_r_d[5] = AHV_Model_X.Integrator_CSTATE_o[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S187>/Cross-flow drag trapezoidal
      // integration' Constant: '<S187>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_b,
                                      AHV_Model_B.nu_r_d[1],
                                      AHV_Model_B.nu_r_d[5], &AHV_Model_B.Sum_l,
                                      &AHV_Model_B.Sum2_g, 82.800003);

      // End of Outputs for SubSystem: '<S187>/Cross-flow drag trapezoidal
      // integration'

      // Product: '<S187>/dx1' incorporates:
      //   Constant: '<S187>/2D drag coefficient '
      //   Constant: '<S187>/Transversal area//Lpp'
      //   Constant: '<S187>/rho'
      //   Gain: '<S187>/Gain1'

      AHV_Model_B.dx1_h =
          -0.5 * AHV_Model_B.Sum_l * 1025.0 * 5.4 * 0.69954741847648816;

      // Product: '<S187>/dx2' incorporates:
      //   Constant: '<S187>/2D drag coefficient '
      //   Constant: '<S187>/Transversal area//Lpp'
      //   Constant: '<S187>/rho'
      //   Gain: '<S187>/Gain'

      AHV_Model_B.dx2_oc =
          -0.5 * AHV_Model_B.Sum2_g * 1025.0 * 5.4 * 0.69954741847648816;
    }

    // Product: '<S181>/G*eta' incorporates:
    //   Constant: '<S181>/Spring stiffness'
    //   Integrator: '<S181>/Integrator1'

    for (i = 0; i < 6; i++) {
      rtb_tau_WF[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_tau_WF[i] += AHV_Model_ConstP.pooled48[6 * i_0 + i] *
                         AHV_Model_X.Integrator1_CSTATE_lj[i_0];
      }
    }

    // End of Product: '<S181>/G*eta'

    // Integrator: '<S231>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init3'
    //   SignalConversion generated from: '<S231>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_n != 0) {
      AHV_Model_X.Integrator3_CSTATE_g[0] = Vessel_init3[0];
      AHV_Model_X.Integrator3_CSTATE_g[1] = Vessel_init3[1];
      AHV_Model_X.Integrator3_CSTATE_g[2] = Vessel_init3[5];
    }

    // Saturate: '<S239>/x_Saturation' incorporates:
    //   Integrator: '<S231>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_g[2] > 1.0E+10) {
      rtb_T33 = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_g[2] < -1.0E+10) {
      rtb_T33 = -1.0E+10;
    } else {
      rtb_T33 = AHV_Model_X.Integrator3_CSTATE_g[2];
    }

    // End of Saturate: '<S239>/x_Saturation'

    // Signum: '<S239>/x_Sign'
    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S239>/x_Sign'

    // Gain: '<S239>/pi'
    rtb_psi_dot *= 3.1415926535897931;

    // Sum: '<S239>/Sum1'
    rtb_T33 += rtb_psi_dot;

    // Math: '<S239>/Math Function' incorporates:
    //   Constant: '<S239>/Constant'

    rtb_T33 = rt_remd_snf(rtb_T33, 6.2831853071795862);

    // Sum: '<S239>/Sum'
    rtb_T33 -= rtb_psi_dot;

    // Signum: '<S248>/x_Sign' incorporates:
    //   Fcn: '<S233>/yaw_angle'

    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S248>/x_Sign'

    // Gain: '<S248>/pi'
    rtb_Switch_a = 3.1415926535897931 * rtb_psi_dot;

    // Sum: '<S248>/Sum' incorporates:
    //   Constant: '<S248>/Constant'
    //   Fcn: '<S233>/yaw_angle'
    //   Math: '<S248>/Math Function'
    //   Sum: '<S248>/Sum1'

    rtb_psi_dot =
        rt_remd_snf(rtb_T33 + rtb_Switch_a, 6.2831853071795862) - rtb_Switch_a;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S232>/Chart' incorporates:
      //   Inport: '<Root>/Hold_Position3'
      //   Inport: '<Root>/Vessel_X_Ref3'
      //   Inport: '<Root>/Vessel_Y_Ref3'
      //   Integrator: '<S181>/Integrator1'

      AHV_Model_Chart(Hold_Position3, Vessel_X_Ref3, Vessel_Y_Ref3,
                      AHV_Model_X.Integrator1_CSTATE_lj[0],
                      AHV_Model_X.Integrator1_CSTATE_lj[1],
                      &AHV_Model_B.x_ref_rel_p, &AHV_Model_B.y_ref_rel_g,
                      &AHV_Model_DW.sf_Chart_k);
    }

    // Switch: '<S232>/Switch' incorporates:
    //   Gain: '<S235>/Gain'
    //   Inport: '<Root>/heading_angle_ref3'
    //   Inport: '<Root>/heading_mode3'
    //   Integrator: '<S181>/Integrator1'
    //   Trigonometry: '<S230>/Trigonometric Function'

    if (heading_mode3 > 0.5) {
      AHV_Model_B.RateLimiter_l =
          57.295779513082323 *
          rt_atan2d_snf(AHV_Model_X.Integrator1_CSTATE_lj[1],
                        AHV_Model_X.Integrator1_CSTATE_lj[0]);
    } else {
      AHV_Model_B.RateLimiter_l = heading_angle_ref3;
    }

    // End of Switch: '<S232>/Switch'

    // RateLimiter: '<S232>/Rate Limiter'
    if (!(AHV_Model_DW.LastMajorTime_a == (rtInf))) {
      rtb_Row3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_a;
      rtb_Row1 = rtb_Row3 * 10.0;
      rtb_Switch1 = AHV_Model_B.RateLimiter_l - AHV_Model_DW.PrevY_l;
      if (rtb_Switch1 > rtb_Row1) {
        AHV_Model_B.RateLimiter_l = AHV_Model_DW.PrevY_l + rtb_Row1;
      } else {
        rtb_Row3 *= -10.0;
        if (rtb_Switch1 < rtb_Row3) {
          AHV_Model_B.RateLimiter_l = AHV_Model_DW.PrevY_l + rtb_Row3;
        }
      }
    }

    // End of RateLimiter: '<S232>/Rate Limiter'

    // Gain: '<S242>/Gain1'
    rtb_Switch_a = 0.017453292519943295 * AHV_Model_B.RateLimiter_l;

    // Saturate: '<S243>/x_Saturation'
    if (rtb_Switch_a > 1.0E+10) {
      rtb_Switch_a = 1.0E+10;
    } else {
      if (rtb_Switch_a < -1.0E+10) {
        rtb_Switch_a = -1.0E+10;
      }
    }

    // End of Saturate: '<S243>/x_Saturation'

    // Signum: '<S243>/x_Sign'
    if (rtb_Switch_a < 0.0) {
      psi_mean_m = -1.0;
    } else if (rtb_Switch_a > 0.0) {
      psi_mean_m = 1.0;
    } else if (rtb_Switch_a == 0.0) {
      psi_mean_m = 0.0;
    } else {
      psi_mean_m = (rtNaN);
    }

    // End of Signum: '<S243>/x_Sign'

    // Gain: '<S243>/pi'
    rtb_Switch1 = 3.1415926535897931 * psi_mean_m;

    // Sum: '<S243>/Sum1'
    rtb_Switch_a += rtb_Switch1;

    // Math: '<S243>/Math Function' incorporates:
    //   Constant: '<S243>/Constant'

    rtb_Switch_a = rt_remd_snf(rtb_Switch_a, 6.2831853071795862);

    // Sum: '<S243>/Sum'
    rtb_Switch_a -= rtb_Switch1;

    // Sum: '<S233>/Sum2' incorporates:
    //   Integrator: '<S231>/Integrator3'

    rtb_Gain2[0] =
        AHV_Model_B.x_ref_rel_p - AHV_Model_X.Integrator3_CSTATE_g[0];
    rtb_Gain2[1] =
        AHV_Model_B.y_ref_rel_g - AHV_Model_X.Integrator3_CSTATE_g[1];
    rtb_Gain2[2] = rtb_Switch_a - rtb_T33;

    // Signum: '<S247>/x_Sign' incorporates:
    //   Saturate: '<S247>/x_Saturation'

    if (rtb_Gain2[2] < 0.0) {
      rtb_Switch_a = -1.0;
    } else if (rtb_Gain2[2] > 0.0) {
      rtb_Switch_a = 1.0;
    } else if (rtb_Gain2[2] == 0.0) {
      rtb_Switch_a = 0.0;
    } else {
      rtb_Switch_a = (rtNaN);
    }

    // End of Signum: '<S247>/x_Sign'

    // Gain: '<S247>/pi'
    rtb_Switch1 = 3.1415926535897931 * rtb_Switch_a;

    // Fcn: '<S246>/Row3' incorporates:
    //   Constant: '<S247>/Constant'
    //   Math: '<S247>/Math Function'
    //   Saturate: '<S247>/x_Saturation'
    //   Sum: '<S247>/Sum'
    //   Sum: '<S247>/Sum1'

    AHV_Model_B.Row3_o =
        rt_remd_snf(rtb_Gain2[2] + rtb_Switch1, 6.2831853071795862) -
        rtb_Switch1;

    // Fcn: '<S246>/Row1' incorporates:
    //   Fcn: '<S246>/Row2'

    psi_mean_m = std::sin(rtb_psi_dot);
    rtb_psi_dot = std::cos(rtb_psi_dot);
    rtb_Row1_e = rtb_psi_dot * rtb_Gain2[0] + psi_mean_m * rtb_Gain2[1];

    // If: '<S251>/If' incorporates:
    //   Abs: '<S251>/Abs'
    //   Constant: '<Root>/Constant41'
    //   Inport: '<S254>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_e =
          static_cast<int8_T>(!(std::abs(rtb_Row1_e) <= Featured));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_e) {
      case 0:
        // Outputs for IfAction SubSystem: '<S251>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S254>/Action Port'

        AHV_Model_B.Merge_p = rtb_Row1_e;

        // End of Outputs for SubSystem: '<S251>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S251>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S253>/Action Port'

        AHV_Model_IfActionSubsystem_l(1.0, rtb_Row1_e, &AHV_Model_B.Merge_p,
                                      &AHV_Model_ConstB.IfActionSubsystem_c);

        // End of Outputs for SubSystem: '<S251>/If Action Subsystem'
        break;
    }

    // End of If: '<S251>/If'

    // Fcn: '<S246>/Row2'
    rtb_psi_dot = -psi_mean_m * rtb_Gain2[0] + rtb_psi_dot * rtb_Gain2[1];

    // If: '<S252>/If' incorporates:
    //   Abs: '<S252>/Abs'
    //   Constant: '<Root>/Constant42'
    //   Inport: '<S256>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_k =
          static_cast<int8_T>(!(std::abs(rtb_psi_dot) <= SidePush));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_k) {
      case 0:
        // Outputs for IfAction SubSystem: '<S252>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S256>/Action Port'

        AHV_Model_B.Merge_i = rtb_psi_dot;

        // End of Outputs for SubSystem: '<S252>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S252>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S255>/Action Port'

        AHV_Model_IfActionSubsystem_l(10.0, rtb_psi_dot, &AHV_Model_B.Merge_i,
                                      &AHV_Model_ConstB.IfActionSubsystem_e);

        // End of Outputs for SubSystem: '<S252>/If Action Subsystem'
        break;
    }

    // End of If: '<S252>/If'

    // If: '<S244>/If' incorporates:
    //   Inport: '<Root>/Hold_Position3'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_ef = static_cast<int8_T>(!Hold_Position3);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_ef) {
      case 0:
        // Outputs for IfAction SubSystem: '<S244>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S249>/Action Port'

        AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_p, AHV_Model_B.Merge_i,
                                    AHV_Model_B.Row3_o, rtb_Gain2,
                                    AHV_Model_ConstP.pooled2);

        // End of Outputs for SubSystem: '<S244>/If Action Subsystem'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S244>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S250>/Action Port'

        AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_p, AHV_Model_B.Merge_i,
                                     AHV_Model_B.Row3_o, rtb_Gain2,
                                     AHV_Model_ConstP.pooled1);

        // End of Outputs for SubSystem: '<S244>/If Action Subsystem1'
        break;
    }

    // End of If: '<S244>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S233>/Ki' incorporates:
      //   DiscreteIntegrator: '<S233>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki_o[i] = 0.0;
        AHV_Model_B.Ki_o[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_i[0];
        AHV_Model_B.Ki_o[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_i[1];
        AHV_Model_B.Ki_o[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_i[2];
      }

      // End of Gain: '<S233>/Ki'
    }

    // Sum: '<S233>/Sum1' incorporates:
    //   Gain: '<S233>/Kd'
    //   Integrator: '<S231>/Integrator4'
    //   Sum: '<S233>/Sum3'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = (rtb_Gain2[i] + AHV_Model_B.Ki_o[i]) -
                    ((AHV_Model_ConstP.pooled64[i + 3] *
                          AHV_Model_X.Integrator4_CSTATE_d[1] +
                      AHV_Model_ConstP.pooled64[i] *
                          AHV_Model_X.Integrator4_CSTATE_d[0]) +
                     AHV_Model_ConstP.pooled64[i + 6] *
                         AHV_Model_X.Integrator4_CSTATE_d[2]);
    }

    // End of Sum: '<S233>/Sum1'

    // Saturate: '<S234>/Surge Force Saturation'
    if (rtb_K2_u[0] > 2.0E+6) {
      rtb_psi_dot = 2.0E+6;
    } else if (rtb_K2_u[0] < -2.0E+6) {
      rtb_psi_dot = -2.0E+6;
    } else {
      rtb_psi_dot = rtb_K2_u[0] * MainThrust;
    }

    // End of Saturate: '<S234>/Surge Force Saturation'

    // Saturate: '<S234>/Sway Force Saturation'
    if (rtb_K2_u[1] > 1.5E+6) {
      rtb_Switch_a = 1.5E+6;
    } else if (rtb_K2_u[1] < -1.5E+6) {
      rtb_Switch_a = -1.5E+6;
    } else {
      rtb_Switch_a = rtb_K2_u[1] * FuThrust;
    }

    // End of Saturate: '<S234>/Sway Force Saturation'

    // Saturate: '<S234>/Yaw Moment Saturation'
    if (rtb_K2_u[2] > 2.0E+7) {
      rtb_Switch1 = 2.0E+7;
    } else if (rtb_K2_u[2] < -2.0E+7) {
      rtb_Switch1 = -2.0E+7;
    } else {
      rtb_Switch1 = rtb_K2_u[2] * AzimuthThrust;
    }

    // End of Saturate: '<S234>/Yaw Moment Saturation'

    // Sum: '<S181>/Sum1' incorporates:
    //   Constant: '<S178>/Constant'

    rtb_Row1 = rtb_tau_WF[0];
    tmp = rtb_tau_WF[1];
    tmp_0 = rtb_tau_WF[2];
    tmp_1 = rtb_tau_WF[3];
    tmp_2 = rtb_tau_WF[4];
    tmp_3 = rtb_tau_WF[5];
    rtb_tau_WF[0] = rtb_psi_dot - rtb_Row1;
    rtb_tau_WF[1] = rtb_Switch_a - tmp;
    rtb_tau_WF[2] = 0.0 - tmp_0;
    rtb_tau_WF[3] = 0.0 - tmp_1;
    rtb_tau_WF[4] = 0.0 - tmp_2;
    rtb_tau_WF[5] = rtb_Switch1 - tmp_3;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S176>/Subsystem3'
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

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread,
                depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff,
                rand_freq, rand_dir, AHV_Model_B.Zeta_a_o, AHV_Model_B.Omega_o,
                AHV_Model_B.Phase_k, AHV_Model_B.Wavenum_e, AHV_Model_B.Psi_e,
                &psi_mean_m, &AHV_Model_B.Subsystem3_d,
                &AHV_Model_DW.Subsystem3_d, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S176>/Subsystem3'
    }

    // Outputs for Atomic SubSystem: '<S176>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_X.Integrator1_CSTATE_lj, AHV_Model_B.Psi_e,
                   AHV_Model_B.Wavenum_e, AHV_Model_B.Omega_o,
                   AHV_Model_B.Phase_k, AHV_Model_B.Zeta_a_o, rtb_Sum1_di,
                   rtb_tau_WD, &AHV_Model_B.WaveloadsU0_d,
                   &AHV_Model_ConstB.WaveloadsU0_d,
                   &AHV_Model_DW.WaveloadsU0_d);

    // End of Outputs for SubSystem: '<S176>/Wave loads (U=0)'

    // Sum: '<S181>/Sum7' incorporates:
    //   Sum: '<S176>/Sum'
    //   Sum: '<S181>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Integrator[i] += (rtb_Sum1_di[i] + rtb_tau_WD[i]) + rtb_tau_WF[i];
    }

    // End of Sum: '<S181>/Sum7'

    // Sum: '<S181>/Sum5' incorporates:
    //   Abs: '<S187>/Abs1'
    //   Gain: '<S187>/Gain3'
    //   Product: '<S187>/Product7'
    //   Product: '<S187>/Product8'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] =
        -(std::abs(AHV_Model_B.nu_r_d[0]) * AHV_Model_B.nu_r_d[0] *
          AHV_Model_ConstB.Product5_or) +
        rtb_Row1;
    rtb_Integrator[1] = tmp + AHV_Model_B.dx1_h;
    rtb_Integrator[2] = tmp_0;
    rtb_Integrator[3] = tmp_1;
    rtb_Integrator[4] = tmp_2;
    rtb_Integrator[5] = tmp_3 + AHV_Model_B.dx2_oc;

    // Sum: '<S181>/Sum4' incorporates:
    //   Constant: '<S181>/damping'
    //   Product: '<S181>/D*eta '

    for (i = 0; i < 6; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Row1 +=
            AHV_Model_ConstP.pooled51[6 * i_0 + i] * AHV_Model_B.nu_r_d[i_0];
      }

      rtb_Integrator[i] -= rtb_Row1;
    }

    // End of Sum: '<S181>/Sum4'

    // Product: '<S200>/Product2' incorporates:
    //   Integrator: '<S200>/Integrator'

    rtb_Sum5_0 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Sum5_0 += AHV_Model_ConstB.MathFunction_d[i] *
                    AHV_Model_X.Integrator_CSTATE_h[i];
    }

    // Sum: '<S181>/Sum' incorporates:
    //   Constant: '<S200>/B44_inf'
    //   Constant: '<S200>/D44'
    //   Product: '<S200>/Product2'
    //   Product: '<S200>/Product3'
    //   Product: '<S200>/Product4'
    //   StateSpace: '<S188>/Dp(1,1)'
    //   StateSpace: '<S188>/Dp(2,2)'
    //   StateSpace: '<S188>/Dp(2,4)'
    //   StateSpace: '<S188>/Dp(2,6)'
    //   StateSpace: '<S188>/Dp(3,3)'
    //   StateSpace: '<S188>/Dp(3,5)'
    //   StateSpace: '<S188>/Dp(4,2)'
    //   StateSpace: '<S188>/Dp(4,6)'
    //   StateSpace: '<S188>/Dp(5,3)'
    //   StateSpace: '<S188>/Dp(5,5)'
    //   StateSpace: '<S188>/Dp(6,2)'
    //   StateSpace: '<S188>/Dp(6,4)'
    //   StateSpace: '<S188>/Dp(6,6)'
    //   Sum: '<S188>/Sum'
    //   Sum: '<S188>/Sum1'
    //   Sum: '<S188>/Sum3'
    //   Sum: '<S188>/Sum4'
    //   Sum: '<S188>/Sum5'
    //   Sum: '<S200>/Sum1'
    //   Sum: '<S200>/Sum2'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] =
        rtb_Row1 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE_j[0] +
                        -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE_j[1]) +
                       -26196.671588656758 * AHV_Model_X.Dp11_CSTATE_j[2]) +
                      -69550.629589629229 * AHV_Model_X.Dp11_CSTATE_j[3]) +
                     -6667.9011063139569 * AHV_Model_X.Dp11_CSTATE_j[4]) +
                    28360.700724836512 * AHV_Model_B.nu_r_d[0]);
    rtb_Integrator[1] =
        tmp - (((((((-1.9219252802461751E+6 * AHV_Model_X.Dp22_CSTATE_g[0] +
                     -43506.445868266834 * AHV_Model_X.Dp22_CSTATE_g[1]) +
                    -261966.71588657616 * AHV_Model_X.Dp22_CSTATE_g[2]) +
                   -695506.29589629406 * AHV_Model_X.Dp22_CSTATE_g[3]) +
                  -66679.0110631567 * AHV_Model_X.Dp22_CSTATE_g[4]) +
                 283607.00724836387 * AHV_Model_B.nu_r_d[1]) +
                (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_f[0] +
                     -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_f[1]) +
                    782719.9224317756 * AHV_Model_X.Dp24_CSTATE_f[2]) +
                   906439.02406890341 * AHV_Model_X.Dp24_CSTATE_f[3]) +
                  48702.05056677026 * AHV_Model_X.Dp24_CSTATE_f[4]) +
                 342290.49224698928 * AHV_Model_B.nu_r_d[3])) +
               (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE_i[0] +
                    9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE_i[1]) +
                   1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE_i[2]) +
                  3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE_i[3]) +
                 -523575.51382008579 * AHV_Model_X.Dp26_CSTATE_i[4]) +
                1.4324880870688609E+6 * AHV_Model_B.nu_r_d[5]));
    rtb_Integrator[2] =
        tmp_0 - ((((((-2.8230572248621574E+6 * AHV_Model_X.Dp33_CSTATE_o[0] +
                      2971.0437052436932 * AHV_Model_X.Dp33_CSTATE_o[1]) +
                     -112975.07052200216 * AHV_Model_X.Dp33_CSTATE_o[2]) +
                    -459867.38057091518 * AHV_Model_X.Dp33_CSTATE_o[3]) +
                   35755.997807249405 * AHV_Model_X.Dp33_CSTATE_o[4]) +
                  425986.43085235963 * AHV_Model_B.nu_r_d[2]) +
                 (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE_i[0] +
                      -189169.67892356269 * AHV_Model_X.Dp35_CSTATE_i[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE_i[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE_i[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE_i[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_d[4]));
    rtb_Integrator[3] =
        tmp_1 -
        (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE_f[0] +
               -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE_f[1]) +
              782719.9224317756 * AHV_Model_X.Dp42_CSTATE_f[2]) +
             906439.02406890341 * AHV_Model_X.Dp42_CSTATE_f[3]) +
            48702.05056677026 * AHV_Model_X.Dp42_CSTATE_f[4]) +
           342290.49224698928 * AHV_Model_B.nu_r_d[1]) +
          ((2.3074903324953854E+7 * AHV_Model_B.nu_r_d[3] + rtb_Sum5_0) +
           4.97124674352808E+7 * AHV_Model_B.nu_r_d[3])) +
         (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE_h[0] +
              2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE_h[1]) +
             -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE_h[2]) +
            1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE_h[3]) +
           -70143.531785937739 * AHV_Model_X.Dp46_CSTATE_h[4]) +
          -8.0522169282609783E+6 * AHV_Model_B.nu_r_d[5]));
    rtb_Integrator[4] =
        tmp_2 - ((((((-4.3928932046187744E+7 * AHV_Model_X.Dp53_CSTATE_c[0] +
                      -189169.67892356269 * AHV_Model_X.Dp53_CSTATE_c[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp53_CSTATE_c[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp53_CSTATE_c[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp53_CSTATE_c[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_d[2]) +
                 (((((-1.3510265271416564E+7 * AHV_Model_X.Dp55_CSTATE_f[0] +
                      2.1371969925573766E+9 * AHV_Model_X.Dp55_CSTATE_f[1]) +
                     9.8710379288756877E+7 * AHV_Model_X.Dp55_CSTATE_f[2]) +
                    6.1487736935078049E+8 * AHV_Model_X.Dp55_CSTATE_f[3]) +
                   5.7220438213245332E+8 * AHV_Model_X.Dp55_CSTATE_f[4]) +
                  3.0813057877131975E+8 * AHV_Model_B.nu_r_d[4]));
    rtb_Integrator[5] =
        tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE_n[0] +
                       9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE_n[1]) +
                      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE_n[2]) +
                     3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE_n[3]) +
                    -523575.51382008579 * AHV_Model_X.Dp62_CSTATE_n[4]) +
                   1.4324880870688609E+6 * AHV_Model_B.nu_r_d[1]) +
                  (((((-1.5746048421550707E+6 * AHV_Model_X.Dp64_CSTATE_b[0] +
                       2.9884305406230576E+7 * AHV_Model_X.Dp64_CSTATE_b[1]) +
                      -4.3953040127794016E+6 * AHV_Model_X.Dp64_CSTATE_b[2]) +
                     1.205872403100216E+7 * AHV_Model_X.Dp64_CSTATE_b[3]) +
                    -70143.531785937739 * AHV_Model_X.Dp64_CSTATE_b[4]) +
                   -8.0522169282609783E+6 * AHV_Model_B.nu_r_d[3])) +
                 (((((-8.07757267677547E+8 * AHV_Model_X.Dp66_CSTATE_n[0] +
                      -1.0798841515040772E+7 * AHV_Model_X.Dp66_CSTATE_n[1]) +
                     -1.27285250509213E+8 * AHV_Model_X.Dp66_CSTATE_n[2]) +
                    -3.0027175358715254E+8 * AHV_Model_X.Dp66_CSTATE_n[3]) +
                   1.7969070654437996E+7 * AHV_Model_X.Dp66_CSTATE_n[4]) +
                  1.1943182502908507E+8 * AHV_Model_B.nu_r_d[5]));

    // Product: '<S181>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_l, rtb_Integrator,
                                    AHV_Model_B.Minvtau_a);
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S181>/Integrator1' incorporates:
      //   Integrator: '<S181>/Integrator'
      //   Product: '<S186>/Product'
      //   Product: '<S186>/Product1'

      AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i] =
          rtb_TmpSignalConversionAtProd_k[i + 6] *
              AHV_Model_X.Integrator_CSTATE_o[2] +
          (rtb_TmpSignalConversionAtProd_k[i + 3] *
               AHV_Model_X.Integrator_CSTATE_o[1] +
           rtb_TmpSignalConversionAtProd_k[i] *
               AHV_Model_X.Integrator_CSTATE_o[0]);
      AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i + 3] =
          rtb_TmpSignalConversionAtProduc[i + 6] *
              AHV_Model_X.Integrator_CSTATE_o[5] +
          (rtb_TmpSignalConversionAtProduc[i + 3] *
               AHV_Model_X.Integrator_CSTATE_o[4] +
           rtb_TmpSignalConversionAtProduc[i] *
               AHV_Model_X.Integrator_CSTATE_o[3]);
    }

    // Sum: '<S200>/Sum' incorporates:
    //   Constant: '<S200>/A44'
    //   Constant: '<S200>/B44'
    //   Integrator: '<S200>/Integrator'
    //   Product: '<S200>/Product'
    //   Product: '<S200>/Product1'

    for (i = 0; i < 5; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        rtb_Row1 += AHV_Model_ConstP.pooled52[5 * i_0 + i] *
                    AHV_Model_X.Integrator_CSTATE_h[i_0];
      }

      AHV_Model_B.Sum_e[i] =
          AHV_Model_ConstP.pooled53[i] * AHV_Model_B.nu_r_d[3] + rtb_Row1;
    }

    // End of Sum: '<S200>/Sum'

    // Fcn: '<S184>/Fcn' incorporates:
    //   Integrator: '<S181>/Integrator'

    rtb_Row1 =
        AHV_Model_X.Integrator_CSTATE_o[0] *
            AHV_Model_X.Integrator_CSTATE_o[0] +
        AHV_Model_X.Integrator_CSTATE_o[1] * AHV_Model_X.Integrator_CSTATE_o[1];
    if (rtb_Row1 < 0.0) {
      AHV_Model_B.Fcn_d = -std::sqrt(-rtb_Row1);
    } else {
      AHV_Model_B.Fcn_d = std::sqrt(rtb_Row1);
    }

    // End of Fcn: '<S184>/Fcn'

    // Sum: '<S231>/Sum2' incorporates:
    //   Integrator: '<S181>/Integrator1'
    //   Integrator: '<S231>/Integrator1'
    //   Integrator: '<S231>/Integrator3'
    //   Sum: '<S231>/Sum4'

    psi_mean_m = AHV_Model_X.Integrator1_CSTATE_lj[0] -
                 (AHV_Model_X.Integrator1_CSTATE_m[0] +
                  AHV_Model_X.Integrator3_CSTATE_g[0]);
    rtb_Row1_e = AHV_Model_X.Integrator1_CSTATE_lj[1] -
                 (AHV_Model_X.Integrator1_CSTATE_m[1] +
                  AHV_Model_X.Integrator3_CSTATE_g[1]);
    rtb_Gain2[2] = AHV_Model_X.Integrator1_CSTATE_lj[5] -
                   (AHV_Model_X.Integrator1_CSTATE_m[2] + rtb_T33);

    // Saturate: '<S240>/x_Saturation'
    if (rtb_Gain2[2] > 1.0E+10) {
      rtb_x_dot = 1.0E+10;
    } else if (rtb_Gain2[2] < -1.0E+10) {
      rtb_x_dot = -1.0E+10;
    } else {
      rtb_x_dot = rtb_Gain2[2];
    }

    // End of Saturate: '<S240>/x_Saturation'

    // Signum: '<S240>/x_Sign'
    if (rtb_x_dot < 0.0) {
      rtb_T33 = -1.0;
    } else if (rtb_x_dot > 0.0) {
      rtb_T33 = 1.0;
    } else if (rtb_x_dot == 0.0) {
      rtb_T33 = 0.0;
    } else {
      rtb_T33 = (rtNaN);
    }

    // End of Signum: '<S240>/x_Sign'

    // Gain: '<S240>/pi'
    rtb_T33 *= 3.1415926535897931;

    // Sum: '<S240>/Sum1'
    rtb_x_dot += rtb_T33;

    // Math: '<S240>/Math Function' incorporates:
    //   Constant: '<S240>/Constant'

    rtb_x_dot = rt_remd_snf(rtb_x_dot, 6.2831853071795862);

    // Sum: '<S240>/Sum'
    rtb_Row1 = rtb_x_dot - rtb_T33;

    // Gain: '<S231>/K4' incorporates:
    //   Sum: '<S240>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = AHV_Model_ConstP.pooled109[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled109[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled109[i] * psi_mean_m);
    }

    // End of Gain: '<S231>/K4'

    // Fcn: '<S237>/Row1'
    rtb_x_dot = rtb_K2_u_tmp * rtb_K2_u[0] + rtb_K3_u_tmp * rtb_K2_u[1];

    // Fcn: '<S237>/Row2' incorporates:
    //   Integrator: '<S181>/Integrator1'

    rtb_T33 = -std::sin(AHV_Model_X.Integrator1_CSTATE_lj[5]) * rtb_K2_u[0] +
              rtb_K2_u_tmp * rtb_K2_u[1];

    // Fcn: '<S237>/Row3'
    rtb_Row3 = rtb_K2_u[2];

    // Gain: '<S231>/Gain6' incorporates:
    //   Integrator: '<S231>/Integrator4'

    for (i = 0; i < 3; i++) {
      rtb_Gain2[i] =
          AHV_Model_ConstP.pooled110[i + 6] *
              AHV_Model_X.Integrator4_CSTATE_d[2] +
          (AHV_Model_ConstP.pooled110[i + 3] *
               AHV_Model_X.Integrator4_CSTATE_d[1] +
           AHV_Model_ConstP.pooled110[i] * AHV_Model_X.Integrator4_CSTATE_d[0]);
    }

    // End of Gain: '<S231>/Gain6'

    // Sum: '<S231>/Sum8' incorporates:
    //   Fcn: '<S238>/Row1'
    //   Fcn: '<S238>/Row2'
    //   Integrator: '<S181>/Integrator1'
    //   Integrator: '<S231>/Integrator6'
    //   Sum: '<S231>/Sum1'

    rtb_K2_u[0] = (((rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_e[0] +
                     rtb_K3_u_tmp * AHV_Model_X.Integrator6_CSTATE_e[1]) +
                    rtb_psi_dot) +
                   rtb_x_dot) -
                  rtb_Gain2[0];
    rtb_K2_u[1] = (((-std::sin(AHV_Model_X.Integrator1_CSTATE_lj[5]) *
                         AHV_Model_X.Integrator6_CSTATE_e[0] +
                     rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_e[1]) +
                    rtb_Switch_a) +
                   rtb_T33) -
                  rtb_Gain2[1];
    rtb_K2_u[2] =
        ((AHV_Model_X.Integrator6_CSTATE_e[2] + rtb_Switch1) + rtb_Row3) -
        rtb_Gain2[2];

    // Gain: '<S231>/Gain3'
    for (i = 0; i < 3; i++) {
      AHV_Model_B.M_u_i[i] = 0.0;
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled111[i] * rtb_K2_u[0];
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled111[i + 3] * rtb_K2_u[1];
      AHV_Model_B.M_u_i[i] += AHV_Model_ConstP.pooled111[i + 6] * rtb_K2_u[2];
    }

    // End of Gain: '<S231>/Gain3'
    for (i = 0; i < 3; i++) {
      // Sum: '<S231>/Sum5' incorporates:
      //   Gain: '<S231>/K11'
      //   Integrator: '<S231>/Integrator1'

      AHV_Model_B.psi_WF_c[i] =
          ((AHV_Model_ConstP.pooled112[i + 3] * rtb_Row1_e +
            AHV_Model_ConstP.pooled112[i] * psi_mean_m) +
           AHV_Model_ConstP.pooled112[i + 6] * rtb_Row1) +
          AHV_Model_X.Integrator1_CSTATE_m[i];

      // Sum: '<S231>/Sum6' incorporates:
      //   Gain: '<S231>/Gain1'
      //   Gain: '<S231>/Gain2'
      //   Gain: '<S231>/K12'
      //   Integrator: '<S231>/Integrator1'
      //   Integrator: '<S231>/Integrator2'
      //   Sum: '<S240>/Sum'

      AHV_Model_B.Sum6_b[i] =
          ((AHV_Model_ConstP.pooled113[i + 6] * rtb_Row1 +
            (AHV_Model_ConstP.pooled113[i + 3] * rtb_Row1_e +
             AHV_Model_ConstP.pooled113[i] * psi_mean_m)) -
           (AHV_Model_ConstP.pooled108[i + 6] *
                AHV_Model_X.Integrator2_CSTATE_j[2] +
            (AHV_Model_ConstP.pooled108[i + 3] *
                 AHV_Model_X.Integrator2_CSTATE_j[1] +
             AHV_Model_ConstP.pooled108[i] *
                 AHV_Model_X.Integrator2_CSTATE_j[0]))) -
          ((AHV_Model_ConstP.pooled107[i + 3] *
                AHV_Model_X.Integrator1_CSTATE_m[1] +
            AHV_Model_ConstP.pooled107[i] *
                AHV_Model_X.Integrator1_CSTATE_m[0]) +
           AHV_Model_ConstP.pooled107[i + 6] *
               AHV_Model_X.Integrator1_CSTATE_m[2]);

      // Sum: '<S231>/Sum7' incorporates:
      //   Gain: '<S231>/K3'
      //   Gain: '<S231>/inv(T_b)'
      //   Integrator: '<S231>/Integrator6'
      //   Sum: '<S240>/Sum'

      AHV_Model_B.Sum7_e[i] = (AHV_Model_ConstP.pooled115[i + 6] * rtb_Row1 +
                               (AHV_Model_ConstP.pooled115[i + 3] * rtb_Row1_e +
                                AHV_Model_ConstP.pooled115[i] * psi_mean_m)) -
                              (AHV_Model_ConstP.pooled116[i + 6] *
                                   AHV_Model_X.Integrator6_CSTATE_e[2] +
                               (AHV_Model_ConstP.pooled116[i + 3] *
                                    AHV_Model_X.Integrator6_CSTATE_e[1] +
                                AHV_Model_ConstP.pooled116[i] *
                                    AHV_Model_X.Integrator6_CSTATE_e[0]));

      // Gain: '<S231>/K2' incorporates:
      //   Sum: '<S231>/Sum2'

      rtb_K2_u[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled114[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled114[i] * psi_mean_m);
    }

    // Sum: '<S231>/Sum3' incorporates:
    //   Fcn: '<S236>/Fcn'
    //   Fcn: '<S236>/Fcn1'
    //   Fcn: '<S236>/Fcn2'
    //   Integrator: '<S231>/Integrator4'

    AHV_Model_B.sun_k2_n[0] =
        (rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_d[0] -
         rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_d[1]) +
        rtb_K2_u[0];
    AHV_Model_B.sun_k2_n[1] =
        (rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_d[0] +
         rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_d[1]) +
        rtb_K2_u[1];
    AHV_Model_B.sun_k2_n[2] = rtb_K2_u[2] + AHV_Model_X.Integrator4_CSTATE_d[2];
    for (i = 0; i < 6; i++) {
      // Integrator: '<S265>/Integrator1' incorporates:
      //   Inport: '<Root>/Vessel_init4'

      if (AHV_Model_DW.Integrator1_IWORK_f != 0) {
        AHV_Model_X.Integrator1_CSTATE_b[i] = Vessel_init4[i];
      }

      // Outport: '<Root>/eta_AHV4' incorporates:
      //   Integrator: '<S265>/Integrator1'

      eta_AHV4[i] = AHV_Model_X.Integrator1_CSTATE_b[i];

      // Outport: '<Root>/nu4' incorporates:
      //   Integrator: '<S265>/Integrator'

      nu4[i] = AHV_Model_X.Integrator_CSTATE_a[i];
    }

    // Outport: '<Root>/AHV_speed4' incorporates:
    //   Gain: '<S268>/Gain'
    //   Gain: '<S268>/Gain1'
    //   Rounding: '<S268>/Rounding Function'
    //   TransferFcn: '<S268>/Transfer Fcn'

    AHV_speed4 =
        std::floor(0.2 * AHV_Model_X.TransferFcn_CSTATE_m * 10.0) * 0.1;

    // Gain: '<S264>/Gain' incorporates:
    //   Inport: '<Root>/Current_direction'

    rtb_x_dot = 0.017453292519943295 * Current_direction;

    // Trigonometry: '<S313>/sin(theta)' incorporates:
    //   Fcn: '<S274>/T23'
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/sin(theta)'

    rtb_T33 = std::sin(AHV_Model_X.Integrator1_CSTATE_b[3]);

    // Trigonometry: '<S313>/cos(theta)' incorporates:
    //   Fcn: '<S274>/T31 '
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/cos(theta)'

    psi_mean_m = std::cos(AHV_Model_X.Integrator1_CSTATE_b[3]);

    // Trigonometry: '<S313>/sin(theta)' incorporates:
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/sin(theta)'

    rtb_psi_dot = std::sin(AHV_Model_X.Integrator1_CSTATE_b[4]);

    // Trigonometry: '<S313>/cos(theta)' incorporates:
    //   Fcn: '<S274>/T23'
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/cos(theta)'

    rtb_Row1_e = std::cos(AHV_Model_X.Integrator1_CSTATE_b[4]);

    // Trigonometry: '<S313>/sin(theta)' incorporates:
    //   Fcn: '<S320>/Fcn'
    //   Fcn: '<S320>/Fcn1'
    //   Fcn: '<S321>/Row1'
    //   Fcn: '<S322>/Row1'
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/sin(theta)'

    rtb_K3_u_tmp = std::sin(AHV_Model_X.Integrator1_CSTATE_b[5]);

    // Trigonometry: '<S313>/cos(theta)' incorporates:
    //   Fcn: '<S320>/Fcn'
    //   Fcn: '<S320>/Fcn1'
    //   Fcn: '<S321>/Row1'
    //   Fcn: '<S321>/Row2'
    //   Fcn: '<S322>/Row1'
    //   Fcn: '<S322>/Row2'
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S273>/cos(theta)'

    rtb_K2_u_tmp = std::cos(AHV_Model_X.Integrator1_CSTATE_b[5]);

    // Product: '<S261>/Product1' incorporates:
    //   Fcn: '<S313>/R11'
    //   Fcn: '<S313>/R12'
    //   Fcn: '<S313>/R13'
    //   Fcn: '<S313>/R21 '
    //   Fcn: '<S313>/R22'
    //   Fcn: '<S313>/R23'
    //   Fcn: '<S313>/R31 '
    //   Fcn: '<S313>/R32'
    //   Fcn: '<S313>/R33'
    //   Inport: '<Root>/ahv_fairlead4'
    //   Trigonometry: '<S313>/cos(theta)'
    //   Trigonometry: '<S313>/sin(theta)'

    rtb_K3_u_0[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[3] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_K3_u_0[6] = -rtb_psi_dot;
    rtb_K3_u_0[1] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 + -rtb_K3_u_tmp * psi_mean_m;
    rtb_K3_u_0[4] =
        rtb_K3_u_tmp * rtb_psi_dot * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_K3_u_0[7] = rtb_Row1_e * rtb_T33;
    rtb_K3_u_0[2] =
        rtb_K2_u_tmp * rtb_psi_dot * psi_mean_m + rtb_K3_u_tmp * rtb_T33;
    rtb_K3_u_0[5] =
        rtb_K3_u_tmp * rtb_psi_dot * psi_mean_m + -rtb_K2_u_tmp * rtb_T33;
    rtb_K3_u_0[8] = rtb_Row1_e * psi_mean_m;
    for (i = 0; i < 3; i++) {
      rtb_Switch_a = rtb_K3_u_0[i + 6] * ahv_fairlead4[2] +
                     (rtb_K3_u_0[i + 3] * ahv_fairlead4[1] +
                      rtb_K3_u_0[i] * ahv_fairlead4[0]);

      // Gain: '<S260>/Gain1' incorporates:
      //   Inport: '<Root>/ahv_fairlead4'

      rtb_Integrator[i] = AHV_Model_ConstP.pooled60[i] * rtb_Switch_a;
      rtb_Gain2[i] = rtb_Switch_a;
    }

    // End of Product: '<S261>/Product1'

    // Gain: '<S260>/Gain1' incorporates:
    //   Inport: '<Root>/tau_cable4'
    //   Product: '<S312>/ae'
    //   Product: '<S312>/af'
    //   Product: '<S312>/bd'
    //   Product: '<S312>/bf'
    //   Product: '<S312>/cd'
    //   Product: '<S312>/ce'
    //   Sum: '<S312>/Sum'
    //   Sum: '<S312>/Sum1'
    //   Sum: '<S312>/Sum2'

    rtb_Integrator[3] =
        tau_cable4[1] * rtb_Gain2[2] - tau_cable4[2] * rtb_Gain2[1];
    rtb_Integrator[4] =
        (tau_cable4[2] * rtb_Gain2[0] - tau_cable4[0] * rtb_Gain2[2]) * 0.5;
    rtb_Integrator[5] =
        tau_cable4[0] * rtb_Gain2[1] - tau_cable4[1] * rtb_Gain2[0];

    // Fcn: '<S274>/T21 ' incorporates:
    //   Fcn: '<S274>/T31 '
    //   Integrator: '<S265>/Integrator1'

    rtb_Switch_a = std::tan(AHV_Model_X.Integrator1_CSTATE_b[4]);

    // Reshape: '<S274>/Reshape 9x1->3x3' incorporates:
    //   Constant: '<S274>/Constant'
    //   Constant: '<S274>/Constant '
    //   Fcn: '<S274>/T21 '
    //   Fcn: '<S274>/T23'
    //   Fcn: '<S274>/T31 '
    //   Fcn: '<S274>/T32'
    //   Fcn: '<S274>/T33'

    rtb_TmpSignalConversionAtProduc[0] = 1.0;
    rtb_TmpSignalConversionAtProduc[1] = 0.0;
    rtb_TmpSignalConversionAtProduc[2] = 0.0;
    rtb_TmpSignalConversionAtProduc[3] = rtb_T33 * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[4] = psi_mean_m;
    rtb_TmpSignalConversionAtProduc[5] = rtb_T33 / rtb_Row1_e;
    rtb_TmpSignalConversionAtProduc[6] = psi_mean_m * rtb_Switch_a;
    rtb_TmpSignalConversionAtProduc[7] = -rtb_T33;
    rtb_TmpSignalConversionAtProduc[8] = psi_mean_m / rtb_Row1_e;

    // SignalConversion generated from: '<S270>/Product' incorporates:
    //   Fcn: '<S273>/R11'
    //   Fcn: '<S273>/R12'
    //   Fcn: '<S273>/R21 '
    //   Fcn: '<S273>/R31 '

    rtb_TmpSignalConversionAtProd_k[0] = rtb_K2_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[1] = rtb_K3_u_tmp * rtb_Row1_e;
    rtb_TmpSignalConversionAtProd_k[2] = -rtb_psi_dot;
    rtb_TmpSignalConversionAtProd_k[3] =
        rtb_K2_u_tmp * rtb_psi_dot * rtb_T33 - rtb_K3_u_tmp * psi_mean_m;

    // Fcn: '<S273>/R22' incorporates:
    //   Fcn: '<S273>/R23'

    rtb_Switch_a = rtb_K3_u_tmp * rtb_psi_dot;

    // SignalConversion generated from: '<S270>/Product' incorporates:
    //   Fcn: '<S273>/R13'
    //   Fcn: '<S273>/R22'
    //   Fcn: '<S273>/R23'
    //   Fcn: '<S273>/R32'
    //   Fcn: '<S273>/R33'

    rtb_TmpSignalConversionAtProd_k[4] =
        rtb_Switch_a * rtb_T33 + rtb_K2_u_tmp * psi_mean_m;
    rtb_TmpSignalConversionAtProd_k[5] = rtb_Row1_e * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[6] =
        psi_mean_m * rtb_psi_dot * rtb_K2_u_tmp + rtb_K3_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[7] =
        rtb_Switch_a * psi_mean_m - rtb_K2_u_tmp * rtb_T33;
    rtb_TmpSignalConversionAtProd_k[8] = rtb_Row1_e * psi_mean_m;

    // Sum: '<S265>/Sum6' incorporates:
    //   Inport: '<Root>/Current_speed'
    //   Integrator: '<S265>/Integrator'
    //   Product: '<S264>/Product'
    //   Product: '<S264>/Product1'
    //   Trigonometry: '<S264>/cos'
    //   Trigonometry: '<S264>/sin'

    AHV_Model_B.nu_r_f[0] = AHV_Model_X.Integrator_CSTATE_a[0] -
                            std::cos(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_f[1] = AHV_Model_X.Integrator_CSTATE_a[1] -
                            std::sin(rtb_x_dot) * Current_speed;
    AHV_Model_B.nu_r_f[2] = AHV_Model_X.Integrator_CSTATE_a[2];
    AHV_Model_B.nu_r_f[3] = AHV_Model_X.Integrator_CSTATE_a[3];
    AHV_Model_B.nu_r_f[4] = AHV_Model_X.Integrator_CSTATE_a[4];
    AHV_Model_B.nu_r_f[5] = AHV_Model_X.Integrator_CSTATE_a[5];
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Iterator SubSystem: '<S271>/Cross-flow drag trapezoidal
      // integration' Constant: '<S271>/N'
      Crossflowdragtrapezoidalintegra(20.0, AHV_Model_ConstB.dx_h,
                                      AHV_Model_B.nu_r_f[1],
                                      AHV_Model_B.nu_r_f[5], &AHV_Model_B.Sum_f,
                                      &AHV_Model_B.Sum2, 82.800003);

      // End of Outputs for SubSystem: '<S271>/Cross-flow drag trapezoidal
      // integration'

      // Product: '<S271>/dx1' incorporates:
      //   Constant: '<S271>/2D drag coefficient '
      //   Constant: '<S271>/Transversal area//Lpp'
      //   Constant: '<S271>/rho'
      //   Gain: '<S271>/Gain1'

      AHV_Model_B.dx1_f =
          -0.5 * AHV_Model_B.Sum_f * 1025.0 * 5.4 * 0.69954741847648816;

      // Product: '<S271>/dx2' incorporates:
      //   Constant: '<S271>/2D drag coefficient '
      //   Constant: '<S271>/Transversal area//Lpp'
      //   Constant: '<S271>/rho'
      //   Gain: '<S271>/Gain'

      AHV_Model_B.dx2_j =
          -0.5 * AHV_Model_B.Sum2 * 1025.0 * 5.4 * 0.69954741847648816;
    }

    // Product: '<S265>/G*eta' incorporates:
    //   Constant: '<S265>/Spring stiffness'
    //   Integrator: '<S265>/Integrator1'

    for (i = 0; i < 6; i++) {
      rtb_tau_WF[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_tau_WF[i] += AHV_Model_ConstP.pooled48[6 * i_0 + i] *
                         AHV_Model_X.Integrator1_CSTATE_b[i_0];
      }
    }

    // End of Product: '<S265>/G*eta'

    // Integrator: '<S315>/Integrator3' incorporates:
    //   Inport: '<Root>/Vessel_init4'
    //   SignalConversion generated from: '<S315>/Integrator3'

    if (AHV_Model_DW.Integrator3_IWORK_c != 0) {
      AHV_Model_X.Integrator3_CSTATE_d[0] = Vessel_init4[0];
      AHV_Model_X.Integrator3_CSTATE_d[1] = Vessel_init4[1];
      AHV_Model_X.Integrator3_CSTATE_d[2] = Vessel_init4[5];
    }

    // Saturate: '<S323>/x_Saturation' incorporates:
    //   Integrator: '<S315>/Integrator3'

    if (AHV_Model_X.Integrator3_CSTATE_d[2] > 1.0E+10) {
      rtb_T33 = 1.0E+10;
    } else if (AHV_Model_X.Integrator3_CSTATE_d[2] < -1.0E+10) {
      rtb_T33 = -1.0E+10;
    } else {
      rtb_T33 = AHV_Model_X.Integrator3_CSTATE_d[2];
    }

    // End of Saturate: '<S323>/x_Saturation'

    // Signum: '<S323>/x_Sign'
    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S323>/x_Sign'

    // Gain: '<S323>/pi'
    rtb_psi_dot *= 3.1415926535897931;

    // Sum: '<S323>/Sum1'
    rtb_T33 += rtb_psi_dot;

    // Math: '<S323>/Math Function' incorporates:
    //   Constant: '<S323>/Constant'

    rtb_T33 = rt_remd_snf(rtb_T33, 6.2831853071795862);

    // Sum: '<S323>/Sum'
    rtb_T33 -= rtb_psi_dot;

    // Signum: '<S332>/x_Sign' incorporates:
    //   Fcn: '<S317>/yaw_angle'

    if (rtb_T33 < 0.0) {
      rtb_psi_dot = -1.0;
    } else if (rtb_T33 > 0.0) {
      rtb_psi_dot = 1.0;
    } else if (rtb_T33 == 0.0) {
      rtb_psi_dot = 0.0;
    } else {
      rtb_psi_dot = (rtNaN);
    }

    // End of Signum: '<S332>/x_Sign'

    // Gain: '<S332>/pi'
    rtb_Switch_a = 3.1415926535897931 * rtb_psi_dot;

    // Sum: '<S332>/Sum' incorporates:
    //   Constant: '<S332>/Constant'
    //   Fcn: '<S317>/yaw_angle'
    //   Math: '<S332>/Math Function'
    //   Sum: '<S332>/Sum1'

    rtb_psi_dot =
        rt_remd_snf(rtb_T33 + rtb_Switch_a, 6.2831853071795862) - rtb_Switch_a;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Chart: '<S316>/Chart' incorporates:
      //   Inport: '<Root>/Hold_Position4'
      //   Inport: '<Root>/Vessel_X_Ref4'
      //   Inport: '<Root>/Vessel_Y_Ref4'
      //   Integrator: '<S265>/Integrator1'

      AHV_Model_Chart(Hold_Position4, Vessel_X_Ref4, Vessel_Y_Ref4,
                      AHV_Model_X.Integrator1_CSTATE_b[0],
                      AHV_Model_X.Integrator1_CSTATE_b[1],
                      &AHV_Model_B.x_ref_rel, &AHV_Model_B.y_ref_rel,
                      &AHV_Model_DW.sf_Chart_j);
    }

    // Switch: '<S316>/Switch' incorporates:
    //   Gain: '<S319>/Gain'
    //   Inport: '<Root>/heading_angle_ref4'
    //   Inport: '<Root>/heading_mode4'
    //   Integrator: '<S265>/Integrator1'
    //   Trigonometry: '<S314>/Trigonometric Function'

    if (heading_mode4 > 0.5) {
      AHV_Model_B.RateLimiter_lx =
          57.295779513082323 *
          rt_atan2d_snf(AHV_Model_X.Integrator1_CSTATE_b[1],
                        AHV_Model_X.Integrator1_CSTATE_b[0]);
    } else {
      AHV_Model_B.RateLimiter_lx = heading_angle_ref4;
    }

    // End of Switch: '<S316>/Switch'

    // RateLimiter: '<S316>/Rate Limiter'
    if (!(AHV_Model_DW.LastMajorTime_m == (rtInf))) {
      rtb_Row3 = (&AHV_Model_M)->Timing.t[0] - AHV_Model_DW.LastMajorTime_m;
      rtb_Row1 = rtb_Row3 * 10.0;
      rtb_Switch1 = AHV_Model_B.RateLimiter_lx - AHV_Model_DW.PrevY_d;
      if (rtb_Switch1 > rtb_Row1) {
        AHV_Model_B.RateLimiter_lx = AHV_Model_DW.PrevY_d + rtb_Row1;
      } else {
        rtb_Row3 *= -10.0;
        if (rtb_Switch1 < rtb_Row3) {
          AHV_Model_B.RateLimiter_lx = AHV_Model_DW.PrevY_d + rtb_Row3;
        }
      }
    }

    // End of RateLimiter: '<S316>/Rate Limiter'

    // Gain: '<S326>/Gain1'
    rtb_Switch_a = 0.017453292519943295 * AHV_Model_B.RateLimiter_lx;

    // Saturate: '<S327>/x_Saturation'
    if (rtb_Switch_a > 1.0E+10) {
      rtb_Switch_a = 1.0E+10;
    } else {
      if (rtb_Switch_a < -1.0E+10) {
        rtb_Switch_a = -1.0E+10;
      }
    }

    // End of Saturate: '<S327>/x_Saturation'

    // Signum: '<S327>/x_Sign'
    if (rtb_Switch_a < 0.0) {
      psi_mean_m = -1.0;
    } else if (rtb_Switch_a > 0.0) {
      psi_mean_m = 1.0;
    } else if (rtb_Switch_a == 0.0) {
      psi_mean_m = 0.0;
    } else {
      psi_mean_m = (rtNaN);
    }

    // End of Signum: '<S327>/x_Sign'

    // Gain: '<S327>/pi'
    rtb_Switch1 = 3.1415926535897931 * psi_mean_m;

    // Sum: '<S327>/Sum1'
    rtb_Switch_a += rtb_Switch1;

    // Math: '<S327>/Math Function' incorporates:
    //   Constant: '<S327>/Constant'

    rtb_Switch_a = rt_remd_snf(rtb_Switch_a, 6.2831853071795862);

    // Sum: '<S327>/Sum'
    rtb_Switch_a -= rtb_Switch1;

    // Sum: '<S317>/Sum2' incorporates:
    //   Integrator: '<S315>/Integrator3'

    rtb_Gain2[0] = AHV_Model_B.x_ref_rel - AHV_Model_X.Integrator3_CSTATE_d[0];
    rtb_Gain2[1] = AHV_Model_B.y_ref_rel - AHV_Model_X.Integrator3_CSTATE_d[1];
    rtb_Gain2[2] = rtb_Switch_a - rtb_T33;

    // Signum: '<S331>/x_Sign' incorporates:
    //   Saturate: '<S331>/x_Saturation'

    if (rtb_Gain2[2] < 0.0) {
      rtb_Switch_a = -1.0;
    } else if (rtb_Gain2[2] > 0.0) {
      rtb_Switch_a = 1.0;
    } else if (rtb_Gain2[2] == 0.0) {
      rtb_Switch_a = 0.0;
    } else {
      rtb_Switch_a = (rtNaN);
    }

    // End of Signum: '<S331>/x_Sign'

    // Gain: '<S331>/pi'
    rtb_Switch1 = 3.1415926535897931 * rtb_Switch_a;

    // Fcn: '<S330>/Row3' incorporates:
    //   Constant: '<S331>/Constant'
    //   Math: '<S331>/Math Function'
    //   Saturate: '<S331>/x_Saturation'
    //   Sum: '<S331>/Sum'
    //   Sum: '<S331>/Sum1'

    AHV_Model_B.Row3_l =
        rt_remd_snf(rtb_Gain2[2] + rtb_Switch1, 6.2831853071795862) -
        rtb_Switch1;

    // Fcn: '<S330>/Row1' incorporates:
    //   Fcn: '<S330>/Row2'

    psi_mean_m = std::sin(rtb_psi_dot);
    rtb_psi_dot = std::cos(rtb_psi_dot);
    rtb_Row1_e = rtb_psi_dot * rtb_Gain2[0] + psi_mean_m * rtb_Gain2[1];

    // If: '<S335>/If' incorporates:
    //   Abs: '<S335>/Abs'
    //   Constant: '<Root>/Constant41'
    //   Inport: '<S338>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_f =
          static_cast<int8_T>(!(std::abs(rtb_Row1_e) <= Featured));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_f) {
      case 0:
        // Outputs for IfAction SubSystem: '<S335>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S338>/Action Port'

        AHV_Model_B.Merge_k = rtb_Row1_e;

        // End of Outputs for SubSystem: '<S335>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S335>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S337>/Action Port'

        AHV_Model_IfActionSubsystem_l(1.0, rtb_Row1_e, &AHV_Model_B.Merge_k,
                                      &AHV_Model_ConstB.IfActionSubsystem_fv);

        // End of Outputs for SubSystem: '<S335>/If Action Subsystem'
        break;
    }

    // End of If: '<S335>/If'

    // Fcn: '<S330>/Row2'
    rtb_psi_dot = -psi_mean_m * rtb_Gain2[0] + rtb_psi_dot * rtb_Gain2[1];

    // If: '<S336>/If' incorporates:
    //   Abs: '<S336>/Abs'
    //   Constant: '<Root>/Constant42'
    //   Inport: '<S340>/main_T'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_hp =
          static_cast<int8_T>(!(std::abs(rtb_psi_dot) <= SidePush));
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_hp) {
      case 0:
        // Outputs for IfAction SubSystem: '<S336>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S340>/Action Port'

        AHV_Model_B.Merge_d = rtb_psi_dot;

        // End of Outputs for SubSystem: '<S336>/If Action Subsystem1'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S336>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S339>/Action Port'

        AHV_Model_IfActionSubsystem_l(10.0, rtb_psi_dot, &AHV_Model_B.Merge_d,
                                      &AHV_Model_ConstB.IfActionSubsystem_i);

        // End of Outputs for SubSystem: '<S336>/If Action Subsystem'
        break;
    }

    // End of If: '<S336>/If'

    // If: '<S328>/If' incorporates:
    //   Inport: '<Root>/Hold_Position4'

    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      AHV_Model_DW.If_ActiveSubsystem_i4 = static_cast<int8_T>(!Hold_Position4);
    }

    switch (AHV_Model_DW.If_ActiveSubsystem_i4) {
      case 0:
        // Outputs for IfAction SubSystem: '<S328>/If Action Subsystem'
        // incorporates:
        //   ActionPort: '<S333>/Action Port'

        AHV_Model_IfActionSubsystem(AHV_Model_B.Merge_k, AHV_Model_B.Merge_d,
                                    AHV_Model_B.Row3_l, rtb_Gain2,
                                    AHV_Model_ConstP.pooled2);

        // End of Outputs for SubSystem: '<S328>/If Action Subsystem'
        break;

      case 1:
        // Outputs for IfAction SubSystem: '<S328>/If Action Subsystem1'
        // incorporates:
        //   ActionPort: '<S334>/Action Port'

        AHV_Model_IfActionSubsystem1(AHV_Model_B.Merge_k, AHV_Model_B.Merge_d,
                                     AHV_Model_B.Row3_l, rtb_Gain2,
                                     AHV_Model_ConstP.pooled1);

        // End of Outputs for SubSystem: '<S328>/If Action Subsystem1'
        break;
    }

    // End of If: '<S328>/If'
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Gain: '<S317>/Ki' incorporates:
      //   DiscreteIntegrator: '<S317>/Integrator1'

      for (i = 0; i < 3; i++) {
        AHV_Model_B.Ki_f[i] = 0.0;
        AHV_Model_B.Ki_f[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_in[0];
        AHV_Model_B.Ki_f[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_in[1];
        AHV_Model_B.Ki_f[i] += 0.0 * AHV_Model_DW.Integrator1_DSTATE_in[2];
      }

      // End of Gain: '<S317>/Ki'
    }

    // Sum: '<S317>/Sum1' incorporates:
    //   Gain: '<S317>/Kd'
    //   Integrator: '<S315>/Integrator4'
    //   Sum: '<S317>/Sum3'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = (rtb_Gain2[i] + AHV_Model_B.Ki_f[i]) -
                    ((AHV_Model_ConstP.pooled64[i + 3] *
                          AHV_Model_X.Integrator4_CSTATE_o[1] +
                      AHV_Model_ConstP.pooled64[i] *
                          AHV_Model_X.Integrator4_CSTATE_o[0]) +
                     AHV_Model_ConstP.pooled64[i + 6] *
                         AHV_Model_X.Integrator4_CSTATE_o[2]);
    }

    // End of Sum: '<S317>/Sum1'

    // Saturate: '<S318>/Surge Force Saturation'
    if (rtb_K2_u[0] > 2.0E+6) {
      rtb_psi_dot = 2.0E+6;
    } else if (rtb_K2_u[0] < -2.0E+6) {
      rtb_psi_dot = -2.0E+6;
    } else {
      rtb_psi_dot = rtb_K2_u[0] * MainThrust;
    }

    // End of Saturate: '<S318>/Surge Force Saturation'

    // Saturate: '<S318>/Sway Force Saturation'
    if (rtb_K2_u[1] > 1.5E+6) {
      rtb_Switch_a = 1.5E+6;
    } else if (rtb_K2_u[1] < -1.5E+6) {
      rtb_Switch_a = -1.5E+6;
    } else {
      rtb_Switch_a = rtb_K2_u[1] * FuThrust;
    }

    // End of Saturate: '<S318>/Sway Force Saturation'

    // Saturate: '<S318>/Yaw Moment Saturation'
    if (rtb_K2_u[2] > 2.0E+7) {
      rtb_Switch1 = 2.0E+7;
    } else if (rtb_K2_u[2] < -2.0E+7) {
      rtb_Switch1 = -2.0E+7;
    } else {
      rtb_Switch1 = rtb_K2_u[2] * AzimuthThrust;
    }

    // End of Saturate: '<S318>/Yaw Moment Saturation'

    // Sum: '<S265>/Sum1' incorporates:
    //   Constant: '<S262>/Constant'

    rtb_Row1 = rtb_tau_WF[0];
    tmp = rtb_tau_WF[1];
    tmp_0 = rtb_tau_WF[2];
    tmp_1 = rtb_tau_WF[3];
    tmp_2 = rtb_tau_WF[4];
    tmp_3 = rtb_tau_WF[5];
    rtb_tau_WF[0] = rtb_psi_dot - rtb_Row1;
    rtb_tau_WF[1] = rtb_Switch_a - tmp;
    rtb_tau_WF[2] = 0.0 - tmp_0;
    rtb_tau_WF[3] = 0.0 - tmp_1;
    rtb_tau_WF[4] = 0.0 - tmp_2;
    rtb_tau_WF[5] = rtb_Switch1 - tmp_3;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Outputs for Atomic SubSystem: '<S260>/Subsystem3'
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

      Wave_init(spectrum_type, hs, omega_peak, psi_mean, gamma_value, spread,
                depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff,
                rand_freq, rand_dir, AHV_Model_B.Zeta_a, AHV_Model_B.Omega,
                AHV_Model_B.Phase, AHV_Model_B.Wavenum, AHV_Model_B.Psi,
                &psi_mean_m, &AHV_Model_B.Subsystem3_m,
                &AHV_Model_DW.Subsystem3_m, 20.0, 10.0);

      // End of Outputs for SubSystem: '<S260>/Subsystem3'
    }

    // Outputs for Atomic SubSystem: '<S260>/Wave loads (U=0)'
    Wave_loads_fun(AHV_Model_X.Integrator1_CSTATE_b, AHV_Model_B.Psi,
                   AHV_Model_B.Wavenum, AHV_Model_B.Omega, AHV_Model_B.Phase,
                   AHV_Model_B.Zeta_a, rtb_Sum1_di, rtb_tau_WD,
                   &AHV_Model_B.WaveloadsU0_e, &AHV_Model_ConstB.WaveloadsU0_e,
                   &AHV_Model_DW.WaveloadsU0_e);

    // End of Outputs for SubSystem: '<S260>/Wave loads (U=0)'

    // Sum: '<S265>/Sum7' incorporates:
    //   Sum: '<S260>/Sum'
    //   Sum: '<S265>/Sum3'

    for (i = 0; i < 6; i++) {
      rtb_Integrator[i] += (rtb_Sum1_di[i] + rtb_tau_WD[i]) + rtb_tau_WF[i];
    }

    // End of Sum: '<S265>/Sum7'

    // Sum: '<S265>/Sum5' incorporates:
    //   Abs: '<S271>/Abs1'
    //   Gain: '<S271>/Gain3'
    //   Product: '<S271>/Product7'
    //   Product: '<S271>/Product8'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] = -(std::abs(AHV_Model_B.nu_r_f[0]) *
                          AHV_Model_B.nu_r_f[0] * AHV_Model_ConstB.Product5_f) +
                        rtb_Row1;
    rtb_Integrator[1] = tmp + AHV_Model_B.dx1_f;
    rtb_Integrator[2] = tmp_0;
    rtb_Integrator[3] = tmp_1;
    rtb_Integrator[4] = tmp_2;
    rtb_Integrator[5] = tmp_3 + AHV_Model_B.dx2_j;

    // Sum: '<S265>/Sum4' incorporates:
    //   Constant: '<S265>/damping'
    //   Product: '<S265>/D*eta '

    for (i = 0; i < 6; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Row1 +=
            AHV_Model_ConstP.pooled51[6 * i_0 + i] * AHV_Model_B.nu_r_f[i_0];
      }

      rtb_Integrator[i] -= rtb_Row1;
    }

    // End of Sum: '<S265>/Sum4'

    // Product: '<S284>/Product2' incorporates:
    //   Integrator: '<S284>/Integrator'

    rtb_Sum5_0 = 0.0;
    for (i = 0; i < 5; i++) {
      rtb_Sum5_0 += AHV_Model_ConstB.MathFunction_a[i] *
                    AHV_Model_X.Integrator_CSTATE_bp[i];
    }

    // Sum: '<S265>/Sum' incorporates:
    //   Constant: '<S284>/B44_inf'
    //   Constant: '<S284>/D44'
    //   Product: '<S284>/Product2'
    //   Product: '<S284>/Product3'
    //   Product: '<S284>/Product4'
    //   StateSpace: '<S272>/Dp(1,1)'
    //   StateSpace: '<S272>/Dp(2,2)'
    //   StateSpace: '<S272>/Dp(2,4)'
    //   StateSpace: '<S272>/Dp(2,6)'
    //   StateSpace: '<S272>/Dp(3,3)'
    //   StateSpace: '<S272>/Dp(3,5)'
    //   StateSpace: '<S272>/Dp(4,2)'
    //   StateSpace: '<S272>/Dp(4,6)'
    //   StateSpace: '<S272>/Dp(5,3)'
    //   StateSpace: '<S272>/Dp(5,5)'
    //   StateSpace: '<S272>/Dp(6,2)'
    //   StateSpace: '<S272>/Dp(6,4)'
    //   StateSpace: '<S272>/Dp(6,6)'
    //   Sum: '<S272>/Sum'
    //   Sum: '<S272>/Sum1'
    //   Sum: '<S272>/Sum3'
    //   Sum: '<S272>/Sum4'
    //   Sum: '<S272>/Sum5'
    //   Sum: '<S284>/Sum1'
    //   Sum: '<S284>/Sum2'

    rtb_Row1 = rtb_Integrator[0];
    tmp = rtb_Integrator[1];
    tmp_0 = rtb_Integrator[2];
    tmp_1 = rtb_Integrator[3];
    tmp_2 = rtb_Integrator[4];
    tmp_3 = rtb_Integrator[5];
    rtb_Integrator[0] =
        rtb_Row1 - (((((-192192.528024618 * AHV_Model_X.Dp11_CSTATE_e[0] +
                        -4350.6445868200408 * AHV_Model_X.Dp11_CSTATE_e[1]) +
                       -26196.671588656758 * AHV_Model_X.Dp11_CSTATE_e[2]) +
                      -69550.629589629229 * AHV_Model_X.Dp11_CSTATE_e[3]) +
                     -6667.9011063139569 * AHV_Model_X.Dp11_CSTATE_e[4]) +
                    28360.700724836512 * AHV_Model_B.nu_r_f[0]);
    rtb_Integrator[1] =
        tmp - (((((((-1.9219252802461751E+6 * AHV_Model_X.Dp22_CSTATE_a[0] +
                     -43506.445868266834 * AHV_Model_X.Dp22_CSTATE_a[1]) +
                    -261966.71588657616 * AHV_Model_X.Dp22_CSTATE_a[2]) +
                   -695506.29589629406 * AHV_Model_X.Dp22_CSTATE_a[3]) +
                  -66679.0110631567 * AHV_Model_X.Dp22_CSTATE_a[4]) +
                 283607.00724836387 * AHV_Model_B.nu_r_f[1]) +
                (((((163377.73376110662 * AHV_Model_X.Dp24_CSTATE_d[0] +
                     -2.24222697524951E+6 * AHV_Model_X.Dp24_CSTATE_d[1]) +
                    782719.9224317756 * AHV_Model_X.Dp24_CSTATE_d[2]) +
                   906439.02406890341 * AHV_Model_X.Dp24_CSTATE_d[3]) +
                  48702.05056677026 * AHV_Model_X.Dp24_CSTATE_d[4]) +
                 342290.49224698928 * AHV_Model_B.nu_r_f[3])) +
               (((((-377217.71116477274 * AHV_Model_X.Dp26_CSTATE_f[0] +
                    9.7759770791847613E+6 * AHV_Model_X.Dp26_CSTATE_f[1]) +
                   1.2719330237464791E+6 * AHV_Model_X.Dp26_CSTATE_f[2]) +
                  3.67278239005485E+6 * AHV_Model_X.Dp26_CSTATE_f[3]) +
                 -523575.51382008579 * AHV_Model_X.Dp26_CSTATE_f[4]) +
                1.4324880870688609E+6 * AHV_Model_B.nu_r_f[5]));
    rtb_Integrator[2] =
        tmp_0 - ((((((-2.8230572248621574E+6 * AHV_Model_X.Dp33_CSTATE_a[0] +
                      2971.0437052436932 * AHV_Model_X.Dp33_CSTATE_a[1]) +
                     -112975.07052200216 * AHV_Model_X.Dp33_CSTATE_a[2]) +
                    -459867.38057091518 * AHV_Model_X.Dp33_CSTATE_a[3]) +
                   35755.997807249405 * AHV_Model_X.Dp33_CSTATE_a[4]) +
                  425986.43085235963 * AHV_Model_B.nu_r_f[2]) +
                 (((((-4.3928932046187744E+7 * AHV_Model_X.Dp35_CSTATE_l[0] +
                      -189169.67892356269 * AHV_Model_X.Dp35_CSTATE_l[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp35_CSTATE_l[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp35_CSTATE_l[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp35_CSTATE_l[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_f[4]));
    rtb_Integrator[3] =
        tmp_1 -
        (((((((163377.73376110662 * AHV_Model_X.Dp42_CSTATE_h[0] +
               -2.24222697524951E+6 * AHV_Model_X.Dp42_CSTATE_h[1]) +
              782719.9224317756 * AHV_Model_X.Dp42_CSTATE_h[2]) +
             906439.02406890341 * AHV_Model_X.Dp42_CSTATE_h[3]) +
            48702.05056677026 * AHV_Model_X.Dp42_CSTATE_h[4]) +
           342290.49224698928 * AHV_Model_B.nu_r_f[1]) +
          ((2.3074903324953854E+7 * AHV_Model_B.nu_r_f[3] + rtb_Sum5_0) +
           4.97124674352808E+7 * AHV_Model_B.nu_r_f[3])) +
         (((((-1.5746048421550707E+6 * AHV_Model_X.Dp46_CSTATE_i[0] +
              2.9884305406230576E+7 * AHV_Model_X.Dp46_CSTATE_i[1]) +
             -4.3953040127794016E+6 * AHV_Model_X.Dp46_CSTATE_i[2]) +
            1.205872403100216E+7 * AHV_Model_X.Dp46_CSTATE_i[3]) +
           -70143.531785937739 * AHV_Model_X.Dp46_CSTATE_i[4]) +
          -8.0522169282609783E+6 * AHV_Model_B.nu_r_f[5]));
    rtb_Integrator[4] =
        tmp_2 - ((((((-4.3928932046187744E+7 * AHV_Model_X.Dp53_CSTATE_pz[0] +
                      -189169.67892356269 * AHV_Model_X.Dp53_CSTATE_pz[1]) +
                     1.5146159421991596E+7 * AHV_Model_X.Dp53_CSTATE_pz[2]) +
                    1.2040006787420809E+7 * AHV_Model_X.Dp53_CSTATE_pz[3]) +
                   2.1294858953942191E+6 * AHV_Model_X.Dp53_CSTATE_pz[4]) +
                  6.2116453321941122E+6 * AHV_Model_B.nu_r_f[2]) +
                 (((((-1.3510265271416564E+7 * AHV_Model_X.Dp55_CSTATE_a[0] +
                      2.1371969925573766E+9 * AHV_Model_X.Dp55_CSTATE_a[1]) +
                     9.8710379288756877E+7 * AHV_Model_X.Dp55_CSTATE_a[2]) +
                    6.1487736935078049E+8 * AHV_Model_X.Dp55_CSTATE_a[3]) +
                   5.7220438213245332E+8 * AHV_Model_X.Dp55_CSTATE_a[4]) +
                  3.0813057877131975E+8 * AHV_Model_B.nu_r_f[4]));
    rtb_Integrator[5] =
        tmp_3 - (((((((-377217.71116477274 * AHV_Model_X.Dp62_CSTATE_i[0] +
                       9.7759770791847613E+6 * AHV_Model_X.Dp62_CSTATE_i[1]) +
                      1.2719330237464791E+6 * AHV_Model_X.Dp62_CSTATE_i[2]) +
                     3.67278239005485E+6 * AHV_Model_X.Dp62_CSTATE_i[3]) +
                    -523575.51382008579 * AHV_Model_X.Dp62_CSTATE_i[4]) +
                   1.4324880870688609E+6 * AHV_Model_B.nu_r_f[1]) +
                  (((((-1.5746048421550707E+6 * AHV_Model_X.Dp64_CSTATE_e[0] +
                       2.9884305406230576E+7 * AHV_Model_X.Dp64_CSTATE_e[1]) +
                      -4.3953040127794016E+6 * AHV_Model_X.Dp64_CSTATE_e[2]) +
                     1.205872403100216E+7 * AHV_Model_X.Dp64_CSTATE_e[3]) +
                    -70143.531785937739 * AHV_Model_X.Dp64_CSTATE_e[4]) +
                   -8.0522169282609783E+6 * AHV_Model_B.nu_r_f[3])) +
                 (((((-8.07757267677547E+8 * AHV_Model_X.Dp66_CSTATE_g3[0] +
                      -1.0798841515040772E+7 * AHV_Model_X.Dp66_CSTATE_g3[1]) +
                     -1.27285250509213E+8 * AHV_Model_X.Dp66_CSTATE_g3[2]) +
                    -3.0027175358715254E+8 * AHV_Model_X.Dp66_CSTATE_g3[3]) +
                   1.7969070654437996E+7 * AHV_Model_X.Dp66_CSTATE_g3[4]) +
                  1.1943182502908507E+8 * AHV_Model_B.nu_r_f[5]));

    // Product: '<S265>/Minv*tau'
    rt_mldivide_U1d6x6_U2d_4sw8yi9v(AHV_Model_ConstB.Sum2_e, rtb_Integrator,
                                    AHV_Model_B.Minvtau_l);
    for (i = 0; i < 3; i++) {
      // SignalConversion generated from: '<S265>/Integrator1' incorporates:
      //   Integrator: '<S265>/Integrator'
      //   Product: '<S270>/Product'
      //   Product: '<S270>/Product1'

      AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i] =
          rtb_TmpSignalConversionAtProd_k[i + 6] *
              AHV_Model_X.Integrator_CSTATE_a[2] +
          (rtb_TmpSignalConversionAtProd_k[i + 3] *
               AHV_Model_X.Integrator_CSTATE_a[1] +
           rtb_TmpSignalConversionAtProd_k[i] *
               AHV_Model_X.Integrator_CSTATE_a[0]);
      AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i + 3] =
          rtb_TmpSignalConversionAtProduc[i + 6] *
              AHV_Model_X.Integrator_CSTATE_a[5] +
          (rtb_TmpSignalConversionAtProduc[i + 3] *
               AHV_Model_X.Integrator_CSTATE_a[4] +
           rtb_TmpSignalConversionAtProduc[i] *
               AHV_Model_X.Integrator_CSTATE_a[3]);
    }

    // Sum: '<S284>/Sum' incorporates:
    //   Constant: '<S284>/A44'
    //   Constant: '<S284>/B44'
    //   Integrator: '<S284>/Integrator'
    //   Product: '<S284>/Product'
    //   Product: '<S284>/Product1'

    for (i = 0; i < 5; i++) {
      rtb_Row1 = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        rtb_Row1 += AHV_Model_ConstP.pooled52[5 * i_0 + i] *
                    AHV_Model_X.Integrator_CSTATE_bp[i_0];
      }

      AHV_Model_B.Sum_c[i] =
          AHV_Model_ConstP.pooled53[i] * AHV_Model_B.nu_r_f[3] + rtb_Row1;
    }

    // End of Sum: '<S284>/Sum'

    // Fcn: '<S268>/Fcn' incorporates:
    //   Integrator: '<S265>/Integrator'

    rtb_Row1 =
        AHV_Model_X.Integrator_CSTATE_a[0] *
            AHV_Model_X.Integrator_CSTATE_a[0] +
        AHV_Model_X.Integrator_CSTATE_a[1] * AHV_Model_X.Integrator_CSTATE_a[1];
    if (rtb_Row1 < 0.0) {
      AHV_Model_B.Fcn_h = -std::sqrt(-rtb_Row1);
    } else {
      AHV_Model_B.Fcn_h = std::sqrt(rtb_Row1);
    }

    // End of Fcn: '<S268>/Fcn'

    // Sum: '<S315>/Sum2' incorporates:
    //   Integrator: '<S265>/Integrator1'
    //   Integrator: '<S315>/Integrator1'
    //   Integrator: '<S315>/Integrator3'
    //   Sum: '<S315>/Sum4'

    psi_mean_m = AHV_Model_X.Integrator1_CSTATE_b[0] -
                 (AHV_Model_X.Integrator1_CSTATE_me[0] +
                  AHV_Model_X.Integrator3_CSTATE_d[0]);
    rtb_Row1_e = AHV_Model_X.Integrator1_CSTATE_b[1] -
                 (AHV_Model_X.Integrator1_CSTATE_me[1] +
                  AHV_Model_X.Integrator3_CSTATE_d[1]);
    rtb_Gain2[2] = AHV_Model_X.Integrator1_CSTATE_b[5] -
                   (AHV_Model_X.Integrator1_CSTATE_me[2] + rtb_T33);

    // Saturate: '<S324>/x_Saturation'
    if (rtb_Gain2[2] > 1.0E+10) {
      rtb_x_dot = 1.0E+10;
    } else if (rtb_Gain2[2] < -1.0E+10) {
      rtb_x_dot = -1.0E+10;
    } else {
      rtb_x_dot = rtb_Gain2[2];
    }

    // End of Saturate: '<S324>/x_Saturation'

    // Signum: '<S324>/x_Sign'
    if (rtb_x_dot < 0.0) {
      rtb_T33 = -1.0;
    } else if (rtb_x_dot > 0.0) {
      rtb_T33 = 1.0;
    } else if (rtb_x_dot == 0.0) {
      rtb_T33 = 0.0;
    } else {
      rtb_T33 = (rtNaN);
    }

    // End of Signum: '<S324>/x_Sign'

    // Gain: '<S324>/pi'
    rtb_T33 *= 3.1415926535897931;

    // Sum: '<S324>/Sum1'
    rtb_x_dot += rtb_T33;

    // Math: '<S324>/Math Function' incorporates:
    //   Constant: '<S324>/Constant'

    rtb_x_dot = rt_remd_snf(rtb_x_dot, 6.2831853071795862);

    // Sum: '<S324>/Sum'
    rtb_Row1 = rtb_x_dot - rtb_T33;

    // Gain: '<S315>/K4' incorporates:
    //   Sum: '<S324>/Sum'

    for (i = 0; i < 3; i++) {
      rtb_K2_u[i] = AHV_Model_ConstP.pooled109[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled109[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled109[i] * psi_mean_m);
    }

    // End of Gain: '<S315>/K4'

    // Fcn: '<S321>/Row1'
    rtb_x_dot = rtb_K2_u_tmp * rtb_K2_u[0] + rtb_K3_u_tmp * rtb_K2_u[1];

    // Fcn: '<S321>/Row2' incorporates:
    //   Integrator: '<S265>/Integrator1'

    rtb_T33 = -std::sin(AHV_Model_X.Integrator1_CSTATE_b[5]) * rtb_K2_u[0] +
              rtb_K2_u_tmp * rtb_K2_u[1];

    // Fcn: '<S321>/Row3'
    rtb_Row3 = rtb_K2_u[2];

    // Gain: '<S315>/Gain6' incorporates:
    //   Integrator: '<S315>/Integrator4'

    for (i = 0; i < 3; i++) {
      rtb_Gain2[i] =
          AHV_Model_ConstP.pooled110[i + 6] *
              AHV_Model_X.Integrator4_CSTATE_o[2] +
          (AHV_Model_ConstP.pooled110[i + 3] *
               AHV_Model_X.Integrator4_CSTATE_o[1] +
           AHV_Model_ConstP.pooled110[i] * AHV_Model_X.Integrator4_CSTATE_o[0]);
    }

    // End of Gain: '<S315>/Gain6'

    // Sum: '<S315>/Sum8' incorporates:
    //   Fcn: '<S322>/Row1'
    //   Fcn: '<S322>/Row2'
    //   Integrator: '<S265>/Integrator1'
    //   Integrator: '<S315>/Integrator6'
    //   Sum: '<S315>/Sum1'

    rtb_K2_u[0] = (((rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_f[0] +
                     rtb_K3_u_tmp * AHV_Model_X.Integrator6_CSTATE_f[1]) +
                    rtb_psi_dot) +
                   rtb_x_dot) -
                  rtb_Gain2[0];
    rtb_K2_u[1] = (((-std::sin(AHV_Model_X.Integrator1_CSTATE_b[5]) *
                         AHV_Model_X.Integrator6_CSTATE_f[0] +
                     rtb_K2_u_tmp * AHV_Model_X.Integrator6_CSTATE_f[1]) +
                    rtb_Switch_a) +
                   rtb_T33) -
                  rtb_Gain2[1];
    rtb_K2_u[2] =
        ((AHV_Model_X.Integrator6_CSTATE_f[2] + rtb_Switch1) + rtb_Row3) -
        rtb_Gain2[2];

    // Gain: '<S315>/Gain3'
    for (i = 0; i < 3; i++) {
      AHV_Model_B.M_u_j[i] = 0.0;
      AHV_Model_B.M_u_j[i] += AHV_Model_ConstP.pooled111[i] * rtb_K2_u[0];
      AHV_Model_B.M_u_j[i] += AHV_Model_ConstP.pooled111[i + 3] * rtb_K2_u[1];
      AHV_Model_B.M_u_j[i] += AHV_Model_ConstP.pooled111[i + 6] * rtb_K2_u[2];
    }

    // End of Gain: '<S315>/Gain3'
    for (i = 0; i < 3; i++) {
      // Sum: '<S315>/Sum5' incorporates:
      //   Gain: '<S315>/K11'
      //   Integrator: '<S315>/Integrator1'

      AHV_Model_B.psi_WF_d[i] =
          ((AHV_Model_ConstP.pooled112[i + 3] * rtb_Row1_e +
            AHV_Model_ConstP.pooled112[i] * psi_mean_m) +
           AHV_Model_ConstP.pooled112[i + 6] * rtb_Row1) +
          AHV_Model_X.Integrator1_CSTATE_me[i];

      // Gain: '<S315>/K2' incorporates:
      //   Sum: '<S315>/Sum2'

      rtb_K2_u[i] = AHV_Model_ConstP.pooled114[i + 6] * rtb_Row1 +
                    (AHV_Model_ConstP.pooled114[i + 3] * rtb_Row1_e +
                     AHV_Model_ConstP.pooled114[i] * psi_mean_m);
    }

    // Sum: '<S315>/Sum3' incorporates:
    //   Fcn: '<S320>/Fcn'
    //   Fcn: '<S320>/Fcn1'
    //   Fcn: '<S320>/Fcn2'
    //   Integrator: '<S315>/Integrator4'

    AHV_Model_B.sun_k2_k[0] =
        (rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_o[0] -
         rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_o[1]) +
        rtb_K2_u[0];
    AHV_Model_B.sun_k2_k[1] =
        (rtb_K3_u_tmp * AHV_Model_X.Integrator4_CSTATE_o[0] +
         rtb_K2_u_tmp * AHV_Model_X.Integrator4_CSTATE_o[1]) +
        rtb_K2_u[1];
    AHV_Model_B.sun_k2_k[2] = rtb_K2_u[2] + AHV_Model_X.Integrator4_CSTATE_o[2];
    for (i = 0; i < 3; i++) {
      // Sum: '<S315>/Sum6' incorporates:
      //   Gain: '<S315>/Gain1'
      //   Gain: '<S315>/Gain2'
      //   Gain: '<S315>/K12'
      //   Integrator: '<S315>/Integrator1'
      //   Integrator: '<S315>/Integrator2'
      //   Sum: '<S324>/Sum'

      AHV_Model_B.Sum6_m[i] =
          ((AHV_Model_ConstP.pooled113[i + 6] * rtb_Row1 +
            (AHV_Model_ConstP.pooled113[i + 3] * rtb_Row1_e +
             AHV_Model_ConstP.pooled113[i] * psi_mean_m)) -
           (AHV_Model_ConstP.pooled108[i + 6] *
                AHV_Model_X.Integrator2_CSTATE_h[2] +
            (AHV_Model_ConstP.pooled108[i + 3] *
                 AHV_Model_X.Integrator2_CSTATE_h[1] +
             AHV_Model_ConstP.pooled108[i] *
                 AHV_Model_X.Integrator2_CSTATE_h[0]))) -
          ((AHV_Model_ConstP.pooled107[i + 3] *
                AHV_Model_X.Integrator1_CSTATE_me[1] +
            AHV_Model_ConstP.pooled107[i] *
                AHV_Model_X.Integrator1_CSTATE_me[0]) +
           AHV_Model_ConstP.pooled107[i + 6] *
               AHV_Model_X.Integrator1_CSTATE_me[2]);

      // Sum: '<S315>/Sum7' incorporates:
      //   Gain: '<S315>/K3'
      //   Gain: '<S315>/inv(T_b)'
      //   Integrator: '<S315>/Integrator6'
      //   Sum: '<S324>/Sum'

      AHV_Model_B.Sum7_g[i] = (AHV_Model_ConstP.pooled115[i + 6] * rtb_Row1 +
                               (AHV_Model_ConstP.pooled115[i + 3] * rtb_Row1_e +
                                AHV_Model_ConstP.pooled115[i] * psi_mean_m)) -
                              (AHV_Model_ConstP.pooled116[i + 6] *
                                   AHV_Model_X.Integrator6_CSTATE_f[2] +
                               (AHV_Model_ConstP.pooled116[i + 3] *
                                    AHV_Model_X.Integrator6_CSTATE_f[1] +
                                AHV_Model_ConstP.pooled116[i] *
                                    AHV_Model_X.Integrator6_CSTATE_f[0]));
    }
  }

  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    real_T tmp;
    real_T tmp_0;
    real_T tmp_1;

    // Update for Integrator: '<S13>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK = 0;

    // Update for Integrator: '<S63>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK = 0;

    // Update for RateLimiter: '<S64>/Rate Limiter' incorporates:
    //   RateLimiter: '<S148>/Rate Limiter'
    //   RateLimiter: '<S232>/Rate Limiter'
    //   RateLimiter: '<S316>/Rate Limiter'

    AHV_Model_DW.PrevY = AHV_Model_B.RateLimiter;
    AHV_Model_DW.LastMajorTime = (&AHV_Model_M)->Timing.t[0];

    // Update for Integrator: '<S97>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_a = 0;

    // Update for Integrator: '<S147>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_p = 0;

    // Update for RateLimiter: '<S148>/Rate Limiter'
    AHV_Model_DW.PrevY_p = AHV_Model_B.RateLimiter_j;
    AHV_Model_DW.LastMajorTime_n = AHV_Model_DW.LastMajorTime;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for DiscreteIntegrator: '<S149>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE[2];
      AHV_Model_DW.Integrator1_DSTATE[0] = 0.02 * AHV_Model_B.Merge_a + tmp;
      AHV_Model_DW.Integrator1_DSTATE[1] = 0.02 * AHV_Model_B.Merge_g + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE[2] = 0.02 * AHV_Model_B.Row3_j + tmp_1;

      // Update for DiscreteIntegrator: '<S233>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE_i[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE_i[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE_i[2];
      AHV_Model_DW.Integrator1_DSTATE_i[0] = 0.02 * AHV_Model_B.Merge_p + tmp;
      AHV_Model_DW.Integrator1_DSTATE_i[1] = 0.02 * AHV_Model_B.Merge_i + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE_i[2] = 0.02 * AHV_Model_B.Row3_o + tmp_1;
    }

    // Update for Integrator: '<S181>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_p = 0;

    // Update for Integrator: '<S231>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_n = 0;

    // Update for RateLimiter: '<S232>/Rate Limiter'
    AHV_Model_DW.PrevY_l = AHV_Model_B.RateLimiter_l;
    AHV_Model_DW.LastMajorTime_a = AHV_Model_DW.LastMajorTime;

    // Update for Integrator: '<S265>/Integrator1'
    AHV_Model_DW.Integrator1_IWORK_f = 0;

    // Update for Integrator: '<S315>/Integrator3'
    AHV_Model_DW.Integrator3_IWORK_c = 0;

    // Update for RateLimiter: '<S316>/Rate Limiter'
    AHV_Model_DW.PrevY_d = AHV_Model_B.RateLimiter_lx;
    AHV_Model_DW.LastMajorTime_m = AHV_Model_DW.LastMajorTime;
    if (rtmIsMajorTimeStep((&AHV_Model_M))) {
      // Update for DiscreteIntegrator: '<S317>/Integrator1'
      tmp = AHV_Model_DW.Integrator1_DSTATE_in[0];
      tmp_0 = AHV_Model_DW.Integrator1_DSTATE_in[1];
      tmp_1 = AHV_Model_DW.Integrator1_DSTATE_in[2];
      AHV_Model_DW.Integrator1_DSTATE_in[0] = 0.02 * AHV_Model_B.Merge_k + tmp;
      AHV_Model_DW.Integrator1_DSTATE_in[1] =
          0.02 * AHV_Model_B.Merge_d + tmp_0;
      AHV_Model_DW.Integrator1_DSTATE_in[2] = 0.02 * AHV_Model_B.Row3_l + tmp_1;
    }
  }  // end MajorTimeStep

  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    rt_ertODEUpdateContinuousStates(&(&AHV_Model_M)->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++(&AHV_Model_M)->Timing.clockTick0;
    (&AHV_Model_M)->Timing.t[0] =
        rtsiGetSolverStopTime(&(&AHV_Model_M)->solverInfo);

    {
      // Update absolute timer for sample time: [0.02s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.02, which is
      //  the step size of the task. Size of "clockTick1" ensures timer will not
      //  overflow during the application lifespan selected.

      (&AHV_Model_M)->Timing.clockTick1++;
    }
  }  // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void AH_Model_v1ModelClass::AHV_Model_derivatives() {
  int32_T i;
  XDot_AHV_Model_T* _rtXdot;
  _rtXdot = ((XDot_AHV_Model_T*)(&AHV_Model_M)->derivs);
  for (i = 0; i < 6; i++) {
    // Derivatives for Integrator: '<S13>/Integrator1'
    _rtXdot->Integrator1_CSTATE[i] =
        AHV_Model_B.TmpSignalConversionAtIntegrat_o[i];

    // Derivatives for Integrator: '<S13>/Integrator'
    _rtXdot->Integrator_CSTATE[i] = AHV_Model_B.Minvtau[i];
  }

  // Derivatives for TransferFcn: '<S16>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += -0.2 * AHV_Model_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += AHV_Model_B.Fcn;

  // Derivatives for Integrator: '<S63>/Integrator3'
  _rtXdot->Integrator3_CSTATE[0] = AHV_Model_B.sun_k2[0];

  // Derivatives for Integrator: '<S65>/Integrator'
  _rtXdot->Integrator_CSTATE_b[0] =
      AHV_Model_B.TmpSignalConversionAtIntegrat_p[0];

  // Derivatives for Integrator: '<S63>/Integrator4'
  _rtXdot->Integrator4_CSTATE[0] = AHV_Model_B.M_u[0];

  // Derivatives for Integrator: '<S63>/Integrator3'
  _rtXdot->Integrator3_CSTATE[1] = AHV_Model_B.sun_k2[1];

  // Derivatives for Integrator: '<S65>/Integrator'
  _rtXdot->Integrator_CSTATE_b[1] =
      AHV_Model_B.TmpSignalConversionAtIntegrat_p[1];

  // Derivatives for Integrator: '<S63>/Integrator4'
  _rtXdot->Integrator4_CSTATE[1] = AHV_Model_B.M_u[1];

  // Derivatives for Integrator: '<S63>/Integrator3'
  _rtXdot->Integrator3_CSTATE[2] = AHV_Model_B.sun_k2[2];

  // Derivatives for Integrator: '<S65>/Integrator'
  _rtXdot->Integrator_CSTATE_b[2] =
      AHV_Model_B.TmpSignalConversionAtIntegrat_p[2];

  // Derivatives for Integrator: '<S63>/Integrator4'
  _rtXdot->Integrator4_CSTATE[2] = AHV_Model_B.M_u[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S20>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE[i] = 0.0;

    // Derivatives for Integrator: '<S32>/Integrator'
    _rtXdot->Integrator_CSTATE_l[i] = AHV_Model_B.Sum[i];

    // Derivatives for StateSpace: '<S20>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE[i] = 0.0;

    // Derivatives for StateSpace: '<S20>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S20>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[0] += -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[0] += -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[0] += -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[0] += -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE[4];
  _rtXdot->Dp11_CSTATE[1] += 1.2054566685310679 * AHV_Model_X.Dp11_CSTATE[0];
  _rtXdot->Dp11_CSTATE[1] +=
      -0.00075311499999405624 * AHV_Model_X.Dp11_CSTATE[1];
  _rtXdot->Dp11_CSTATE[1] += -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE[2];
  _rtXdot->Dp11_CSTATE[1] += -0.023103995048734959 * AHV_Model_X.Dp11_CSTATE[3];
  _rtXdot->Dp11_CSTATE[1] +=
      -0.0023400735275004216 * AHV_Model_X.Dp11_CSTATE[4];
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

  // Derivatives for StateSpace: '<S20>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE[0] += -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[0] += -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE[1];
  _rtXdot->Dp22_CSTATE[0] += -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE[2];
  _rtXdot->Dp22_CSTATE[0] += -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE[3];
  _rtXdot->Dp22_CSTATE[0] += -0.095423424399736376 * AHV_Model_X.Dp22_CSTATE[4];
  _rtXdot->Dp22_CSTATE[1] += 1.2054566685310135 * AHV_Model_X.Dp22_CSTATE[0];
  _rtXdot->Dp22_CSTATE[1] +=
      -0.00075311499999501316 * AHV_Model_X.Dp22_CSTATE[1];
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

  // Derivatives for StateSpace: '<S20>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp24_CSTATE[0];
  _rtXdot->Dp24_CSTATE[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE[1];
  _rtXdot->Dp24_CSTATE[0] += -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE[2];
  _rtXdot->Dp24_CSTATE[0] += -0.074692039015985229 * AHV_Model_X.Dp24_CSTATE[3];
  _rtXdot->Dp24_CSTATE[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp24_CSTATE[4];
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

  // Derivatives for StateSpace: '<S20>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp26_CSTATE[0];
  _rtXdot->Dp26_CSTATE[0] += 1.2361995937266759 * AHV_Model_X.Dp26_CSTATE[1];
  _rtXdot->Dp26_CSTATE[0] += 0.012808448456780379 * AHV_Model_X.Dp26_CSTATE[2];
  _rtXdot->Dp26_CSTATE[0] += 0.031758809083141569 * AHV_Model_X.Dp26_CSTATE[3];
  _rtXdot->Dp26_CSTATE[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp26_CSTATE[4];
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

  // Derivatives for StateSpace: '<S20>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE[0] += -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[0] += -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[0] += -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[0] += 0.033960047754095488 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[1] += -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[1] += -1.538297528844471E-6 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[1] +=
      0.00010972576839633958 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[1] +=
      0.00050341083172827721 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[1] +=
      -3.8193652830025795E-5 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[2] += 0.11835430912143848 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[2] += 0.0001097257683747581 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[2] += -0.033603468126584782 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[2] += -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[2] += 0.039518459830251443 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[3] += -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[3] +=
      -0.00050341083159401254 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[3] += 1.4905346203698526 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[3] += -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[3] += 0.067736905936386635 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[4] += 0.033960047754138974 * AHV_Model_X.Dp33_CSTATE[0];
  _rtXdot->Dp33_CSTATE[4] += 3.8193652822306167E-5 * AHV_Model_X.Dp33_CSTATE[1];
  _rtXdot->Dp33_CSTATE[4] += -0.039518459830266792 * AHV_Model_X.Dp33_CSTATE[2];
  _rtXdot->Dp33_CSTATE[4] += 0.067736905936572681 * AHV_Model_X.Dp33_CSTATE[3];
  _rtXdot->Dp33_CSTATE[4] +=
      -0.0072896299361453719 * AHV_Model_X.Dp33_CSTATE[4];
  _rtXdot->Dp33_CSTATE[0] += -3.3182263091979736 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[1] += -0.0034921698713633515 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[2] += 0.13279109186599056 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[3] += -0.54052890870000092 * AHV_Model_B.nu_r[2];
  _rtXdot->Dp33_CSTATE[4] += 0.042027661214546957 * AHV_Model_B.nu_r[2];

  // Derivatives for StateSpace: '<S20>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE[0] += -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[0] += -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE[1];
  _rtXdot->Dp35_CSTATE[0] += 2.1365278323686749 * AHV_Model_X.Dp35_CSTATE[2];
  _rtXdot->Dp35_CSTATE[0] += 1.6362305234027079 * AHV_Model_X.Dp35_CSTATE[3];
  _rtXdot->Dp35_CSTATE[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE[4];
  _rtXdot->Dp35_CSTATE[1] += 0.9998076276438177 * AHV_Model_X.Dp35_CSTATE[0];
  _rtXdot->Dp35_CSTATE[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp35_CSTATE[1];
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

  // Derivatives for StateSpace: '<S20>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp42_CSTATE[0];
  _rtXdot->Dp42_CSTATE[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE[1];
  _rtXdot->Dp42_CSTATE[0] += -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE[2];
  _rtXdot->Dp42_CSTATE[0] += -0.074692039015985229 * AHV_Model_X.Dp42_CSTATE[3];
  _rtXdot->Dp42_CSTATE[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp42_CSTATE[4];
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

  // Derivatives for StateSpace: '<S20>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE[0] += -0.004396104715914141 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[0] += 1.3927195985501903 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[0] += -0.027511960605097124 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[0] += 0.065100301994703638 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[1] += -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[1] += -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[1] += 0.47138330474887147 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[1] += -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[1] += 0.008356128979928646 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[2] += 0.02751196060511972 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[2] += 0.47138330474889245 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[2] += -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[2] += 2.5521376035497743 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[3] += 0.06510030199471134 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[3] += 1.5107242348458856 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[3] += -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[3] += -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[3] += 0.12040468647842413 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp46_CSTATE[0];
  _rtXdot->Dp46_CSTATE[4] += 0.0083561289799246388 * AHV_Model_X.Dp46_CSTATE[1];
  _rtXdot->Dp46_CSTATE[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp46_CSTATE[2];
  _rtXdot->Dp46_CSTATE[4] += -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE[3];
  _rtXdot->Dp46_CSTATE[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp46_CSTATE[4];
  _rtXdot->Dp46_CSTATE[0] += -0.33524288945878539 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[1] += -6.3625492730866693 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[2] += 0.93578679415036736 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[3] += 2.5673752417819116 * AHV_Model_B.nu_r[5];
  _rtXdot->Dp46_CSTATE[4] += 0.014933981938330865 * AHV_Model_B.nu_r[5];

  // Derivatives for StateSpace: '<S20>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE[0] += -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[0] += -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE[1];
  _rtXdot->Dp53_CSTATE[0] += 2.1365278323686749 * AHV_Model_X.Dp53_CSTATE[2];
  _rtXdot->Dp53_CSTATE[0] += 1.6362305234027079 * AHV_Model_X.Dp53_CSTATE[3];
  _rtXdot->Dp53_CSTATE[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE[4];
  _rtXdot->Dp53_CSTATE[1] += 0.9998076276438177 * AHV_Model_X.Dp53_CSTATE[0];
  _rtXdot->Dp53_CSTATE[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp53_CSTATE[1];
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

  // Derivatives for StateSpace: '<S20>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE[0] +=
      -8.7466634083809418E-5 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[0] += 0.70495353490101365 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[0] += 0.0013720743421809071 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[0] += 0.0075482812368147384 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[0] += 0.00709599540037472 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[1] += -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE[0];
  _rtXdot->Dp55_CSTATE[1] += -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE[1];
  _rtXdot->Dp55_CSTATE[1] += -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE[2];
  _rtXdot->Dp55_CSTATE[1] += -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE[3];
  _rtXdot->Dp55_CSTATE[1] += -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE[4];
  _rtXdot->Dp55_CSTATE[2] +=
      -0.0013720743421854747 * AHV_Model_X.Dp55_CSTATE[0];
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

  // Derivatives for StateSpace: '<S20>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp62_CSTATE[0];
  _rtXdot->Dp62_CSTATE[0] += 1.2361995937266759 * AHV_Model_X.Dp62_CSTATE[1];
  _rtXdot->Dp62_CSTATE[0] += 0.012808448456780379 * AHV_Model_X.Dp62_CSTATE[2];
  _rtXdot->Dp62_CSTATE[0] += 0.031758809083141569 * AHV_Model_X.Dp62_CSTATE[3];
  _rtXdot->Dp62_CSTATE[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp62_CSTATE[4];
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

  // Derivatives for StateSpace: '<S20>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE[0] += -0.004396104715914141 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[0] += 1.3927195985501903 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[0] += -0.027511960605097124 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[0] += 0.065100301994703638 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[1] += -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[1] += -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[1] += 0.47138330474887147 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[1] += -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[1] += 0.008356128979928646 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[2] += 0.02751196060511972 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[2] += 0.47138330474889245 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[2] += -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[2] += 2.5521376035497743 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[3] += 0.06510030199471134 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[3] += 1.5107242348458856 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[3] += -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[3] += -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[3] += 0.12040468647842413 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp64_CSTATE[0];
  _rtXdot->Dp64_CSTATE[4] += 0.0083561289799246388 * AHV_Model_X.Dp64_CSTATE[1];
  _rtXdot->Dp64_CSTATE[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp64_CSTATE[2];
  _rtXdot->Dp64_CSTATE[4] += -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE[3];
  _rtXdot->Dp64_CSTATE[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp64_CSTATE[4];
  _rtXdot->Dp64_CSTATE[0] += -0.33524288945878539 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[1] += -6.3625492730866693 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[2] += 0.93578679415036736 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[3] += 2.5673752417819116 * AHV_Model_B.nu_r[3];
  _rtXdot->Dp64_CSTATE[4] += 0.014933981938330865 * AHV_Model_B.nu_r[3];

  // Derivatives for StateSpace: '<S20>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE[0] += -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[0] += -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[0] += -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[0] += -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[0] += 0.069477290559771143 * AHV_Model_X.Dp66_CSTATE[4];
  _rtXdot->Dp66_CSTATE[1] += 1.3217886535783951 * AHV_Model_X.Dp66_CSTATE[0];
  _rtXdot->Dp66_CSTATE[1] +=
      -0.00029199832472599611 * AHV_Model_X.Dp66_CSTATE[1];
  _rtXdot->Dp66_CSTATE[1] +=
      -0.0080078570280257711 * AHV_Model_X.Dp66_CSTATE[2];
  _rtXdot->Dp66_CSTATE[1] += -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE[3];
  _rtXdot->Dp66_CSTATE[1] +=
      0.00098472174210678556 * AHV_Model_X.Dp66_CSTATE[4];
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

  // Derivatives for Integrator: '<S63>/Integrator1'
  _rtXdot->Integrator1_CSTATE_c[0] = AHV_Model_B.Sum6[0];

  // Derivatives for Integrator: '<S63>/Integrator2'
  _rtXdot->Integrator2_CSTATE[0] = AHV_Model_B.psi_WF[0];

  // Derivatives for Integrator: '<S63>/Integrator6'
  _rtXdot->Integrator6_CSTATE[0] = AHV_Model_B.Sum7[0];

  // Derivatives for Integrator: '<S63>/Integrator1'
  _rtXdot->Integrator1_CSTATE_c[1] = AHV_Model_B.Sum6[1];

  // Derivatives for Integrator: '<S63>/Integrator2'
  _rtXdot->Integrator2_CSTATE[1] = AHV_Model_B.psi_WF[1];

  // Derivatives for Integrator: '<S63>/Integrator6'
  _rtXdot->Integrator6_CSTATE[1] = AHV_Model_B.Sum7[1];

  // Derivatives for Integrator: '<S63>/Integrator1'
  _rtXdot->Integrator1_CSTATE_c[2] = AHV_Model_B.Sum6[2];

  // Derivatives for Integrator: '<S63>/Integrator2'
  _rtXdot->Integrator2_CSTATE[2] = AHV_Model_B.psi_WF[2];

  // Derivatives for Integrator: '<S63>/Integrator6'
  _rtXdot->Integrator6_CSTATE[2] = AHV_Model_B.Sum7[2];
  for (i = 0; i < 6; i++) {
    // Derivatives for Integrator: '<S97>/Integrator1'
    _rtXdot->Integrator1_CSTATE_l[i] =
        AHV_Model_B.TmpSignalConversionAtIntegra_o3[i];

    // Derivatives for Integrator: '<S97>/Integrator'
    _rtXdot->Integrator_CSTATE_n[i] = AHV_Model_B.Minvtau_n[i];
  }

  // Derivatives for TransferFcn: '<S100>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_e = 0.0;
  _rtXdot->TransferFcn_CSTATE_e += -0.2 * AHV_Model_X.TransferFcn_CSTATE_e;
  _rtXdot->TransferFcn_CSTATE_e += AHV_Model_B.Fcn_j;

  // Derivatives for Integrator: '<S147>/Integrator3'
  _rtXdot->Integrator3_CSTATE_i[0] = AHV_Model_B.sun_k2_d[0];

  // Derivatives for Integrator: '<S147>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[0] = AHV_Model_B.M_u_f[0];

  // Derivatives for Integrator: '<S147>/Integrator3'
  _rtXdot->Integrator3_CSTATE_i[1] = AHV_Model_B.sun_k2_d[1];

  // Derivatives for Integrator: '<S147>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[1] = AHV_Model_B.M_u_f[1];

  // Derivatives for Integrator: '<S147>/Integrator3'
  _rtXdot->Integrator3_CSTATE_i[2] = AHV_Model_B.sun_k2_d[2];

  // Derivatives for Integrator: '<S147>/Integrator4'
  _rtXdot->Integrator4_CSTATE_p[2] = AHV_Model_B.M_u_f[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S104>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_j[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_m[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_n[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_c[i] = 0.0;

    // Derivatives for Integrator: '<S116>/Integrator'
    _rtXdot->Integrator_CSTATE_m[i] = AHV_Model_B.Sum_j[i];

    // Derivatives for StateSpace: '<S104>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_o[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_c[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S104>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_g[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S104>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_k[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_k[0];
  _rtXdot->Dp11_CSTATE_k[0] +=
      -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_k[1];
  _rtXdot->Dp11_CSTATE_k[0] +=
      -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_k[2];
  _rtXdot->Dp11_CSTATE_k[0] +=
      -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_k[3];
  _rtXdot->Dp11_CSTATE_k[0] +=
      -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_k[4];
  _rtXdot->Dp11_CSTATE_k[1] +=
      1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_k[0];
  _rtXdot->Dp11_CSTATE_k[1] +=
      -0.00075311499999405624 * AHV_Model_X.Dp11_CSTATE_k[1];
  _rtXdot->Dp11_CSTATE_k[1] +=
      -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_k[2];
  _rtXdot->Dp11_CSTATE_k[1] +=
      -0.023103995048734959 * AHV_Model_X.Dp11_CSTATE_k[3];
  _rtXdot->Dp11_CSTATE_k[1] +=
      -0.0023400735275004216 * AHV_Model_X.Dp11_CSTATE_k[4];
  _rtXdot->Dp11_CSTATE_k[2] +=
      -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_k[0];
  _rtXdot->Dp11_CSTATE_k[2] +=
      0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_k[1];
  _rtXdot->Dp11_CSTATE_k[2] +=
      -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_k[2];
  _rtXdot->Dp11_CSTATE_k[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_k[3];
  _rtXdot->Dp11_CSTATE_k[2] +=
      -0.089746811088073059 * AHV_Model_X.Dp11_CSTATE_k[4];
  _rtXdot->Dp11_CSTATE_k[3] +=
      1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_k[0];
  _rtXdot->Dp11_CSTATE_k[3] +=
      -0.023103995048766656 * AHV_Model_X.Dp11_CSTATE_k[1];
  _rtXdot->Dp11_CSTATE_k[3] +=
      1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_k[2];
  _rtXdot->Dp11_CSTATE_k[3] +=
      -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_k[3];
  _rtXdot->Dp11_CSTATE_k[3] +=
      -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_k[4];
  _rtXdot->Dp11_CSTATE_k[4] +=
      -0.095423424399719625 * AHV_Model_X.Dp11_CSTATE_k[0];
  _rtXdot->Dp11_CSTATE_k[4] +=
      0.0023400735275035814 * AHV_Model_X.Dp11_CSTATE_k[1];
  _rtXdot->Dp11_CSTATE_k[4] +=
      -0.089746811088075445 * AHV_Model_X.Dp11_CSTATE_k[2];
  _rtXdot->Dp11_CSTATE_k[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_k[3];
  _rtXdot->Dp11_CSTATE_k[4] +=
      -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_k[4];
  _rtXdot->Dp11_CSTATE_k[0] += -3.3881001186476949 * AHV_Model_B.nu_r_b[0];
  _rtXdot->Dp11_CSTATE_k[1] += 0.07669611088588818 * AHV_Model_B.nu_r_b[0];
  _rtXdot->Dp11_CSTATE_k[2] += -0.46181267830731043 * AHV_Model_B.nu_r_b[0];
  _rtXdot->Dp11_CSTATE_k[3] += 1.226085627712131 * AHV_Model_B.nu_r_b[0];
  _rtXdot->Dp11_CSTATE_k[4] += -0.11754627904442222 * AHV_Model_B.nu_r_b[0];

  // Derivatives for StateSpace: '<S104>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_l[0] +=
      -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE_l[0];
  _rtXdot->Dp22_CSTATE_l[0] +=
      -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE_l[1];
  _rtXdot->Dp22_CSTATE_l[0] +=
      -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE_l[2];
  _rtXdot->Dp22_CSTATE_l[0] +=
      -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE_l[3];
  _rtXdot->Dp22_CSTATE_l[0] +=
      -0.095423424399736376 * AHV_Model_X.Dp22_CSTATE_l[4];
  _rtXdot->Dp22_CSTATE_l[1] +=
      1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_l[0];
  _rtXdot->Dp22_CSTATE_l[1] +=
      -0.00075311499999501316 * AHV_Model_X.Dp22_CSTATE_l[1];
  _rtXdot->Dp22_CSTATE_l[1] +=
      -0.010562919761942575 * AHV_Model_X.Dp22_CSTATE_l[2];
  _rtXdot->Dp22_CSTATE_l[1] +=
      -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE_l[3];
  _rtXdot->Dp22_CSTATE_l[1] +=
      -0.002340073527504722 * AHV_Model_X.Dp22_CSTATE_l[4];
  _rtXdot->Dp22_CSTATE_l[2] +=
      -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE_l[0];
  _rtXdot->Dp22_CSTATE_l[2] +=
      0.010562919761928094 * AHV_Model_X.Dp22_CSTATE_l[1];
  _rtXdot->Dp22_CSTATE_l[2] +=
      -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE_l[2];
  _rtXdot->Dp22_CSTATE_l[2] +=
      -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE_l[3];
  _rtXdot->Dp22_CSTATE_l[2] +=
      -0.089746811088091308 * AHV_Model_X.Dp22_CSTATE_l[4];
  _rtXdot->Dp22_CSTATE_l[3] +=
      1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_l[0];
  _rtXdot->Dp22_CSTATE_l[3] +=
      -0.023103995048724405 * AHV_Model_X.Dp22_CSTATE_l[1];
  _rtXdot->Dp22_CSTATE_l[3] +=
      1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_l[2];
  _rtXdot->Dp22_CSTATE_l[3] +=
      -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE_l[3];
  _rtXdot->Dp22_CSTATE_l[3] +=
      -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE_l[4];
  _rtXdot->Dp22_CSTATE_l[4] +=
      -0.095423424399735959 * AHV_Model_X.Dp22_CSTATE_l[0];
  _rtXdot->Dp22_CSTATE_l[4] +=
      0.00234007352750089 * AHV_Model_X.Dp22_CSTATE_l[1];
  _rtXdot->Dp22_CSTATE_l[4] +=
      -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE_l[2];
  _rtXdot->Dp22_CSTATE_l[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_l[3];
  _rtXdot->Dp22_CSTATE_l[4] +=
      -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE_l[4];
  _rtXdot->Dp22_CSTATE_l[0] += -3.3881001186476967 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp22_CSTATE_l[1] += 0.076696110885761712 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp22_CSTATE_l[2] += -0.46181267830733025 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp22_CSTATE_l[3] += 1.2260856277121315 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp22_CSTATE_l[4] += -0.11754627904444465 * AHV_Model_B.nu_r_b[1];

  // Derivatives for StateSpace: '<S104>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_i[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp24_CSTATE_i[0];
  _rtXdot->Dp24_CSTATE_i[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_i[1];
  _rtXdot->Dp24_CSTATE_i[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_i[2];
  _rtXdot->Dp24_CSTATE_i[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp24_CSTATE_i[3];
  _rtXdot->Dp24_CSTATE_i[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp24_CSTATE_i[4];
  _rtXdot->Dp24_CSTATE_i[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_i[0];
  _rtXdot->Dp24_CSTATE_i[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_i[1];
  _rtXdot->Dp24_CSTATE_i[1] +=
      0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_i[2];
  _rtXdot->Dp24_CSTATE_i[1] +=
      1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_i[3];
  _rtXdot->Dp24_CSTATE_i[1] +=
      0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_i[4];
  _rtXdot->Dp24_CSTATE_i[2] +=
      0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_i[0];
  _rtXdot->Dp24_CSTATE_i[2] +=
      0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_i[1];
  _rtXdot->Dp24_CSTATE_i[2] +=
      -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_i[2];
  _rtXdot->Dp24_CSTATE_i[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_i[3];
  _rtXdot->Dp24_CSTATE_i[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_i[4];
  _rtXdot->Dp24_CSTATE_i[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp24_CSTATE_i[0];
  _rtXdot->Dp24_CSTATE_i[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_i[1];
  _rtXdot->Dp24_CSTATE_i[3] +=
      3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_i[2];
  _rtXdot->Dp24_CSTATE_i[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_i[3];
  _rtXdot->Dp24_CSTATE_i[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_i[4];
  _rtXdot->Dp24_CSTATE_i[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp24_CSTATE_i[0];
  _rtXdot->Dp24_CSTATE_i[4] +=
      0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_i[1];
  _rtXdot->Dp24_CSTATE_i[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_i[2];
  _rtXdot->Dp24_CSTATE_i[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_i[3];
  _rtXdot->Dp24_CSTATE_i[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp24_CSTATE_i[4];
  _rtXdot->Dp24_CSTATE_i[0] += -0.23872786278308805 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp24_CSTATE_i[1] += -3.2763464234276056 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp24_CSTATE_i[2] += 1.1437118751635387 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp24_CSTATE_i[3] += -1.3244904674438165 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp24_CSTATE_i[4] += 0.071163531145327891 * AHV_Model_B.nu_r_b[3];

  // Derivatives for StateSpace: '<S104>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_p[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp26_CSTATE_p[0];
  _rtXdot->Dp26_CSTATE_p[0] +=
      1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_p[1];
  _rtXdot->Dp26_CSTATE_p[0] +=
      0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_p[2];
  _rtXdot->Dp26_CSTATE_p[0] +=
      0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_p[3];
  _rtXdot->Dp26_CSTATE_p[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp26_CSTATE_p[4];
  _rtXdot->Dp26_CSTATE_p[1] +=
      -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_p[0];
  _rtXdot->Dp26_CSTATE_p[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_p[1];
  _rtXdot->Dp26_CSTATE_p[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_p[2];
  _rtXdot->Dp26_CSTATE_p[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_p[3];
  _rtXdot->Dp26_CSTATE_p[1] +=
      0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_p[4];
  _rtXdot->Dp26_CSTATE_p[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp26_CSTATE_p[0];
  _rtXdot->Dp26_CSTATE_p[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_p[1];
  _rtXdot->Dp26_CSTATE_p[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_p[2];
  _rtXdot->Dp26_CSTATE_p[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_p[3];
  _rtXdot->Dp26_CSTATE_p[2] +=
      0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_p[4];
  _rtXdot->Dp26_CSTATE_p[3] +=
      0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_p[0];
  _rtXdot->Dp26_CSTATE_p[3] +=
      0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_p[1];
  _rtXdot->Dp26_CSTATE_p[3] +=
      1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_p[2];
  _rtXdot->Dp26_CSTATE_p[3] +=
      -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_p[3];
  _rtXdot->Dp26_CSTATE_p[3] +=
      1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_p[4];
  _rtXdot->Dp26_CSTATE_p[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp26_CSTATE_p[0];
  _rtXdot->Dp26_CSTATE_p[4] +=
      0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_p[1];
  _rtXdot->Dp26_CSTATE_p[4] +=
      0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_p[2];
  _rtXdot->Dp26_CSTATE_p[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_p[3];
  _rtXdot->Dp26_CSTATE_p[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_p[4];
  _rtXdot->Dp26_CSTATE_p[0] += 0.13167647316600378 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp26_CSTATE_p[1] += 3.4125284827245985 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp26_CSTATE_p[2] += 0.44399732492152644 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp26_CSTATE_p[3] += -1.2820687298457427 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp26_CSTATE_p[4] += -0.18276601298221137 * AHV_Model_B.nu_r_b[5];

  // Derivatives for StateSpace: '<S104>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_m[0] +=
      -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[0] +=
      -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[0] +=
      -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[0] +=
      0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[1] +=
      -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[1] +=
      -1.538297528844471E-6 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[1] +=
      0.00010972576839633958 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[1] +=
      0.00050341083172827721 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[1] +=
      -3.8193652830025795E-5 * AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[2] +=
      0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[2] +=
      0.0001097257683747581 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[2] +=
      -0.033603468126584782 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[2] +=
      -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[2] +=
      0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[3] +=
      -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[3] +=
      -0.00050341083159401254 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[3] +=
      1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[3] +=
      -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[3] +=
      0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[4] +=
      0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_m[0];
  _rtXdot->Dp33_CSTATE_m[4] +=
      3.8193652822306167E-5 * AHV_Model_X.Dp33_CSTATE_m[1];
  _rtXdot->Dp33_CSTATE_m[4] +=
      -0.039518459830266792 * AHV_Model_X.Dp33_CSTATE_m[2];
  _rtXdot->Dp33_CSTATE_m[4] +=
      0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_m[3];
  _rtXdot->Dp33_CSTATE_m[4] +=
      -0.0072896299361453719 * AHV_Model_X.Dp33_CSTATE_m[4];
  _rtXdot->Dp33_CSTATE_m[0] += -3.3182263091979736 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp33_CSTATE_m[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp33_CSTATE_m[2] += 0.13279109186599056 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp33_CSTATE_m[3] += -0.54052890870000092 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp33_CSTATE_m[4] += 0.042027661214546957 * AHV_Model_B.nu_r_b[2];

  // Derivatives for StateSpace: '<S104>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_n[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_n[0];
  _rtXdot->Dp35_CSTATE_n[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_n[1];
  _rtXdot->Dp35_CSTATE_n[0] +=
      2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_n[2];
  _rtXdot->Dp35_CSTATE_n[0] +=
      1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_n[3];
  _rtXdot->Dp35_CSTATE_n[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_n[4];
  _rtXdot->Dp35_CSTATE_n[1] +=
      0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_n[0];
  _rtXdot->Dp35_CSTATE_n[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp35_CSTATE_n[1];
  _rtXdot->Dp35_CSTATE_n[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp35_CSTATE_n[2];
  _rtXdot->Dp35_CSTATE_n[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp35_CSTATE_n[3];
  _rtXdot->Dp35_CSTATE_n[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp35_CSTATE_n[4];
  _rtXdot->Dp35_CSTATE_n[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_n[0];
  _rtXdot->Dp35_CSTATE_n[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp35_CSTATE_n[1];
  _rtXdot->Dp35_CSTATE_n[2] +=
      -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_n[2];
  _rtXdot->Dp35_CSTATE_n[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_n[3];
  _rtXdot->Dp35_CSTATE_n[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_n[4];
  _rtXdot->Dp35_CSTATE_n[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_n[0];
  _rtXdot->Dp35_CSTATE_n[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp35_CSTATE_n[1];
  _rtXdot->Dp35_CSTATE_n[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_n[2];
  _rtXdot->Dp35_CSTATE_n[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_n[3];
  _rtXdot->Dp35_CSTATE_n[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_n[4];
  _rtXdot->Dp35_CSTATE_n[4] +=
      0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_n[0];
  _rtXdot->Dp35_CSTATE_n[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_n[1];
  _rtXdot->Dp35_CSTATE_n[4] +=
      1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_n[2];
  _rtXdot->Dp35_CSTATE_n[4] +=
      3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_n[3];
  _rtXdot->Dp35_CSTATE_n[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_n[4];
  _rtXdot->Dp35_CSTATE_n[0] += -3.5374929885435322 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp35_CSTATE_n[1] += 0.015233386783081164 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp35_CSTATE_n[2] += -1.2196843916515756 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp35_CSTATE_n[3] += -0.96955326725759594 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp35_CSTATE_n[4] += 0.17148246208757431 * AHV_Model_B.nu_r_b[4];

  // Derivatives for StateSpace: '<S104>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_c[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp42_CSTATE_c[0];
  _rtXdot->Dp42_CSTATE_c[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_c[1];
  _rtXdot->Dp42_CSTATE_c[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_c[2];
  _rtXdot->Dp42_CSTATE_c[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp42_CSTATE_c[3];
  _rtXdot->Dp42_CSTATE_c[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp42_CSTATE_c[4];
  _rtXdot->Dp42_CSTATE_c[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_c[0];
  _rtXdot->Dp42_CSTATE_c[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_c[1];
  _rtXdot->Dp42_CSTATE_c[1] +=
      0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_c[2];
  _rtXdot->Dp42_CSTATE_c[1] +=
      1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_c[3];
  _rtXdot->Dp42_CSTATE_c[1] +=
      0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_c[4];
  _rtXdot->Dp42_CSTATE_c[2] +=
      0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_c[0];
  _rtXdot->Dp42_CSTATE_c[2] +=
      0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_c[1];
  _rtXdot->Dp42_CSTATE_c[2] +=
      -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_c[2];
  _rtXdot->Dp42_CSTATE_c[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_c[3];
  _rtXdot->Dp42_CSTATE_c[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_c[4];
  _rtXdot->Dp42_CSTATE_c[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp42_CSTATE_c[0];
  _rtXdot->Dp42_CSTATE_c[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_c[1];
  _rtXdot->Dp42_CSTATE_c[3] +=
      3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_c[2];
  _rtXdot->Dp42_CSTATE_c[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_c[3];
  _rtXdot->Dp42_CSTATE_c[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_c[4];
  _rtXdot->Dp42_CSTATE_c[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp42_CSTATE_c[0];
  _rtXdot->Dp42_CSTATE_c[4] +=
      0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_c[1];
  _rtXdot->Dp42_CSTATE_c[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_c[2];
  _rtXdot->Dp42_CSTATE_c[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_c[3];
  _rtXdot->Dp42_CSTATE_c[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp42_CSTATE_c[4];
  _rtXdot->Dp42_CSTATE_c[0] += -0.23872786278308805 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp42_CSTATE_c[1] += -3.2763464234276056 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp42_CSTATE_c[2] += 1.1437118751635387 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp42_CSTATE_c[3] += -1.3244904674438165 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp42_CSTATE_c[4] += 0.071163531145327891 * AHV_Model_B.nu_r_b[1];

  // Derivatives for StateSpace: '<S104>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_o[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[0] +=
      1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[0] +=
      0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[1] +=
      0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[1] +=
      0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[2] +=
      0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[2] +=
      0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[2] +=
      2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[3] +=
      0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[3] +=
      1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[3] +=
      0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp46_CSTATE_o[0];
  _rtXdot->Dp46_CSTATE_o[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp46_CSTATE_o[1];
  _rtXdot->Dp46_CSTATE_o[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp46_CSTATE_o[2];
  _rtXdot->Dp46_CSTATE_o[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_o[3];
  _rtXdot->Dp46_CSTATE_o[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp46_CSTATE_o[4];
  _rtXdot->Dp46_CSTATE_o[0] += -0.33524288945878539 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp46_CSTATE_o[1] += -6.3625492730866693 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp46_CSTATE_o[2] += 0.93578679415036736 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp46_CSTATE_o[3] += 2.5673752417819116 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp46_CSTATE_o[4] += 0.014933981938330865 * AHV_Model_B.nu_r_b[5];

  // Derivatives for StateSpace: '<S104>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_p[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_p[0];
  _rtXdot->Dp53_CSTATE_p[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_p[1];
  _rtXdot->Dp53_CSTATE_p[0] +=
      2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_p[2];
  _rtXdot->Dp53_CSTATE_p[0] +=
      1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_p[3];
  _rtXdot->Dp53_CSTATE_p[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE_p[4];
  _rtXdot->Dp53_CSTATE_p[1] +=
      0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_p[0];
  _rtXdot->Dp53_CSTATE_p[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp53_CSTATE_p[1];
  _rtXdot->Dp53_CSTATE_p[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp53_CSTATE_p[2];
  _rtXdot->Dp53_CSTATE_p[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp53_CSTATE_p[3];
  _rtXdot->Dp53_CSTATE_p[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp53_CSTATE_p[4];
  _rtXdot->Dp53_CSTATE_p[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_p[0];
  _rtXdot->Dp53_CSTATE_p[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp53_CSTATE_p[1];
  _rtXdot->Dp53_CSTATE_p[2] +=
      -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_p[2];
  _rtXdot->Dp53_CSTATE_p[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_p[3];
  _rtXdot->Dp53_CSTATE_p[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_p[4];
  _rtXdot->Dp53_CSTATE_p[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_p[0];
  _rtXdot->Dp53_CSTATE_p[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp53_CSTATE_p[1];
  _rtXdot->Dp53_CSTATE_p[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_p[2];
  _rtXdot->Dp53_CSTATE_p[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_p[3];
  _rtXdot->Dp53_CSTATE_p[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_p[4];
  _rtXdot->Dp53_CSTATE_p[4] +=
      0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_p[0];
  _rtXdot->Dp53_CSTATE_p[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_p[1];
  _rtXdot->Dp53_CSTATE_p[4] +=
      1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_p[2];
  _rtXdot->Dp53_CSTATE_p[4] +=
      3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_p[3];
  _rtXdot->Dp53_CSTATE_p[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_p[4];
  _rtXdot->Dp53_CSTATE_p[0] += -3.5374929885435322 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp53_CSTATE_p[1] += 0.015233386783081164 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp53_CSTATE_p[2] += -1.2196843916515756 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp53_CSTATE_p[3] += -0.96955326725759594 * AHV_Model_B.nu_r_b[2];
  _rtXdot->Dp53_CSTATE_p[4] += 0.17148246208757431 * AHV_Model_B.nu_r_b[2];

  // Derivatives for StateSpace: '<S104>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_e[0] +=
      -8.7466634083809418E-5 * AHV_Model_X.Dp55_CSTATE_e[0];
  _rtXdot->Dp55_CSTATE_e[0] +=
      0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_e[1];
  _rtXdot->Dp55_CSTATE_e[0] +=
      0.0013720743421809071 * AHV_Model_X.Dp55_CSTATE_e[2];
  _rtXdot->Dp55_CSTATE_e[0] +=
      0.0075482812368147384 * AHV_Model_X.Dp55_CSTATE_e[3];
  _rtXdot->Dp55_CSTATE_e[0] +=
      0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_e[4];
  _rtXdot->Dp55_CSTATE_e[1] +=
      -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_e[0];
  _rtXdot->Dp55_CSTATE_e[1] +=
      -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_e[1];
  _rtXdot->Dp55_CSTATE_e[1] +=
      -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_e[2];
  _rtXdot->Dp55_CSTATE_e[1] +=
      -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_e[3];
  _rtXdot->Dp55_CSTATE_e[1] +=
      -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_e[4];
  _rtXdot->Dp55_CSTATE_e[2] +=
      -0.0013720743421854747 * AHV_Model_X.Dp55_CSTATE_e[0];
  _rtXdot->Dp55_CSTATE_e[2] +=
      -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_e[1];
  _rtXdot->Dp55_CSTATE_e[2] +=
      -0.068186331750790974 * AHV_Model_X.Dp55_CSTATE_e[2];
  _rtXdot->Dp55_CSTATE_e[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_e[3];
  _rtXdot->Dp55_CSTATE_e[2] +=
      -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_e[4];
  _rtXdot->Dp55_CSTATE_e[3] +=
      0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_e[0];
  _rtXdot->Dp55_CSTATE_e[3] +=
      1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_e[1];
  _rtXdot->Dp55_CSTATE_e[3] +=
      4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_e[2];
  _rtXdot->Dp55_CSTATE_e[3] +=
      -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_e[3];
  _rtXdot->Dp55_CSTATE_e[3] +=
      -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_e[4];
  _rtXdot->Dp55_CSTATE_e[4] +=
      0.0070959954003254081 * AHV_Model_X.Dp55_CSTATE_e[0];
  _rtXdot->Dp55_CSTATE_e[4] +=
      1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_e[1];
  _rtXdot->Dp55_CSTATE_e[4] +=
      2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_e[2];
  _rtXdot->Dp55_CSTATE_e[4] +=
      -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_e[3];
  _rtXdot->Dp55_CSTATE_e[4] +=
      -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_e[4];
  _rtXdot->Dp55_CSTATE_e[0] += 0.021904981880023978 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp55_CSTATE_e[1] += 3.4651622640804733 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp55_CSTATE_e[2] += 0.16004490113712252 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp55_CSTATE_e[3] += -0.9969365784863482 * AHV_Model_B.nu_r_b[4];
  _rtXdot->Dp55_CSTATE_e[4] += -0.92774837285087253 * AHV_Model_B.nu_r_b[4];

  // Derivatives for StateSpace: '<S104>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_k[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp62_CSTATE_k[0];
  _rtXdot->Dp62_CSTATE_k[0] +=
      1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_k[1];
  _rtXdot->Dp62_CSTATE_k[0] +=
      0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_k[2];
  _rtXdot->Dp62_CSTATE_k[0] +=
      0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_k[3];
  _rtXdot->Dp62_CSTATE_k[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp62_CSTATE_k[4];
  _rtXdot->Dp62_CSTATE_k[1] +=
      -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_k[0];
  _rtXdot->Dp62_CSTATE_k[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_k[1];
  _rtXdot->Dp62_CSTATE_k[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_k[2];
  _rtXdot->Dp62_CSTATE_k[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_k[3];
  _rtXdot->Dp62_CSTATE_k[1] +=
      0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_k[4];
  _rtXdot->Dp62_CSTATE_k[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp62_CSTATE_k[0];
  _rtXdot->Dp62_CSTATE_k[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_k[1];
  _rtXdot->Dp62_CSTATE_k[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_k[2];
  _rtXdot->Dp62_CSTATE_k[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_k[3];
  _rtXdot->Dp62_CSTATE_k[2] +=
      0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_k[4];
  _rtXdot->Dp62_CSTATE_k[3] +=
      0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_k[0];
  _rtXdot->Dp62_CSTATE_k[3] +=
      0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_k[1];
  _rtXdot->Dp62_CSTATE_k[3] +=
      1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_k[2];
  _rtXdot->Dp62_CSTATE_k[3] +=
      -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_k[3];
  _rtXdot->Dp62_CSTATE_k[3] +=
      1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_k[4];
  _rtXdot->Dp62_CSTATE_k[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp62_CSTATE_k[0];
  _rtXdot->Dp62_CSTATE_k[4] +=
      0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_k[1];
  _rtXdot->Dp62_CSTATE_k[4] +=
      0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_k[2];
  _rtXdot->Dp62_CSTATE_k[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_k[3];
  _rtXdot->Dp62_CSTATE_k[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_k[4];
  _rtXdot->Dp62_CSTATE_k[0] += 0.13167647316600378 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp62_CSTATE_k[1] += 3.4125284827245985 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp62_CSTATE_k[2] += 0.44399732492152644 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp62_CSTATE_k[3] += -1.2820687298457427 * AHV_Model_B.nu_r_b[1];
  _rtXdot->Dp62_CSTATE_k[4] += -0.18276601298221137 * AHV_Model_B.nu_r_b[1];

  // Derivatives for StateSpace: '<S104>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_l[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp64_CSTATE_l[0];
  _rtXdot->Dp64_CSTATE_l[0] +=
      1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_l[1];
  _rtXdot->Dp64_CSTATE_l[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp64_CSTATE_l[2];
  _rtXdot->Dp64_CSTATE_l[0] +=
      0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_l[3];
  _rtXdot->Dp64_CSTATE_l[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp64_CSTATE_l[4];
  _rtXdot->Dp64_CSTATE_l[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_l[0];
  _rtXdot->Dp64_CSTATE_l[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_l[1];
  _rtXdot->Dp64_CSTATE_l[1] +=
      0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_l[2];
  _rtXdot->Dp64_CSTATE_l[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_l[3];
  _rtXdot->Dp64_CSTATE_l[1] +=
      0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_l[4];
  _rtXdot->Dp64_CSTATE_l[2] +=
      0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_l[0];
  _rtXdot->Dp64_CSTATE_l[2] +=
      0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_l[1];
  _rtXdot->Dp64_CSTATE_l[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_l[2];
  _rtXdot->Dp64_CSTATE_l[2] +=
      2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_l[3];
  _rtXdot->Dp64_CSTATE_l[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp64_CSTATE_l[4];
  _rtXdot->Dp64_CSTATE_l[3] +=
      0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_l[0];
  _rtXdot->Dp64_CSTATE_l[3] +=
      1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_l[1];
  _rtXdot->Dp64_CSTATE_l[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_l[2];
  _rtXdot->Dp64_CSTATE_l[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_l[3];
  _rtXdot->Dp64_CSTATE_l[3] +=
      0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_l[4];
  _rtXdot->Dp64_CSTATE_l[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp64_CSTATE_l[0];
  _rtXdot->Dp64_CSTATE_l[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp64_CSTATE_l[1];
  _rtXdot->Dp64_CSTATE_l[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp64_CSTATE_l[2];
  _rtXdot->Dp64_CSTATE_l[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_l[3];
  _rtXdot->Dp64_CSTATE_l[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp64_CSTATE_l[4];
  _rtXdot->Dp64_CSTATE_l[0] += -0.33524288945878539 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp64_CSTATE_l[1] += -6.3625492730866693 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp64_CSTATE_l[2] += 0.93578679415036736 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp64_CSTATE_l[3] += 2.5673752417819116 * AHV_Model_B.nu_r_b[3];
  _rtXdot->Dp64_CSTATE_l[4] += 0.014933981938330865 * AHV_Model_B.nu_r_b[3];

  // Derivatives for StateSpace: '<S104>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_g[0] +=
      -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE_g[0];
  _rtXdot->Dp66_CSTATE_g[0] +=
      -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE_g[1];
  _rtXdot->Dp66_CSTATE_g[0] +=
      -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE_g[2];
  _rtXdot->Dp66_CSTATE_g[0] +=
      -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE_g[3];
  _rtXdot->Dp66_CSTATE_g[0] +=
      0.069477290559771143 * AHV_Model_X.Dp66_CSTATE_g[4];
  _rtXdot->Dp66_CSTATE_g[1] +=
      1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_g[0];
  _rtXdot->Dp66_CSTATE_g[1] +=
      -0.00029199832472599611 * AHV_Model_X.Dp66_CSTATE_g[1];
  _rtXdot->Dp66_CSTATE_g[1] +=
      -0.0080078570280257711 * AHV_Model_X.Dp66_CSTATE_g[2];
  _rtXdot->Dp66_CSTATE_g[1] +=
      -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE_g[3];
  _rtXdot->Dp66_CSTATE_g[1] +=
      0.00098472174210678556 * AHV_Model_X.Dp66_CSTATE_g[4];
  _rtXdot->Dp66_CSTATE_g[2] +=
      -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE_g[0];
  _rtXdot->Dp66_CSTATE_g[2] +=
      0.00800785702802341 * AHV_Model_X.Dp66_CSTATE_g[1];
  _rtXdot->Dp66_CSTATE_g[2] +=
      -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE_g[2];
  _rtXdot->Dp66_CSTATE_g[2] +=
      -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE_g[3];
  _rtXdot->Dp66_CSTATE_g[2] +=
      0.074587161707531782 * AHV_Model_X.Dp66_CSTATE_g[4];
  _rtXdot->Dp66_CSTATE_g[3] +=
      1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_g[0];
  _rtXdot->Dp66_CSTATE_g[3] +=
      -0.015684008067190565 * AHV_Model_X.Dp66_CSTATE_g[1];
  _rtXdot->Dp66_CSTATE_g[3] +=
      1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_g[2];
  _rtXdot->Dp66_CSTATE_g[3] +=
      -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE_g[3];
  _rtXdot->Dp66_CSTATE_g[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE_g[4];
  _rtXdot->Dp66_CSTATE_g[4] +=
      0.069477290559780386 * AHV_Model_X.Dp66_CSTATE_g[0];
  _rtXdot->Dp66_CSTATE_g[4] +=
      -0.000984721742106249 * AHV_Model_X.Dp66_CSTATE_g[1];
  _rtXdot->Dp66_CSTATE_g[4] +=
      0.074587161707547089 * AHV_Model_X.Dp66_CSTATE_g[2];
  _rtXdot->Dp66_CSTATE_g[4] +=
      -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE_g[3];
  _rtXdot->Dp66_CSTATE_g[4] +=
      -0.061420217921035136 * AHV_Model_X.Dp66_CSTATE_g[4];
  _rtXdot->Dp66_CSTATE_g[0] += -3.3844821061206916 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp66_CSTATE_g[1] += 0.045246867266885989 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp66_CSTATE_g[2] += -0.53332191483729374 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp66_CSTATE_g[3] += 1.2581308985447515 * AHV_Model_B.nu_r_b[5];
  _rtXdot->Dp66_CSTATE_g[4] += 0.075289942321940723 * AHV_Model_B.nu_r_b[5];

  // Derivatives for Integrator: '<S147>/Integrator1'
  _rtXdot->Integrator1_CSTATE_d[0] = AHV_Model_B.Sum6_i[0];

  // Derivatives for Integrator: '<S147>/Integrator2'
  _rtXdot->Integrator2_CSTATE_i[0] = AHV_Model_B.psi_WF_g[0];

  // Derivatives for Integrator: '<S147>/Integrator6'
  _rtXdot->Integrator6_CSTATE_b[0] = AHV_Model_B.Sum7_k[0];

  // Derivatives for Integrator: '<S147>/Integrator1'
  _rtXdot->Integrator1_CSTATE_d[1] = AHV_Model_B.Sum6_i[1];

  // Derivatives for Integrator: '<S147>/Integrator2'
  _rtXdot->Integrator2_CSTATE_i[1] = AHV_Model_B.psi_WF_g[1];

  // Derivatives for Integrator: '<S147>/Integrator6'
  _rtXdot->Integrator6_CSTATE_b[1] = AHV_Model_B.Sum7_k[1];

  // Derivatives for Integrator: '<S147>/Integrator1'
  _rtXdot->Integrator1_CSTATE_d[2] = AHV_Model_B.Sum6_i[2];

  // Derivatives for Integrator: '<S147>/Integrator2'
  _rtXdot->Integrator2_CSTATE_i[2] = AHV_Model_B.psi_WF_g[2];

  // Derivatives for Integrator: '<S147>/Integrator6'
  _rtXdot->Integrator6_CSTATE_b[2] = AHV_Model_B.Sum7_k[2];
  for (i = 0; i < 6; i++) {
    // Derivatives for Integrator: '<S181>/Integrator1'
    _rtXdot->Integrator1_CSTATE_lj[i] =
        AHV_Model_B.TmpSignalConversionAtIntegr_o3s[i];

    // Derivatives for Integrator: '<S181>/Integrator'
    _rtXdot->Integrator_CSTATE_o[i] = AHV_Model_B.Minvtau_a[i];
  }

  // Derivatives for TransferFcn: '<S184>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_l = 0.0;
  _rtXdot->TransferFcn_CSTATE_l += -0.2 * AHV_Model_X.TransferFcn_CSTATE_l;
  _rtXdot->TransferFcn_CSTATE_l += AHV_Model_B.Fcn_d;

  // Derivatives for Integrator: '<S231>/Integrator3'
  _rtXdot->Integrator3_CSTATE_g[0] = AHV_Model_B.sun_k2_n[0];

  // Derivatives for Integrator: '<S231>/Integrator4'
  _rtXdot->Integrator4_CSTATE_d[0] = AHV_Model_B.M_u_i[0];

  // Derivatives for Integrator: '<S231>/Integrator3'
  _rtXdot->Integrator3_CSTATE_g[1] = AHV_Model_B.sun_k2_n[1];

  // Derivatives for Integrator: '<S231>/Integrator4'
  _rtXdot->Integrator4_CSTATE_d[1] = AHV_Model_B.M_u_i[1];

  // Derivatives for Integrator: '<S231>/Integrator3'
  _rtXdot->Integrator3_CSTATE_g[2] = AHV_Model_B.sun_k2_n[2];

  // Derivatives for Integrator: '<S231>/Integrator4'
  _rtXdot->Integrator4_CSTATE_d[2] = AHV_Model_B.M_u_i[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S188>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_j[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_p[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_g[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_o[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_f[i] = 0.0;

    // Derivatives for Integrator: '<S200>/Integrator'
    _rtXdot->Integrator_CSTATE_h[i] = AHV_Model_B.Sum_e[i];

    // Derivatives for StateSpace: '<S188>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_h[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_c[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_n[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_b[i] = 0.0;

    // Derivatives for StateSpace: '<S188>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_n[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S188>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_j[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_j[0];
  _rtXdot->Dp11_CSTATE_j[0] +=
      -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_j[1];
  _rtXdot->Dp11_CSTATE_j[0] +=
      -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_j[2];
  _rtXdot->Dp11_CSTATE_j[0] +=
      -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_j[3];
  _rtXdot->Dp11_CSTATE_j[0] +=
      -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_j[4];
  _rtXdot->Dp11_CSTATE_j[1] +=
      1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_j[0];
  _rtXdot->Dp11_CSTATE_j[1] +=
      -0.00075311499999405624 * AHV_Model_X.Dp11_CSTATE_j[1];
  _rtXdot->Dp11_CSTATE_j[1] +=
      -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_j[2];
  _rtXdot->Dp11_CSTATE_j[1] +=
      -0.023103995048734959 * AHV_Model_X.Dp11_CSTATE_j[3];
  _rtXdot->Dp11_CSTATE_j[1] +=
      -0.0023400735275004216 * AHV_Model_X.Dp11_CSTATE_j[4];
  _rtXdot->Dp11_CSTATE_j[2] +=
      -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_j[0];
  _rtXdot->Dp11_CSTATE_j[2] +=
      0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_j[1];
  _rtXdot->Dp11_CSTATE_j[2] +=
      -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_j[2];
  _rtXdot->Dp11_CSTATE_j[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_j[3];
  _rtXdot->Dp11_CSTATE_j[2] +=
      -0.089746811088073059 * AHV_Model_X.Dp11_CSTATE_j[4];
  _rtXdot->Dp11_CSTATE_j[3] +=
      1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_j[0];
  _rtXdot->Dp11_CSTATE_j[3] +=
      -0.023103995048766656 * AHV_Model_X.Dp11_CSTATE_j[1];
  _rtXdot->Dp11_CSTATE_j[3] +=
      1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_j[2];
  _rtXdot->Dp11_CSTATE_j[3] +=
      -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_j[3];
  _rtXdot->Dp11_CSTATE_j[3] +=
      -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_j[4];
  _rtXdot->Dp11_CSTATE_j[4] +=
      -0.095423424399719625 * AHV_Model_X.Dp11_CSTATE_j[0];
  _rtXdot->Dp11_CSTATE_j[4] +=
      0.0023400735275035814 * AHV_Model_X.Dp11_CSTATE_j[1];
  _rtXdot->Dp11_CSTATE_j[4] +=
      -0.089746811088075445 * AHV_Model_X.Dp11_CSTATE_j[2];
  _rtXdot->Dp11_CSTATE_j[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_j[3];
  _rtXdot->Dp11_CSTATE_j[4] +=
      -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_j[4];
  _rtXdot->Dp11_CSTATE_j[0] += -3.3881001186476949 * AHV_Model_B.nu_r_d[0];
  _rtXdot->Dp11_CSTATE_j[1] += 0.07669611088588818 * AHV_Model_B.nu_r_d[0];
  _rtXdot->Dp11_CSTATE_j[2] += -0.46181267830731043 * AHV_Model_B.nu_r_d[0];
  _rtXdot->Dp11_CSTATE_j[3] += 1.226085627712131 * AHV_Model_B.nu_r_d[0];
  _rtXdot->Dp11_CSTATE_j[4] += -0.11754627904442222 * AHV_Model_B.nu_r_d[0];

  // Derivatives for StateSpace: '<S188>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_g[0] +=
      -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE_g[0];
  _rtXdot->Dp22_CSTATE_g[0] +=
      -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE_g[1];
  _rtXdot->Dp22_CSTATE_g[0] +=
      -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE_g[2];
  _rtXdot->Dp22_CSTATE_g[0] +=
      -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE_g[3];
  _rtXdot->Dp22_CSTATE_g[0] +=
      -0.095423424399736376 * AHV_Model_X.Dp22_CSTATE_g[4];
  _rtXdot->Dp22_CSTATE_g[1] +=
      1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_g[0];
  _rtXdot->Dp22_CSTATE_g[1] +=
      -0.00075311499999501316 * AHV_Model_X.Dp22_CSTATE_g[1];
  _rtXdot->Dp22_CSTATE_g[1] +=
      -0.010562919761942575 * AHV_Model_X.Dp22_CSTATE_g[2];
  _rtXdot->Dp22_CSTATE_g[1] +=
      -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE_g[3];
  _rtXdot->Dp22_CSTATE_g[1] +=
      -0.002340073527504722 * AHV_Model_X.Dp22_CSTATE_g[4];
  _rtXdot->Dp22_CSTATE_g[2] +=
      -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE_g[0];
  _rtXdot->Dp22_CSTATE_g[2] +=
      0.010562919761928094 * AHV_Model_X.Dp22_CSTATE_g[1];
  _rtXdot->Dp22_CSTATE_g[2] +=
      -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE_g[2];
  _rtXdot->Dp22_CSTATE_g[2] +=
      -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE_g[3];
  _rtXdot->Dp22_CSTATE_g[2] +=
      -0.089746811088091308 * AHV_Model_X.Dp22_CSTATE_g[4];
  _rtXdot->Dp22_CSTATE_g[3] +=
      1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_g[0];
  _rtXdot->Dp22_CSTATE_g[3] +=
      -0.023103995048724405 * AHV_Model_X.Dp22_CSTATE_g[1];
  _rtXdot->Dp22_CSTATE_g[3] +=
      1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_g[2];
  _rtXdot->Dp22_CSTATE_g[3] +=
      -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE_g[3];
  _rtXdot->Dp22_CSTATE_g[3] +=
      -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE_g[4];
  _rtXdot->Dp22_CSTATE_g[4] +=
      -0.095423424399735959 * AHV_Model_X.Dp22_CSTATE_g[0];
  _rtXdot->Dp22_CSTATE_g[4] +=
      0.00234007352750089 * AHV_Model_X.Dp22_CSTATE_g[1];
  _rtXdot->Dp22_CSTATE_g[4] +=
      -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE_g[2];
  _rtXdot->Dp22_CSTATE_g[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_g[3];
  _rtXdot->Dp22_CSTATE_g[4] +=
      -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE_g[4];
  _rtXdot->Dp22_CSTATE_g[0] += -3.3881001186476967 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp22_CSTATE_g[1] += 0.076696110885761712 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp22_CSTATE_g[2] += -0.46181267830733025 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp22_CSTATE_g[3] += 1.2260856277121315 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp22_CSTATE_g[4] += -0.11754627904444465 * AHV_Model_B.nu_r_d[1];

  // Derivatives for StateSpace: '<S188>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_f[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp24_CSTATE_f[0];
  _rtXdot->Dp24_CSTATE_f[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_f[1];
  _rtXdot->Dp24_CSTATE_f[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_f[2];
  _rtXdot->Dp24_CSTATE_f[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp24_CSTATE_f[3];
  _rtXdot->Dp24_CSTATE_f[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp24_CSTATE_f[4];
  _rtXdot->Dp24_CSTATE_f[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_f[0];
  _rtXdot->Dp24_CSTATE_f[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_f[1];
  _rtXdot->Dp24_CSTATE_f[1] +=
      0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_f[2];
  _rtXdot->Dp24_CSTATE_f[1] +=
      1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_f[3];
  _rtXdot->Dp24_CSTATE_f[1] +=
      0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_f[4];
  _rtXdot->Dp24_CSTATE_f[2] +=
      0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_f[0];
  _rtXdot->Dp24_CSTATE_f[2] +=
      0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_f[1];
  _rtXdot->Dp24_CSTATE_f[2] +=
      -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_f[2];
  _rtXdot->Dp24_CSTATE_f[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_f[3];
  _rtXdot->Dp24_CSTATE_f[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_f[4];
  _rtXdot->Dp24_CSTATE_f[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp24_CSTATE_f[0];
  _rtXdot->Dp24_CSTATE_f[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_f[1];
  _rtXdot->Dp24_CSTATE_f[3] +=
      3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_f[2];
  _rtXdot->Dp24_CSTATE_f[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_f[3];
  _rtXdot->Dp24_CSTATE_f[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_f[4];
  _rtXdot->Dp24_CSTATE_f[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp24_CSTATE_f[0];
  _rtXdot->Dp24_CSTATE_f[4] +=
      0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_f[1];
  _rtXdot->Dp24_CSTATE_f[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_f[2];
  _rtXdot->Dp24_CSTATE_f[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_f[3];
  _rtXdot->Dp24_CSTATE_f[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp24_CSTATE_f[4];
  _rtXdot->Dp24_CSTATE_f[0] += -0.23872786278308805 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp24_CSTATE_f[1] += -3.2763464234276056 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp24_CSTATE_f[2] += 1.1437118751635387 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp24_CSTATE_f[3] += -1.3244904674438165 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp24_CSTATE_f[4] += 0.071163531145327891 * AHV_Model_B.nu_r_d[3];

  // Derivatives for StateSpace: '<S188>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_i[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp26_CSTATE_i[0];
  _rtXdot->Dp26_CSTATE_i[0] +=
      1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_i[1];
  _rtXdot->Dp26_CSTATE_i[0] +=
      0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_i[2];
  _rtXdot->Dp26_CSTATE_i[0] +=
      0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_i[3];
  _rtXdot->Dp26_CSTATE_i[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp26_CSTATE_i[4];
  _rtXdot->Dp26_CSTATE_i[1] +=
      -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_i[0];
  _rtXdot->Dp26_CSTATE_i[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_i[1];
  _rtXdot->Dp26_CSTATE_i[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_i[2];
  _rtXdot->Dp26_CSTATE_i[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_i[3];
  _rtXdot->Dp26_CSTATE_i[1] +=
      0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_i[4];
  _rtXdot->Dp26_CSTATE_i[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp26_CSTATE_i[0];
  _rtXdot->Dp26_CSTATE_i[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_i[1];
  _rtXdot->Dp26_CSTATE_i[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_i[2];
  _rtXdot->Dp26_CSTATE_i[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_i[3];
  _rtXdot->Dp26_CSTATE_i[2] +=
      0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_i[4];
  _rtXdot->Dp26_CSTATE_i[3] +=
      0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_i[0];
  _rtXdot->Dp26_CSTATE_i[3] +=
      0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_i[1];
  _rtXdot->Dp26_CSTATE_i[3] +=
      1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_i[2];
  _rtXdot->Dp26_CSTATE_i[3] +=
      -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_i[3];
  _rtXdot->Dp26_CSTATE_i[3] +=
      1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_i[4];
  _rtXdot->Dp26_CSTATE_i[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp26_CSTATE_i[0];
  _rtXdot->Dp26_CSTATE_i[4] +=
      0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_i[1];
  _rtXdot->Dp26_CSTATE_i[4] +=
      0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_i[2];
  _rtXdot->Dp26_CSTATE_i[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_i[3];
  _rtXdot->Dp26_CSTATE_i[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_i[4];
  _rtXdot->Dp26_CSTATE_i[0] += 0.13167647316600378 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp26_CSTATE_i[1] += 3.4125284827245985 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp26_CSTATE_i[2] += 0.44399732492152644 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp26_CSTATE_i[3] += -1.2820687298457427 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp26_CSTATE_i[4] += -0.18276601298221137 * AHV_Model_B.nu_r_d[5];

  // Derivatives for StateSpace: '<S188>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_o[0] +=
      -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[0] +=
      -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[0] +=
      -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[0] +=
      0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[1] +=
      -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[1] +=
      -1.538297528844471E-6 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[1] +=
      0.00010972576839633958 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[1] +=
      0.00050341083172827721 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[1] +=
      -3.8193652830025795E-5 * AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[2] +=
      0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[2] +=
      0.0001097257683747581 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[2] +=
      -0.033603468126584782 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[2] +=
      -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[2] +=
      0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[3] +=
      -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[3] +=
      -0.00050341083159401254 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[3] +=
      1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[3] +=
      -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[3] +=
      0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[4] +=
      0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_o[0];
  _rtXdot->Dp33_CSTATE_o[4] +=
      3.8193652822306167E-5 * AHV_Model_X.Dp33_CSTATE_o[1];
  _rtXdot->Dp33_CSTATE_o[4] +=
      -0.039518459830266792 * AHV_Model_X.Dp33_CSTATE_o[2];
  _rtXdot->Dp33_CSTATE_o[4] +=
      0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_o[3];
  _rtXdot->Dp33_CSTATE_o[4] +=
      -0.0072896299361453719 * AHV_Model_X.Dp33_CSTATE_o[4];
  _rtXdot->Dp33_CSTATE_o[0] += -3.3182263091979736 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp33_CSTATE_o[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp33_CSTATE_o[2] += 0.13279109186599056 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp33_CSTATE_o[3] += -0.54052890870000092 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp33_CSTATE_o[4] += 0.042027661214546957 * AHV_Model_B.nu_r_d[2];

  // Derivatives for StateSpace: '<S188>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_i[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_i[0];
  _rtXdot->Dp35_CSTATE_i[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_i[1];
  _rtXdot->Dp35_CSTATE_i[0] +=
      2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_i[2];
  _rtXdot->Dp35_CSTATE_i[0] +=
      1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_i[3];
  _rtXdot->Dp35_CSTATE_i[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_i[4];
  _rtXdot->Dp35_CSTATE_i[1] +=
      0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_i[0];
  _rtXdot->Dp35_CSTATE_i[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp35_CSTATE_i[1];
  _rtXdot->Dp35_CSTATE_i[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp35_CSTATE_i[2];
  _rtXdot->Dp35_CSTATE_i[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp35_CSTATE_i[3];
  _rtXdot->Dp35_CSTATE_i[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp35_CSTATE_i[4];
  _rtXdot->Dp35_CSTATE_i[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_i[0];
  _rtXdot->Dp35_CSTATE_i[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp35_CSTATE_i[1];
  _rtXdot->Dp35_CSTATE_i[2] +=
      -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_i[2];
  _rtXdot->Dp35_CSTATE_i[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_i[3];
  _rtXdot->Dp35_CSTATE_i[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_i[4];
  _rtXdot->Dp35_CSTATE_i[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_i[0];
  _rtXdot->Dp35_CSTATE_i[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp35_CSTATE_i[1];
  _rtXdot->Dp35_CSTATE_i[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_i[2];
  _rtXdot->Dp35_CSTATE_i[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_i[3];
  _rtXdot->Dp35_CSTATE_i[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_i[4];
  _rtXdot->Dp35_CSTATE_i[4] +=
      0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_i[0];
  _rtXdot->Dp35_CSTATE_i[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_i[1];
  _rtXdot->Dp35_CSTATE_i[4] +=
      1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_i[2];
  _rtXdot->Dp35_CSTATE_i[4] +=
      3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_i[3];
  _rtXdot->Dp35_CSTATE_i[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_i[4];
  _rtXdot->Dp35_CSTATE_i[0] += -3.5374929885435322 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp35_CSTATE_i[1] += 0.015233386783081164 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp35_CSTATE_i[2] += -1.2196843916515756 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp35_CSTATE_i[3] += -0.96955326725759594 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp35_CSTATE_i[4] += 0.17148246208757431 * AHV_Model_B.nu_r_d[4];

  // Derivatives for StateSpace: '<S188>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_f[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp42_CSTATE_f[0];
  _rtXdot->Dp42_CSTATE_f[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_f[1];
  _rtXdot->Dp42_CSTATE_f[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_f[2];
  _rtXdot->Dp42_CSTATE_f[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp42_CSTATE_f[3];
  _rtXdot->Dp42_CSTATE_f[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp42_CSTATE_f[4];
  _rtXdot->Dp42_CSTATE_f[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_f[0];
  _rtXdot->Dp42_CSTATE_f[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_f[1];
  _rtXdot->Dp42_CSTATE_f[1] +=
      0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_f[2];
  _rtXdot->Dp42_CSTATE_f[1] +=
      1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_f[3];
  _rtXdot->Dp42_CSTATE_f[1] +=
      0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_f[4];
  _rtXdot->Dp42_CSTATE_f[2] +=
      0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_f[0];
  _rtXdot->Dp42_CSTATE_f[2] +=
      0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_f[1];
  _rtXdot->Dp42_CSTATE_f[2] +=
      -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_f[2];
  _rtXdot->Dp42_CSTATE_f[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_f[3];
  _rtXdot->Dp42_CSTATE_f[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_f[4];
  _rtXdot->Dp42_CSTATE_f[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp42_CSTATE_f[0];
  _rtXdot->Dp42_CSTATE_f[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_f[1];
  _rtXdot->Dp42_CSTATE_f[3] +=
      3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_f[2];
  _rtXdot->Dp42_CSTATE_f[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_f[3];
  _rtXdot->Dp42_CSTATE_f[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_f[4];
  _rtXdot->Dp42_CSTATE_f[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp42_CSTATE_f[0];
  _rtXdot->Dp42_CSTATE_f[4] +=
      0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_f[1];
  _rtXdot->Dp42_CSTATE_f[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_f[2];
  _rtXdot->Dp42_CSTATE_f[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_f[3];
  _rtXdot->Dp42_CSTATE_f[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp42_CSTATE_f[4];
  _rtXdot->Dp42_CSTATE_f[0] += -0.23872786278308805 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp42_CSTATE_f[1] += -3.2763464234276056 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp42_CSTATE_f[2] += 1.1437118751635387 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp42_CSTATE_f[3] += -1.3244904674438165 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp42_CSTATE_f[4] += 0.071163531145327891 * AHV_Model_B.nu_r_d[1];

  // Derivatives for StateSpace: '<S188>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_h[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[0] +=
      1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[0] +=
      0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[1] +=
      0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[1] +=
      0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[2] +=
      0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[2] +=
      0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[2] +=
      2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[3] +=
      0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[3] +=
      1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[3] +=
      0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp46_CSTATE_h[0];
  _rtXdot->Dp46_CSTATE_h[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp46_CSTATE_h[1];
  _rtXdot->Dp46_CSTATE_h[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp46_CSTATE_h[2];
  _rtXdot->Dp46_CSTATE_h[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_h[3];
  _rtXdot->Dp46_CSTATE_h[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp46_CSTATE_h[4];
  _rtXdot->Dp46_CSTATE_h[0] += -0.33524288945878539 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp46_CSTATE_h[1] += -6.3625492730866693 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp46_CSTATE_h[2] += 0.93578679415036736 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp46_CSTATE_h[3] += 2.5673752417819116 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp46_CSTATE_h[4] += 0.014933981938330865 * AHV_Model_B.nu_r_d[5];

  // Derivatives for StateSpace: '<S188>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_c[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_c[0];
  _rtXdot->Dp53_CSTATE_c[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_c[1];
  _rtXdot->Dp53_CSTATE_c[0] +=
      2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_c[2];
  _rtXdot->Dp53_CSTATE_c[0] +=
      1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_c[3];
  _rtXdot->Dp53_CSTATE_c[0] += 0.259710121707262 * AHV_Model_X.Dp53_CSTATE_c[4];
  _rtXdot->Dp53_CSTATE_c[1] +=
      0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_c[0];
  _rtXdot->Dp53_CSTATE_c[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp53_CSTATE_c[1];
  _rtXdot->Dp53_CSTATE_c[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp53_CSTATE_c[2];
  _rtXdot->Dp53_CSTATE_c[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp53_CSTATE_c[3];
  _rtXdot->Dp53_CSTATE_c[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp53_CSTATE_c[4];
  _rtXdot->Dp53_CSTATE_c[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_c[0];
  _rtXdot->Dp53_CSTATE_c[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp53_CSTATE_c[1];
  _rtXdot->Dp53_CSTATE_c[2] +=
      -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_c[2];
  _rtXdot->Dp53_CSTATE_c[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_c[3];
  _rtXdot->Dp53_CSTATE_c[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_c[4];
  _rtXdot->Dp53_CSTATE_c[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_c[0];
  _rtXdot->Dp53_CSTATE_c[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp53_CSTATE_c[1];
  _rtXdot->Dp53_CSTATE_c[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_c[2];
  _rtXdot->Dp53_CSTATE_c[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_c[3];
  _rtXdot->Dp53_CSTATE_c[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_c[4];
  _rtXdot->Dp53_CSTATE_c[4] +=
      0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_c[0];
  _rtXdot->Dp53_CSTATE_c[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_c[1];
  _rtXdot->Dp53_CSTATE_c[4] +=
      1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_c[2];
  _rtXdot->Dp53_CSTATE_c[4] +=
      3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_c[3];
  _rtXdot->Dp53_CSTATE_c[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_c[4];
  _rtXdot->Dp53_CSTATE_c[0] += -3.5374929885435322 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp53_CSTATE_c[1] += 0.015233386783081164 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp53_CSTATE_c[2] += -1.2196843916515756 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp53_CSTATE_c[3] += -0.96955326725759594 * AHV_Model_B.nu_r_d[2];
  _rtXdot->Dp53_CSTATE_c[4] += 0.17148246208757431 * AHV_Model_B.nu_r_d[2];

  // Derivatives for StateSpace: '<S188>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_f[0] +=
      -8.7466634083809418E-5 * AHV_Model_X.Dp55_CSTATE_f[0];
  _rtXdot->Dp55_CSTATE_f[0] +=
      0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_f[1];
  _rtXdot->Dp55_CSTATE_f[0] +=
      0.0013720743421809071 * AHV_Model_X.Dp55_CSTATE_f[2];
  _rtXdot->Dp55_CSTATE_f[0] +=
      0.0075482812368147384 * AHV_Model_X.Dp55_CSTATE_f[3];
  _rtXdot->Dp55_CSTATE_f[0] +=
      0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_f[4];
  _rtXdot->Dp55_CSTATE_f[1] +=
      -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_f[0];
  _rtXdot->Dp55_CSTATE_f[1] +=
      -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_f[1];
  _rtXdot->Dp55_CSTATE_f[1] +=
      -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_f[2];
  _rtXdot->Dp55_CSTATE_f[1] +=
      -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_f[3];
  _rtXdot->Dp55_CSTATE_f[1] +=
      -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_f[4];
  _rtXdot->Dp55_CSTATE_f[2] +=
      -0.0013720743421854747 * AHV_Model_X.Dp55_CSTATE_f[0];
  _rtXdot->Dp55_CSTATE_f[2] +=
      -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_f[1];
  _rtXdot->Dp55_CSTATE_f[2] +=
      -0.068186331750790974 * AHV_Model_X.Dp55_CSTATE_f[2];
  _rtXdot->Dp55_CSTATE_f[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_f[3];
  _rtXdot->Dp55_CSTATE_f[2] +=
      -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_f[4];
  _rtXdot->Dp55_CSTATE_f[3] +=
      0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_f[0];
  _rtXdot->Dp55_CSTATE_f[3] +=
      1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_f[1];
  _rtXdot->Dp55_CSTATE_f[3] +=
      4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_f[2];
  _rtXdot->Dp55_CSTATE_f[3] +=
      -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_f[3];
  _rtXdot->Dp55_CSTATE_f[3] +=
      -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_f[4];
  _rtXdot->Dp55_CSTATE_f[4] +=
      0.0070959954003254081 * AHV_Model_X.Dp55_CSTATE_f[0];
  _rtXdot->Dp55_CSTATE_f[4] +=
      1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_f[1];
  _rtXdot->Dp55_CSTATE_f[4] +=
      2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_f[2];
  _rtXdot->Dp55_CSTATE_f[4] +=
      -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_f[3];
  _rtXdot->Dp55_CSTATE_f[4] +=
      -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_f[4];
  _rtXdot->Dp55_CSTATE_f[0] += 0.021904981880023978 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp55_CSTATE_f[1] += 3.4651622640804733 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp55_CSTATE_f[2] += 0.16004490113712252 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp55_CSTATE_f[3] += -0.9969365784863482 * AHV_Model_B.nu_r_d[4];
  _rtXdot->Dp55_CSTATE_f[4] += -0.92774837285087253 * AHV_Model_B.nu_r_d[4];

  // Derivatives for StateSpace: '<S188>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_n[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp62_CSTATE_n[0];
  _rtXdot->Dp62_CSTATE_n[0] +=
      1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_n[1];
  _rtXdot->Dp62_CSTATE_n[0] +=
      0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_n[2];
  _rtXdot->Dp62_CSTATE_n[0] +=
      0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_n[3];
  _rtXdot->Dp62_CSTATE_n[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp62_CSTATE_n[4];
  _rtXdot->Dp62_CSTATE_n[1] +=
      -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_n[0];
  _rtXdot->Dp62_CSTATE_n[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_n[1];
  _rtXdot->Dp62_CSTATE_n[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_n[2];
  _rtXdot->Dp62_CSTATE_n[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_n[3];
  _rtXdot->Dp62_CSTATE_n[1] +=
      0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_n[4];
  _rtXdot->Dp62_CSTATE_n[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp62_CSTATE_n[0];
  _rtXdot->Dp62_CSTATE_n[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_n[1];
  _rtXdot->Dp62_CSTATE_n[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_n[2];
  _rtXdot->Dp62_CSTATE_n[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_n[3];
  _rtXdot->Dp62_CSTATE_n[2] +=
      0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_n[4];
  _rtXdot->Dp62_CSTATE_n[3] +=
      0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_n[0];
  _rtXdot->Dp62_CSTATE_n[3] +=
      0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_n[1];
  _rtXdot->Dp62_CSTATE_n[3] +=
      1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_n[2];
  _rtXdot->Dp62_CSTATE_n[3] +=
      -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_n[3];
  _rtXdot->Dp62_CSTATE_n[3] +=
      1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_n[4];
  _rtXdot->Dp62_CSTATE_n[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp62_CSTATE_n[0];
  _rtXdot->Dp62_CSTATE_n[4] +=
      0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_n[1];
  _rtXdot->Dp62_CSTATE_n[4] +=
      0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_n[2];
  _rtXdot->Dp62_CSTATE_n[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_n[3];
  _rtXdot->Dp62_CSTATE_n[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_n[4];
  _rtXdot->Dp62_CSTATE_n[0] += 0.13167647316600378 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp62_CSTATE_n[1] += 3.4125284827245985 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp62_CSTATE_n[2] += 0.44399732492152644 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp62_CSTATE_n[3] += -1.2820687298457427 * AHV_Model_B.nu_r_d[1];
  _rtXdot->Dp62_CSTATE_n[4] += -0.18276601298221137 * AHV_Model_B.nu_r_d[1];

  // Derivatives for StateSpace: '<S188>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_b[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp64_CSTATE_b[0];
  _rtXdot->Dp64_CSTATE_b[0] +=
      1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_b[1];
  _rtXdot->Dp64_CSTATE_b[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp64_CSTATE_b[2];
  _rtXdot->Dp64_CSTATE_b[0] +=
      0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_b[3];
  _rtXdot->Dp64_CSTATE_b[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp64_CSTATE_b[4];
  _rtXdot->Dp64_CSTATE_b[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_b[0];
  _rtXdot->Dp64_CSTATE_b[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_b[1];
  _rtXdot->Dp64_CSTATE_b[1] +=
      0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_b[2];
  _rtXdot->Dp64_CSTATE_b[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_b[3];
  _rtXdot->Dp64_CSTATE_b[1] +=
      0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_b[4];
  _rtXdot->Dp64_CSTATE_b[2] +=
      0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_b[0];
  _rtXdot->Dp64_CSTATE_b[2] +=
      0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_b[1];
  _rtXdot->Dp64_CSTATE_b[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_b[2];
  _rtXdot->Dp64_CSTATE_b[2] +=
      2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_b[3];
  _rtXdot->Dp64_CSTATE_b[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp64_CSTATE_b[4];
  _rtXdot->Dp64_CSTATE_b[3] +=
      0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_b[0];
  _rtXdot->Dp64_CSTATE_b[3] +=
      1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_b[1];
  _rtXdot->Dp64_CSTATE_b[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_b[2];
  _rtXdot->Dp64_CSTATE_b[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_b[3];
  _rtXdot->Dp64_CSTATE_b[3] +=
      0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_b[4];
  _rtXdot->Dp64_CSTATE_b[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp64_CSTATE_b[0];
  _rtXdot->Dp64_CSTATE_b[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp64_CSTATE_b[1];
  _rtXdot->Dp64_CSTATE_b[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp64_CSTATE_b[2];
  _rtXdot->Dp64_CSTATE_b[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_b[3];
  _rtXdot->Dp64_CSTATE_b[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp64_CSTATE_b[4];
  _rtXdot->Dp64_CSTATE_b[0] += -0.33524288945878539 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp64_CSTATE_b[1] += -6.3625492730866693 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp64_CSTATE_b[2] += 0.93578679415036736 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp64_CSTATE_b[3] += 2.5673752417819116 * AHV_Model_B.nu_r_d[3];
  _rtXdot->Dp64_CSTATE_b[4] += 0.014933981938330865 * AHV_Model_B.nu_r_d[3];

  // Derivatives for StateSpace: '<S188>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_n[0] +=
      -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[0] +=
      -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[0] +=
      -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[0] +=
      -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[0] +=
      0.069477290559771143 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[1] +=
      1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[1] +=
      -0.00029199832472599611 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[1] +=
      -0.0080078570280257711 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[1] +=
      -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[1] +=
      0.00098472174210678556 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[2] +=
      -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[2] +=
      0.00800785702802341 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[2] +=
      -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[2] +=
      -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[2] +=
      0.074587161707531782 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[3] +=
      1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[3] +=
      -0.015684008067190565 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[3] +=
      1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[3] +=
      -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[3] += 1.217430977302715 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[4] +=
      0.069477290559780386 * AHV_Model_X.Dp66_CSTATE_n[0];
  _rtXdot->Dp66_CSTATE_n[4] +=
      -0.000984721742106249 * AHV_Model_X.Dp66_CSTATE_n[1];
  _rtXdot->Dp66_CSTATE_n[4] +=
      0.074587161707547089 * AHV_Model_X.Dp66_CSTATE_n[2];
  _rtXdot->Dp66_CSTATE_n[4] +=
      -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE_n[3];
  _rtXdot->Dp66_CSTATE_n[4] +=
      -0.061420217921035136 * AHV_Model_X.Dp66_CSTATE_n[4];
  _rtXdot->Dp66_CSTATE_n[0] += -3.3844821061206916 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp66_CSTATE_n[1] += 0.045246867266885989 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp66_CSTATE_n[2] += -0.53332191483729374 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp66_CSTATE_n[3] += 1.2581308985447515 * AHV_Model_B.nu_r_d[5];
  _rtXdot->Dp66_CSTATE_n[4] += 0.075289942321940723 * AHV_Model_B.nu_r_d[5];

  // Derivatives for Integrator: '<S231>/Integrator1'
  _rtXdot->Integrator1_CSTATE_m[0] = AHV_Model_B.Sum6_b[0];

  // Derivatives for Integrator: '<S231>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[0] = AHV_Model_B.psi_WF_c[0];

  // Derivatives for Integrator: '<S231>/Integrator6'
  _rtXdot->Integrator6_CSTATE_e[0] = AHV_Model_B.Sum7_e[0];

  // Derivatives for Integrator: '<S231>/Integrator1'
  _rtXdot->Integrator1_CSTATE_m[1] = AHV_Model_B.Sum6_b[1];

  // Derivatives for Integrator: '<S231>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[1] = AHV_Model_B.psi_WF_c[1];

  // Derivatives for Integrator: '<S231>/Integrator6'
  _rtXdot->Integrator6_CSTATE_e[1] = AHV_Model_B.Sum7_e[1];

  // Derivatives for Integrator: '<S231>/Integrator1'
  _rtXdot->Integrator1_CSTATE_m[2] = AHV_Model_B.Sum6_b[2];

  // Derivatives for Integrator: '<S231>/Integrator2'
  _rtXdot->Integrator2_CSTATE_j[2] = AHV_Model_B.psi_WF_c[2];

  // Derivatives for Integrator: '<S231>/Integrator6'
  _rtXdot->Integrator6_CSTATE_e[2] = AHV_Model_B.Sum7_e[2];
  for (i = 0; i < 6; i++) {
    // Derivatives for Integrator: '<S265>/Integrator1'
    _rtXdot->Integrator1_CSTATE_b[i] =
        AHV_Model_B.TmpSignalConversionAtInteg_o3s0[i];

    // Derivatives for Integrator: '<S265>/Integrator'
    _rtXdot->Integrator_CSTATE_a[i] = AHV_Model_B.Minvtau_l[i];
  }

  // Derivatives for TransferFcn: '<S268>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE_m = 0.0;
  _rtXdot->TransferFcn_CSTATE_m += -0.2 * AHV_Model_X.TransferFcn_CSTATE_m;
  _rtXdot->TransferFcn_CSTATE_m += AHV_Model_B.Fcn_h;

  // Derivatives for Integrator: '<S315>/Integrator3'
  _rtXdot->Integrator3_CSTATE_d[0] = AHV_Model_B.sun_k2_k[0];

  // Derivatives for Integrator: '<S315>/Integrator4'
  _rtXdot->Integrator4_CSTATE_o[0] = AHV_Model_B.M_u_j[0];

  // Derivatives for Integrator: '<S315>/Integrator3'
  _rtXdot->Integrator3_CSTATE_d[1] = AHV_Model_B.sun_k2_k[1];

  // Derivatives for Integrator: '<S315>/Integrator4'
  _rtXdot->Integrator4_CSTATE_o[1] = AHV_Model_B.M_u_j[1];

  // Derivatives for Integrator: '<S315>/Integrator3'
  _rtXdot->Integrator3_CSTATE_d[2] = AHV_Model_B.sun_k2_k[2];

  // Derivatives for Integrator: '<S315>/Integrator4'
  _rtXdot->Integrator4_CSTATE_o[2] = AHV_Model_B.M_u_j[2];
  for (i = 0; i < 5; i++) {
    // Derivatives for StateSpace: '<S272>/Dp(1,1)'
    _rtXdot->Dp11_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(1,3)'
    _rtXdot->Dp13_CSTATE_c[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(1,5)'
    _rtXdot->Dp15_CSTATE_d[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(2,2)'
    _rtXdot->Dp22_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(2,4)'
    _rtXdot->Dp24_CSTATE_d[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(2,6)'
    _rtXdot->Dp26_CSTATE_f[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(3,1)'
    _rtXdot->Dp31_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(3,3)'
    _rtXdot->Dp33_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(3,5)'
    _rtXdot->Dp35_CSTATE_l[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(4,2)'
    _rtXdot->Dp42_CSTATE_h[i] = 0.0;

    // Derivatives for Integrator: '<S284>/Integrator'
    _rtXdot->Integrator_CSTATE_bp[i] = AHV_Model_B.Sum_c[i];

    // Derivatives for StateSpace: '<S272>/Dp(4,6)'
    _rtXdot->Dp46_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(5,1)'
    _rtXdot->Dp51_CSTATE_k[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(5,3)'
    _rtXdot->Dp53_CSTATE_pz[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(5,5)'
    _rtXdot->Dp55_CSTATE_a[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(6,2)'
    _rtXdot->Dp62_CSTATE_i[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(6,4)'
    _rtXdot->Dp64_CSTATE_e[i] = 0.0;

    // Derivatives for StateSpace: '<S272>/Dp(6,6)'
    _rtXdot->Dp66_CSTATE_g3[i] = 0.0;
  }

  // Derivatives for StateSpace: '<S272>/Dp(1,1)'
  _rtXdot->Dp11_CSTATE_e[0] += -1.39281413011323 * AHV_Model_X.Dp11_CSTATE_e[0];
  _rtXdot->Dp11_CSTATE_e[0] +=
      -1.2054566685310337 * AHV_Model_X.Dp11_CSTATE_e[1];
  _rtXdot->Dp11_CSTATE_e[0] +=
      -0.33483039198972808 * AHV_Model_X.Dp11_CSTATE_e[2];
  _rtXdot->Dp11_CSTATE_e[0] +=
      -1.0500590385158857 * AHV_Model_X.Dp11_CSTATE_e[3];
  _rtXdot->Dp11_CSTATE_e[0] +=
      -0.0954234243997144 * AHV_Model_X.Dp11_CSTATE_e[4];
  _rtXdot->Dp11_CSTATE_e[1] +=
      1.2054566685310679 * AHV_Model_X.Dp11_CSTATE_e[0];
  _rtXdot->Dp11_CSTATE_e[1] +=
      -0.00075311499999405624 * AHV_Model_X.Dp11_CSTATE_e[1];
  _rtXdot->Dp11_CSTATE_e[1] +=
      -0.01056291976193271 * AHV_Model_X.Dp11_CSTATE_e[2];
  _rtXdot->Dp11_CSTATE_e[1] +=
      -0.023103995048734959 * AHV_Model_X.Dp11_CSTATE_e[3];
  _rtXdot->Dp11_CSTATE_e[1] +=
      -0.0023400735275004216 * AHV_Model_X.Dp11_CSTATE_e[4];
  _rtXdot->Dp11_CSTATE_e[2] +=
      -0.33483039198971015 * AHV_Model_X.Dp11_CSTATE_e[0];
  _rtXdot->Dp11_CSTATE_e[2] +=
      0.010562919761939651 * AHV_Model_X.Dp11_CSTATE_e[1];
  _rtXdot->Dp11_CSTATE_e[2] +=
      -0.19313064085990217 * AHV_Model_X.Dp11_CSTATE_e[2];
  _rtXdot->Dp11_CSTATE_e[2] += -1.46185154439905 * AHV_Model_X.Dp11_CSTATE_e[3];
  _rtXdot->Dp11_CSTATE_e[2] +=
      -0.089746811088073059 * AHV_Model_X.Dp11_CSTATE_e[4];
  _rtXdot->Dp11_CSTATE_e[3] +=
      1.0500590385158888 * AHV_Model_X.Dp11_CSTATE_e[0];
  _rtXdot->Dp11_CSTATE_e[3] +=
      -0.023103995048766656 * AHV_Model_X.Dp11_CSTATE_e[1];
  _rtXdot->Dp11_CSTATE_e[3] +=
      1.4618515443990421 * AHV_Model_X.Dp11_CSTATE_e[2];
  _rtXdot->Dp11_CSTATE_e[3] +=
      -4.560672872596915 * AHV_Model_X.Dp11_CSTATE_e[3];
  _rtXdot->Dp11_CSTATE_e[3] +=
      -1.2857667586957913 * AHV_Model_X.Dp11_CSTATE_e[4];
  _rtXdot->Dp11_CSTATE_e[4] +=
      -0.095423424399719625 * AHV_Model_X.Dp11_CSTATE_e[0];
  _rtXdot->Dp11_CSTATE_e[4] +=
      0.0023400735275035814 * AHV_Model_X.Dp11_CSTATE_e[1];
  _rtXdot->Dp11_CSTATE_e[4] +=
      -0.089746811088075445 * AHV_Model_X.Dp11_CSTATE_e[2];
  _rtXdot->Dp11_CSTATE_e[4] += 1.285766758695809 * AHV_Model_X.Dp11_CSTATE_e[3];
  _rtXdot->Dp11_CSTATE_e[4] +=
      -0.13104378720661222 * AHV_Model_X.Dp11_CSTATE_e[4];
  _rtXdot->Dp11_CSTATE_e[0] += -3.3881001186476949 * AHV_Model_B.nu_r_f[0];
  _rtXdot->Dp11_CSTATE_e[1] += 0.07669611088588818 * AHV_Model_B.nu_r_f[0];
  _rtXdot->Dp11_CSTATE_e[2] += -0.46181267830731043 * AHV_Model_B.nu_r_f[0];
  _rtXdot->Dp11_CSTATE_e[3] += 1.226085627712131 * AHV_Model_B.nu_r_f[0];
  _rtXdot->Dp11_CSTATE_e[4] += -0.11754627904442222 * AHV_Model_B.nu_r_f[0];

  // Derivatives for StateSpace: '<S272>/Dp(2,2)'
  _rtXdot->Dp22_CSTATE_a[0] +=
      -1.3928141301132206 * AHV_Model_X.Dp22_CSTATE_a[0];
  _rtXdot->Dp22_CSTATE_a[0] +=
      -1.2054566685310828 * AHV_Model_X.Dp22_CSTATE_a[1];
  _rtXdot->Dp22_CSTATE_a[0] +=
      -0.33483039198973219 * AHV_Model_X.Dp22_CSTATE_a[2];
  _rtXdot->Dp22_CSTATE_a[0] +=
      -1.0500590385158863 * AHV_Model_X.Dp22_CSTATE_a[3];
  _rtXdot->Dp22_CSTATE_a[0] +=
      -0.095423424399736376 * AHV_Model_X.Dp22_CSTATE_a[4];
  _rtXdot->Dp22_CSTATE_a[1] +=
      1.2054566685310135 * AHV_Model_X.Dp22_CSTATE_a[0];
  _rtXdot->Dp22_CSTATE_a[1] +=
      -0.00075311499999501316 * AHV_Model_X.Dp22_CSTATE_a[1];
  _rtXdot->Dp22_CSTATE_a[1] +=
      -0.010562919761942575 * AHV_Model_X.Dp22_CSTATE_a[2];
  _rtXdot->Dp22_CSTATE_a[1] +=
      -0.0231039950487718 * AHV_Model_X.Dp22_CSTATE_a[3];
  _rtXdot->Dp22_CSTATE_a[1] +=
      -0.002340073527504722 * AHV_Model_X.Dp22_CSTATE_a[4];
  _rtXdot->Dp22_CSTATE_a[2] +=
      -0.33483039198972253 * AHV_Model_X.Dp22_CSTATE_a[0];
  _rtXdot->Dp22_CSTATE_a[2] +=
      0.010562919761928094 * AHV_Model_X.Dp22_CSTATE_a[1];
  _rtXdot->Dp22_CSTATE_a[2] +=
      -0.19313064085990478 * AHV_Model_X.Dp22_CSTATE_a[2];
  _rtXdot->Dp22_CSTATE_a[2] +=
      -1.4618515443989919 * AHV_Model_X.Dp22_CSTATE_a[3];
  _rtXdot->Dp22_CSTATE_a[2] +=
      -0.089746811088091308 * AHV_Model_X.Dp22_CSTATE_a[4];
  _rtXdot->Dp22_CSTATE_a[3] +=
      1.0500590385158797 * AHV_Model_X.Dp22_CSTATE_a[0];
  _rtXdot->Dp22_CSTATE_a[3] +=
      -0.023103995048724405 * AHV_Model_X.Dp22_CSTATE_a[1];
  _rtXdot->Dp22_CSTATE_a[3] +=
      1.4618515443989983 * AHV_Model_X.Dp22_CSTATE_a[2];
  _rtXdot->Dp22_CSTATE_a[3] +=
      -4.5606728725966441 * AHV_Model_X.Dp22_CSTATE_a[3];
  _rtXdot->Dp22_CSTATE_a[3] +=
      -1.2857667586959629 * AHV_Model_X.Dp22_CSTATE_a[4];
  _rtXdot->Dp22_CSTATE_a[4] +=
      -0.095423424399735959 * AHV_Model_X.Dp22_CSTATE_a[0];
  _rtXdot->Dp22_CSTATE_a[4] +=
      0.00234007352750089 * AHV_Model_X.Dp22_CSTATE_a[1];
  _rtXdot->Dp22_CSTATE_a[4] +=
      -0.0897468110880912 * AHV_Model_X.Dp22_CSTATE_a[2];
  _rtXdot->Dp22_CSTATE_a[4] += 1.285766758695962 * AHV_Model_X.Dp22_CSTATE_a[3];
  _rtXdot->Dp22_CSTATE_a[4] +=
      -0.13104378720667384 * AHV_Model_X.Dp22_CSTATE_a[4];
  _rtXdot->Dp22_CSTATE_a[0] += -3.3881001186476967 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp22_CSTATE_a[1] += 0.076696110885761712 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp22_CSTATE_a[2] += -0.46181267830733025 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp22_CSTATE_a[3] += 1.2260856277121315 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp22_CSTATE_a[4] += -0.11754627904444465 * AHV_Model_B.nu_r_f[1];

  // Derivatives for StateSpace: '<S272>/Dp(2,4)'
  _rtXdot->Dp24_CSTATE_d[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp24_CSTATE_d[0];
  _rtXdot->Dp24_CSTATE_d[0] += 1.046159104970257 * AHV_Model_X.Dp24_CSTATE_d[1];
  _rtXdot->Dp24_CSTATE_d[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp24_CSTATE_d[2];
  _rtXdot->Dp24_CSTATE_d[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp24_CSTATE_d[3];
  _rtXdot->Dp24_CSTATE_d[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp24_CSTATE_d[4];
  _rtXdot->Dp24_CSTATE_d[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp24_CSTATE_d[0];
  _rtXdot->Dp24_CSTATE_d[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp24_CSTATE_d[1];
  _rtXdot->Dp24_CSTATE_d[1] +=
      0.95129605479096324 * AHV_Model_X.Dp24_CSTATE_d[2];
  _rtXdot->Dp24_CSTATE_d[1] +=
      1.3710797263241543 * AHV_Model_X.Dp24_CSTATE_d[3];
  _rtXdot->Dp24_CSTATE_d[1] +=
      0.067632454414034648 * AHV_Model_X.Dp24_CSTATE_d[4];
  _rtXdot->Dp24_CSTATE_d[2] +=
      0.078929006579498431 * AHV_Model_X.Dp24_CSTATE_d[0];
  _rtXdot->Dp24_CSTATE_d[2] +=
      0.9512960547909608 * AHV_Model_X.Dp24_CSTATE_d[1];
  _rtXdot->Dp24_CSTATE_d[2] +=
      -1.065705011247746 * AHV_Model_X.Dp24_CSTATE_d[2];
  _rtXdot->Dp24_CSTATE_d[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp24_CSTATE_d[3];
  _rtXdot->Dp24_CSTATE_d[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp24_CSTATE_d[4];
  _rtXdot->Dp24_CSTATE_d[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp24_CSTATE_d[0];
  _rtXdot->Dp24_CSTATE_d[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp24_CSTATE_d[1];
  _rtXdot->Dp24_CSTATE_d[3] +=
      3.3410455288652465 * AHV_Model_X.Dp24_CSTATE_d[2];
  _rtXdot->Dp24_CSTATE_d[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp24_CSTATE_d[3];
  _rtXdot->Dp24_CSTATE_d[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp24_CSTATE_d[4];
  _rtXdot->Dp24_CSTATE_d[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp24_CSTATE_d[0];
  _rtXdot->Dp24_CSTATE_d[4] +=
      0.067632454414049734 * AHV_Model_X.Dp24_CSTATE_d[1];
  _rtXdot->Dp24_CSTATE_d[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp24_CSTATE_d[2];
  _rtXdot->Dp24_CSTATE_d[4] += 2.46497845458085 * AHV_Model_X.Dp24_CSTATE_d[3];
  _rtXdot->Dp24_CSTATE_d[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp24_CSTATE_d[4];
  _rtXdot->Dp24_CSTATE_d[0] += -0.23872786278308805 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp24_CSTATE_d[1] += -3.2763464234276056 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp24_CSTATE_d[2] += 1.1437118751635387 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp24_CSTATE_d[3] += -1.3244904674438165 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp24_CSTATE_d[4] += 0.071163531145327891 * AHV_Model_B.nu_r_f[3];

  // Derivatives for StateSpace: '<S272>/Dp(2,6)'
  _rtXdot->Dp26_CSTATE_f[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp26_CSTATE_f[0];
  _rtXdot->Dp26_CSTATE_f[0] +=
      1.2361995937266759 * AHV_Model_X.Dp26_CSTATE_f[1];
  _rtXdot->Dp26_CSTATE_f[0] +=
      0.012808448456780379 * AHV_Model_X.Dp26_CSTATE_f[2];
  _rtXdot->Dp26_CSTATE_f[0] +=
      0.031758809083141569 * AHV_Model_X.Dp26_CSTATE_f[3];
  _rtXdot->Dp26_CSTATE_f[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp26_CSTATE_f[4];
  _rtXdot->Dp26_CSTATE_f[1] +=
      -1.236199593726705 * AHV_Model_X.Dp26_CSTATE_f[0];
  _rtXdot->Dp26_CSTATE_f[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp26_CSTATE_f[1];
  _rtXdot->Dp26_CSTATE_f[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp26_CSTATE_f[2];
  _rtXdot->Dp26_CSTATE_f[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp26_CSTATE_f[3];
  _rtXdot->Dp26_CSTATE_f[1] +=
      0.12995794495601007 * AHV_Model_X.Dp26_CSTATE_f[4];
  _rtXdot->Dp26_CSTATE_f[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp26_CSTATE_f[0];
  _rtXdot->Dp26_CSTATE_f[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp26_CSTATE_f[1];
  _rtXdot->Dp26_CSTATE_f[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp26_CSTATE_f[2];
  _rtXdot->Dp26_CSTATE_f[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp26_CSTATE_f[3];
  _rtXdot->Dp26_CSTATE_f[2] +=
      0.1356516649888091 * AHV_Model_X.Dp26_CSTATE_f[4];
  _rtXdot->Dp26_CSTATE_f[3] +=
      0.031758809083164766 * AHV_Model_X.Dp26_CSTATE_f[0];
  _rtXdot->Dp26_CSTATE_f[3] +=
      0.96895750705376427 * AHV_Model_X.Dp26_CSTATE_f[1];
  _rtXdot->Dp26_CSTATE_f[3] +=
      1.8113135858263727 * AHV_Model_X.Dp26_CSTATE_f[2];
  _rtXdot->Dp26_CSTATE_f[3] +=
      -3.762337151306161 * AHV_Model_X.Dp26_CSTATE_f[3];
  _rtXdot->Dp26_CSTATE_f[3] +=
      1.5320613622870083 * AHV_Model_X.Dp26_CSTATE_f[4];
  _rtXdot->Dp26_CSTATE_f[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp26_CSTATE_f[0];
  _rtXdot->Dp26_CSTATE_f[4] +=
      0.12995794495600341 * AHV_Model_X.Dp26_CSTATE_f[1];
  _rtXdot->Dp26_CSTATE_f[4] +=
      0.13565166498879488 * AHV_Model_X.Dp26_CSTATE_f[2];
  _rtXdot->Dp26_CSTATE_f[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp26_CSTATE_f[3];
  _rtXdot->Dp26_CSTATE_f[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp26_CSTATE_f[4];
  _rtXdot->Dp26_CSTATE_f[0] += 0.13167647316600378 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp26_CSTATE_f[1] += 3.4125284827245985 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp26_CSTATE_f[2] += 0.44399732492152644 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp26_CSTATE_f[3] += -1.2820687298457427 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp26_CSTATE_f[4] += -0.18276601298221137 * AHV_Model_B.nu_r_f[5];

  // Derivatives for StateSpace: '<S272>/Dp(3,3)'
  _rtXdot->Dp33_CSTATE_a[0] +=
      -1.3813850547222439 * AHV_Model_X.Dp33_CSTATE_a[0];
  _rtXdot->Dp33_CSTATE_a[0] += 0.539512803048723 * AHV_Model_X.Dp33_CSTATE_a[1];
  _rtXdot->Dp33_CSTATE_a[0] +=
      -0.11835430912143184 * AHV_Model_X.Dp33_CSTATE_a[2];
  _rtXdot->Dp33_CSTATE_a[0] +=
      -0.42708992140225116 * AHV_Model_X.Dp33_CSTATE_a[3];
  _rtXdot->Dp33_CSTATE_a[0] +=
      0.033960047754095488 * AHV_Model_X.Dp33_CSTATE_a[4];
  _rtXdot->Dp33_CSTATE_a[1] +=
      -0.53951280304832772 * AHV_Model_X.Dp33_CSTATE_a[0];
  _rtXdot->Dp33_CSTATE_a[1] +=
      -1.538297528844471E-6 * AHV_Model_X.Dp33_CSTATE_a[1];
  _rtXdot->Dp33_CSTATE_a[1] +=
      0.00010972576839633958 * AHV_Model_X.Dp33_CSTATE_a[2];
  _rtXdot->Dp33_CSTATE_a[1] +=
      0.00050341083172827721 * AHV_Model_X.Dp33_CSTATE_a[3];
  _rtXdot->Dp33_CSTATE_a[1] +=
      -3.8193652830025795E-5 * AHV_Model_X.Dp33_CSTATE_a[4];
  _rtXdot->Dp33_CSTATE_a[2] +=
      0.11835430912143848 * AHV_Model_X.Dp33_CSTATE_a[0];
  _rtXdot->Dp33_CSTATE_a[2] +=
      0.0001097257683747581 * AHV_Model_X.Dp33_CSTATE_a[1];
  _rtXdot->Dp33_CSTATE_a[2] +=
      -0.033603468126584782 * AHV_Model_X.Dp33_CSTATE_a[2];
  _rtXdot->Dp33_CSTATE_a[2] +=
      -1.4905346203698435 * AHV_Model_X.Dp33_CSTATE_a[3];
  _rtXdot->Dp33_CSTATE_a[2] +=
      0.039518459830251443 * AHV_Model_X.Dp33_CSTATE_a[4];
  _rtXdot->Dp33_CSTATE_a[3] +=
      -0.42708992140225777 * AHV_Model_X.Dp33_CSTATE_a[0];
  _rtXdot->Dp33_CSTATE_a[3] +=
      -0.00050341083159401254 * AHV_Model_X.Dp33_CSTATE_a[1];
  _rtXdot->Dp33_CSTATE_a[3] +=
      1.4905346203698526 * AHV_Model_X.Dp33_CSTATE_a[2];
  _rtXdot->Dp33_CSTATE_a[3] +=
      -0.68194162307127393 * AHV_Model_X.Dp33_CSTATE_a[3];
  _rtXdot->Dp33_CSTATE_a[3] +=
      0.067736905936386635 * AHV_Model_X.Dp33_CSTATE_a[4];
  _rtXdot->Dp33_CSTATE_a[4] +=
      0.033960047754138974 * AHV_Model_X.Dp33_CSTATE_a[0];
  _rtXdot->Dp33_CSTATE_a[4] +=
      3.8193652822306167E-5 * AHV_Model_X.Dp33_CSTATE_a[1];
  _rtXdot->Dp33_CSTATE_a[4] +=
      -0.039518459830266792 * AHV_Model_X.Dp33_CSTATE_a[2];
  _rtXdot->Dp33_CSTATE_a[4] +=
      0.067736905936572681 * AHV_Model_X.Dp33_CSTATE_a[3];
  _rtXdot->Dp33_CSTATE_a[4] +=
      -0.0072896299361453719 * AHV_Model_X.Dp33_CSTATE_a[4];
  _rtXdot->Dp33_CSTATE_a[0] += -3.3182263091979736 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp33_CSTATE_a[1] += -0.0034921698713633515 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp33_CSTATE_a[2] += 0.13279109186599056 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp33_CSTATE_a[3] += -0.54052890870000092 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp33_CSTATE_a[4] += 0.042027661214546957 * AHV_Model_B.nu_r_f[2];

  // Derivatives for StateSpace: '<S272>/Dp(3,5)'
  _rtXdot->Dp35_CSTATE_l[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp35_CSTATE_l[0];
  _rtXdot->Dp35_CSTATE_l[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp35_CSTATE_l[1];
  _rtXdot->Dp35_CSTATE_l[0] +=
      2.1365278323686749 * AHV_Model_X.Dp35_CSTATE_l[2];
  _rtXdot->Dp35_CSTATE_l[0] +=
      1.6362305234027079 * AHV_Model_X.Dp35_CSTATE_l[3];
  _rtXdot->Dp35_CSTATE_l[0] += 0.259710121707262 * AHV_Model_X.Dp35_CSTATE_l[4];
  _rtXdot->Dp35_CSTATE_l[1] +=
      0.9998076276438177 * AHV_Model_X.Dp35_CSTATE_l[0];
  _rtXdot->Dp35_CSTATE_l[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp35_CSTATE_l[1];
  _rtXdot->Dp35_CSTATE_l[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp35_CSTATE_l[2];
  _rtXdot->Dp35_CSTATE_l[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp35_CSTATE_l[3];
  _rtXdot->Dp35_CSTATE_l[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp35_CSTATE_l[4];
  _rtXdot->Dp35_CSTATE_l[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp35_CSTATE_l[0];
  _rtXdot->Dp35_CSTATE_l[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp35_CSTATE_l[1];
  _rtXdot->Dp35_CSTATE_l[2] +=
      -3.350088004888828 * AHV_Model_X.Dp35_CSTATE_l[2];
  _rtXdot->Dp35_CSTATE_l[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp35_CSTATE_l[3];
  _rtXdot->Dp35_CSTATE_l[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp35_CSTATE_l[4];
  _rtXdot->Dp35_CSTATE_l[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp35_CSTATE_l[0];
  _rtXdot->Dp35_CSTATE_l[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp35_CSTATE_l[1];
  _rtXdot->Dp35_CSTATE_l[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp35_CSTATE_l[2];
  _rtXdot->Dp35_CSTATE_l[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp35_CSTATE_l[3];
  _rtXdot->Dp35_CSTATE_l[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp35_CSTATE_l[4];
  _rtXdot->Dp35_CSTATE_l[4] +=
      0.2597101217072279 * AHV_Model_X.Dp35_CSTATE_l[0];
  _rtXdot->Dp35_CSTATE_l[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp35_CSTATE_l[1];
  _rtXdot->Dp35_CSTATE_l[4] +=
      1.6370488166464041 * AHV_Model_X.Dp35_CSTATE_l[2];
  _rtXdot->Dp35_CSTATE_l[4] +=
      3.2557107662584994 * AHV_Model_X.Dp35_CSTATE_l[3];
  _rtXdot->Dp35_CSTATE_l[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp35_CSTATE_l[4];
  _rtXdot->Dp35_CSTATE_l[0] += -3.5374929885435322 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp35_CSTATE_l[1] += 0.015233386783081164 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp35_CSTATE_l[2] += -1.2196843916515756 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp35_CSTATE_l[3] += -0.96955326725759594 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp35_CSTATE_l[4] += 0.17148246208757431 * AHV_Model_B.nu_r_f[4];

  // Derivatives for StateSpace: '<S272>/Dp(4,2)'
  _rtXdot->Dp42_CSTATE_h[0] +=
      -0.0069962390232986621 * AHV_Model_X.Dp42_CSTATE_h[0];
  _rtXdot->Dp42_CSTATE_h[0] += 1.046159104970257 * AHV_Model_X.Dp42_CSTATE_h[1];
  _rtXdot->Dp42_CSTATE_h[0] +=
      -0.07892900657949925 * AHV_Model_X.Dp42_CSTATE_h[2];
  _rtXdot->Dp42_CSTATE_h[0] +=
      -0.074692039015985229 * AHV_Model_X.Dp42_CSTATE_h[3];
  _rtXdot->Dp42_CSTATE_h[0] +=
      -0.0042999629253483171 * AHV_Model_X.Dp42_CSTATE_h[4];
  _rtXdot->Dp42_CSTATE_h[1] +=
      -1.0461591049702557 * AHV_Model_X.Dp42_CSTATE_h[0];
  _rtXdot->Dp42_CSTATE_h[1] +=
      -1.6140427135911086 * AHV_Model_X.Dp42_CSTATE_h[1];
  _rtXdot->Dp42_CSTATE_h[1] +=
      0.95129605479096324 * AHV_Model_X.Dp42_CSTATE_h[2];
  _rtXdot->Dp42_CSTATE_h[1] +=
      1.3710797263241543 * AHV_Model_X.Dp42_CSTATE_h[3];
  _rtXdot->Dp42_CSTATE_h[1] +=
      0.067632454414034648 * AHV_Model_X.Dp42_CSTATE_h[4];
  _rtXdot->Dp42_CSTATE_h[2] +=
      0.078929006579498431 * AHV_Model_X.Dp42_CSTATE_h[0];
  _rtXdot->Dp42_CSTATE_h[2] +=
      0.9512960547909608 * AHV_Model_X.Dp42_CSTATE_h[1];
  _rtXdot->Dp42_CSTATE_h[2] +=
      -1.065705011247746 * AHV_Model_X.Dp42_CSTATE_h[2];
  _rtXdot->Dp42_CSTATE_h[2] +=
      -3.3410455288652421 * AHV_Model_X.Dp42_CSTATE_h[3];
  _rtXdot->Dp42_CSTATE_h[2] +=
      -0.11061670514386444 * AHV_Model_X.Dp42_CSTATE_h[4];
  _rtXdot->Dp42_CSTATE_h[3] +=
      -0.074692039015982259 * AHV_Model_X.Dp42_CSTATE_h[0];
  _rtXdot->Dp42_CSTATE_h[3] +=
      -1.3710797263241621 * AHV_Model_X.Dp42_CSTATE_h[1];
  _rtXdot->Dp42_CSTATE_h[3] +=
      3.3410455288652465 * AHV_Model_X.Dp42_CSTATE_h[2];
  _rtXdot->Dp42_CSTATE_h[3] +=
      -5.4714163714089334 * AHV_Model_X.Dp42_CSTATE_h[3];
  _rtXdot->Dp42_CSTATE_h[3] +=
      -2.4649784545808426 * AHV_Model_X.Dp42_CSTATE_h[4];
  _rtXdot->Dp42_CSTATE_h[4] +=
      0.0042999629253487083 * AHV_Model_X.Dp42_CSTATE_h[0];
  _rtXdot->Dp42_CSTATE_h[4] +=
      0.067632454414049734 * AHV_Model_X.Dp42_CSTATE_h[1];
  _rtXdot->Dp42_CSTATE_h[4] +=
      -0.11061670514388669 * AHV_Model_X.Dp42_CSTATE_h[2];
  _rtXdot->Dp42_CSTATE_h[4] += 2.46497845458085 * AHV_Model_X.Dp42_CSTATE_h[3];
  _rtXdot->Dp42_CSTATE_h[4] +=
      -0.020742386481004484 * AHV_Model_X.Dp42_CSTATE_h[4];
  _rtXdot->Dp42_CSTATE_h[0] += -0.23872786278308805 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp42_CSTATE_h[1] += -3.2763464234276056 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp42_CSTATE_h[2] += 1.1437118751635387 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp42_CSTATE_h[3] += -1.3244904674438165 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp42_CSTATE_h[4] += 0.071163531145327891 * AHV_Model_B.nu_r_f[1];

  // Derivatives for StateSpace: '<S272>/Dp(4,6)'
  _rtXdot->Dp46_CSTATE_i[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp46_CSTATE_i[0];
  _rtXdot->Dp46_CSTATE_i[0] +=
      1.3927195985501903 * AHV_Model_X.Dp46_CSTATE_i[1];
  _rtXdot->Dp46_CSTATE_i[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp46_CSTATE_i[2];
  _rtXdot->Dp46_CSTATE_i[0] +=
      0.065100301994703638 * AHV_Model_X.Dp46_CSTATE_i[3];
  _rtXdot->Dp46_CSTATE_i[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp46_CSTATE_i[4];
  _rtXdot->Dp46_CSTATE_i[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp46_CSTATE_i[0];
  _rtXdot->Dp46_CSTATE_i[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp46_CSTATE_i[1];
  _rtXdot->Dp46_CSTATE_i[1] +=
      0.47138330474887147 * AHV_Model_X.Dp46_CSTATE_i[2];
  _rtXdot->Dp46_CSTATE_i[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp46_CSTATE_i[3];
  _rtXdot->Dp46_CSTATE_i[1] +=
      0.008356128979928646 * AHV_Model_X.Dp46_CSTATE_i[4];
  _rtXdot->Dp46_CSTATE_i[2] +=
      0.02751196060511972 * AHV_Model_X.Dp46_CSTATE_i[0];
  _rtXdot->Dp46_CSTATE_i[2] +=
      0.47138330474889245 * AHV_Model_X.Dp46_CSTATE_i[1];
  _rtXdot->Dp46_CSTATE_i[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp46_CSTATE_i[2];
  _rtXdot->Dp46_CSTATE_i[2] +=
      2.5521376035497743 * AHV_Model_X.Dp46_CSTATE_i[3];
  _rtXdot->Dp46_CSTATE_i[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp46_CSTATE_i[4];
  _rtXdot->Dp46_CSTATE_i[3] +=
      0.06510030199471134 * AHV_Model_X.Dp46_CSTATE_i[0];
  _rtXdot->Dp46_CSTATE_i[3] +=
      1.5107242348458856 * AHV_Model_X.Dp46_CSTATE_i[1];
  _rtXdot->Dp46_CSTATE_i[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp46_CSTATE_i[2];
  _rtXdot->Dp46_CSTATE_i[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp46_CSTATE_i[3];
  _rtXdot->Dp46_CSTATE_i[3] +=
      0.12040468647842413 * AHV_Model_X.Dp46_CSTATE_i[4];
  _rtXdot->Dp46_CSTATE_i[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp46_CSTATE_i[0];
  _rtXdot->Dp46_CSTATE_i[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp46_CSTATE_i[1];
  _rtXdot->Dp46_CSTATE_i[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp46_CSTATE_i[2];
  _rtXdot->Dp46_CSTATE_i[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp46_CSTATE_i[3];
  _rtXdot->Dp46_CSTATE_i[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp46_CSTATE_i[4];
  _rtXdot->Dp46_CSTATE_i[0] += -0.33524288945878539 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp46_CSTATE_i[1] += -6.3625492730866693 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp46_CSTATE_i[2] += 0.93578679415036736 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp46_CSTATE_i[3] += 2.5673752417819116 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp46_CSTATE_i[4] += 0.014933981938330865 * AHV_Model_B.nu_r_f[5];

  // Derivatives for StateSpace: '<S272>/Dp(5,3)'
  _rtXdot->Dp53_CSTATE_pz[0] +=
      -2.7914210326106108 * AHV_Model_X.Dp53_CSTATE_pz[0];
  _rtXdot->Dp53_CSTATE_pz[0] +=
      -0.99980762764400322 * AHV_Model_X.Dp53_CSTATE_pz[1];
  _rtXdot->Dp53_CSTATE_pz[0] +=
      2.1365278323686749 * AHV_Model_X.Dp53_CSTATE_pz[2];
  _rtXdot->Dp53_CSTATE_pz[0] +=
      1.6362305234027079 * AHV_Model_X.Dp53_CSTATE_pz[3];
  _rtXdot->Dp53_CSTATE_pz[0] +=
      0.259710121707262 * AHV_Model_X.Dp53_CSTATE_pz[4];
  _rtXdot->Dp53_CSTATE_pz[1] +=
      0.9998076276438177 * AHV_Model_X.Dp53_CSTATE_pz[0];
  _rtXdot->Dp53_CSTATE_pz[1] +=
      -5.3039287993722015E-5 * AHV_Model_X.Dp53_CSTATE_pz[1];
  _rtXdot->Dp53_CSTATE_pz[1] +=
      0.0077107435119831529 * AHV_Model_X.Dp53_CSTATE_pz[2];
  _rtXdot->Dp53_CSTATE_pz[1] +=
      0.0063309425482545112 * AHV_Model_X.Dp53_CSTATE_pz[3];
  _rtXdot->Dp53_CSTATE_pz[1] +=
      0.0012479014747342309 * AHV_Model_X.Dp53_CSTATE_pz[4];
  _rtXdot->Dp53_CSTATE_pz[2] +=
      -2.1365278323686576 * AHV_Model_X.Dp53_CSTATE_pz[0];
  _rtXdot->Dp53_CSTATE_pz[2] +=
      0.0077107435118432491 * AHV_Model_X.Dp53_CSTATE_pz[1];
  _rtXdot->Dp53_CSTATE_pz[2] +=
      -3.350088004888828 * AHV_Model_X.Dp53_CSTATE_pz[2];
  _rtXdot->Dp53_CSTATE_pz[2] +=
      -3.2190337160141729 * AHV_Model_X.Dp53_CSTATE_pz[3];
  _rtXdot->Dp53_CSTATE_pz[2] +=
      -1.6370488166464789 * AHV_Model_X.Dp53_CSTATE_pz[4];
  _rtXdot->Dp53_CSTATE_pz[3] +=
      -1.6362305234027152 * AHV_Model_X.Dp53_CSTATE_pz[0];
  _rtXdot->Dp53_CSTATE_pz[3] +=
      0.0063309425481496455 * AHV_Model_X.Dp53_CSTATE_pz[1];
  _rtXdot->Dp53_CSTATE_pz[3] +=
      -3.2190337160142368 * AHV_Model_X.Dp53_CSTATE_pz[2];
  _rtXdot->Dp53_CSTATE_pz[3] +=
      -3.2340681753232938 * AHV_Model_X.Dp53_CSTATE_pz[3];
  _rtXdot->Dp53_CSTATE_pz[3] +=
      -3.2557107662585683 * AHV_Model_X.Dp53_CSTATE_pz[4];
  _rtXdot->Dp53_CSTATE_pz[4] +=
      0.2597101217072279 * AHV_Model_X.Dp53_CSTATE_pz[0];
  _rtXdot->Dp53_CSTATE_pz[4] +=
      -0.00124790147473145 * AHV_Model_X.Dp53_CSTATE_pz[1];
  _rtXdot->Dp53_CSTATE_pz[4] +=
      1.6370488166464041 * AHV_Model_X.Dp53_CSTATE_pz[2];
  _rtXdot->Dp53_CSTATE_pz[4] +=
      3.2557107662584994 * AHV_Model_X.Dp53_CSTATE_pz[3];
  _rtXdot->Dp53_CSTATE_pz[4] +=
      -0.15597581472272126 * AHV_Model_X.Dp53_CSTATE_pz[4];
  _rtXdot->Dp53_CSTATE_pz[0] += -3.5374929885435322 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp53_CSTATE_pz[1] += 0.015233386783081164 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp53_CSTATE_pz[2] += -1.2196843916515756 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp53_CSTATE_pz[3] += -0.96955326725759594 * AHV_Model_B.nu_r_f[2];
  _rtXdot->Dp53_CSTATE_pz[4] += 0.17148246208757431 * AHV_Model_B.nu_r_f[2];

  // Derivatives for StateSpace: '<S272>/Dp(5,5)'
  _rtXdot->Dp55_CSTATE_a[0] +=
      -8.7466634083809418E-5 * AHV_Model_X.Dp55_CSTATE_a[0];
  _rtXdot->Dp55_CSTATE_a[0] +=
      0.70495353490101365 * AHV_Model_X.Dp55_CSTATE_a[1];
  _rtXdot->Dp55_CSTATE_a[0] +=
      0.0013720743421809071 * AHV_Model_X.Dp55_CSTATE_a[2];
  _rtXdot->Dp55_CSTATE_a[0] +=
      0.0075482812368147384 * AHV_Model_X.Dp55_CSTATE_a[3];
  _rtXdot->Dp55_CSTATE_a[0] +=
      0.00709599540037472 * AHV_Model_X.Dp55_CSTATE_a[4];
  _rtXdot->Dp55_CSTATE_a[1] +=
      -0.70495353490093393 * AHV_Model_X.Dp55_CSTATE_a[0];
  _rtXdot->Dp55_CSTATE_a[1] +=
      -2.2782193684765151 * AHV_Model_X.Dp55_CSTATE_a[1];
  _rtXdot->Dp55_CSTATE_a[1] +=
      -0.19644591246833848 * AHV_Model_X.Dp55_CSTATE_a[2];
  _rtXdot->Dp55_CSTATE_a[1] +=
      -1.3901163749547494 * AHV_Model_X.Dp55_CSTATE_a[3];
  _rtXdot->Dp55_CSTATE_a[1] +=
      -1.2786273904776069 * AHV_Model_X.Dp55_CSTATE_a[4];
  _rtXdot->Dp55_CSTATE_a[2] +=
      -0.0013720743421854747 * AHV_Model_X.Dp55_CSTATE_a[0];
  _rtXdot->Dp55_CSTATE_a[2] +=
      -0.19644591246830914 * AHV_Model_X.Dp55_CSTATE_a[1];
  _rtXdot->Dp55_CSTATE_a[2] +=
      -0.068186331750790974 * AHV_Model_X.Dp55_CSTATE_a[2];
  _rtXdot->Dp55_CSTATE_a[2] += -4.2372658827384 * AHV_Model_X.Dp55_CSTATE_a[3];
  _rtXdot->Dp55_CSTATE_a[2] +=
      -2.2216067374576487 * AHV_Model_X.Dp55_CSTATE_a[4];
  _rtXdot->Dp55_CSTATE_a[3] +=
      0.007548281236770463 * AHV_Model_X.Dp55_CSTATE_a[0];
  _rtXdot->Dp55_CSTATE_a[3] +=
      1.3901163749547758 * AHV_Model_X.Dp55_CSTATE_a[1];
  _rtXdot->Dp55_CSTATE_a[3] +=
      4.2372658827383063 * AHV_Model_X.Dp55_CSTATE_a[2];
  _rtXdot->Dp55_CSTATE_a[3] +=
      -3.3091650540438455 * AHV_Model_X.Dp55_CSTATE_a[3];
  _rtXdot->Dp55_CSTATE_a[3] +=
      -3.4108951000309653 * AHV_Model_X.Dp55_CSTATE_a[4];
  _rtXdot->Dp55_CSTATE_a[4] +=
      0.0070959954003254081 * AHV_Model_X.Dp55_CSTATE_a[0];
  _rtXdot->Dp55_CSTATE_a[4] +=
      1.2786273904776095 * AHV_Model_X.Dp55_CSTATE_a[1];
  _rtXdot->Dp55_CSTATE_a[4] +=
      2.2216067374575723 * AHV_Model_X.Dp55_CSTATE_a[2];
  _rtXdot->Dp55_CSTATE_a[4] +=
      -3.4108951000310745 * AHV_Model_X.Dp55_CSTATE_a[3];
  _rtXdot->Dp55_CSTATE_a[4] +=
      -3.5569423734797287 * AHV_Model_X.Dp55_CSTATE_a[4];
  _rtXdot->Dp55_CSTATE_a[0] += 0.021904981880023978 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp55_CSTATE_a[1] += 3.4651622640804733 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp55_CSTATE_a[2] += 0.16004490113712252 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp55_CSTATE_a[3] += -0.9969365784863482 * AHV_Model_B.nu_r_f[4];
  _rtXdot->Dp55_CSTATE_a[4] += -0.92774837285087253 * AHV_Model_B.nu_r_f[4];

  // Derivatives for StateSpace: '<S272>/Dp(6,2)'
  _rtXdot->Dp62_CSTATE_i[0] +=
      -0.0017008072822443973 * AHV_Model_X.Dp62_CSTATE_i[0];
  _rtXdot->Dp62_CSTATE_i[0] +=
      1.2361995937266759 * AHV_Model_X.Dp62_CSTATE_i[1];
  _rtXdot->Dp62_CSTATE_i[0] +=
      0.012808448456780379 * AHV_Model_X.Dp62_CSTATE_i[2];
  _rtXdot->Dp62_CSTATE_i[0] +=
      0.031758809083141569 * AHV_Model_X.Dp62_CSTATE_i[3];
  _rtXdot->Dp62_CSTATE_i[0] +=
      -0.0047828760333047419 * AHV_Model_X.Dp62_CSTATE_i[4];
  _rtXdot->Dp62_CSTATE_i[1] +=
      -1.236199593726705 * AHV_Model_X.Dp62_CSTATE_i[0];
  _rtXdot->Dp62_CSTATE_i[1] +=
      -1.2300467663140022 * AHV_Model_X.Dp62_CSTATE_i[1];
  _rtXdot->Dp62_CSTATE_i[1] +=
      -0.28770138559075209 * AHV_Model_X.Dp62_CSTATE_i[2];
  _rtXdot->Dp62_CSTATE_i[1] +=
      -0.96895750705376871 * AHV_Model_X.Dp62_CSTATE_i[3];
  _rtXdot->Dp62_CSTATE_i[1] +=
      0.12995794495601007 * AHV_Model_X.Dp62_CSTATE_i[4];
  _rtXdot->Dp62_CSTATE_i[2] +=
      -0.012808448456765344 * AHV_Model_X.Dp62_CSTATE_i[0];
  _rtXdot->Dp62_CSTATE_i[2] +=
      -0.28770138559072356 * AHV_Model_X.Dp62_CSTATE_i[1];
  _rtXdot->Dp62_CSTATE_i[2] +=
      -0.18502978066497214 * AHV_Model_X.Dp62_CSTATE_i[2];
  _rtXdot->Dp62_CSTATE_i[2] +=
      -1.8113135858263947 * AHV_Model_X.Dp62_CSTATE_i[3];
  _rtXdot->Dp62_CSTATE_i[2] +=
      0.1356516649888091 * AHV_Model_X.Dp62_CSTATE_i[4];
  _rtXdot->Dp62_CSTATE_i[3] +=
      0.031758809083164766 * AHV_Model_X.Dp62_CSTATE_i[0];
  _rtXdot->Dp62_CSTATE_i[3] +=
      0.96895750705376427 * AHV_Model_X.Dp62_CSTATE_i[1];
  _rtXdot->Dp62_CSTATE_i[3] +=
      1.8113135858263727 * AHV_Model_X.Dp62_CSTATE_i[2];
  _rtXdot->Dp62_CSTATE_i[3] +=
      -3.762337151306161 * AHV_Model_X.Dp62_CSTATE_i[3];
  _rtXdot->Dp62_CSTATE_i[3] +=
      1.5320613622870083 * AHV_Model_X.Dp62_CSTATE_i[4];
  _rtXdot->Dp62_CSTATE_i[4] +=
      0.0047828760333030167 * AHV_Model_X.Dp62_CSTATE_i[0];
  _rtXdot->Dp62_CSTATE_i[4] +=
      0.12995794495600341 * AHV_Model_X.Dp62_CSTATE_i[1];
  _rtXdot->Dp62_CSTATE_i[4] +=
      0.13565166498879488 * AHV_Model_X.Dp62_CSTATE_i[2];
  _rtXdot->Dp62_CSTATE_i[4] +=
      -1.5320613622870141 * AHV_Model_X.Dp62_CSTATE_i[3];
  _rtXdot->Dp62_CSTATE_i[4] +=
      -0.25499556157548836 * AHV_Model_X.Dp62_CSTATE_i[4];
  _rtXdot->Dp62_CSTATE_i[0] += 0.13167647316600378 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp62_CSTATE_i[1] += 3.4125284827245985 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp62_CSTATE_i[2] += 0.44399732492152644 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp62_CSTATE_i[3] += -1.2820687298457427 * AHV_Model_B.nu_r_f[1];
  _rtXdot->Dp62_CSTATE_i[4] += -0.18276601298221137 * AHV_Model_B.nu_r_f[1];

  // Derivatives for StateSpace: '<S272>/Dp(6,4)'
  _rtXdot->Dp64_CSTATE_e[0] +=
      -0.004396104715914141 * AHV_Model_X.Dp64_CSTATE_e[0];
  _rtXdot->Dp64_CSTATE_e[0] +=
      1.3927195985501903 * AHV_Model_X.Dp64_CSTATE_e[1];
  _rtXdot->Dp64_CSTATE_e[0] +=
      -0.027511960605097124 * AHV_Model_X.Dp64_CSTATE_e[2];
  _rtXdot->Dp64_CSTATE_e[0] +=
      0.065100301994703638 * AHV_Model_X.Dp64_CSTATE_e[3];
  _rtXdot->Dp64_CSTATE_e[0] +=
      -0.0003953747088775357 * AHV_Model_X.Dp64_CSTATE_e[4];
  _rtXdot->Dp64_CSTATE_e[1] +=
      -1.3927195985502012 * AHV_Model_X.Dp64_CSTATE_e[0];
  _rtXdot->Dp64_CSTATE_e[1] +=
      -1.7990221690034174 * AHV_Model_X.Dp64_CSTATE_e[1];
  _rtXdot->Dp64_CSTATE_e[1] +=
      0.47138330474887147 * AHV_Model_X.Dp64_CSTATE_e[2];
  _rtXdot->Dp64_CSTATE_e[1] +=
      -1.5107242348458776 * AHV_Model_X.Dp64_CSTATE_e[3];
  _rtXdot->Dp64_CSTATE_e[1] +=
      0.008356128979928646 * AHV_Model_X.Dp64_CSTATE_e[4];
  _rtXdot->Dp64_CSTATE_e[2] +=
      0.02751196060511972 * AHV_Model_X.Dp64_CSTATE_e[0];
  _rtXdot->Dp64_CSTATE_e[2] +=
      0.47138330474889245 * AHV_Model_X.Dp64_CSTATE_e[1];
  _rtXdot->Dp64_CSTATE_e[2] +=
      -0.31733520843172491 * AHV_Model_X.Dp64_CSTATE_e[2];
  _rtXdot->Dp64_CSTATE_e[2] +=
      2.5521376035497743 * AHV_Model_X.Dp64_CSTATE_e[3];
  _rtXdot->Dp64_CSTATE_e[2] +=
      -0.0093184009624265578 * AHV_Model_X.Dp64_CSTATE_e[4];
  _rtXdot->Dp64_CSTATE_e[3] +=
      0.06510030199471134 * AHV_Model_X.Dp64_CSTATE_e[0];
  _rtXdot->Dp64_CSTATE_e[3] +=
      1.5107242348458856 * AHV_Model_X.Dp64_CSTATE_e[1];
  _rtXdot->Dp64_CSTATE_e[3] +=
      -2.5521376035497325 * AHV_Model_X.Dp64_CSTATE_e[2];
  _rtXdot->Dp64_CSTATE_e[3] +=
      -7.5177108865705646 * AHV_Model_X.Dp64_CSTATE_e[3];
  _rtXdot->Dp64_CSTATE_e[3] +=
      0.12040468647842413 * AHV_Model_X.Dp64_CSTATE_e[4];
  _rtXdot->Dp64_CSTATE_e[4] +=
      0.00039537470888004991 * AHV_Model_X.Dp64_CSTATE_e[0];
  _rtXdot->Dp64_CSTATE_e[4] +=
      0.0083561289799246388 * AHV_Model_X.Dp64_CSTATE_e[1];
  _rtXdot->Dp64_CSTATE_e[4] +=
      -0.0093184009624202816 * AHV_Model_X.Dp64_CSTATE_e[2];
  _rtXdot->Dp64_CSTATE_e[4] +=
      -0.12040468647838372 * AHV_Model_X.Dp64_CSTATE_e[3];
  _rtXdot->Dp64_CSTATE_e[4] +=
      -0.00092959984608804514 * AHV_Model_X.Dp64_CSTATE_e[4];
  _rtXdot->Dp64_CSTATE_e[0] += -0.33524288945878539 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp64_CSTATE_e[1] += -6.3625492730866693 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp64_CSTATE_e[2] += 0.93578679415036736 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp64_CSTATE_e[3] += 2.5673752417819116 * AHV_Model_B.nu_r_f[3];
  _rtXdot->Dp64_CSTATE_e[4] += 0.014933981938330865 * AHV_Model_B.nu_r_f[3];

  // Derivatives for StateSpace: '<S272>/Dp(6,6)'
  _rtXdot->Dp66_CSTATE_g3[0] +=
      -1.5814921936539046 * AHV_Model_X.Dp66_CSTATE_g3[0];
  _rtXdot->Dp66_CSTATE_g3[0] +=
      -1.3217886535784036 * AHV_Model_X.Dp66_CSTATE_g3[1];
  _rtXdot->Dp66_CSTATE_g3[0] +=
      -0.43878308112001396 * AHV_Model_X.Dp66_CSTATE_g3[2];
  _rtXdot->Dp66_CSTATE_g3[0] +=
      -1.2174625158903036 * AHV_Model_X.Dp66_CSTATE_g3[3];
  _rtXdot->Dp66_CSTATE_g3[0] +=
      0.069477290559771143 * AHV_Model_X.Dp66_CSTATE_g3[4];
  _rtXdot->Dp66_CSTATE_g3[1] +=
      1.3217886535783951 * AHV_Model_X.Dp66_CSTATE_g3[0];
  _rtXdot->Dp66_CSTATE_g3[1] +=
      -0.00029199832472599611 * AHV_Model_X.Dp66_CSTATE_g3[1];
  _rtXdot->Dp66_CSTATE_g3[1] +=
      -0.0080078570280257711 * AHV_Model_X.Dp66_CSTATE_g3[2];
  _rtXdot->Dp66_CSTATE_g3[1] +=
      -0.01568400806720225 * AHV_Model_X.Dp66_CSTATE_g3[3];
  _rtXdot->Dp66_CSTATE_g3[1] +=
      0.00098472174210678556 * AHV_Model_X.Dp66_CSTATE_g3[4];
  _rtXdot->Dp66_CSTATE_g3[2] +=
      -0.43878308111999881 * AHV_Model_X.Dp66_CSTATE_g3[0];
  _rtXdot->Dp66_CSTATE_g3[2] +=
      0.00800785702802341 * AHV_Model_X.Dp66_CSTATE_g3[1];
  _rtXdot->Dp66_CSTATE_g3[2] +=
      -0.28893903358444806 * AHV_Model_X.Dp66_CSTATE_g3[2];
  _rtXdot->Dp66_CSTATE_g3[2] +=
      -1.8221092762782123 * AHV_Model_X.Dp66_CSTATE_g3[3];
  _rtXdot->Dp66_CSTATE_g3[2] +=
      0.074587161707531782 * AHV_Model_X.Dp66_CSTATE_g3[4];
  _rtXdot->Dp66_CSTATE_g3[3] +=
      1.2174625158903034 * AHV_Model_X.Dp66_CSTATE_g3[0];
  _rtXdot->Dp66_CSTATE_g3[3] +=
      -0.015684008067190565 * AHV_Model_X.Dp66_CSTATE_g3[1];
  _rtXdot->Dp66_CSTATE_g3[3] +=
      1.8221092762782234 * AHV_Model_X.Dp66_CSTATE_g3[2];
  _rtXdot->Dp66_CSTATE_g3[3] +=
      -6.3850648427943675 * AHV_Model_X.Dp66_CSTATE_g3[3];
  _rtXdot->Dp66_CSTATE_g3[3] +=
      1.217430977302715 * AHV_Model_X.Dp66_CSTATE_g3[4];
  _rtXdot->Dp66_CSTATE_g3[4] +=
      0.069477290559780386 * AHV_Model_X.Dp66_CSTATE_g3[0];
  _rtXdot->Dp66_CSTATE_g3[4] +=
      -0.000984721742106249 * AHV_Model_X.Dp66_CSTATE_g3[1];
  _rtXdot->Dp66_CSTATE_g3[4] +=
      0.074587161707547089 * AHV_Model_X.Dp66_CSTATE_g3[2];
  _rtXdot->Dp66_CSTATE_g3[4] +=
      -1.2174309773027607 * AHV_Model_X.Dp66_CSTATE_g3[3];
  _rtXdot->Dp66_CSTATE_g3[4] +=
      -0.061420217921035136 * AHV_Model_X.Dp66_CSTATE_g3[4];
  _rtXdot->Dp66_CSTATE_g3[0] += -3.3844821061206916 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp66_CSTATE_g3[1] += 0.045246867266885989 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp66_CSTATE_g3[2] += -0.53332191483729374 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp66_CSTATE_g3[3] += 1.2581308985447515 * AHV_Model_B.nu_r_f[5];
  _rtXdot->Dp66_CSTATE_g3[4] += 0.075289942321940723 * AHV_Model_B.nu_r_f[5];

  // Derivatives for Integrator: '<S315>/Integrator1'
  _rtXdot->Integrator1_CSTATE_me[0] = AHV_Model_B.Sum6_m[0];

  // Derivatives for Integrator: '<S315>/Integrator2'
  _rtXdot->Integrator2_CSTATE_h[0] = AHV_Model_B.psi_WF_d[0];

  // Derivatives for Integrator: '<S315>/Integrator6'
  _rtXdot->Integrator6_CSTATE_f[0] = AHV_Model_B.Sum7_g[0];

  // Derivatives for Integrator: '<S315>/Integrator1'
  _rtXdot->Integrator1_CSTATE_me[1] = AHV_Model_B.Sum6_m[1];

  // Derivatives for Integrator: '<S315>/Integrator2'
  _rtXdot->Integrator2_CSTATE_h[1] = AHV_Model_B.psi_WF_d[1];

  // Derivatives for Integrator: '<S315>/Integrator6'
  _rtXdot->Integrator6_CSTATE_f[1] = AHV_Model_B.Sum7_g[1];

  // Derivatives for Integrator: '<S315>/Integrator1'
  _rtXdot->Integrator1_CSTATE_me[2] = AHV_Model_B.Sum6_m[2];

  // Derivatives for Integrator: '<S315>/Integrator2'
  _rtXdot->Integrator2_CSTATE_h[2] = AHV_Model_B.psi_WF_d[2];

  // Derivatives for Integrator: '<S315>/Integrator6'
  _rtXdot->Integrator6_CSTATE_f[2] = AHV_Model_B.Sum7_g[2];
}

// Model initialize function
void AH_Model_v1ModelClass::initialize(double dtime) {
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&(&AHV_Model_M)->solverInfo,
                          &(&AHV_Model_M)->Timing.simTimeStep);
    rtsiSetTPtr(&(&AHV_Model_M)->solverInfo, &rtmGetTPtr((&AHV_Model_M)));
    rtsiSetStepSizePtr(&(&AHV_Model_M)->solverInfo,
                       &(&AHV_Model_M)->Timing.stepSize0);
    rtsiSetdXPtr(&(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)->derivs);
    rtsiSetContStatesPtr(&(&AHV_Model_M)->solverInfo,
                         (real_T**)&(&AHV_Model_M)->contStates);
    rtsiSetNumContStatesPtr(&(&AHV_Model_M)->solverInfo,
                            &(&AHV_Model_M)->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(
        &(&AHV_Model_M)->solverInfo,
        &(&AHV_Model_M)->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(
        &(&AHV_Model_M)->solverInfo, &(&AHV_Model_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&AHV_Model_M)->solverInfo,
                                      &(&AHV_Model_M)->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&(&AHV_Model_M)->solverInfo,
                          (&rtmGetErrorStatus((&AHV_Model_M))));
    rtsiSetRTModelPtr(&(&AHV_Model_M)->solverInfo, (&AHV_Model_M));
  }

  rtsiSetSimTimeStep(&(&AHV_Model_M)->solverInfo, MAJOR_TIME_STEP);
  (&AHV_Model_M)->intgData.y = (&AHV_Model_M)->odeY;
  (&AHV_Model_M)->intgData.f[0] = (&AHV_Model_M)->odeF[0];
  (&AHV_Model_M)->intgData.f[1] = (&AHV_Model_M)->odeF[1];
  (&AHV_Model_M)->intgData.f[2] = (&AHV_Model_M)->odeF[2];
  (&AHV_Model_M)->intgData.f[3] = (&AHV_Model_M)->odeF[3];
  (&AHV_Model_M)->contStates = ((X_AHV_Model_T*)&AHV_Model_X);
  rtsiSetSolverData(&(&AHV_Model_M)->solverInfo,
                    static_cast<void*>(&(&AHV_Model_M)->intgData));
  rtsiSetSolverName(&(&AHV_Model_M)->solverInfo, "ode4");
  rtmSetTPtr((&AHV_Model_M), &(&AHV_Model_M)->Timing.tArray[0]);
  (&AHV_Model_M)->Timing.stepSize0 = dtime;
  rtmSetFirstInitCond((&AHV_Model_M), 1);

  {
    int32_T i;

    // Start for If: '<S83>/If'
    AHV_Model_DW.If_ActiveSubsystem = -1;

    // Start for If: '<S84>/If'
    AHV_Model_DW.If_ActiveSubsystem_a = -1;

    // Start for If: '<S76>/If'
    AHV_Model_DW.If_ActiveSubsystem_i = -1;

    // Start for If: '<S167>/If'
    AHV_Model_DW.If_ActiveSubsystem_h = -1;

    // Start for If: '<S168>/If'
    AHV_Model_DW.If_ActiveSubsystem_ig = -1;

    // Start for If: '<S160>/If'
    AHV_Model_DW.If_ActiveSubsystem_ab = -1;

    // Start for If: '<S251>/If'
    AHV_Model_DW.If_ActiveSubsystem_e = -1;

    // Start for If: '<S252>/If'
    AHV_Model_DW.If_ActiveSubsystem_k = -1;

    // Start for If: '<S244>/If'
    AHV_Model_DW.If_ActiveSubsystem_ef = -1;

    // Start for If: '<S335>/If'
    AHV_Model_DW.If_ActiveSubsystem_f = -1;

    // Start for If: '<S336>/If'
    AHV_Model_DW.If_ActiveSubsystem_hp = -1;

    // Start for If: '<S328>/If'
    AHV_Model_DW.If_ActiveSubsystem_i4 = -1;

    // InitializeConditions for Integrator: '<S13>/Integrator1'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK = 1;

    // End of InitializeConditions for Integrator: '<S13>/Integrator1'

    // InitializeConditions for Integrator: '<S13>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S13>/Integrator'

    // InitializeConditions for TransferFcn: '<S16>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK = 1;

    // End of InitializeConditions for Integrator: '<S63>/Integrator3'

    // InitializeConditions for RateLimiter: '<S64>/Rate Limiter'
    AHV_Model_DW.LastMajorTime = (rtInf);

    // InitializeConditions for Integrator: '<S65>/Integrator'
    AHV_Model_X.Integrator_CSTATE_b[0] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE[0] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator'
    AHV_Model_X.Integrator_CSTATE_b[1] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE[1] = 0.0;

    // InitializeConditions for Integrator: '<S65>/Integrator'
    AHV_Model_X.Integrator_CSTATE_b[2] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S20>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE[i] = 0.0;

      // InitializeConditions for Integrator: '<S32>/Integrator'
      AHV_Model_X.Integrator_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE[i] = 0.0;

      // InitializeConditions for StateSpace: '<S20>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S63>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_c[0] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[0] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[0] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_c[1] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[1] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[1] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_c[2] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE[2] = 0.0;

    // InitializeConditions for Integrator: '<S63>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE[2] = 0.0;

    // InitializeConditions for Integrator: '<S97>/Integrator1'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE_l[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_l[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_l[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_l[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_l[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_l[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK_a = 1;

    // End of InitializeConditions for Integrator: '<S97>/Integrator1'

    // InitializeConditions for Integrator: '<S97>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_n[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S97>/Integrator'

    // InitializeConditions for TransferFcn: '<S100>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_e = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_i[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_i[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_i[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_p = 1;

    // End of InitializeConditions for Integrator: '<S147>/Integrator3'

    // InitializeConditions for RateLimiter: '<S148>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_n = (rtInf);

    // InitializeConditions for Integrator: '<S147>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_p[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_p[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_p[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S104>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_j[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_m[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_n[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_c[i] = 0.0;

      // InitializeConditions for Integrator: '<S116>/Integrator'
      AHV_Model_X.Integrator_CSTATE_m[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_o[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S104>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_g[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S147>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_d[0] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_i[0] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_b[0] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_d[1] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_i[1] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_b[1] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_d[2] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_i[2] = 0.0;

    // InitializeConditions for Integrator: '<S147>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_b[2] = 0.0;

    // InitializeConditions for Integrator: '<S181>/Integrator1'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE_lj[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_lj[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_lj[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_lj[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_lj[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_lj[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK_p = 1;

    // End of InitializeConditions for Integrator: '<S181>/Integrator1'

    // InitializeConditions for Integrator: '<S181>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_o[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S181>/Integrator'

    // InitializeConditions for TransferFcn: '<S184>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_l = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_g[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_g[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_g[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_n = 1;

    // End of InitializeConditions for Integrator: '<S231>/Integrator3'

    // InitializeConditions for RateLimiter: '<S232>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_a = (rtInf);

    // InitializeConditions for Integrator: '<S231>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_d[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_d[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_d[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S188>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_j[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_p[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_g[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_o[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_f[i] = 0.0;

      // InitializeConditions for Integrator: '<S200>/Integrator'
      AHV_Model_X.Integrator_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_h[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_n[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_b[i] = 0.0;

      // InitializeConditions for StateSpace: '<S188>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_n[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S231>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_m[0] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[0] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_e[0] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_m[1] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[1] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_e[1] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_m[2] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_j[2] = 0.0;

    // InitializeConditions for Integrator: '<S231>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_e[2] = 0.0;

    // InitializeConditions for Integrator: '<S265>/Integrator1'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator1_CSTATE_b[0] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_b[1] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_b[2] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_b[3] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_b[4] = 0.0;
      AHV_Model_X.Integrator1_CSTATE_b[5] = 0.0;
    }

    AHV_Model_DW.Integrator1_IWORK_f = 1;

    // End of InitializeConditions for Integrator: '<S265>/Integrator1'

    // InitializeConditions for Integrator: '<S265>/Integrator'
    for (i = 0; i < 6; i++) {
      AHV_Model_X.Integrator_CSTATE_a[i] = 0.0;
    }

    // End of InitializeConditions for Integrator: '<S265>/Integrator'

    // InitializeConditions for TransferFcn: '<S268>/Transfer Fcn'
    AHV_Model_X.TransferFcn_CSTATE_m = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator3'
    if (rtmIsFirstInitCond((&AHV_Model_M))) {
      AHV_Model_X.Integrator3_CSTATE_d[0] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_d[1] = 0.0;
      AHV_Model_X.Integrator3_CSTATE_d[2] = 0.0;
    }

    AHV_Model_DW.Integrator3_IWORK_c = 1;

    // End of InitializeConditions for Integrator: '<S315>/Integrator3'

    // InitializeConditions for RateLimiter: '<S316>/Rate Limiter'
    AHV_Model_DW.LastMajorTime_m = (rtInf);

    // InitializeConditions for Integrator: '<S315>/Integrator4'
    AHV_Model_X.Integrator4_CSTATE_o[0] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_o[1] = 0.0;
    AHV_Model_X.Integrator4_CSTATE_o[2] = 0.0;
    for (i = 0; i < 5; i++) {
      // InitializeConditions for StateSpace: '<S272>/Dp(1,1)'
      AHV_Model_X.Dp11_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(1,3)'
      AHV_Model_X.Dp13_CSTATE_c[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(1,5)'
      AHV_Model_X.Dp15_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(2,2)'
      AHV_Model_X.Dp22_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(2,4)'
      AHV_Model_X.Dp24_CSTATE_d[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(2,6)'
      AHV_Model_X.Dp26_CSTATE_f[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(3,1)'
      AHV_Model_X.Dp31_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(3,3)'
      AHV_Model_X.Dp33_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(3,5)'
      AHV_Model_X.Dp35_CSTATE_l[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(4,2)'
      AHV_Model_X.Dp42_CSTATE_h[i] = 0.0;

      // InitializeConditions for Integrator: '<S284>/Integrator'
      AHV_Model_X.Integrator_CSTATE_bp[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(4,6)'
      AHV_Model_X.Dp46_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(5,1)'
      AHV_Model_X.Dp51_CSTATE_k[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(5,3)'
      AHV_Model_X.Dp53_CSTATE_pz[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(5,5)'
      AHV_Model_X.Dp55_CSTATE_a[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(6,2)'
      AHV_Model_X.Dp62_CSTATE_i[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(6,4)'
      AHV_Model_X.Dp64_CSTATE_e[i] = 0.0;

      // InitializeConditions for StateSpace: '<S272>/Dp(6,6)'
      AHV_Model_X.Dp66_CSTATE_g3[i] = 0.0;
    }

    // InitializeConditions for Integrator: '<S315>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_me[0] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_h[0] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_f[0] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_me[1] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_h[1] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_f[1] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator1'
    AHV_Model_X.Integrator1_CSTATE_me[2] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator2'
    AHV_Model_X.Integrator2_CSTATE_h[2] = 0.0;

    // InitializeConditions for Integrator: '<S315>/Integrator6'
    AHV_Model_X.Integrator6_CSTATE_f[2] = 0.0;

    // SystemInitialize for Atomic SubSystem: '<S8>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3);

    // End of SystemInitialize for SubSystem: '<S8>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S92>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_b);

    // End of SystemInitialize for SubSystem: '<S92>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S176>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_d);

    // End of SystemInitialize for SubSystem: '<S176>/Subsystem3'

    // SystemInitialize for Atomic SubSystem: '<S260>/Subsystem3'
    Wave_init_Init(&AHV_Model_DW.Subsystem3_m);

    // End of SystemInitialize for SubSystem: '<S260>/Subsystem3'
  }

  // set "at time zero" to false
  if (rtmIsFirstInitCond((&AHV_Model_M))) {
    rtmSetFirstInitCond((&AHV_Model_M), 0);
  }
}

// Model terminate function
void AH_Model_v1ModelClass::terminate() {
  // (no terminate code required)
}

// Constructor
AH_Model_v1ModelClass::AH_Model_v1ModelClass()
    : AHV_Model_B(), AHV_Model_DW(), AHV_Model_X(), AHV_Model_M() {
  // Currently there is no constructor body generated.
}

// Destructor
AH_Model_v1ModelClass::~AH_Model_v1ModelClass() {
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
AH_Model_v1ModelClass::RT_MODEL_AHV_Model_T* AH_Model_v1ModelClass::getRTM() {
  return (&AHV_Model_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
