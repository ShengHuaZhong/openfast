//
// File: AHV_Model.h
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
#ifndef RTW_HEADER_AHV_Model_h_
#define RTW_HEADER_AHV_Model_h_
#include <cfloat>
#include <cmath>
#include <cstring>
#include <math.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "AHV_Model_types.h"

// Child system includes
#include "Wave.h"
#include "Wave_loads.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "rtGetNaN.h"
#include "rt_defines.h"

// Macros for accessing real-time model data structure
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef ODE4_INTG
#define ODE4_INTG

// ODE4 Integration Data
typedef struct {
  real_T *y;                           // output
  real_T *f[4];                        // derivatives
} ODE4_IntgData;

#endif

//
//  Exported Global Signals
//
//  Note: Exported global signals are block signals with an exported global
//  storage class designation.  Code generation will declare the memory for
//  these signals and export their symbols.
//

#ifdef __cplusplus

extern "C" {

#endif

  extern real_T spectrum_type;         // '<Root>/spectrum_type'
  extern real_T hs;                    // '<Root>/hs'
  extern real_T omega_peak;            // '<Root>/omega_peak'
  extern real_T psi_mean;              // '<Root>/psi_mean'
  extern real_T gamma_value;           // '<Root>/gamma_value'
  extern real_T spread;                // '<Root>/spread'
  extern real_T depth;                 // '<Root>/depth'
  extern real_T nfreq;                 // '<Root>/nfreq'
  extern real_T ndir;                  // '<Root>/ndir'
  extern real_T energylim;             // '<Root>/energylim'
  extern real_T freq_cutoff;           // '<Root>/freq_cutoff'
  extern real_T dir_cutoff;            // '<Root>/dir_cutoff'
  extern real_T rand_freq;             // '<Root>/rand_freq'
  extern real_T rand_dir;              // '<Root>/rand_dir'
  extern real_T Current_direction;     // '<Root>/Current_direction'
  extern real_T Current_speed;         // '<Root>/Current_speed'
  extern real_T Vessel_init1[6];       // '<Root>/Vessel_init1'
  extern boolean_T Hold_Position1;     // '<Root>/Hold_Position1'
  extern real_T heading_mode1;         // '<Root>/heading_mode1'
  extern real_T heading_angle_ref1;    // '<Root>/heading_angle_ref1'
  extern real_T Vessel_X_Ref1;         // '<Root>/Vessel_X_Ref1'
  extern real_T Vessel_Y_Ref1;         // '<Root>/Vessel_Y_Ref1'
  extern real_T ahv_fairlead1[3];      // '<Root>/ahv_fairlead1'
  extern real_T tau_cable1[3];         // '<Root>/tau_cable1'
  extern real_T Vessel_init2[6];       // '<Root>/Vessel_init2'
  extern boolean_T Hold_Position2;     // '<Root>/Hold_Position2'
  extern real_T heading_mode2;         // '<Root>/heading_mode2'
  extern real_T heading_angle_ref2;    // '<Root>/heading_angle_ref2'
  extern real_T Vessel_X_Ref2;         // '<Root>/Vessel_X_Ref2'
  extern real_T Vessel_Y_Ref2;         // '<Root>/Vessel_Y_Ref2'
  extern real_T ahv_fairlead2[3];      // '<Root>/ahv_fairlead2'
  extern real_T tau_cable2[3];         // '<Root>/tau_cable2'
  extern real_T Vessel_init3[6];       // '<Root>/Vessel_init3'
  extern boolean_T Hold_Position3;     // '<Root>/Hold_Position3'
  extern real_T heading_mode3;         // '<Root>/heading_mode3'
  extern real_T heading_angle_ref3;    // '<Root>/heading_angle_ref3'
  extern real_T Vessel_X_Ref3;         // '<Root>/Vessel_X_Ref3'
  extern real_T Vessel_Y_Ref3;         // '<Root>/Vessel_Y_Ref3'
  extern real_T ahv_fairlead3[3];      // '<Root>/ahv_fairlead3'
  extern real_T tau_cable3[3];         // '<Root>/tau_cable3'
  extern real_T Vessel_init4[6];       // '<Root>/Vessel_init4'
  extern boolean_T Hold_Position4;     // '<Root>/Hold_Position4'
  extern real_T heading_mode4;         // '<Root>/heading_mode4'
  extern real_T heading_angle_ref4;    // '<Root>/heading_angle_ref4'
  extern real_T Vessel_X_Ref4;         // '<Root>/Vessel_X_Ref4'
  extern real_T Vessel_Y_Ref4;         // '<Root>/Vessel_Y_Ref4'
  extern real_T ahv_fairlead4[3];      // '<Root>/ahv_fairlead4'
  extern real_T tau_cable4[3];         // '<Root>/tau_cable4'
  extern real_T eta_AHV1[6];           // '<Root>/eta_AHV1'
  extern real_T nu1[6];                // '<Root>/nu1'
  extern real_T AHV_speed1;            // '<Root>/AHV_speed1'
  extern real_T eta_AHV2[6];           // '<Root>/eta_AHV2'
  extern real_T nu2[6];                // '<Root>/nu2'
  extern real_T AHV_speed2;            // '<Root>/AHV_speed2'
  extern real_T eta_AHV3[6];           // '<Root>/eta_AHV3'
  extern real_T nu3[6];                // '<Root>/nu3'
  extern real_T AHV_speed3;            // '<Root>/AHV_speed3'
  extern real_T eta_AHV4[6];           // '<Root>/eta_AHV4'
  extern real_T nu4[6];                // '<Root>/nu4'
  extern real_T AHV_speed4;            // '<Root>/AHV_speed4'

#ifdef __cplusplus

}
#endif

// Class declaration for model AHV_Model
class AH_Model_v1ModelClass {
  // public data and function members
 public:
  // Block signals for system '<S8>/Subsystem3'
  typedef struct {
    real_T Psi[900];
    real_T Wavenum[900];
    real_T Phase[900];
    real_T Omega[900];
    real_T Zeta_a[900];
  } B_Wave_init_T;

  // Block states (default storage) for system '<S8>/Subsystem3'
  typedef struct {
    uint32_T state[625];               // '<S14>/Wave'
  } DW_Wave_init_T;

  // Block signals for system '<S15>/Wave_loads_for_heading1'
  typedef struct {
    real_T sin_h[900];                 // '<S40>/sin'
    real_T cos_d[900];                 // '<S40>/cos'
    real_T clock_Omega[900];           // '<S35>/Product1'
    real_T Phase_tot[900];             // '<S35>/Sum'
    real_T psi_r[900];                 // '<S35>/Sum1'
  } B_Wave_loads_for_heading1_AHV_T;

  // Block states (default storage) for system '<S15>/Wave_loads_for_heading1'
  typedef struct {
    real_T Delay1_DSTATE;              // '<S41>/Delay1'
    real_T m_bpLambda[3];            // '<S44>/Lookup Table (n-D)  Wavedrift 1'
    real_T m_bpLambda_c[3];          // '<S44>/Lookup Table (n-D)  Wavedrift 2'
    real_T m_bpLambda_k[3];          // '<S44>/Lookup Table (n-D)  Wavedrift 3'
    real_T m_bpLambda_o[3];            // '<S43>/Lookup Table (n-D) Amp1'
    real_T m_bpLambda_p[3];            // '<S43>/Lookup Table (n-D) Amp2'
    real_T m_bpLambda_b[3];            // '<S43>/Lookup Table (n-D) Amp3'
    real_T m_bpLambda_d[3];            // '<S43>/Lookup Table (n-D) Amp4'
    real_T m_bpLambda_n[3];            // '<S43>/Lookup Table (n-D) Amp5'
    real_T m_bpLambda_d1[3];           // '<S43>/Lookup Table (n-D) Amp6'
    real_T m_bpLambda_j[3];            // '<S43>/Lookup Table (n-D) Phase1'
    real_T m_bpLambda_l[3];            // '<S43>/Lookup Table (n-D) Phase2'
    real_T m_bpLambda_dj[3];           // '<S43>/Lookup Table (n-D) Phase'
    real_T m_bpLambda_e[3];            // '<S43>/Lookup Table (n-D) Phase4'
    real_T m_bpLambda_j3[3];           // '<S43>/Lookup Table (n-D) Phase5'
    real_T m_bpLambda_g[3];            // '<S43>/Lookup Table (n-D) Phase6'
    uint32_T Prelookup2_DWORK1;        // '<S42>/Prelookup2'
    uint32_T m_bpIndex[3];           // '<S44>/Lookup Table (n-D)  Wavedrift 1'
    uint32_T m_bpIndex_n[3];         // '<S44>/Lookup Table (n-D)  Wavedrift 2'
    uint32_T m_bpIndex_l[3];         // '<S44>/Lookup Table (n-D)  Wavedrift 3'
    uint32_T m_bpIndex_h[3];           // '<S43>/Lookup Table (n-D) Amp1'
    uint32_T m_bpIndex_b[3];           // '<S43>/Lookup Table (n-D) Amp2'
    uint32_T m_bpIndex_hm[3];          // '<S43>/Lookup Table (n-D) Amp3'
    uint32_T m_bpIndex_c[3];           // '<S43>/Lookup Table (n-D) Amp4'
    uint32_T m_bpIndex_ni[3];          // '<S43>/Lookup Table (n-D) Amp5'
    uint32_T m_bpIndex_nf[3];          // '<S43>/Lookup Table (n-D) Amp6'
    uint32_T m_bpIndex_nq[3];          // '<S43>/Lookup Table (n-D) Phase1'
    uint32_T m_bpIndex_cl[3];          // '<S43>/Lookup Table (n-D) Phase2'
    uint32_T m_bpIndex_j[3];           // '<S43>/Lookup Table (n-D) Phase'
    uint32_T m_bpIndex_a[3];           // '<S43>/Lookup Table (n-D) Phase4'
    uint32_T m_bpIndex_d[3];           // '<S43>/Lookup Table (n-D) Phase5'
    uint32_T m_bpIndex_f[3];           // '<S43>/Lookup Table (n-D) Phase6'
  } DW_Wave_loads_for_heading1_AH_T;

  // Block signals for system '<S8>/Wave loads (U=0)'
  typedef struct {
    real_T tau_WD[6];                  // '<S49>/Sum1'
    real_T tau_WF[6];                  // '<S49>/Sum4'
    real_T tau_WD_k[6];                // '<S42>/Sum1'
    real_T tau_WF_j[6];                // '<S42>/Sum4'
    B_Wave_loads_for_heading1_AHV_T Wave_loads_for_heading2;// '<S15>/Wave_loads_for_heading2' 
    B_Wave_loads_for_heading1_AHV_T Wave_loads_for_heading1;// '<S15>/Wave_loads_for_heading1' 
  } B_Wave_loads_fun_T;

  // Block states (default storage) for system '<S8>/Wave loads (U=0)'
  typedef struct {
    DW_Wave_loads_for_heading1_AH_T Wave_loads_for_heading2;// '<S15>/Wave_loads_for_heading2' 
    DW_Wave_loads_for_heading1_AH_T Wave_loads_for_heading1;// '<S15>/Wave_loads_for_heading1' 
  } DW_Wave_loads_fun_T;

  // Block states (default storage) for system '<S64>/Chart'
  typedef struct {
    real_T x_local;                    // '<S64>/Chart'
    real_T y_local;                    // '<S64>/Chart'
    uint8_T is_active_c129_AHV_Model;  // '<S64>/Chart'
    uint8_T is_c129_AHV_Model;         // '<S64>/Chart'
  } DW_Chart_AHV_Model_T;

  // Block signals (default storage)
  typedef struct {
    real_T nu_r[6];                    // '<S13>/Sum6'
    real_T dx1;                        // '<S19>/dx1'
    real_T dx2;                        // '<S19>/dx2'
    real_T RateLimiter;                // '<S64>/Rate Limiter'
    real_T Merge;                      // '<S83>/Merge'
    real_T Merge_m;                    // '<S84>/Merge'
    real_T Kp[3];                      // '<S76>/Merge'
    real_T Minvtau[6];                 // '<S13>/Minv*tau'
    real_T TmpSignalConversionAtIntegrat_o[6];// '<S13>/6 DOF transformation'
    real_T Sum[5];                     // '<S32>/Sum'
    real_T Fcn;                        // '<S16>/Fcn'
    real_T M_u[3];                     // '<S63>/Gain3'
    real_T sun_k2[3];                  // '<S63>/Sum3'
    real_T psi_WF[3];                  // '<S63>/Sum5'
    real_T Sum6[3];                    // '<S63>/Sum6'
    real_T Sum7[3];                    // '<S63>/Sum7'
    real_T TmpSignalConversionAtIntegrat_p[3];// '<S65>/Subsystem'
    real_T nu_r_b[6];                  // '<S97>/Sum6'
    real_T dx1_g;                      // '<S103>/dx1'
    real_T dx2_o;                      // '<S103>/dx2'
    real_T RateLimiter_j;              // '<S148>/Rate Limiter'
    real_T Row3_j;                     // '<S162>/Row3'
    real_T Merge_a;                    // '<S167>/Merge'
    real_T Merge_g;                    // '<S168>/Merge'
    real_T Ki[3];                      // '<S149>/Ki'
    real_T Minvtau_n[6];               // '<S97>/Minv*tau'
    real_T TmpSignalConversionAtIntegra_o3[6];// '<S97>/6 DOF transformation'
    real_T Sum_j[5];                   // '<S116>/Sum'
    real_T Fcn_j;                      // '<S100>/Fcn'
    real_T M_u_f[3];                   // '<S147>/Gain3'
    real_T sun_k2_d[3];                // '<S147>/Sum3'
    real_T psi_WF_g[3];                // '<S147>/Sum5'
    real_T Sum6_i[3];                  // '<S147>/Sum6'
    real_T Sum7_k[3];                  // '<S147>/Sum7'
    real_T nu_r_d[6];                  // '<S181>/Sum6'
    real_T dx1_h;                      // '<S187>/dx1'
    real_T dx2_oc;                     // '<S187>/dx2'
    real_T RateLimiter_l;              // '<S232>/Rate Limiter'
    real_T Row3_o;                     // '<S246>/Row3'
    real_T Merge_p;                    // '<S251>/Merge'
    real_T Merge_i;                    // '<S252>/Merge'
    real_T Ki_o[3];                    // '<S233>/Ki'
    real_T Minvtau_a[6];               // '<S181>/Minv*tau'
    real_T TmpSignalConversionAtIntegr_o3s[6];// '<S181>/6 DOF transformation'
    real_T Sum_e[5];                   // '<S200>/Sum'
    real_T Fcn_d;                      // '<S184>/Fcn'
    real_T M_u_i[3];                   // '<S231>/Gain3'
    real_T sun_k2_n[3];                // '<S231>/Sum3'
    real_T psi_WF_c[3];                // '<S231>/Sum5'
    real_T Sum6_b[3];                  // '<S231>/Sum6'
    real_T Sum7_e[3];                  // '<S231>/Sum7'
    real_T nu_r_f[6];                  // '<S265>/Sum6'
    real_T dx1_f;                      // '<S271>/dx1'
    real_T dx2_j;                      // '<S271>/dx2'
    real_T RateLimiter_lx;             // '<S316>/Rate Limiter'
    real_T Row3_l;                     // '<S330>/Row3'
    real_T Merge_k;                    // '<S335>/Merge'
    real_T Merge_d;                    // '<S336>/Merge'
    real_T Ki_f[3];                    // '<S317>/Ki'
    real_T Minvtau_l[6];               // '<S265>/Minv*tau'
    real_T TmpSignalConversionAtInteg_o3s0[6];// '<S265>/6 DOF transformation'
    real_T Sum_c[5];                   // '<S284>/Sum'
    real_T Fcn_h;                      // '<S268>/Fcn'
    real_T M_u_j[3];                   // '<S315>/Gain3'
    real_T sun_k2_k[3];                // '<S315>/Sum3'
    real_T psi_WF_d[3];                // '<S315>/Sum5'
    real_T Sum6_m[3];                  // '<S315>/Sum6'
    real_T Sum7_g[3];                  // '<S315>/Sum7'
    real_T x_ref_rel;                  // '<S316>/Chart'
    real_T y_ref_rel;                  // '<S316>/Chart'
    real_T Zeta_a[900];                // '<S266>/Wave'
    real_T Omega[900];                 // '<S266>/Wave'
    real_T Phase[900];                 // '<S266>/Wave'
    real_T Wavenum[900];               // '<S266>/Wave'
    real_T Psi[900];                   // '<S266>/Wave'
    real_T Sum2;                       // '<S275>/Sum2'
    real_T Sum_f;                      // '<S275>/Sum'
    real_T x_ref_rel_p;                // '<S232>/Chart'
    real_T y_ref_rel_g;                // '<S232>/Chart'
    real_T Zeta_a_o[900];              // '<S182>/Wave'
    real_T Omega_o[900];               // '<S182>/Wave'
    real_T Phase_k[900];               // '<S182>/Wave'
    real_T Wavenum_e[900];             // '<S182>/Wave'
    real_T Psi_e[900];                 // '<S182>/Wave'
    real_T Sum2_g;                     // '<S191>/Sum2'
    real_T Sum_l;                      // '<S191>/Sum'
    real_T x_ref_rel_e;                // '<S148>/Chart'
    real_T y_ref_rel_n;                // '<S148>/Chart'
    real_T Zeta_a_d[900];              // '<S98>/Wave'
    real_T Omega_e[900];               // '<S98>/Wave'
    real_T Phase_j[900];               // '<S98>/Wave'
    real_T Wavenum_m[900];             // '<S98>/Wave'
    real_T Psi_a[900];                 // '<S98>/Wave'
    real_T Sum2_p;                     // '<S107>/Sum2'
    real_T Sum_k;                      // '<S107>/Sum'
    real_T x_ref_rel_eg;               // '<S64>/Chart'
    real_T y_ref_rel_no;               // '<S64>/Chart'
    real_T Zeta_a_k[900];              // '<S14>/Wave'
    real_T Omega_b[900];               // '<S14>/Wave'
    real_T Phase_kf[900];              // '<S14>/Wave'
    real_T Wavenum_d[900];             // '<S14>/Wave'
    real_T Psi_l[900];                 // '<S14>/Wave'
    real_T Sum2_pj;                    // '<S23>/Sum2'
    real_T Sum_n;                      // '<S23>/Sum'
    B_Wave_loads_fun_T WaveloadsU0_e;  // '<S260>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_m;        // '<S260>/Subsystem3'
    B_Wave_loads_fun_T WaveloadsU0_d;  // '<S176>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_d;        // '<S176>/Subsystem3'
    B_Wave_loads_fun_T WaveloadsU0_b;  // '<S92>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_b;        // '<S92>/Subsystem3'
    B_Wave_loads_fun_T WaveloadsU0;    // '<S8>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3;          // '<S8>/Subsystem3'
  } B_AHV_Model_T;

  // Block states (default storage) for system '<Root>'
  typedef struct {
    real_T Integrator1_DSTATE[3];      // '<S149>/Integrator1'
    real_T Integrator1_DSTATE_i[3];    // '<S233>/Integrator1'
    real_T Integrator1_DSTATE_in[3];   // '<S317>/Integrator1'
    real_T PrevY;                      // '<S64>/Rate Limiter'
    real_T LastMajorTime;              // '<S64>/Rate Limiter'
    real_T PrevY_p;                    // '<S148>/Rate Limiter'
    real_T LastMajorTime_n;            // '<S148>/Rate Limiter'
    real_T PrevY_l;                    // '<S232>/Rate Limiter'
    real_T LastMajorTime_a;            // '<S232>/Rate Limiter'
    real_T PrevY_d;                    // '<S316>/Rate Limiter'
    real_T LastMajorTime_m;            // '<S316>/Rate Limiter'
    int_T Integrator1_IWORK;           // '<S13>/Integrator1'
    int_T Integrator3_IWORK;           // '<S63>/Integrator3'
    int_T Integrator1_IWORK_a;         // '<S97>/Integrator1'
    int_T Integrator3_IWORK_p;         // '<S147>/Integrator3'
    int_T Integrator1_IWORK_p;         // '<S181>/Integrator1'
    int_T Integrator3_IWORK_n;         // '<S231>/Integrator3'
    int_T Integrator1_IWORK_f;         // '<S265>/Integrator1'
    int_T Integrator3_IWORK_c;         // '<S315>/Integrator3'
    int8_T If_ActiveSubsystem;         // '<S83>/If'
    int8_T If_ActiveSubsystem_a;       // '<S84>/If'
    int8_T If_ActiveSubsystem_i;       // '<S76>/If'
    int8_T If_ActiveSubsystem_h;       // '<S167>/If'
    int8_T If_ActiveSubsystem_ig;      // '<S168>/If'
    int8_T If_ActiveSubsystem_ab;      // '<S160>/If'
    int8_T If_ActiveSubsystem_e;       // '<S251>/If'
    int8_T If_ActiveSubsystem_k;       // '<S252>/If'
    int8_T If_ActiveSubsystem_ef;      // '<S244>/If'
    int8_T If_ActiveSubsystem_f;       // '<S335>/If'
    int8_T If_ActiveSubsystem_hp;      // '<S336>/If'
    int8_T If_ActiveSubsystem_i4;      // '<S328>/If'
    DW_Chart_AHV_Model_T sf_Chart_j;   // '<S316>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_e; // '<S260>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_m;       // '<S260>/Subsystem3'
    DW_Chart_AHV_Model_T sf_Chart_k;   // '<S232>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_d; // '<S176>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_d;       // '<S176>/Subsystem3'
    DW_Chart_AHV_Model_T sf_Chart_o;   // '<S148>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_b; // '<S92>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_b;       // '<S92>/Subsystem3'
    DW_Chart_AHV_Model_T sf_Chart;     // '<S64>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0;   // '<S8>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3;         // '<S8>/Subsystem3'
  } DW_AHV_Model_T;

  // Continuous states (default storage)
  typedef struct {
    real_T Integrator1_CSTATE[6];      // '<S13>/Integrator1'
    real_T Integrator_CSTATE[6];       // '<S13>/Integrator'
    real_T TransferFcn_CSTATE;         // '<S16>/Transfer Fcn'
    real_T Integrator3_CSTATE[3];      // '<S63>/Integrator3'
    real_T Integrator_CSTATE_b[3];     // '<S65>/Integrator'
    real_T Integrator4_CSTATE[3];      // '<S63>/Integrator4'
    real_T Dp11_CSTATE[5];             // '<S20>/Dp(1,1)'
    real_T Dp13_CSTATE[5];             // '<S20>/Dp(1,3)'
    real_T Dp15_CSTATE[5];             // '<S20>/Dp(1,5)'
    real_T Dp22_CSTATE[5];             // '<S20>/Dp(2,2)'
    real_T Dp24_CSTATE[5];             // '<S20>/Dp(2,4)'
    real_T Dp26_CSTATE[5];             // '<S20>/Dp(2,6)'
    real_T Dp31_CSTATE[5];             // '<S20>/Dp(3,1)'
    real_T Dp33_CSTATE[5];             // '<S20>/Dp(3,3)'
    real_T Dp35_CSTATE[5];             // '<S20>/Dp(3,5)'
    real_T Dp42_CSTATE[5];             // '<S20>/Dp(4,2)'
    real_T Integrator_CSTATE_l[5];     // '<S32>/Integrator'
    real_T Dp46_CSTATE[5];             // '<S20>/Dp(4,6)'
    real_T Dp51_CSTATE[5];             // '<S20>/Dp(5,1)'
    real_T Dp53_CSTATE[5];             // '<S20>/Dp(5,3)'
    real_T Dp55_CSTATE[5];             // '<S20>/Dp(5,5)'
    real_T Dp62_CSTATE[5];             // '<S20>/Dp(6,2)'
    real_T Dp64_CSTATE[5];             // '<S20>/Dp(6,4)'
    real_T Dp66_CSTATE[5];             // '<S20>/Dp(6,6)'
    real_T Integrator1_CSTATE_c[3];    // '<S63>/Integrator1'
    real_T Integrator2_CSTATE[3];      // '<S63>/Integrator2'
    real_T Integrator6_CSTATE[3];      // '<S63>/Integrator6'
    real_T Integrator1_CSTATE_l[6];    // '<S97>/Integrator1'
    real_T Integrator_CSTATE_n[6];     // '<S97>/Integrator'
    real_T TransferFcn_CSTATE_e;       // '<S100>/Transfer Fcn'
    real_T Integrator3_CSTATE_i[3];    // '<S147>/Integrator3'
    real_T Integrator4_CSTATE_p[3];    // '<S147>/Integrator4'
    real_T Dp11_CSTATE_k[5];           // '<S104>/Dp(1,1)'
    real_T Dp13_CSTATE_k[5];           // '<S104>/Dp(1,3)'
    real_T Dp15_CSTATE_j[5];           // '<S104>/Dp(1,5)'
    real_T Dp22_CSTATE_l[5];           // '<S104>/Dp(2,2)'
    real_T Dp24_CSTATE_i[5];           // '<S104>/Dp(2,4)'
    real_T Dp26_CSTATE_p[5];           // '<S104>/Dp(2,6)'
    real_T Dp31_CSTATE_l[5];           // '<S104>/Dp(3,1)'
    real_T Dp33_CSTATE_m[5];           // '<S104>/Dp(3,3)'
    real_T Dp35_CSTATE_n[5];           // '<S104>/Dp(3,5)'
    real_T Dp42_CSTATE_c[5];           // '<S104>/Dp(4,2)'
    real_T Integrator_CSTATE_m[5];     // '<S116>/Integrator'
    real_T Dp46_CSTATE_o[5];           // '<S104>/Dp(4,6)'
    real_T Dp51_CSTATE_c[5];           // '<S104>/Dp(5,1)'
    real_T Dp53_CSTATE_p[5];           // '<S104>/Dp(5,3)'
    real_T Dp55_CSTATE_e[5];           // '<S104>/Dp(5,5)'
    real_T Dp62_CSTATE_k[5];           // '<S104>/Dp(6,2)'
    real_T Dp64_CSTATE_l[5];           // '<S104>/Dp(6,4)'
    real_T Dp66_CSTATE_g[5];           // '<S104>/Dp(6,6)'
    real_T Integrator1_CSTATE_d[3];    // '<S147>/Integrator1'
    real_T Integrator2_CSTATE_i[3];    // '<S147>/Integrator2'
    real_T Integrator6_CSTATE_b[3];    // '<S147>/Integrator6'
    real_T Integrator1_CSTATE_lj[6];   // '<S181>/Integrator1'
    real_T Integrator_CSTATE_o[6];     // '<S181>/Integrator'
    real_T TransferFcn_CSTATE_l;       // '<S184>/Transfer Fcn'
    real_T Integrator3_CSTATE_g[3];    // '<S231>/Integrator3'
    real_T Integrator4_CSTATE_d[3];    // '<S231>/Integrator4'
    real_T Dp11_CSTATE_j[5];           // '<S188>/Dp(1,1)'
    real_T Dp13_CSTATE_a[5];           // '<S188>/Dp(1,3)'
    real_T Dp15_CSTATE_p[5];           // '<S188>/Dp(1,5)'
    real_T Dp22_CSTATE_g[5];           // '<S188>/Dp(2,2)'
    real_T Dp24_CSTATE_f[5];           // '<S188>/Dp(2,4)'
    real_T Dp26_CSTATE_i[5];           // '<S188>/Dp(2,6)'
    real_T Dp31_CSTATE_e[5];           // '<S188>/Dp(3,1)'
    real_T Dp33_CSTATE_o[5];           // '<S188>/Dp(3,3)'
    real_T Dp35_CSTATE_i[5];           // '<S188>/Dp(3,5)'
    real_T Dp42_CSTATE_f[5];           // '<S188>/Dp(4,2)'
    real_T Integrator_CSTATE_h[5];     // '<S200>/Integrator'
    real_T Dp46_CSTATE_h[5];           // '<S188>/Dp(4,6)'
    real_T Dp51_CSTATE_f[5];           // '<S188>/Dp(5,1)'
    real_T Dp53_CSTATE_c[5];           // '<S188>/Dp(5,3)'
    real_T Dp55_CSTATE_f[5];           // '<S188>/Dp(5,5)'
    real_T Dp62_CSTATE_n[5];           // '<S188>/Dp(6,2)'
    real_T Dp64_CSTATE_b[5];           // '<S188>/Dp(6,4)'
    real_T Dp66_CSTATE_n[5];           // '<S188>/Dp(6,6)'
    real_T Integrator1_CSTATE_m[3];    // '<S231>/Integrator1'
    real_T Integrator2_CSTATE_j[3];    // '<S231>/Integrator2'
    real_T Integrator6_CSTATE_e[3];    // '<S231>/Integrator6'
    real_T Integrator1_CSTATE_b[6];    // '<S265>/Integrator1'
    real_T Integrator_CSTATE_a[6];     // '<S265>/Integrator'
    real_T TransferFcn_CSTATE_m;       // '<S268>/Transfer Fcn'
    real_T Integrator3_CSTATE_d[3];    // '<S315>/Integrator3'
    real_T Integrator4_CSTATE_o[3];    // '<S315>/Integrator4'
    real_T Dp11_CSTATE_e[5];           // '<S272>/Dp(1,1)'
    real_T Dp13_CSTATE_c[5];           // '<S272>/Dp(1,3)'
    real_T Dp15_CSTATE_d[5];           // '<S272>/Dp(1,5)'
    real_T Dp22_CSTATE_a[5];           // '<S272>/Dp(2,2)'
    real_T Dp24_CSTATE_d[5];           // '<S272>/Dp(2,4)'
    real_T Dp26_CSTATE_f[5];           // '<S272>/Dp(2,6)'
    real_T Dp31_CSTATE_a[5];           // '<S272>/Dp(3,1)'
    real_T Dp33_CSTATE_a[5];           // '<S272>/Dp(3,3)'
    real_T Dp35_CSTATE_l[5];           // '<S272>/Dp(3,5)'
    real_T Dp42_CSTATE_h[5];           // '<S272>/Dp(4,2)'
    real_T Integrator_CSTATE_bp[5];    // '<S284>/Integrator'
    real_T Dp46_CSTATE_i[5];           // '<S272>/Dp(4,6)'
    real_T Dp51_CSTATE_k[5];           // '<S272>/Dp(5,1)'
    real_T Dp53_CSTATE_pz[5];          // '<S272>/Dp(5,3)'
    real_T Dp55_CSTATE_a[5];           // '<S272>/Dp(5,5)'
    real_T Dp62_CSTATE_i[5];           // '<S272>/Dp(6,2)'
    real_T Dp64_CSTATE_e[5];           // '<S272>/Dp(6,4)'
    real_T Dp66_CSTATE_g3[5];          // '<S272>/Dp(6,6)'
    real_T Integrator1_CSTATE_me[3];   // '<S315>/Integrator1'
    real_T Integrator2_CSTATE_h[3];    // '<S315>/Integrator2'
    real_T Integrator6_CSTATE_f[3];    // '<S315>/Integrator6'
  } X_AHV_Model_T;

  // State derivatives (default storage)
  typedef struct {
    real_T Integrator1_CSTATE[6];      // '<S13>/Integrator1'
    real_T Integrator_CSTATE[6];       // '<S13>/Integrator'
    real_T TransferFcn_CSTATE;         // '<S16>/Transfer Fcn'
    real_T Integrator3_CSTATE[3];      // '<S63>/Integrator3'
    real_T Integrator_CSTATE_b[3];     // '<S65>/Integrator'
    real_T Integrator4_CSTATE[3];      // '<S63>/Integrator4'
    real_T Dp11_CSTATE[5];             // '<S20>/Dp(1,1)'
    real_T Dp13_CSTATE[5];             // '<S20>/Dp(1,3)'
    real_T Dp15_CSTATE[5];             // '<S20>/Dp(1,5)'
    real_T Dp22_CSTATE[5];             // '<S20>/Dp(2,2)'
    real_T Dp24_CSTATE[5];             // '<S20>/Dp(2,4)'
    real_T Dp26_CSTATE[5];             // '<S20>/Dp(2,6)'
    real_T Dp31_CSTATE[5];             // '<S20>/Dp(3,1)'
    real_T Dp33_CSTATE[5];             // '<S20>/Dp(3,3)'
    real_T Dp35_CSTATE[5];             // '<S20>/Dp(3,5)'
    real_T Dp42_CSTATE[5];             // '<S20>/Dp(4,2)'
    real_T Integrator_CSTATE_l[5];     // '<S32>/Integrator'
    real_T Dp46_CSTATE[5];             // '<S20>/Dp(4,6)'
    real_T Dp51_CSTATE[5];             // '<S20>/Dp(5,1)'
    real_T Dp53_CSTATE[5];             // '<S20>/Dp(5,3)'
    real_T Dp55_CSTATE[5];             // '<S20>/Dp(5,5)'
    real_T Dp62_CSTATE[5];             // '<S20>/Dp(6,2)'
    real_T Dp64_CSTATE[5];             // '<S20>/Dp(6,4)'
    real_T Dp66_CSTATE[5];             // '<S20>/Dp(6,6)'
    real_T Integrator1_CSTATE_c[3];    // '<S63>/Integrator1'
    real_T Integrator2_CSTATE[3];      // '<S63>/Integrator2'
    real_T Integrator6_CSTATE[3];      // '<S63>/Integrator6'
    real_T Integrator1_CSTATE_l[6];    // '<S97>/Integrator1'
    real_T Integrator_CSTATE_n[6];     // '<S97>/Integrator'
    real_T TransferFcn_CSTATE_e;       // '<S100>/Transfer Fcn'
    real_T Integrator3_CSTATE_i[3];    // '<S147>/Integrator3'
    real_T Integrator4_CSTATE_p[3];    // '<S147>/Integrator4'
    real_T Dp11_CSTATE_k[5];           // '<S104>/Dp(1,1)'
    real_T Dp13_CSTATE_k[5];           // '<S104>/Dp(1,3)'
    real_T Dp15_CSTATE_j[5];           // '<S104>/Dp(1,5)'
    real_T Dp22_CSTATE_l[5];           // '<S104>/Dp(2,2)'
    real_T Dp24_CSTATE_i[5];           // '<S104>/Dp(2,4)'
    real_T Dp26_CSTATE_p[5];           // '<S104>/Dp(2,6)'
    real_T Dp31_CSTATE_l[5];           // '<S104>/Dp(3,1)'
    real_T Dp33_CSTATE_m[5];           // '<S104>/Dp(3,3)'
    real_T Dp35_CSTATE_n[5];           // '<S104>/Dp(3,5)'
    real_T Dp42_CSTATE_c[5];           // '<S104>/Dp(4,2)'
    real_T Integrator_CSTATE_m[5];     // '<S116>/Integrator'
    real_T Dp46_CSTATE_o[5];           // '<S104>/Dp(4,6)'
    real_T Dp51_CSTATE_c[5];           // '<S104>/Dp(5,1)'
    real_T Dp53_CSTATE_p[5];           // '<S104>/Dp(5,3)'
    real_T Dp55_CSTATE_e[5];           // '<S104>/Dp(5,5)'
    real_T Dp62_CSTATE_k[5];           // '<S104>/Dp(6,2)'
    real_T Dp64_CSTATE_l[5];           // '<S104>/Dp(6,4)'
    real_T Dp66_CSTATE_g[5];           // '<S104>/Dp(6,6)'
    real_T Integrator1_CSTATE_d[3];    // '<S147>/Integrator1'
    real_T Integrator2_CSTATE_i[3];    // '<S147>/Integrator2'
    real_T Integrator6_CSTATE_b[3];    // '<S147>/Integrator6'
    real_T Integrator1_CSTATE_lj[6];   // '<S181>/Integrator1'
    real_T Integrator_CSTATE_o[6];     // '<S181>/Integrator'
    real_T TransferFcn_CSTATE_l;       // '<S184>/Transfer Fcn'
    real_T Integrator3_CSTATE_g[3];    // '<S231>/Integrator3'
    real_T Integrator4_CSTATE_d[3];    // '<S231>/Integrator4'
    real_T Dp11_CSTATE_j[5];           // '<S188>/Dp(1,1)'
    real_T Dp13_CSTATE_a[5];           // '<S188>/Dp(1,3)'
    real_T Dp15_CSTATE_p[5];           // '<S188>/Dp(1,5)'
    real_T Dp22_CSTATE_g[5];           // '<S188>/Dp(2,2)'
    real_T Dp24_CSTATE_f[5];           // '<S188>/Dp(2,4)'
    real_T Dp26_CSTATE_i[5];           // '<S188>/Dp(2,6)'
    real_T Dp31_CSTATE_e[5];           // '<S188>/Dp(3,1)'
    real_T Dp33_CSTATE_o[5];           // '<S188>/Dp(3,3)'
    real_T Dp35_CSTATE_i[5];           // '<S188>/Dp(3,5)'
    real_T Dp42_CSTATE_f[5];           // '<S188>/Dp(4,2)'
    real_T Integrator_CSTATE_h[5];     // '<S200>/Integrator'
    real_T Dp46_CSTATE_h[5];           // '<S188>/Dp(4,6)'
    real_T Dp51_CSTATE_f[5];           // '<S188>/Dp(5,1)'
    real_T Dp53_CSTATE_c[5];           // '<S188>/Dp(5,3)'
    real_T Dp55_CSTATE_f[5];           // '<S188>/Dp(5,5)'
    real_T Dp62_CSTATE_n[5];           // '<S188>/Dp(6,2)'
    real_T Dp64_CSTATE_b[5];           // '<S188>/Dp(6,4)'
    real_T Dp66_CSTATE_n[5];           // '<S188>/Dp(6,6)'
    real_T Integrator1_CSTATE_m[3];    // '<S231>/Integrator1'
    real_T Integrator2_CSTATE_j[3];    // '<S231>/Integrator2'
    real_T Integrator6_CSTATE_e[3];    // '<S231>/Integrator6'
    real_T Integrator1_CSTATE_b[6];    // '<S265>/Integrator1'
    real_T Integrator_CSTATE_a[6];     // '<S265>/Integrator'
    real_T TransferFcn_CSTATE_m;       // '<S268>/Transfer Fcn'
    real_T Integrator3_CSTATE_d[3];    // '<S315>/Integrator3'
    real_T Integrator4_CSTATE_o[3];    // '<S315>/Integrator4'
    real_T Dp11_CSTATE_e[5];           // '<S272>/Dp(1,1)'
    real_T Dp13_CSTATE_c[5];           // '<S272>/Dp(1,3)'
    real_T Dp15_CSTATE_d[5];           // '<S272>/Dp(1,5)'
    real_T Dp22_CSTATE_a[5];           // '<S272>/Dp(2,2)'
    real_T Dp24_CSTATE_d[5];           // '<S272>/Dp(2,4)'
    real_T Dp26_CSTATE_f[5];           // '<S272>/Dp(2,6)'
    real_T Dp31_CSTATE_a[5];           // '<S272>/Dp(3,1)'
    real_T Dp33_CSTATE_a[5];           // '<S272>/Dp(3,3)'
    real_T Dp35_CSTATE_l[5];           // '<S272>/Dp(3,5)'
    real_T Dp42_CSTATE_h[5];           // '<S272>/Dp(4,2)'
    real_T Integrator_CSTATE_bp[5];    // '<S284>/Integrator'
    real_T Dp46_CSTATE_i[5];           // '<S272>/Dp(4,6)'
    real_T Dp51_CSTATE_k[5];           // '<S272>/Dp(5,1)'
    real_T Dp53_CSTATE_pz[5];          // '<S272>/Dp(5,3)'
    real_T Dp55_CSTATE_a[5];           // '<S272>/Dp(5,5)'
    real_T Dp62_CSTATE_i[5];           // '<S272>/Dp(6,2)'
    real_T Dp64_CSTATE_e[5];           // '<S272>/Dp(6,4)'
    real_T Dp66_CSTATE_g3[5];          // '<S272>/Dp(6,6)'
    real_T Integrator1_CSTATE_me[3];   // '<S315>/Integrator1'
    real_T Integrator2_CSTATE_h[3];    // '<S315>/Integrator2'
    real_T Integrator6_CSTATE_f[3];    // '<S315>/Integrator6'
  } XDot_AHV_Model_T;

  // State disabled
  typedef struct {
    boolean_T Integrator1_CSTATE[6];   // '<S13>/Integrator1'
    boolean_T Integrator_CSTATE[6];    // '<S13>/Integrator'
    boolean_T TransferFcn_CSTATE;      // '<S16>/Transfer Fcn'
    boolean_T Integrator3_CSTATE[3];   // '<S63>/Integrator3'
    boolean_T Integrator_CSTATE_b[3];  // '<S65>/Integrator'
    boolean_T Integrator4_CSTATE[3];   // '<S63>/Integrator4'
    boolean_T Dp11_CSTATE[5];          // '<S20>/Dp(1,1)'
    boolean_T Dp13_CSTATE[5];          // '<S20>/Dp(1,3)'
    boolean_T Dp15_CSTATE[5];          // '<S20>/Dp(1,5)'
    boolean_T Dp22_CSTATE[5];          // '<S20>/Dp(2,2)'
    boolean_T Dp24_CSTATE[5];          // '<S20>/Dp(2,4)'
    boolean_T Dp26_CSTATE[5];          // '<S20>/Dp(2,6)'
    boolean_T Dp31_CSTATE[5];          // '<S20>/Dp(3,1)'
    boolean_T Dp33_CSTATE[5];          // '<S20>/Dp(3,3)'
    boolean_T Dp35_CSTATE[5];          // '<S20>/Dp(3,5)'
    boolean_T Dp42_CSTATE[5];          // '<S20>/Dp(4,2)'
    boolean_T Integrator_CSTATE_l[5];  // '<S32>/Integrator'
    boolean_T Dp46_CSTATE[5];          // '<S20>/Dp(4,6)'
    boolean_T Dp51_CSTATE[5];          // '<S20>/Dp(5,1)'
    boolean_T Dp53_CSTATE[5];          // '<S20>/Dp(5,3)'
    boolean_T Dp55_CSTATE[5];          // '<S20>/Dp(5,5)'
    boolean_T Dp62_CSTATE[5];          // '<S20>/Dp(6,2)'
    boolean_T Dp64_CSTATE[5];          // '<S20>/Dp(6,4)'
    boolean_T Dp66_CSTATE[5];          // '<S20>/Dp(6,6)'
    boolean_T Integrator1_CSTATE_c[3]; // '<S63>/Integrator1'
    boolean_T Integrator2_CSTATE[3];   // '<S63>/Integrator2'
    boolean_T Integrator6_CSTATE[3];   // '<S63>/Integrator6'
    boolean_T Integrator1_CSTATE_l[6]; // '<S97>/Integrator1'
    boolean_T Integrator_CSTATE_n[6];  // '<S97>/Integrator'
    boolean_T TransferFcn_CSTATE_e;    // '<S100>/Transfer Fcn'
    boolean_T Integrator3_CSTATE_i[3]; // '<S147>/Integrator3'
    boolean_T Integrator4_CSTATE_p[3]; // '<S147>/Integrator4'
    boolean_T Dp11_CSTATE_k[5];        // '<S104>/Dp(1,1)'
    boolean_T Dp13_CSTATE_k[5];        // '<S104>/Dp(1,3)'
    boolean_T Dp15_CSTATE_j[5];        // '<S104>/Dp(1,5)'
    boolean_T Dp22_CSTATE_l[5];        // '<S104>/Dp(2,2)'
    boolean_T Dp24_CSTATE_i[5];        // '<S104>/Dp(2,4)'
    boolean_T Dp26_CSTATE_p[5];        // '<S104>/Dp(2,6)'
    boolean_T Dp31_CSTATE_l[5];        // '<S104>/Dp(3,1)'
    boolean_T Dp33_CSTATE_m[5];        // '<S104>/Dp(3,3)'
    boolean_T Dp35_CSTATE_n[5];        // '<S104>/Dp(3,5)'
    boolean_T Dp42_CSTATE_c[5];        // '<S104>/Dp(4,2)'
    boolean_T Integrator_CSTATE_m[5];  // '<S116>/Integrator'
    boolean_T Dp46_CSTATE_o[5];        // '<S104>/Dp(4,6)'
    boolean_T Dp51_CSTATE_c[5];        // '<S104>/Dp(5,1)'
    boolean_T Dp53_CSTATE_p[5];        // '<S104>/Dp(5,3)'
    boolean_T Dp55_CSTATE_e[5];        // '<S104>/Dp(5,5)'
    boolean_T Dp62_CSTATE_k[5];        // '<S104>/Dp(6,2)'
    boolean_T Dp64_CSTATE_l[5];        // '<S104>/Dp(6,4)'
    boolean_T Dp66_CSTATE_g[5];        // '<S104>/Dp(6,6)'
    boolean_T Integrator1_CSTATE_d[3]; // '<S147>/Integrator1'
    boolean_T Integrator2_CSTATE_i[3]; // '<S147>/Integrator2'
    boolean_T Integrator6_CSTATE_b[3]; // '<S147>/Integrator6'
    boolean_T Integrator1_CSTATE_lj[6];// '<S181>/Integrator1'
    boolean_T Integrator_CSTATE_o[6];  // '<S181>/Integrator'
    boolean_T TransferFcn_CSTATE_l;    // '<S184>/Transfer Fcn'
    boolean_T Integrator3_CSTATE_g[3]; // '<S231>/Integrator3'
    boolean_T Integrator4_CSTATE_d[3]; // '<S231>/Integrator4'
    boolean_T Dp11_CSTATE_j[5];        // '<S188>/Dp(1,1)'
    boolean_T Dp13_CSTATE_a[5];        // '<S188>/Dp(1,3)'
    boolean_T Dp15_CSTATE_p[5];        // '<S188>/Dp(1,5)'
    boolean_T Dp22_CSTATE_g[5];        // '<S188>/Dp(2,2)'
    boolean_T Dp24_CSTATE_f[5];        // '<S188>/Dp(2,4)'
    boolean_T Dp26_CSTATE_i[5];        // '<S188>/Dp(2,6)'
    boolean_T Dp31_CSTATE_e[5];        // '<S188>/Dp(3,1)'
    boolean_T Dp33_CSTATE_o[5];        // '<S188>/Dp(3,3)'
    boolean_T Dp35_CSTATE_i[5];        // '<S188>/Dp(3,5)'
    boolean_T Dp42_CSTATE_f[5];        // '<S188>/Dp(4,2)'
    boolean_T Integrator_CSTATE_h[5];  // '<S200>/Integrator'
    boolean_T Dp46_CSTATE_h[5];        // '<S188>/Dp(4,6)'
    boolean_T Dp51_CSTATE_f[5];        // '<S188>/Dp(5,1)'
    boolean_T Dp53_CSTATE_c[5];        // '<S188>/Dp(5,3)'
    boolean_T Dp55_CSTATE_f[5];        // '<S188>/Dp(5,5)'
    boolean_T Dp62_CSTATE_n[5];        // '<S188>/Dp(6,2)'
    boolean_T Dp64_CSTATE_b[5];        // '<S188>/Dp(6,4)'
    boolean_T Dp66_CSTATE_n[5];        // '<S188>/Dp(6,6)'
    boolean_T Integrator1_CSTATE_m[3]; // '<S231>/Integrator1'
    boolean_T Integrator2_CSTATE_j[3]; // '<S231>/Integrator2'
    boolean_T Integrator6_CSTATE_e[3]; // '<S231>/Integrator6'
    boolean_T Integrator1_CSTATE_b[6]; // '<S265>/Integrator1'
    boolean_T Integrator_CSTATE_a[6];  // '<S265>/Integrator'
    boolean_T TransferFcn_CSTATE_m;    // '<S268>/Transfer Fcn'
    boolean_T Integrator3_CSTATE_d[3]; // '<S315>/Integrator3'
    boolean_T Integrator4_CSTATE_o[3]; // '<S315>/Integrator4'
    boolean_T Dp11_CSTATE_e[5];        // '<S272>/Dp(1,1)'
    boolean_T Dp13_CSTATE_c[5];        // '<S272>/Dp(1,3)'
    boolean_T Dp15_CSTATE_d[5];        // '<S272>/Dp(1,5)'
    boolean_T Dp22_CSTATE_a[5];        // '<S272>/Dp(2,2)'
    boolean_T Dp24_CSTATE_d[5];        // '<S272>/Dp(2,4)'
    boolean_T Dp26_CSTATE_f[5];        // '<S272>/Dp(2,6)'
    boolean_T Dp31_CSTATE_a[5];        // '<S272>/Dp(3,1)'
    boolean_T Dp33_CSTATE_a[5];        // '<S272>/Dp(3,3)'
    boolean_T Dp35_CSTATE_l[5];        // '<S272>/Dp(3,5)'
    boolean_T Dp42_CSTATE_h[5];        // '<S272>/Dp(4,2)'
    boolean_T Integrator_CSTATE_bp[5]; // '<S284>/Integrator'
    boolean_T Dp46_CSTATE_i[5];        // '<S272>/Dp(4,6)'
    boolean_T Dp51_CSTATE_k[5];        // '<S272>/Dp(5,1)'
    boolean_T Dp53_CSTATE_pz[5];       // '<S272>/Dp(5,3)'
    boolean_T Dp55_CSTATE_a[5];        // '<S272>/Dp(5,5)'
    boolean_T Dp62_CSTATE_i[5];        // '<S272>/Dp(6,2)'
    boolean_T Dp64_CSTATE_e[5];        // '<S272>/Dp(6,4)'
    boolean_T Dp66_CSTATE_g3[5];       // '<S272>/Dp(6,6)'
    boolean_T Integrator1_CSTATE_me[3];// '<S315>/Integrator1'
    boolean_T Integrator2_CSTATE_h[3]; // '<S315>/Integrator2'
    boolean_T Integrator6_CSTATE_f[3]; // '<S315>/Integrator6'
  } XDis_AHV_Model_T;

  // Invariant block signals for system '<S15>/Wave_loads_for_heading1'
  typedef const struct tag_ConstB_Wave_loads_for_hea_T {
    real_T number_of_iterations;       // '<S35>/Width '
    real_T U;                     // '<S42>/Direct Lookup Table (n-D) Velocity'
    uint32_T Prelookup2;               // '<S42>/Prelookup2'
  } ConstB_Wave_loads_for_heading_T;

  // Invariant block signals for system '<S8>/Wave loads (U=0)'
  typedef const struct tag_ConstB_Wave_loads_fun_T {
    ConstB_Wave_loads_for_heading_T Wave_loads_for_heading2;// '<S15>/Wave_loads_for_heading2' 
    ConstB_Wave_loads_for_heading_T Wave_loads_for_heading1;// '<S15>/Wave_loads_for_heading1' 
  } ConstB_Wave_loads_fun_T;

  // Invariant block signals for system '<S83>/If Action Subsystem'
  typedef const struct tag_ConstB_IfActionSubsyste_b_T {
    real_T Gain;                       // '<S85>/Gain'
  } ConstB_IfActionSubsystem_AH_p_T;

  // Invariant block signals (default storage)
  typedef const struct tag_ConstB_AHV_Model_T {
    real_T N1;                         // '<S19>/Sum'
    real_T dx;                         // '<S19>/dx'
    real_T Gain2;                      // '<S19>/Gain2'
    real_T Product5;                   // '<S19>/Product5'
    real_T Sum2[36];                   // '<S13>/Sum2'
    real_T MathFunction[5];            // '<S32>/Math Function'
    real_T N1_p;                       // '<S103>/Sum'
    real_T dx_d;                       // '<S103>/dx'
    real_T Gain2_f;                    // '<S103>/Gain2'
    real_T Product5_o;                 // '<S103>/Product5'
    real_T Sum2_o[36];                 // '<S97>/Sum2'
    real_T MathFunction_i[5];          // '<S116>/Math Function'
    real_T N1_b;                       // '<S187>/Sum'
    real_T dx_b;                       // '<S187>/dx'
    real_T Gain2_k;                    // '<S187>/Gain2'
    real_T Product5_or;                // '<S187>/Product5'
    real_T Sum2_l[36];                 // '<S181>/Sum2'
    real_T MathFunction_d[5];          // '<S200>/Math Function'
    real_T N1_l;                       // '<S271>/Sum'
    real_T dx_h;                       // '<S271>/dx'
    real_T Gain2_l;                    // '<S271>/Gain2'
    real_T Product5_f;                 // '<S271>/Product5'
    real_T Sum2_e[36];                 // '<S265>/Sum2'
    real_T MathFunction_a[5];          // '<S284>/Math Function'
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_i;// '<S336>/If Action Subsystem' 
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_fv;// '<S335>/If Action Subsystem' 
    ConstB_Wave_loads_fun_T WaveloadsU0_e;// '<S260>/Wave loads (U=0)'
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_e;// '<S252>/If Action Subsystem' 
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_c;// '<S251>/If Action Subsystem' 
    ConstB_Wave_loads_fun_T WaveloadsU0_d;// '<S176>/Wave loads (U=0)'
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_f;// '<S168>/If Action Subsystem' 
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_d;// '<S167>/If Action Subsystem' 
    ConstB_Wave_loads_fun_T WaveloadsU0_b;// '<S92>/Wave loads (U=0)'
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_m;// '<S84>/If Action Subsystem' 
    ConstB_IfActionSubsystem_AH_p_T IfActionSubsystem_l;// '<S83>/If Action Subsystem' 
    ConstB_Wave_loads_fun_T WaveloadsU0;// '<S8>/Wave loads (U=0)'
  } ConstB_AHV_Model_T;

  // Constant parameters (default storage)
  typedef struct {
    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S76>/If Action Subsystem1'
    //    '<S160>/If Action Subsystem1'
    //    '<S244>/If Action Subsystem1'
    //    '<S328>/If Action Subsystem1'

    real_T pooled1[9];

    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S76>/If Action Subsystem'
    //    '<S160>/If Action Subsystem'
    //    '<S244>/If Action Subsystem'
    //    '<S328>/If Action Subsystem'

    real_T pooled2[9];

    // Pooled Parameter (Expression: phase_vector)
    //  Referenced by:
    //    '<S14>/Constant'
    //    '<S98>/Constant'
    //    '<S182>/Constant'
    //    '<S266>/Constant'

    real_T pooled7[200];

    // Pooled Parameter (Expression: psi_vector)
    //  Referenced by:
    //    '<S14>/Constant1'
    //    '<S98>/Constant1'
    //    '<S182>/Constant1'
    //    '<S266>/Constant1'

    real_T pooled8[200];

    // Pooled Parameter (Expression: vessel_dp.headings)
    //  Referenced by:
    //    '<S42>/Prelookup1'
    //    '<S42>/Direct Lookup Table (n-D) Beta'
    //    '<S49>/Prelookup1'
    //    '<S49>/Direct Lookup Table (n-D) Beta'
    //    '<S126>/Prelookup1'
    //    '<S126>/Direct Lookup Table (n-D) Beta'
    //    '<S133>/Prelookup1'
    //    '<S133>/Direct Lookup Table (n-D) Beta'
    //    '<S210>/Prelookup1'
    //    '<S210>/Direct Lookup Table (n-D) Beta'
    //    '<S217>/Prelookup1'
    //    '<S217>/Direct Lookup Table (n-D) Beta'
    //    '<S294>/Prelookup1'
    //    '<S294>/Direct Lookup Table (n-D) Beta'
    //    '<S301>/Prelookup1'
    //    '<S301>/Direct Lookup Table (n-D) Beta'
    //    '<S43>/Lookup Table (n-D) Amp1'
    //    '<S43>/Lookup Table (n-D) Amp2'
    //    '<S43>/Lookup Table (n-D) Amp3'
    //    '<S43>/Lookup Table (n-D) Amp4'
    //    '<S43>/Lookup Table (n-D) Amp5'
    //    '<S43>/Lookup Table (n-D) Amp6'
    //    '<S43>/Lookup Table (n-D) Phase'
    //    '<S43>/Lookup Table (n-D) Phase1'
    //    '<S43>/Lookup Table (n-D) Phase2'
    //    '<S43>/Lookup Table (n-D) Phase4'
    //    '<S43>/Lookup Table (n-D) Phase5'
    //    '<S43>/Lookup Table (n-D) Phase6'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S50>/Lookup Table (n-D) Amp1'
    //    '<S50>/Lookup Table (n-D) Amp2'
    //    '<S50>/Lookup Table (n-D) Amp3'
    //    '<S50>/Lookup Table (n-D) Amp4'
    //    '<S50>/Lookup Table (n-D) Amp5'
    //    '<S50>/Lookup Table (n-D) Amp6'
    //    '<S50>/Lookup Table (n-D) Phase'
    //    '<S50>/Lookup Table (n-D) Phase1'
    //    '<S50>/Lookup Table (n-D) Phase2'
    //    '<S50>/Lookup Table (n-D) Phase4'
    //    '<S50>/Lookup Table (n-D) Phase5'
    //    '<S50>/Lookup Table (n-D) Phase6'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S127>/Lookup Table (n-D) Amp1'
    //    '<S127>/Lookup Table (n-D) Amp2'
    //    '<S127>/Lookup Table (n-D) Amp3'
    //    '<S127>/Lookup Table (n-D) Amp4'
    //    '<S127>/Lookup Table (n-D) Amp5'
    //    '<S127>/Lookup Table (n-D) Amp6'
    //    '<S127>/Lookup Table (n-D) Phase'
    //    '<S127>/Lookup Table (n-D) Phase1'
    //    '<S127>/Lookup Table (n-D) Phase2'
    //    '<S127>/Lookup Table (n-D) Phase4'
    //    '<S127>/Lookup Table (n-D) Phase5'
    //    '<S127>/Lookup Table (n-D) Phase6'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S134>/Lookup Table (n-D) Amp1'
    //    '<S134>/Lookup Table (n-D) Amp2'
    //    '<S134>/Lookup Table (n-D) Amp3'
    //    '<S134>/Lookup Table (n-D) Amp4'
    //    '<S134>/Lookup Table (n-D) Amp5'
    //    '<S134>/Lookup Table (n-D) Amp6'
    //    '<S134>/Lookup Table (n-D) Phase'
    //    '<S134>/Lookup Table (n-D) Phase1'
    //    '<S134>/Lookup Table (n-D) Phase2'
    //    '<S134>/Lookup Table (n-D) Phase4'
    //    '<S134>/Lookup Table (n-D) Phase5'
    //    '<S134>/Lookup Table (n-D) Phase6'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S211>/Lookup Table (n-D) Amp1'
    //    '<S211>/Lookup Table (n-D) Amp2'
    //    '<S211>/Lookup Table (n-D) Amp3'
    //    '<S211>/Lookup Table (n-D) Amp4'
    //    '<S211>/Lookup Table (n-D) Amp5'
    //    '<S211>/Lookup Table (n-D) Amp6'
    //    '<S211>/Lookup Table (n-D) Phase'
    //    '<S211>/Lookup Table (n-D) Phase1'
    //    '<S211>/Lookup Table (n-D) Phase2'
    //    '<S211>/Lookup Table (n-D) Phase4'
    //    '<S211>/Lookup Table (n-D) Phase5'
    //    '<S211>/Lookup Table (n-D) Phase6'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S218>/Lookup Table (n-D) Amp1'
    //    '<S218>/Lookup Table (n-D) Amp2'
    //    '<S218>/Lookup Table (n-D) Amp3'
    //    '<S218>/Lookup Table (n-D) Amp4'
    //    '<S218>/Lookup Table (n-D) Amp5'
    //    '<S218>/Lookup Table (n-D) Amp6'
    //    '<S218>/Lookup Table (n-D) Phase'
    //    '<S218>/Lookup Table (n-D) Phase1'
    //    '<S218>/Lookup Table (n-D) Phase2'
    //    '<S218>/Lookup Table (n-D) Phase4'
    //    '<S218>/Lookup Table (n-D) Phase5'
    //    '<S218>/Lookup Table (n-D) Phase6'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S295>/Lookup Table (n-D) Amp1'
    //    '<S295>/Lookup Table (n-D) Amp2'
    //    '<S295>/Lookup Table (n-D) Amp3'
    //    '<S295>/Lookup Table (n-D) Amp4'
    //    '<S295>/Lookup Table (n-D) Amp5'
    //    '<S295>/Lookup Table (n-D) Amp6'
    //    '<S295>/Lookup Table (n-D) Phase'
    //    '<S295>/Lookup Table (n-D) Phase1'
    //    '<S295>/Lookup Table (n-D) Phase2'
    //    '<S295>/Lookup Table (n-D) Phase4'
    //    '<S295>/Lookup Table (n-D) Phase5'
    //    '<S295>/Lookup Table (n-D) Phase6'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S302>/Lookup Table (n-D) Amp1'
    //    '<S302>/Lookup Table (n-D) Amp2'
    //    '<S302>/Lookup Table (n-D) Amp3'
    //    '<S302>/Lookup Table (n-D) Amp4'
    //    '<S302>/Lookup Table (n-D) Amp5'
    //    '<S302>/Lookup Table (n-D) Amp6'
    //    '<S302>/Lookup Table (n-D) Phase'
    //    '<S302>/Lookup Table (n-D) Phase1'
    //    '<S302>/Lookup Table (n-D) Phase2'
    //    '<S302>/Lookup Table (n-D) Phase4'
    //    '<S302>/Lookup Table (n-D) Phase5'
    //    '<S302>/Lookup Table (n-D) Phase6'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled14[36];

    // Pooled Parameter (Expression: WD_Amp{1})
    //  Referenced by:
    //    '<S44>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 1'

    real_T pooled15[2808];

    // Pooled Parameter (Expression: WD_Freq)
    //  Referenced by:
    //    '<S44>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled16[39];

    // Pooled Parameter (Expression: RAO_Vel)
    //  Referenced by:
    //    '<S42>/Prelookup2'
    //    '<S42>/Direct Lookup Table (n-D) Velocity'
    //    '<S49>/Prelookup2'
    //    '<S49>/Direct Lookup Table (n-D) Velocity'
    //    '<S126>/Prelookup2'
    //    '<S126>/Direct Lookup Table (n-D) Velocity'
    //    '<S133>/Prelookup2'
    //    '<S133>/Direct Lookup Table (n-D) Velocity'
    //    '<S210>/Prelookup2'
    //    '<S210>/Direct Lookup Table (n-D) Velocity'
    //    '<S217>/Prelookup2'
    //    '<S217>/Direct Lookup Table (n-D) Velocity'
    //    '<S294>/Prelookup2'
    //    '<S294>/Direct Lookup Table (n-D) Velocity'
    //    '<S301>/Prelookup2'
    //    '<S301>/Direct Lookup Table (n-D) Velocity'
    //    '<S43>/Lookup Table (n-D) Amp1'
    //    '<S43>/Lookup Table (n-D) Amp2'
    //    '<S43>/Lookup Table (n-D) Amp3'
    //    '<S43>/Lookup Table (n-D) Amp4'
    //    '<S43>/Lookup Table (n-D) Amp5'
    //    '<S43>/Lookup Table (n-D) Amp6'
    //    '<S43>/Lookup Table (n-D) Phase'
    //    '<S43>/Lookup Table (n-D) Phase1'
    //    '<S43>/Lookup Table (n-D) Phase2'
    //    '<S43>/Lookup Table (n-D) Phase4'
    //    '<S43>/Lookup Table (n-D) Phase5'
    //    '<S43>/Lookup Table (n-D) Phase6'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S50>/Lookup Table (n-D) Amp1'
    //    '<S50>/Lookup Table (n-D) Amp2'
    //    '<S50>/Lookup Table (n-D) Amp3'
    //    '<S50>/Lookup Table (n-D) Amp4'
    //    '<S50>/Lookup Table (n-D) Amp5'
    //    '<S50>/Lookup Table (n-D) Amp6'
    //    '<S50>/Lookup Table (n-D) Phase'
    //    '<S50>/Lookup Table (n-D) Phase1'
    //    '<S50>/Lookup Table (n-D) Phase2'
    //    '<S50>/Lookup Table (n-D) Phase4'
    //    '<S50>/Lookup Table (n-D) Phase5'
    //    '<S50>/Lookup Table (n-D) Phase6'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S127>/Lookup Table (n-D) Amp1'
    //    '<S127>/Lookup Table (n-D) Amp2'
    //    '<S127>/Lookup Table (n-D) Amp3'
    //    '<S127>/Lookup Table (n-D) Amp4'
    //    '<S127>/Lookup Table (n-D) Amp5'
    //    '<S127>/Lookup Table (n-D) Amp6'
    //    '<S127>/Lookup Table (n-D) Phase'
    //    '<S127>/Lookup Table (n-D) Phase1'
    //    '<S127>/Lookup Table (n-D) Phase2'
    //    '<S127>/Lookup Table (n-D) Phase4'
    //    '<S127>/Lookup Table (n-D) Phase5'
    //    '<S127>/Lookup Table (n-D) Phase6'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S134>/Lookup Table (n-D) Amp1'
    //    '<S134>/Lookup Table (n-D) Amp2'
    //    '<S134>/Lookup Table (n-D) Amp3'
    //    '<S134>/Lookup Table (n-D) Amp4'
    //    '<S134>/Lookup Table (n-D) Amp5'
    //    '<S134>/Lookup Table (n-D) Amp6'
    //    '<S134>/Lookup Table (n-D) Phase'
    //    '<S134>/Lookup Table (n-D) Phase1'
    //    '<S134>/Lookup Table (n-D) Phase2'
    //    '<S134>/Lookup Table (n-D) Phase4'
    //    '<S134>/Lookup Table (n-D) Phase5'
    //    '<S134>/Lookup Table (n-D) Phase6'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S211>/Lookup Table (n-D) Amp1'
    //    '<S211>/Lookup Table (n-D) Amp2'
    //    '<S211>/Lookup Table (n-D) Amp3'
    //    '<S211>/Lookup Table (n-D) Amp4'
    //    '<S211>/Lookup Table (n-D) Amp5'
    //    '<S211>/Lookup Table (n-D) Amp6'
    //    '<S211>/Lookup Table (n-D) Phase'
    //    '<S211>/Lookup Table (n-D) Phase1'
    //    '<S211>/Lookup Table (n-D) Phase2'
    //    '<S211>/Lookup Table (n-D) Phase4'
    //    '<S211>/Lookup Table (n-D) Phase5'
    //    '<S211>/Lookup Table (n-D) Phase6'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S218>/Lookup Table (n-D) Amp1'
    //    '<S218>/Lookup Table (n-D) Amp2'
    //    '<S218>/Lookup Table (n-D) Amp3'
    //    '<S218>/Lookup Table (n-D) Amp4'
    //    '<S218>/Lookup Table (n-D) Amp5'
    //    '<S218>/Lookup Table (n-D) Amp6'
    //    '<S218>/Lookup Table (n-D) Phase'
    //    '<S218>/Lookup Table (n-D) Phase1'
    //    '<S218>/Lookup Table (n-D) Phase2'
    //    '<S218>/Lookup Table (n-D) Phase4'
    //    '<S218>/Lookup Table (n-D) Phase5'
    //    '<S218>/Lookup Table (n-D) Phase6'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S295>/Lookup Table (n-D) Amp1'
    //    '<S295>/Lookup Table (n-D) Amp2'
    //    '<S295>/Lookup Table (n-D) Amp3'
    //    '<S295>/Lookup Table (n-D) Amp4'
    //    '<S295>/Lookup Table (n-D) Amp5'
    //    '<S295>/Lookup Table (n-D) Amp6'
    //    '<S295>/Lookup Table (n-D) Phase'
    //    '<S295>/Lookup Table (n-D) Phase1'
    //    '<S295>/Lookup Table (n-D) Phase2'
    //    '<S295>/Lookup Table (n-D) Phase4'
    //    '<S295>/Lookup Table (n-D) Phase5'
    //    '<S295>/Lookup Table (n-D) Phase6'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S302>/Lookup Table (n-D) Amp1'
    //    '<S302>/Lookup Table (n-D) Amp2'
    //    '<S302>/Lookup Table (n-D) Amp3'
    //    '<S302>/Lookup Table (n-D) Amp4'
    //    '<S302>/Lookup Table (n-D) Amp5'
    //    '<S302>/Lookup Table (n-D) Amp6'
    //    '<S302>/Lookup Table (n-D) Phase'
    //    '<S302>/Lookup Table (n-D) Phase1'
    //    '<S302>/Lookup Table (n-D) Phase2'
    //    '<S302>/Lookup Table (n-D) Phase4'
    //    '<S302>/Lookup Table (n-D) Phase5'
    //    '<S302>/Lookup Table (n-D) Phase6'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled17[2];

    // Pooled Parameter (Expression: WD_Amp{2})
    //  Referenced by:
    //    '<S44>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 2'

    real_T pooled18[2808];

    // Pooled Parameter (Expression: WD_Amp{3})
    //  Referenced by:
    //    '<S44>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled19[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{1})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp1'
    //    '<S50>/Lookup Table (n-D) Amp1'
    //    '<S127>/Lookup Table (n-D) Amp1'
    //    '<S134>/Lookup Table (n-D) Amp1'
    //    '<S211>/Lookup Table (n-D) Amp1'
    //    '<S218>/Lookup Table (n-D) Amp1'
    //    '<S295>/Lookup Table (n-D) Amp1'
    //    '<S302>/Lookup Table (n-D) Amp1'

    real_T pooled20[2808];

    // Pooled Parameter (Expression: ForceRAO_Freq)
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp1'
    //    '<S43>/Lookup Table (n-D) Amp2'
    //    '<S43>/Lookup Table (n-D) Amp3'
    //    '<S43>/Lookup Table (n-D) Amp4'
    //    '<S43>/Lookup Table (n-D) Amp5'
    //    '<S43>/Lookup Table (n-D) Amp6'
    //    '<S43>/Lookup Table (n-D) Phase'
    //    '<S43>/Lookup Table (n-D) Phase1'
    //    '<S43>/Lookup Table (n-D) Phase2'
    //    '<S43>/Lookup Table (n-D) Phase4'
    //    '<S43>/Lookup Table (n-D) Phase5'
    //    '<S43>/Lookup Table (n-D) Phase6'
    //    '<S50>/Lookup Table (n-D) Amp1'
    //    '<S50>/Lookup Table (n-D) Amp2'
    //    '<S50>/Lookup Table (n-D) Amp3'
    //    '<S50>/Lookup Table (n-D) Amp4'
    //    '<S50>/Lookup Table (n-D) Amp5'
    //    '<S50>/Lookup Table (n-D) Amp6'
    //    '<S50>/Lookup Table (n-D) Phase'
    //    '<S50>/Lookup Table (n-D) Phase1'
    //    '<S50>/Lookup Table (n-D) Phase2'
    //    '<S50>/Lookup Table (n-D) Phase4'
    //    '<S50>/Lookup Table (n-D) Phase5'
    //    '<S50>/Lookup Table (n-D) Phase6'
    //    '<S127>/Lookup Table (n-D) Amp1'
    //    '<S127>/Lookup Table (n-D) Amp2'
    //    '<S127>/Lookup Table (n-D) Amp3'
    //    '<S127>/Lookup Table (n-D) Amp4'
    //    '<S127>/Lookup Table (n-D) Amp5'
    //    '<S127>/Lookup Table (n-D) Amp6'
    //    '<S127>/Lookup Table (n-D) Phase'
    //    '<S127>/Lookup Table (n-D) Phase1'
    //    '<S127>/Lookup Table (n-D) Phase2'
    //    '<S127>/Lookup Table (n-D) Phase4'
    //    '<S127>/Lookup Table (n-D) Phase5'
    //    '<S127>/Lookup Table (n-D) Phase6'
    //    '<S134>/Lookup Table (n-D) Amp1'
    //    '<S134>/Lookup Table (n-D) Amp2'
    //    '<S134>/Lookup Table (n-D) Amp3'
    //    '<S134>/Lookup Table (n-D) Amp4'
    //    '<S134>/Lookup Table (n-D) Amp5'
    //    '<S134>/Lookup Table (n-D) Amp6'
    //    '<S134>/Lookup Table (n-D) Phase'
    //    '<S134>/Lookup Table (n-D) Phase1'
    //    '<S134>/Lookup Table (n-D) Phase2'
    //    '<S134>/Lookup Table (n-D) Phase4'
    //    '<S134>/Lookup Table (n-D) Phase5'
    //    '<S134>/Lookup Table (n-D) Phase6'
    //    '<S211>/Lookup Table (n-D) Amp1'
    //    '<S211>/Lookup Table (n-D) Amp2'
    //    '<S211>/Lookup Table (n-D) Amp3'
    //    '<S211>/Lookup Table (n-D) Amp4'
    //    '<S211>/Lookup Table (n-D) Amp5'
    //    '<S211>/Lookup Table (n-D) Amp6'
    //    '<S211>/Lookup Table (n-D) Phase'
    //    '<S211>/Lookup Table (n-D) Phase1'
    //    '<S211>/Lookup Table (n-D) Phase2'
    //    '<S211>/Lookup Table (n-D) Phase4'
    //    '<S211>/Lookup Table (n-D) Phase5'
    //    '<S211>/Lookup Table (n-D) Phase6'
    //    '<S218>/Lookup Table (n-D) Amp1'
    //    '<S218>/Lookup Table (n-D) Amp2'
    //    '<S218>/Lookup Table (n-D) Amp3'
    //    '<S218>/Lookup Table (n-D) Amp4'
    //    '<S218>/Lookup Table (n-D) Amp5'
    //    '<S218>/Lookup Table (n-D) Amp6'
    //    '<S218>/Lookup Table (n-D) Phase'
    //    '<S218>/Lookup Table (n-D) Phase1'
    //    '<S218>/Lookup Table (n-D) Phase2'
    //    '<S218>/Lookup Table (n-D) Phase4'
    //    '<S218>/Lookup Table (n-D) Phase5'
    //    '<S218>/Lookup Table (n-D) Phase6'
    //    '<S295>/Lookup Table (n-D) Amp1'
    //    '<S295>/Lookup Table (n-D) Amp2'
    //    '<S295>/Lookup Table (n-D) Amp3'
    //    '<S295>/Lookup Table (n-D) Amp4'
    //    '<S295>/Lookup Table (n-D) Amp5'
    //    '<S295>/Lookup Table (n-D) Amp6'
    //    '<S295>/Lookup Table (n-D) Phase'
    //    '<S295>/Lookup Table (n-D) Phase1'
    //    '<S295>/Lookup Table (n-D) Phase2'
    //    '<S295>/Lookup Table (n-D) Phase4'
    //    '<S295>/Lookup Table (n-D) Phase5'
    //    '<S295>/Lookup Table (n-D) Phase6'
    //    '<S302>/Lookup Table (n-D) Amp1'
    //    '<S302>/Lookup Table (n-D) Amp2'
    //    '<S302>/Lookup Table (n-D) Amp3'
    //    '<S302>/Lookup Table (n-D) Amp4'
    //    '<S302>/Lookup Table (n-D) Amp5'
    //    '<S302>/Lookup Table (n-D) Amp6'
    //    '<S302>/Lookup Table (n-D) Phase'
    //    '<S302>/Lookup Table (n-D) Phase1'
    //    '<S302>/Lookup Table (n-D) Phase2'
    //    '<S302>/Lookup Table (n-D) Phase4'
    //    '<S302>/Lookup Table (n-D) Phase5'
    //    '<S302>/Lookup Table (n-D) Phase6'

    real_T pooled21[39];

    // Pooled Parameter (Expression: ForceRAO_Amp{2})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp2'
    //    '<S50>/Lookup Table (n-D) Amp2'
    //    '<S127>/Lookup Table (n-D) Amp2'
    //    '<S134>/Lookup Table (n-D) Amp2'
    //    '<S211>/Lookup Table (n-D) Amp2'
    //    '<S218>/Lookup Table (n-D) Amp2'
    //    '<S295>/Lookup Table (n-D) Amp2'
    //    '<S302>/Lookup Table (n-D) Amp2'

    real_T pooled22[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{3})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp3'
    //    '<S50>/Lookup Table (n-D) Amp3'
    //    '<S127>/Lookup Table (n-D) Amp3'
    //    '<S134>/Lookup Table (n-D) Amp3'
    //    '<S211>/Lookup Table (n-D) Amp3'
    //    '<S218>/Lookup Table (n-D) Amp3'
    //    '<S295>/Lookup Table (n-D) Amp3'
    //    '<S302>/Lookup Table (n-D) Amp3'

    real_T pooled23[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{4})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp4'
    //    '<S50>/Lookup Table (n-D) Amp4'
    //    '<S127>/Lookup Table (n-D) Amp4'
    //    '<S134>/Lookup Table (n-D) Amp4'
    //    '<S211>/Lookup Table (n-D) Amp4'
    //    '<S218>/Lookup Table (n-D) Amp4'
    //    '<S295>/Lookup Table (n-D) Amp4'
    //    '<S302>/Lookup Table (n-D) Amp4'

    real_T pooled24[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{5})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp5'
    //    '<S50>/Lookup Table (n-D) Amp5'
    //    '<S127>/Lookup Table (n-D) Amp5'
    //    '<S134>/Lookup Table (n-D) Amp5'
    //    '<S211>/Lookup Table (n-D) Amp5'
    //    '<S218>/Lookup Table (n-D) Amp5'
    //    '<S295>/Lookup Table (n-D) Amp5'
    //    '<S302>/Lookup Table (n-D) Amp5'

    real_T pooled25[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{6})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp6'
    //    '<S50>/Lookup Table (n-D) Amp6'
    //    '<S127>/Lookup Table (n-D) Amp6'
    //    '<S134>/Lookup Table (n-D) Amp6'
    //    '<S211>/Lookup Table (n-D) Amp6'
    //    '<S218>/Lookup Table (n-D) Amp6'
    //    '<S295>/Lookup Table (n-D) Amp6'
    //    '<S302>/Lookup Table (n-D) Amp6'

    real_T pooled26[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{1})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase1'
    //    '<S50>/Lookup Table (n-D) Phase1'
    //    '<S127>/Lookup Table (n-D) Phase1'
    //    '<S134>/Lookup Table (n-D) Phase1'
    //    '<S211>/Lookup Table (n-D) Phase1'
    //    '<S218>/Lookup Table (n-D) Phase1'
    //    '<S295>/Lookup Table (n-D) Phase1'
    //    '<S302>/Lookup Table (n-D) Phase1'

    real_T pooled27[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{2})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase2'
    //    '<S50>/Lookup Table (n-D) Phase2'
    //    '<S127>/Lookup Table (n-D) Phase2'
    //    '<S134>/Lookup Table (n-D) Phase2'
    //    '<S211>/Lookup Table (n-D) Phase2'
    //    '<S218>/Lookup Table (n-D) Phase2'
    //    '<S295>/Lookup Table (n-D) Phase2'
    //    '<S302>/Lookup Table (n-D) Phase2'

    real_T pooled28[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{3})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase'
    //    '<S50>/Lookup Table (n-D) Phase'
    //    '<S127>/Lookup Table (n-D) Phase'
    //    '<S134>/Lookup Table (n-D) Phase'
    //    '<S211>/Lookup Table (n-D) Phase'
    //    '<S218>/Lookup Table (n-D) Phase'
    //    '<S295>/Lookup Table (n-D) Phase'
    //    '<S302>/Lookup Table (n-D) Phase'

    real_T pooled29[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{4})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase4'
    //    '<S50>/Lookup Table (n-D) Phase4'
    //    '<S127>/Lookup Table (n-D) Phase4'
    //    '<S134>/Lookup Table (n-D) Phase4'
    //    '<S211>/Lookup Table (n-D) Phase4'
    //    '<S218>/Lookup Table (n-D) Phase4'
    //    '<S295>/Lookup Table (n-D) Phase4'
    //    '<S302>/Lookup Table (n-D) Phase4'

    real_T pooled30[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{5})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase5'
    //    '<S50>/Lookup Table (n-D) Phase5'
    //    '<S127>/Lookup Table (n-D) Phase5'
    //    '<S134>/Lookup Table (n-D) Phase5'
    //    '<S211>/Lookup Table (n-D) Phase5'
    //    '<S218>/Lookup Table (n-D) Phase5'
    //    '<S295>/Lookup Table (n-D) Phase5'
    //    '<S302>/Lookup Table (n-D) Phase5'

    real_T pooled31[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{6})
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Phase6'
    //    '<S50>/Lookup Table (n-D) Phase6'
    //    '<S127>/Lookup Table (n-D) Phase6'
    //    '<S134>/Lookup Table (n-D) Phase6'
    //    '<S211>/Lookup Table (n-D) Phase6'
    //    '<S218>/Lookup Table (n-D) Phase6'
    //    '<S295>/Lookup Table (n-D) Phase6'
    //    '<S302>/Lookup Table (n-D) Phase6'

    real_T pooled32[2808];

    // Pooled Parameter (Expression: ABC.G)
    //  Referenced by:
    //    '<S13>/Spring stiffness'
    //    '<S97>/Spring stiffness'
    //    '<S181>/Spring stiffness'
    //    '<S265>/Spring stiffness'

    real_T pooled48[36];

    // Pooled Parameter (Expression: ABC.Binf)
    //  Referenced by:
    //    '<S13>/damping'
    //    '<S97>/damping'
    //    '<S181>/damping'
    //    '<S265>/damping'

    real_T pooled51[36];

    // Pooled Parameter (Expression: vesselABC_dp.A44(:,:,1))
    //  Referenced by:
    //    '<S32>/A44'
    //    '<S116>/A44'
    //    '<S200>/A44'
    //    '<S284>/A44'

    real_T pooled52[25];

    // Pooled Parameter (Expression: vesselABC_dp.B44(:,:,1))
    //  Referenced by:
    //    '<S32>/B44'
    //    '<S116>/B44'
    //    '<S200>/B44'
    //    '<S284>/B44'

    real_T pooled53[5];

    // Pooled Parameter (Expression: [0 0 1 1 0.5 1])
    //  Referenced by:
    //    '<S8>/Gain1'
    //    '<S92>/Gain1'
    //    '<S176>/Gain1'
    //    '<S260>/Gain1'

    real_T pooled60[6];

    // Pooled Parameter (Expression: Kd)
    //  Referenced by:
    //    '<S65>/Kd'
    //    '<S149>/Kd'
    //    '<S233>/Kd'
    //    '<S317>/Kd'

    real_T pooled64[9];

    // Pooled Parameter (Expression: 2*lambda*w_o)
    //  Referenced by:
    //    '<S63>/Gain1'
    //    '<S147>/Gain1'
    //    '<S231>/Gain1'
    //    '<S315>/Gain1'

    real_T pooled107[9];

    // Pooled Parameter (Expression: w_o*w_o)
    //  Referenced by:
    //    '<S63>/Gain2'
    //    '<S147>/Gain2'
    //    '<S231>/Gain2'
    //    '<S315>/Gain2'

    real_T pooled108[9];

    // Pooled Parameter (Expression: K4)
    //  Referenced by:
    //    '<S63>/K4'
    //    '<S147>/K4'
    //    '<S231>/K4'
    //    '<S315>/K4'

    real_T pooled109[9];

    // Pooled Parameter (Expression: D)
    //  Referenced by:
    //    '<S63>/Gain6'
    //    '<S147>/Gain6'
    //    '<S231>/Gain6'
    //    '<S315>/Gain6'

    real_T pooled110[9];

    // Pooled Parameter (Expression: inv(M))
    //  Referenced by:
    //    '<S63>/Gain3'
    //    '<S147>/Gain3'
    //    '<S231>/Gain3'
    //    '<S315>/Gain3'

    real_T pooled111[9];

    // Pooled Parameter (Expression: -2*(eye(3)-lambda)*diag([w_c(1,1)/w_o(1,1) w_c(2,2)/w_o(2,2) w_c(3,3)/w_o(3,3)]))
    //  Referenced by:
    //    '<S63>/K11'
    //    '<S147>/K11'
    //    '<S231>/K11'
    //    '<S315>/K11'

    real_T pooled112[9];

    // Pooled Parameter (Expression: 2*w_o*(eye(3)-lambda))
    //  Referenced by:
    //    '<S63>/K12'
    //    '<S147>/K12'
    //    '<S231>/K12'
    //    '<S315>/K12'

    real_T pooled113[9];

    // Pooled Parameter (Expression: w_c)
    //  Referenced by:
    //    '<S63>/K2'
    //    '<S147>/K2'
    //    '<S231>/K2'
    //    '<S315>/K2'

    real_T pooled114[9];

    // Pooled Parameter (Expression: K3)
    //  Referenced by:
    //    '<S63>/K3'
    //    '<S147>/K3'
    //    '<S231>/K3'
    //    '<S315>/K3'

    real_T pooled115[9];

    // Pooled Parameter (Expression: diag([1/T_b(1,1) 1/T_b(2,2) 1/T_b(3,3)]))
    //  Referenced by:
    //    '<S63>/inv(T_b)'
    //    '<S147>/inv(T_b)'
    //    '<S231>/inv(T_b)'
    //    '<S315>/inv(T_b)'

    real_T pooled116[9];

    // Pooled Parameter (Expression: )
    //  Referenced by:
    //    '<S43>/Lookup Table (n-D) Amp1'
    //    '<S43>/Lookup Table (n-D) Amp2'
    //    '<S43>/Lookup Table (n-D) Amp3'
    //    '<S43>/Lookup Table (n-D) Amp4'
    //    '<S43>/Lookup Table (n-D) Amp5'
    //    '<S43>/Lookup Table (n-D) Amp6'
    //    '<S43>/Lookup Table (n-D) Phase'
    //    '<S43>/Lookup Table (n-D) Phase1'
    //    '<S43>/Lookup Table (n-D) Phase2'
    //    '<S43>/Lookup Table (n-D) Phase4'
    //    '<S43>/Lookup Table (n-D) Phase5'
    //    '<S43>/Lookup Table (n-D) Phase6'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S44>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S50>/Lookup Table (n-D) Amp1'
    //    '<S50>/Lookup Table (n-D) Amp2'
    //    '<S50>/Lookup Table (n-D) Amp3'
    //    '<S50>/Lookup Table (n-D) Amp4'
    //    '<S50>/Lookup Table (n-D) Amp5'
    //    '<S50>/Lookup Table (n-D) Amp6'
    //    '<S50>/Lookup Table (n-D) Phase'
    //    '<S50>/Lookup Table (n-D) Phase1'
    //    '<S50>/Lookup Table (n-D) Phase2'
    //    '<S50>/Lookup Table (n-D) Phase4'
    //    '<S50>/Lookup Table (n-D) Phase5'
    //    '<S50>/Lookup Table (n-D) Phase6'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S51>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S127>/Lookup Table (n-D) Amp1'
    //    '<S127>/Lookup Table (n-D) Amp2'
    //    '<S127>/Lookup Table (n-D) Amp3'
    //    '<S127>/Lookup Table (n-D) Amp4'
    //    '<S127>/Lookup Table (n-D) Amp5'
    //    '<S127>/Lookup Table (n-D) Amp6'
    //    '<S127>/Lookup Table (n-D) Phase'
    //    '<S127>/Lookup Table (n-D) Phase1'
    //    '<S127>/Lookup Table (n-D) Phase2'
    //    '<S127>/Lookup Table (n-D) Phase4'
    //    '<S127>/Lookup Table (n-D) Phase5'
    //    '<S127>/Lookup Table (n-D) Phase6'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S128>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S134>/Lookup Table (n-D) Amp1'
    //    '<S134>/Lookup Table (n-D) Amp2'
    //    '<S134>/Lookup Table (n-D) Amp3'
    //    '<S134>/Lookup Table (n-D) Amp4'
    //    '<S134>/Lookup Table (n-D) Amp5'
    //    '<S134>/Lookup Table (n-D) Amp6'
    //    '<S134>/Lookup Table (n-D) Phase'
    //    '<S134>/Lookup Table (n-D) Phase1'
    //    '<S134>/Lookup Table (n-D) Phase2'
    //    '<S134>/Lookup Table (n-D) Phase4'
    //    '<S134>/Lookup Table (n-D) Phase5'
    //    '<S134>/Lookup Table (n-D) Phase6'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S135>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S211>/Lookup Table (n-D) Amp1'
    //    '<S211>/Lookup Table (n-D) Amp2'
    //    '<S211>/Lookup Table (n-D) Amp3'
    //    '<S211>/Lookup Table (n-D) Amp4'
    //    '<S211>/Lookup Table (n-D) Amp5'
    //    '<S211>/Lookup Table (n-D) Amp6'
    //    '<S211>/Lookup Table (n-D) Phase'
    //    '<S211>/Lookup Table (n-D) Phase1'
    //    '<S211>/Lookup Table (n-D) Phase2'
    //    '<S211>/Lookup Table (n-D) Phase4'
    //    '<S211>/Lookup Table (n-D) Phase5'
    //    '<S211>/Lookup Table (n-D) Phase6'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S212>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S218>/Lookup Table (n-D) Amp1'
    //    '<S218>/Lookup Table (n-D) Amp2'
    //    '<S218>/Lookup Table (n-D) Amp3'
    //    '<S218>/Lookup Table (n-D) Amp4'
    //    '<S218>/Lookup Table (n-D) Amp5'
    //    '<S218>/Lookup Table (n-D) Amp6'
    //    '<S218>/Lookup Table (n-D) Phase'
    //    '<S218>/Lookup Table (n-D) Phase1'
    //    '<S218>/Lookup Table (n-D) Phase2'
    //    '<S218>/Lookup Table (n-D) Phase4'
    //    '<S218>/Lookup Table (n-D) Phase5'
    //    '<S218>/Lookup Table (n-D) Phase6'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S219>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S295>/Lookup Table (n-D) Amp1'
    //    '<S295>/Lookup Table (n-D) Amp2'
    //    '<S295>/Lookup Table (n-D) Amp3'
    //    '<S295>/Lookup Table (n-D) Amp4'
    //    '<S295>/Lookup Table (n-D) Amp5'
    //    '<S295>/Lookup Table (n-D) Amp6'
    //    '<S295>/Lookup Table (n-D) Phase'
    //    '<S295>/Lookup Table (n-D) Phase1'
    //    '<S295>/Lookup Table (n-D) Phase2'
    //    '<S295>/Lookup Table (n-D) Phase4'
    //    '<S295>/Lookup Table (n-D) Phase5'
    //    '<S295>/Lookup Table (n-D) Phase6'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S296>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S302>/Lookup Table (n-D) Amp1'
    //    '<S302>/Lookup Table (n-D) Amp2'
    //    '<S302>/Lookup Table (n-D) Amp3'
    //    '<S302>/Lookup Table (n-D) Amp4'
    //    '<S302>/Lookup Table (n-D) Amp5'
    //    '<S302>/Lookup Table (n-D) Amp6'
    //    '<S302>/Lookup Table (n-D) Phase'
    //    '<S302>/Lookup Table (n-D) Phase1'
    //    '<S302>/Lookup Table (n-D) Phase2'
    //    '<S302>/Lookup Table (n-D) Phase4'
    //    '<S302>/Lookup Table (n-D) Phase5'
    //    '<S302>/Lookup Table (n-D) Phase6'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S303>/Lookup Table (n-D)  Wavedrift 3'

    uint32_T pooled119[3];
  } ConstP_AHV_Model_T;

  // Real-time Model Data Structure
  struct RT_MODEL_AHV_Model_T {
    const char_T *errorStatus;
    RTWSolverInfo solverInfo;
    X_AHV_Model_T *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T CTOutputIncnstWithState;
    real_T odeY[475];
    real_T odeF[4][475];
    ODE4_IntgData intgData;

    //
    //  Sizes:
    //  The following substructure contains sizes information
    //  for many of the model attributes such as inputs, outputs,
    //  dwork, sample times, etc.

    struct {
      int_T numContStates;
      int_T numPeriodicContStates;
      int_T numSampTimes;
    } Sizes;

    //
    //  Timing:
    //  The following substructure contains information regarding
    //  the timing information for the model.

    struct {
      uint32_T clockTick0;
      time_T stepSize0;
      uint32_T clockTick1;
      boolean_T firstInitCondFlag;
      SimTimeStep simTimeStep;
      boolean_T stopRequestedFlag;
      time_T *t;
      time_T tArray[2];
    } Timing;
  };

  // model initialize function
  void initialize(double dtime);

  // model step function
  void step();

  // model terminate function
  void terminate();

  // Constructor
  AH_Model_v1ModelClass();

  // Destructor
  ~AH_Model_v1ModelClass();

  // Real-Time Model get method
  AH_Model_v1ModelClass::RT_MODEL_AHV_Model_T * getRTM();

  // private data and function members
 private:
  // Block signals
  B_AHV_Model_T AHV_Model_B;

  // Block states
  DW_AHV_Model_T AHV_Model_DW;
  X_AHV_Model_T AHV_Model_X;           // Block continuous states

  // Real-Time Model
  RT_MODEL_AHV_Model_T AHV_Model_M;

  // private member function(s) for subsystem '<S19>/Cross-flow drag trapezoidal integration'
  void Crossflowdragtrapezoidali_Reset(real_T *memory2_PreviousInput, real_T
    *memory1_PreviousInput);
  void Crossflowdragtrapezoidalintegra(real_T rtu_N, real_T rtu_dx, real_T
    rtu_v_r, real_T rtu_r, real_T *rty_sum1, real_T *rty_sum2, real_T rtp_Lpp);

  // private member function(s) for subsystem '<S8>/Subsystem3'
  void Wave_init_Init(DW_Wave_init_T *localDW);
  void Wave_init(real_T rtu_spectrum_type, real_T rtu_hs, real_T rtu_omega_peak,
                 real_T rtu_psi_mean, real_T rtu_gamma, real_T rtu_spread,
                 real_T rtu_depth, real_T rtu_nfreq, real_T rtu_ndir, real_T
                 rtu_energylim, real_T rtu_freq_cutoff, real_T rtu_dir_cutoff,
                 real_T rtu_rand_freq, real_T rtu_rand_dir, real_T rty_Zeta_a
                 [900], real_T rty_Omega[900], real_T rty_Phase[900], real_T
                 rty_Wavenum[900], real_T rty_Psi[900], real_T
                 *rty_wave_direction, B_Wave_init_T *localB, DW_Wave_init_T
                 *localDW, real_T rtp_nfreq, real_T rtp_ndir);
  real_T AHV_Model_rad2pipi2(real_T angle);
  void AHV_Model_emxInit_real_T(emxArray_real_T_AHV_Model_T **pEmxArray, int32_T
    numDimensions);
  void AHV_Mo_emxEnsureCapacity_real_T(emxArray_real_T_AHV_Model_T *emxArray,
    int32_T oldNumel);
  real_T AHV_Model_sum(const real_T x_data[], const int32_T x_size[2]);
  void AHV_Model_gamma(real_T *x);
  void AHV_Model_power(const real_T a_data[], const int32_T a_size[2], real_T
                       y_data[], int32_T y_size[2]);
  void AHV_Model_exp(real_T x_data[], const int32_T x_size[2]);
  void AHV_Model_power_ky(const real_T a_data[], const int32_T a_size[2], real_T
    b, real_T y_data[], int32_T y_size[2]);
  void AHV_Model_power_k(real_T a, const real_T b_data[], const int32_T b_size[2],
    real_T y_data[], int32_T y_size[2]);
  void AHV_Model_emxFree_real_T(emxArray_real_T_AHV_Model_T **pEmxArray);
  void AHV_Model_torset_spec(real_T Hs, real_T wo, const real_T omg_data[],
    const int32_T omg_size[2], real_T S_data[], int32_T *S_size);
  void AHV_Model_wavespec2(real_T SpecType, const real_T Par_data[], const
    real_T W_data[], const int32_T *W_size, emxArray_real_T_AHV_Model_T *S);
  real_T AHV_Model_eml_rand_mt19937ar(uint32_T state[625]);
  void AHV_Model_rand(real_T varargin_2, emxArray_real_T_AHV_Model_T *r,
                      DW_Wave_init_T *localDW);
  real_T AHV_Model_fzero(real_T FunFcn_tunableEnvironment_f1, const real_T
    FunFcn_tunableEnvironment_f2[900], real_T FunFcn_tunableEnvironment_f3,
    const real_T x[2]);

  // private member function(s) for subsystem '<S15>/Wave_loads_for_heading1'
  void AHV_Mod_Wave_loads_for_heading1(const real_T rtu_eta[6], real_T rtu_psi,
    const real_T rtu_psi_wave[900], const real_T rtu_wavenum[900], const real_T
    rtu_Omega[900], const real_T rtu_Phase[900], const real_T rtu_Zeta_a[900],
    real_T rty_tau_WF[6], real_T rty_tau_WD[6], B_Wave_loads_for_heading1_AHV_T *
    localB, const ConstB_Wave_loads_for_heading_T *localC,
    DW_Wave_loads_for_heading1_AH_T *localDW);

  // private member function(s) for subsystem '<S8>/Wave loads (U=0)'
  void Wave_loads_fun(const real_T rtu_eta[6], const real_T rtu_waves[900],
                      const real_T rtu_waves_p[900], const real_T rtu_waves_c
                      [900], const real_T rtu_waves_m[900], const real_T
                      rtu_waves_mv[900], real_T rty_tau_WF[6], real_T
                      rty_tau_WD[6], B_Wave_loads_fun_T *localB, const
                      ConstB_Wave_loads_fun_T *localC, DW_Wave_loads_fun_T
                      *localDW);

  // private member function(s) for subsystem '<S64>/Chart'
  void AHV_Model_Chart(boolean_T rtu_hold, real_T rtu_x_ref, real_T rtu_y_ref,
                       real_T rtu_x_hold, real_T rtu_y_hold, real_T
                       *rty_x_ref_rel, real_T *rty_y_ref_rel,
                       DW_Chart_AHV_Model_T *localDW);

  // private member function(s) for subsystem '<S76>/If Action Subsystem'
  void AHV_Model_IfActionSubsystem(real_T rtu_In1, real_T rtu_In1_a, real_T
    rtu_In1_d, real_T rty_Kp[3], const real_T rtp_Kp_A[9]);

  // private member function(s) for subsystem '<S76>/If Action Subsystem1'
  void AHV_Model_IfActionSubsystem1(real_T rtu_In1, real_T rtu_In1_g, real_T
    rtu_In1_l, real_T rty_Kp[3], const real_T rtp_Kp[9]);

  // private member function(s) for subsystem '<S83>/If Action Subsystem'
  void AHV_Model_IfActionSubsystem_l(real_T rtu_Increment, real_T rtu_main_T,
    real_T *rty_Out1, const ConstB_IfActionSubsystem_AH_p_T *localC);

  // Continuous states update member function
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  // Derivatives member function
  void AHV_Model_derivatives();
};

extern const AH_Model_v1ModelClass::ConstB_AHV_Model_T AHV_Model_ConstB;// constant block i/o 

// Constant parameters (default storage)
extern const AH_Model_v1ModelClass::ConstP_AHV_Model_T AHV_Model_ConstP;

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S8>/Display' : Unused code path elimination
//  Block '<S8>/Display1' : Unused code path elimination
//  Block '<S8>/Display2' : Unused code path elimination
//  Block '<S8>/Display3' : Unused code path elimination
//  Block '<S8>/Display4' : Unused code path elimination
//  Block '<S8>/Display5' : Unused code path elimination
//  Block '<S54>/Gain' : Unused code path elimination
//  Block '<S55>/Gain' : Unused code path elimination
//  Block '<S56>/Gain' : Unused code path elimination
//  Block '<S57>/Gain' : Unused code path elimination
//  Block '<S58>/Gain' : Unused code path elimination
//  Block '<S59>/Gain' : Unused code path elimination
//  Block '<S17>/Gain' : Unused code path elimination
//  Block '<S8>/wave direction1' : Unused code path elimination
//  Block '<S62>/Display' : Unused code path elimination
//  Block '<S92>/Display' : Unused code path elimination
//  Block '<S92>/Display1' : Unused code path elimination
//  Block '<S92>/Display2' : Unused code path elimination
//  Block '<S92>/Display3' : Unused code path elimination
//  Block '<S92>/Display4' : Unused code path elimination
//  Block '<S92>/Display5' : Unused code path elimination
//  Block '<S138>/Gain' : Unused code path elimination
//  Block '<S139>/Gain' : Unused code path elimination
//  Block '<S140>/Gain' : Unused code path elimination
//  Block '<S141>/Gain' : Unused code path elimination
//  Block '<S142>/Gain' : Unused code path elimination
//  Block '<S143>/Gain' : Unused code path elimination
//  Block '<S101>/Gain' : Unused code path elimination
//  Block '<S92>/wave direction1' : Unused code path elimination
//  Block '<S176>/Display' : Unused code path elimination
//  Block '<S176>/Display1' : Unused code path elimination
//  Block '<S176>/Display2' : Unused code path elimination
//  Block '<S176>/Display3' : Unused code path elimination
//  Block '<S176>/Display4' : Unused code path elimination
//  Block '<S176>/Display5' : Unused code path elimination
//  Block '<S222>/Gain' : Unused code path elimination
//  Block '<S223>/Gain' : Unused code path elimination
//  Block '<S224>/Gain' : Unused code path elimination
//  Block '<S225>/Gain' : Unused code path elimination
//  Block '<S226>/Gain' : Unused code path elimination
//  Block '<S227>/Gain' : Unused code path elimination
//  Block '<S185>/Gain' : Unused code path elimination
//  Block '<S176>/wave direction1' : Unused code path elimination
//  Block '<S260>/Display' : Unused code path elimination
//  Block '<S260>/Display1' : Unused code path elimination
//  Block '<S260>/Display2' : Unused code path elimination
//  Block '<S260>/Display3' : Unused code path elimination
//  Block '<S260>/Display4' : Unused code path elimination
//  Block '<S260>/Display5' : Unused code path elimination
//  Block '<S306>/Gain' : Unused code path elimination
//  Block '<S307>/Gain' : Unused code path elimination
//  Block '<S308>/Gain' : Unused code path elimination
//  Block '<S309>/Gain' : Unused code path elimination
//  Block '<S310>/Gain' : Unused code path elimination
//  Block '<S311>/Gain' : Unused code path elimination
//  Block '<S269>/Gain' : Unused code path elimination
//  Block '<S260>/wave direction1' : Unused code path elimination
//  Block '<Root>/Display' : Unused code path elimination
//  Block '<S21>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S8>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S61>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S63>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S105>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S92>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S145>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S147>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S189>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S176>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S229>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S231>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S273>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S260>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S313>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S315>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S66>/Constant' : Unused code path elimination
//  Block '<S66>/Sway Failure' : Unused code path elimination
//  Block '<S66>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S150>/Constant' : Unused code path elimination
//  Block '<S150>/Sway Failure' : Unused code path elimination
//  Block '<S150>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S234>/Constant' : Unused code path elimination
//  Block '<S234>/Sway Failure' : Unused code path elimination
//  Block '<S234>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S318>/Constant' : Unused code path elimination
//  Block '<S318>/Sway Failure' : Unused code path elimination
//  Block '<S318>/Yaw Moment Failure' : Unused code path elimination


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'AHV_Model'
//  '<S1>'   : 'AHV_Model/AH_Model1'
//  '<S2>'   : 'AHV_Model/AH_Model2'
//  '<S3>'   : 'AHV_Model/AH_Model3'
//  '<S4>'   : 'AHV_Model/AH_Model4'
//  '<S5>'   : 'AHV_Model/AH_Model1/AHV_Model'
//  '<S6>'   : 'AHV_Model/AH_Model1/DP_Thrust'
//  '<S7>'   : 'AHV_Model/AH_Model1/Manual Control'
//  '<S8>'   : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1'
//  '<S9>'   : 'AHV_Model/AH_Model1/AHV_Model/Subsystem'
//  '<S10>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S11>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S12>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S13>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S14>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S15>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S16>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S17>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S18>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S19>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S20>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S21>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S22>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S23>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S24>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S25>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S26>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S27>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S28>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S29>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S30>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S31>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S32>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S33>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S34>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S35>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S36>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S37>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S38>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S39>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S40>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S41>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S42>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S43>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S44>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S45>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S46>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S47>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S48>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S49>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S50>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S51>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S52>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S53>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S54>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S55>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S56>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S57>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S58>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S59>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S60>'  : 'AHV_Model/AH_Model1/AHV_Model/Subsystem/Cross product'
//  '<S61>'  : 'AHV_Model/AH_Model1/AHV_Model/Subsystem/Rbn_gnc'
//  '<S62>'  : 'AHV_Model/AH_Model1/DP_Thrust/Calculate Cable Angles'
//  '<S63>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter'
//  '<S64>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference'
//  '<S65>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem'
//  '<S66>'  : 'AHV_Model/AH_Model1/DP_Thrust/Thrust Limitations'
//  '<S67>'  : 'AHV_Model/AH_Model1/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S68>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S69>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S70>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S71>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S72>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S73>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference/Chart'
//  '<S74>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference/Degrees to Radians'
//  '<S75>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference/[-inf inf] to [-pi pi]'
//  '<S76>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp'
//  '<S77>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem'
//  '<S78>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S79>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S80>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S81>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp/If Action Subsystem'
//  '<S82>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp/If Action Subsystem1'
//  '<S83>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem'
//  '<S84>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem2'
//  '<S85>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem'
//  '<S86>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem1'
//  '<S87>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem2/If Action Subsystem'
//  '<S88>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem/Subsystem2/If Action Subsystem1'
//  '<S89>'  : 'AHV_Model/AH_Model2/AHV_Model'
//  '<S90>'  : 'AHV_Model/AH_Model2/DP_Thrust'
//  '<S91>'  : 'AHV_Model/AH_Model2/Manual Control'
//  '<S92>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1'
//  '<S93>'  : 'AHV_Model/AH_Model2/AHV_Model/Subsystem'
//  '<S94>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S95>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S96>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S97>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S98>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S99>'  : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S100>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S101>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S102>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S103>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S104>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S105>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S106>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S107>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S108>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S109>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S110>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S111>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S112>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S113>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S114>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S115>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S116>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S117>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S118>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S119>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S120>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S121>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S122>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S123>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S124>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S125>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S126>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S127>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S128>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S129>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S130>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S131>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S132>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S133>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S134>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S135>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S136>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S137>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S138>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S139>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S140>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S141>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S142>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S143>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S144>' : 'AHV_Model/AH_Model2/AHV_Model/Subsystem/Cross product'
//  '<S145>' : 'AHV_Model/AH_Model2/AHV_Model/Subsystem/Rbn_gnc'
//  '<S146>' : 'AHV_Model/AH_Model2/DP_Thrust/Calculate Cable Angles'
//  '<S147>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter'
//  '<S148>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference'
//  '<S149>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem'
//  '<S150>' : 'AHV_Model/AH_Model2/DP_Thrust/Thrust Limitations'
//  '<S151>' : 'AHV_Model/AH_Model2/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S152>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S153>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S154>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S155>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S156>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S157>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference/Chart'
//  '<S158>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference/Degrees to Radians'
//  '<S159>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference/[-inf inf] to [-pi pi]'
//  '<S160>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp'
//  '<S161>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem'
//  '<S162>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S163>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S164>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S165>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp/If Action Subsystem'
//  '<S166>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp/If Action Subsystem1'
//  '<S167>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem'
//  '<S168>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem1'
//  '<S169>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem'
//  '<S170>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem1'
//  '<S171>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem'
//  '<S172>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem1'
//  '<S173>' : 'AHV_Model/AH_Model3/AHV_Model'
//  '<S174>' : 'AHV_Model/AH_Model3/DP_Thrust'
//  '<S175>' : 'AHV_Model/AH_Model3/Manual Control'
//  '<S176>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1'
//  '<S177>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem'
//  '<S178>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S179>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S180>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S181>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S182>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S183>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S184>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S185>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S186>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S187>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S188>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S189>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S190>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S191>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S192>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S193>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S194>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S195>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S196>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S197>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S198>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S199>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S200>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S201>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S202>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S203>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S204>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S205>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S206>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S207>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S208>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S209>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S210>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S211>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S212>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S213>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S214>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S215>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S216>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S217>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S218>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S219>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S220>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S221>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S222>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S223>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S224>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S225>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S226>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S227>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S228>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem/Cross product'
//  '<S229>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem/Rbn_gnc'
//  '<S230>' : 'AHV_Model/AH_Model3/DP_Thrust/Calculate Cable Angles'
//  '<S231>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter'
//  '<S232>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference'
//  '<S233>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem'
//  '<S234>' : 'AHV_Model/AH_Model3/DP_Thrust/Thrust Limitations'
//  '<S235>' : 'AHV_Model/AH_Model3/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S236>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S237>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S238>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S239>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S240>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S241>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference/Chart'
//  '<S242>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference/Degrees to Radians'
//  '<S243>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference/[-inf inf] to [-pi pi]'
//  '<S244>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp'
//  '<S245>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem'
//  '<S246>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S247>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S248>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S249>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp/If Action Subsystem'
//  '<S250>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp/If Action Subsystem1'
//  '<S251>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem'
//  '<S252>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem1'
//  '<S253>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem'
//  '<S254>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem1'
//  '<S255>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem'
//  '<S256>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem1'
//  '<S257>' : 'AHV_Model/AH_Model4/AHV_Model'
//  '<S258>' : 'AHV_Model/AH_Model4/DP_Thrust'
//  '<S259>' : 'AHV_Model/AH_Model4/Manual Control'
//  '<S260>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1'
//  '<S261>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem'
//  '<S262>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S263>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S264>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S265>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S266>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S267>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S268>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S269>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S270>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S271>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S272>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S273>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S274>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S275>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S276>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S277>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S278>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S279>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S280>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S281>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S282>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S283>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S284>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S285>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S286>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S287>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S288>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S289>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S290>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S291>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S292>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S293>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S294>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S295>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S296>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S297>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S298>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S299>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S300>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S301>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S302>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S303>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S304>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S305>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S306>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S307>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S308>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S309>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S310>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S311>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S312>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem/Cross product'
//  '<S313>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem/Rbn_gnc'
//  '<S314>' : 'AHV_Model/AH_Model4/DP_Thrust/Calculate Cable Angles'
//  '<S315>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter'
//  '<S316>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference'
//  '<S317>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem'
//  '<S318>' : 'AHV_Model/AH_Model4/DP_Thrust/Thrust Limitations'
//  '<S319>' : 'AHV_Model/AH_Model4/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S320>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S321>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S322>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S323>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S324>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S325>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference/Chart'
//  '<S326>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference/Degrees to Radians'
//  '<S327>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference/[-inf inf] to [-pi pi]'
//  '<S328>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp'
//  '<S329>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem'
//  '<S330>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S331>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S332>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S333>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp/If Action Subsystem'
//  '<S334>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp/If Action Subsystem1'
//  '<S335>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem'
//  '<S336>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem1'
//  '<S337>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem'
//  '<S338>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem/If Action Subsystem1'
//  '<S339>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem'
//  '<S340>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem/Subsystem1/If Action Subsystem1'

#endif                                 // RTW_HEADER_AHV_Model_h_

//
// File trailer for generated code.
//
// [EOF]
//
