//
// File: AHV_Model.h
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
  extern real_T gamma_value;                 // '<Root>/gamma_value'
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
  extern boolean_T Drving_Mode1;          // '<Root>/Drving_Mode1'
  extern real_T heading_mode1;         // '<Root>/heading_mode1'
  extern boolean_T Hold_Position1;     // '<Root>/Hold_Position1'
  extern real_T Rudde_angle1;          // '<Root>/Rudde_ angle1'
  extern real_T Thruster_percentage1;  // '<Root>/Thruster_percentage1'
  extern real_T heading_angle_ref1;    // '<Root>/heading_angle_ref1'
  extern real_T Vessel_X_Ref1;         // '<Root>/Vessel_X_Ref1'
  extern real_T Vessel_Y_Ref1;         // '<Root>/Vessel_Y_Ref1'
  extern real_T traget_speed1;         // '<Root>/traget_speed1'
  extern real_T tau_cable1[3];         // '<Root>/tau_cable1'
  extern real_T ahv_fairlead1[3];      // '<Root>/ahv_fairlead1'
  extern real_T Vessel_init2[6];       // '<Root>/Vessel_init2'
  extern real_T Drving_Mode2;          // '<Root>/Drving_Mode2'
  extern real_T heading_mode2;         // '<Root>/heading_mode2'
  extern boolean_T Hold_Position2;     // '<Root>/Hold_Position2'
  extern real_T Rudde_angle2;          // '<Root>/Rudde_ angle2'
  extern real_T Thruster_percentage2;  // '<Root>/Thruster_percentage2'
  extern real_T heading_angle_ref2;    // '<Root>/heading_angle_ref2'
  extern real_T Vessel_X_Ref2;         // '<Root>/Vessel_X_Ref2'
  extern real_T Vessel_Y_Ref2;         // '<Root>/Vessel_Y_Ref2'
  extern real_T traget_speed2;         // '<Root>/traget_speed2'
  extern real_T tau_cable2[3];         // '<Root>/tau_cable2'
  extern real_T ahv_fairlead2[3];      // '<Root>/ahv_fairlead2'
  extern real_T tau_cable3[3];         // '<Root>/tau_cable3'
  extern real_T ahv_fairlead3[3];      // '<Root>/ahv_fairlead3'
  extern real_T tau_cable4[3];         // '<Root>/tau_cable4'
  extern real_T ahv_fairlead4[3];      // '<Root>/ahv_fairlead4'
  extern real_T Vessel_init3[6];       // '<Root>/Vessel_init3'
  extern real_T Drving_Mode3;          // '<Root>/Drving_Mode3'
  extern real_T heading_mode3;         // '<Root>/heading_mode3'
  extern boolean_T Hold_Position3;     // '<Root>/Hold_Position3'
  extern real_T Rudde_angle3;          // '<Root>/Rudde_ angle3'
  extern real_T Thruster_percentage3;  // '<Root>/Thruster_percentage3'
  extern real_T heading_angle_ref3;    // '<Root>/heading_angle_ref3'
  extern real_T Vessel_X_Ref3;         // '<Root>/Vessel_X_Ref3'
  extern real_T Vessel_Y_Ref3;         // '<Root>/Vessel_Y_Ref3'
  extern real_T traget_speed3;         // '<Root>/traget_speed3'
  extern real_T Vessel_init4[6];       // '<Root>/Vessel_init4'
  extern real_T Drving_Mode4;          // '<Root>/Drving_Mode4'
  extern real_T heading_mode4;         // '<Root>/heading_mode4'
  extern boolean_T Hold_Position4;     // '<Root>/Hold_Position4'
  extern real_T Rudde_angle4;          // '<Root>/Rudde_ angle4'
  extern real_T Thruster_percentage4;  // '<Root>/Thruster_percentage4'
  extern real_T heading_angle_ref4;    // '<Root>/heading_angle_ref4'
  extern real_T Vessel_X_Ref4;         // '<Root>/Vessel_X_Ref4'
  extern real_T Vessel_Y_Ref4;         // '<Root>/Vessel_Y_Ref4'
  extern real_T traget_speed4;         // '<Root>/traget_speed4'
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
  // Block signals for system '<S10>/Subsystem3'
  typedef struct {
    real_T Psi[900];
    real_T Wavenum[900];
    real_T Phase[900];
    real_T Omega[900];
    real_T Zeta_a[900];
  } B_Wave_init_T;

  // Block states (default storage) for system '<S10>/Subsystem3'
  typedef struct {
    uint32_T state[625];               // '<S16>/Wave'
  } DW_Wave_init_T;

  // Block signals for system '<S17>/Wave_loads_for_heading1'
  typedef struct {
    real_T sin_c[900];                 // '<S42>/sin'
    real_T cos_l[900];                 // '<S42>/cos'
    real_T clock_Omega[900];           // '<S37>/Product1'
    real_T Phase_tot[900];             // '<S37>/Sum'
    real_T psi_r[900];                 // '<S37>/Sum1'
  } B_Wave_loads_for_heading1_AHV_T;

  // Block states (default storage) for system '<S17>/Wave_loads_for_heading1'
  typedef struct {
    real_T Delay1_DSTATE;              // '<S43>/Delay1'
    real_T m_bpLambda[3];            // '<S46>/Lookup Table (n-D)  Wavedrift 1'
    real_T m_bpLambda_j[3];          // '<S46>/Lookup Table (n-D)  Wavedrift 2'
    real_T m_bpLambda_jj[3];         // '<S46>/Lookup Table (n-D)  Wavedrift 3'
    real_T m_bpLambda_i[3];            // '<S45>/Lookup Table (n-D) Amp1'
    real_T m_bpLambda_a[3];            // '<S45>/Lookup Table (n-D) Amp2'
    real_T m_bpLambda_h[3];            // '<S45>/Lookup Table (n-D) Amp3'
    real_T m_bpLambda_c[3];            // '<S45>/Lookup Table (n-D) Amp4'
    real_T m_bpLambda_p[3];            // '<S45>/Lookup Table (n-D) Amp5'
    real_T m_bpLambda_k[3];            // '<S45>/Lookup Table (n-D) Amp6'
    real_T m_bpLambda_cl[3];           // '<S45>/Lookup Table (n-D) Phase1'
    real_T m_bpLambda_g[3];            // '<S45>/Lookup Table (n-D) Phase2'
    real_T m_bpLambda_jb[3];           // '<S45>/Lookup Table (n-D) Phase'
    real_T m_bpLambda_d[3];            // '<S45>/Lookup Table (n-D) Phase4'
    real_T m_bpLambda_b[3];            // '<S45>/Lookup Table (n-D) Phase5'
    real_T m_bpLambda_gz[3];           // '<S45>/Lookup Table (n-D) Phase6'
    uint32_T Prelookup2_DWORK1;        // '<S44>/Prelookup2'
    uint32_T m_bpIndex[3];           // '<S46>/Lookup Table (n-D)  Wavedrift 1'
    uint32_T m_bpIndex_n[3];         // '<S46>/Lookup Table (n-D)  Wavedrift 2'
    uint32_T m_bpIndex_a[3];         // '<S46>/Lookup Table (n-D)  Wavedrift 3'
    uint32_T m_bpIndex_c[3];           // '<S45>/Lookup Table (n-D) Amp1'
    uint32_T m_bpIndex_h[3];           // '<S45>/Lookup Table (n-D) Amp2'
    uint32_T m_bpIndex_m[3];           // '<S45>/Lookup Table (n-D) Amp3'
    uint32_T m_bpIndex_d[3];           // '<S45>/Lookup Table (n-D) Amp4'
    uint32_T m_bpIndex_e[3];           // '<S45>/Lookup Table (n-D) Amp5'
    uint32_T m_bpIndex_eu[3];          // '<S45>/Lookup Table (n-D) Amp6'
    uint32_T m_bpIndex_j[3];           // '<S45>/Lookup Table (n-D) Phase1'
    uint32_T m_bpIndex_ct[3];          // '<S45>/Lookup Table (n-D) Phase2'
    uint32_T m_bpIndex_mw[3];          // '<S45>/Lookup Table (n-D) Phase'
    uint32_T m_bpIndex_g[3];           // '<S45>/Lookup Table (n-D) Phase4'
    uint32_T m_bpIndex_nq[3];          // '<S45>/Lookup Table (n-D) Phase5'
    uint32_T m_bpIndex_f[3];           // '<S45>/Lookup Table (n-D) Phase6'
  } DW_Wave_loads_for_heading1_AH_T;

  // Block signals for system '<S10>/Wave loads (U=0)'
  typedef struct {
    real_T tau_WD[6];                  // '<S51>/Sum1'
    real_T tau_WF[6];                  // '<S51>/Sum4'
    real_T tau_WD_k[6];                // '<S44>/Sum1'
    real_T tau_WF_k[6];                // '<S44>/Sum4'
    B_Wave_loads_for_heading1_AHV_T Wave_loads_for_heading2;// '<S17>/Wave_loads_for_heading2' 
    B_Wave_loads_for_heading1_AHV_T Wave_loads_for_heading1;// '<S17>/Wave_loads_for_heading1' 
  } B_Wave_loads_fun_T;

  // Block states (default storage) for system '<S10>/Wave loads (U=0)'
  typedef struct {
    DW_Wave_loads_for_heading1_AH_T Wave_loads_for_heading2;// '<S17>/Wave_loads_for_heading2' 
    DW_Wave_loads_for_heading1_AH_T Wave_loads_for_heading1;// '<S17>/Wave_loads_for_heading1' 
  } DW_Wave_loads_fun_T;

  // Block states (default storage) for system '<S66>/Chart'
  typedef struct {
    real_T x_local;                    // '<S66>/Chart'
    real_T y_local;                    // '<S66>/Chart'
    real_T yaw_local;                  // '<S66>/Chart'
    uint8_T is_active_c129_AHV_Model;  // '<S66>/Chart'
    uint8_T is_c129_AHV_Model;         // '<S66>/Chart'
  } DW_Chart_AHV_Model_T;

  // Block states (default storage) for system '<S94>/Chart1'
  typedef struct {
    real_T local;                      // '<S94>/Chart1'
    real_T yaw_local_deg;              // '<S94>/Chart1'
    uint8_T is_c1_AHV_Model;           // '<S94>/Chart1'
    uint8_T is_hold_angel;             // '<S94>/Chart1'
  } DW_Chart1_AHV_Model_T;

  // Block states (default storage) for system '<S94>/Chart3'
  typedef struct {
    real_T local;                      // '<S94>/Chart3'
    real_T yaw_local_deg;              // '<S94>/Chart3'
    uint8_T is_c3_AHV_Model;           // '<S94>/Chart3'
    uint8_T is_hold_angel;             // '<S94>/Chart3'
  } DW_Chart3_AHV_Model_T;

  // Block signals for system '<S104>/If Action Subsystem'
  typedef struct {
    real_T RateLimiter1;               // '<S190>/Rate Limiter1'
    real_T Gain4;                      // '<S190>/Gain4'
    real_T RateLimiter2;               // '<S190>/Rate Limiter2'
    real_T Product;                    // '<S190>/Product'
    real_T Product1;                   // '<S190>/Product1'
    real_T heading_deg;                // '<S190>/Chart3'
    real_T heading_deg_n;              // '<S190>/Chart2'
  } B_IfActionSubsystem_AHV_Mod_d_T;

  // Block states (default storage) for system '<S104>/If Action Subsystem'
  typedef struct {
    real_T PrevY;                      // '<S190>/Rate Limiter1'
    real_T LastMajorTime;              // '<S190>/Rate Limiter1'
    real_T PrevY_c;                    // '<S190>/Rate Limiter2'
    real_T LastMajorTime_l;            // '<S190>/Rate Limiter2'
    real_T local;                      // '<S190>/Chart1'
    real_T yaw_local_deg;              // '<S190>/Chart1'
    uint8_T is_c9_AHV_Model;           // '<S190>/Chart1'
    uint8_T is_hold_angel;             // '<S190>/Chart1'
    DW_Chart3_AHV_Model_T sf_Chart3;   // '<S190>/Chart3'
    DW_Chart1_AHV_Model_T sf_Chart2;   // '<S190>/Chart2'
  } DW_IfActionSubsystem_AHV_M_ci_T;

  // Block signals (default storage)
  typedef struct {
    real_T Integrator1[6];             // '<S15>/Integrator1'
    real_T Integrator1_b[6];           // '<S111>/Integrator1'
    real_T Integrator1_p[6];           // '<S207>/Integrator1'
    real_T Integrator1_n[6];           // '<S303>/Integrator1'
    real_T nu_r[6];                    // '<S303>/Sum6'
    real_T dx1;                        // '<S309>/dx1'
    real_T dx2;                        // '<S309>/dx2'
    real_T Switch;                     // '<S9>/Switch'
    real_T Merge2;                     // '<S8>/Merge2'
    real_T Merge5;                     // '<S8>/Merge5'
    real_T RateLimiter;                // '<S66>/Rate Limiter'
    real_T Switch_a;                   // '<S105>/Switch'
    real_T Merge2_i;                   // '<S104>/Merge2'
    real_T Merge5_g;                   // '<S104>/Merge5'
    real_T RateLimiter_b;              // '<S162>/Rate Limiter'
    real_T Switch_b;                   // '<S201>/Switch'
    real_T Merge2_h;                   // '<S200>/Merge2'
    real_T Merge5_l;                   // '<S200>/Merge5'
    real_T RateLimiter_c;              // '<S258>/Rate Limiter'
    real_T Switch_p;                   // '<S297>/Switch'
    real_T Merge7;                     // '<S296>/Merge7'
    real_T Merge2_j;                   // '<S296>/Merge2'
    real_T Merge5_o;                   // '<S296>/Merge5'
    real_T RateLimiter_m;              // '<S354>/Rate Limiter'
    real_T regulation_error[3];        // '<S355>/Sum2'
    real_T rxpi_j5;                    // '<S369>/Sum'
    real_T Row3;                       // '<S368>/Row3'
    real_T Merge;                      // '<S373>/Merge'
    real_T Merge_e;                    // '<S374>/Merge'
    real_T Ki[3];                      // '<S355>/Ki'
    real_T Merge8;                     // '<S296>/Merge8'
    real_T Merge9;                     // '<S296>/Merge9'
    real_T Minvtau[6];                 // '<S303>/Minv*tau'
    real_T nu_r_m[6];                  // '<S207>/Sum6'
    real_T dx1_a;                      // '<S213>/dx1'
    real_T dx2_c;                      // '<S213>/dx2'
    real_T Merge7_h;                   // '<S200>/Merge7'
    real_T regulation_error_d[3];      // '<S259>/Sum2'
    real_T rxpi_k;                     // '<S273>/Sum'
    real_T Row3_j;                     // '<S272>/Row3'
    real_T Merge_o;                    // '<S277>/Merge'
    real_T Merge_p;                    // '<S278>/Merge'
    real_T Ki_p[3];                    // '<S259>/Ki'
    real_T Merge8_j;                   // '<S200>/Merge8'
    real_T Merge9_e;                   // '<S200>/Merge9'
    real_T Minvtau_n[6];               // '<S207>/Minv*tau'
    real_T nu_r_j[6];                  // '<S111>/Sum6'
    real_T dx1_a0;                     // '<S117>/dx1'
    real_T dx2_f;                      // '<S117>/dx2'
    real_T Merge7_m;                   // '<S104>/Merge7'
    real_T regulation_error_k[3];      // '<S163>/Sum2'
    real_T rxpi_dyz;                   // '<S177>/Sum'
    real_T Row3_h;                     // '<S176>/Row3'
    real_T Merge_j;                    // '<S181>/Merge'
    real_T Merge_i;                    // '<S182>/Merge'
    real_T Ki_k[3];                    // '<S163>/Ki'
    real_T Merge8_m;                   // '<S104>/Merge8'
    real_T Merge9_l;                   // '<S104>/Merge9'
    real_T Minvtau_f[6];               // '<S111>/Minv*tau'
    real_T nu_r_c[6];                  // '<S15>/Sum6'
    real_T dx1_o;                      // '<S21>/dx1'
    real_T dx2_g;                      // '<S21>/dx2'
    real_T Merge7_e;                   // '<S8>/Merge7'
    real_T regulation_error_i[3];      // '<S67>/Sum2'
    real_T rxpi_e;                     // '<S81>/Sum'
    real_T Row3_l;                     // '<S80>/Row3'
    real_T Merge_k;                    // '<S85>/Merge'
    real_T Merge_b;                    // '<S86>/Merge'
    real_T Ki_m[3];                    // '<S67>/Ki'
    real_T Merge8_mo;                  // '<S8>/Merge8'
    real_T Merge9_p;                   // '<S8>/Merge9'
    real_T Minvtau_o[6];               // '<S15>/Minv*tau'
    real_T TmpSignalConversionAtIntegrat_o[6];// '<S15>/6 DOF transformation'
    real_T TmpSignalConversionAtIntegra_o3[6];// '<S111>/6 DOF transformation'
    real_T TmpSignalConversionAtIntegr_o3s[6];// '<S207>/6 DOF transformation'
    real_T TmpSignalConversionAtInteg_o3s0[6];// '<S303>/6 DOF transformation'
    real_T Sum[5];                     // '<S34>/Sum'
    real_T Sum_h[5];                   // '<S130>/Sum'
    real_T Sum_n[5];                   // '<S226>/Sum'
    real_T Sum_d[5];                   // '<S322>/Sum'
    real_T M_u[3];                     // '<S65>/Gain3'
    real_T psi_WF[3];                  // '<S65>/Sum5'
    real_T Sum6[3];                    // '<S65>/Sum6'
    real_T sun_k2[3];                  // '<S65>/Sum3'
    real_T Sum7[3];                    // '<S65>/Sum7'
    real_T M_u_o[3];                   // '<S161>/Gain3'
    real_T psi_WF_o[3];                // '<S161>/Sum5'
    real_T Sum6_d[3];                  // '<S161>/Sum6'
    real_T sun_k2_e[3];                // '<S161>/Sum3'
    real_T Sum7_o[3];                  // '<S161>/Sum7'
    real_T M_u_i[3];                   // '<S257>/Gain3'
    real_T psi_WF_m[3];                // '<S257>/Sum5'
    real_T Sum6_j[3];                  // '<S257>/Sum6'
    real_T sun_k2_e3[3];               // '<S257>/Sum3'
    real_T Sum7_n[3];                  // '<S257>/Sum7'
    real_T M_u_n[3];                   // '<S353>/Gain3'
    real_T psi_WF_i[3];                // '<S353>/Sum5'
    real_T Sum6_i[3];                  // '<S353>/Sum6'
    real_T sun_k2_a[3];                // '<S353>/Sum3'
    real_T Sum7_n4[3];                 // '<S353>/Sum7'
    real_T Fcn;                        // '<S18>/Fcn'
    real_T Fcn_j;                      // '<S114>/Fcn'
    real_T Fcn_b;                      // '<S210>/Fcn'
    real_T Fcn_e;                      // '<S306>/Fcn'
    real_T x_ref_rel;                  // '<S354>/Chart'
    real_T y_ref_rel;                  // '<S354>/Chart'
    real_T yaw_ref_rel;                // '<S354>/Chart'
    real_T Zeta_a[900];                // '<S304>/Wave'
    real_T Omega[900];                 // '<S304>/Wave'
    real_T Phase[900];                 // '<S304>/Wave'
    real_T Wavenum[900];               // '<S304>/Wave'
    real_T Psi[900];                   // '<S304>/Wave'
    real_T Sum2;                       // '<S313>/Sum2'
    real_T Sum_p;                      // '<S313>/Sum'
    real_T x_ref_rel_p;                // '<S258>/Chart'
    real_T y_ref_rel_c;                // '<S258>/Chart'
    real_T yaw_ref_rel_m;              // '<S258>/Chart'
    real_T Zeta_a_f[900];              // '<S208>/Wave'
    real_T Omega_c[900];               // '<S208>/Wave'
    real_T Phase_f[900];               // '<S208>/Wave'
    real_T Wavenum_d[900];             // '<S208>/Wave'
    real_T Psi_k[900];                 // '<S208>/Wave'
    real_T Sum2_p;                     // '<S217>/Sum2'
    real_T Sum_ps;                     // '<S217>/Sum'
    real_T x_ref_rel_l;                // '<S162>/Chart'
    real_T y_ref_rel_g;                // '<S162>/Chart'
    real_T yaw_ref_rel_g;              // '<S162>/Chart'
    real_T Zeta_a_l[900];              // '<S112>/Wave'
    real_T Omega_f[900];               // '<S112>/Wave'
    real_T Phase_i[900];               // '<S112>/Wave'
    real_T Wavenum_do[900];            // '<S112>/Wave'
    real_T Psi_b[900];                 // '<S112>/Wave'
    real_T Sum2_h;                     // '<S121>/Sum2'
    real_T Sum_nu;                     // '<S121>/Sum'
    real_T RateLimiter1;               // '<S94>/Rate Limiter1'
    real_T Gain4;                      // '<S94>/Gain4'
    real_T RateLimiter2;               // '<S94>/Rate Limiter2'
    real_T Product;                    // '<S94>/Product'
    real_T Product1;                   // '<S94>/Product1'
    real_T heading_deg_g;              // '<S94>/Chart2'
    real_T x_ref_rel_m;                // '<S66>/Chart'
    real_T y_ref_rel_d;                // '<S66>/Chart'
    real_T yaw_ref_rel_c;              // '<S66>/Chart'
    real_T Zeta_a_c[900];              // '<S16>/Wave'
    real_T Omega_d[900];               // '<S16>/Wave'
    real_T Phase_a[900];               // '<S16>/Wave'
    real_T Wavenum_b[900];             // '<S16>/Wave'
    real_T Psi_g[900];                 // '<S16>/Wave'
    real_T Sum2_n;                     // '<S25>/Sum2'
    real_T Sum_f;                      // '<S25>/Sum'
    boolean_T Merge_pf;                // '<S285>/Merge'
    boolean_T Merge_l;                 // '<S189>/Merge'
    boolean_T Merge_g;                 // '<S93>/Merge'
    boolean_T Merge_ik;                // '<S381>/Merge'
    B_IfActionSubsystem_AHV_Mod_d_T IfActionSubsystem_eo;// '<S296>/If Action Subsystem' 
    B_Wave_loads_fun_T WaveloadsU0_j;  // '<S298>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_o;        // '<S298>/Subsystem3'
    B_IfActionSubsystem_AHV_Mod_d_T IfActionSubsystem_e;// '<S200>/If Action Subsystem' 
    B_Wave_loads_fun_T WaveloadsU0_p;  // '<S202>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_h;        // '<S202>/Subsystem3'
    B_IfActionSubsystem_AHV_Mod_d_T IfActionSubsystem_c;// '<S104>/If Action Subsystem' 
    B_Wave_loads_fun_T WaveloadsU0_h;  // '<S106>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3_d;        // '<S106>/Subsystem3'
    B_Wave_loads_fun_T WaveloadsU0;    // '<S10>/Wave loads (U=0)'
    B_Wave_init_T Subsystem3;          // '<S10>/Subsystem3'
  } B_AHV_Model_T;

  // Block states (default storage) for system '<Root>'
  typedef struct {
    real_T UnitDelay_DSTATE[3];        // '<S9>/Unit Delay'
    real_T UnitDelay_DSTATE_h[3];      // '<S105>/Unit Delay'
    real_T UnitDelay_DSTATE_c[3];      // '<S201>/Unit Delay'
    real_T UnitDelay_DSTATE_a[3];      // '<S297>/Unit Delay'
    real_T Integrator1_DSTATE[3];      // '<S355>/Integrator1'
    real_T Integrator1_DSTATE_d[3];    // '<S259>/Integrator1'
    real_T Integrator1_DSTATE_p[3];    // '<S163>/Integrator1'
    real_T Integrator1_DSTATE_a[3];    // '<S67>/Integrator1'
    real_T PrevY;                      // '<S66>/Rate Limiter'
    real_T LastMajorTime;              // '<S66>/Rate Limiter'
    real_T PrevY_p;                    // '<S162>/Rate Limiter'
    real_T LastMajorTime_g;            // '<S162>/Rate Limiter'
    real_T PrevY_i;                    // '<S258>/Rate Limiter'
    real_T LastMajorTime_a;            // '<S258>/Rate Limiter'
    real_T PrevY_o;                    // '<S354>/Rate Limiter'
    real_T LastMajorTime_c;            // '<S354>/Rate Limiter'
    real_T PrevY_a;                    // '<S94>/Rate Limiter1'
    real_T LastMajorTime_e;            // '<S94>/Rate Limiter1'
    real_T PrevY_ic;                   // '<S94>/Rate Limiter2'
    real_T LastMajorTime_h;            // '<S94>/Rate Limiter2'
    int_T Integrator1_IWORK;           // '<S15>/Integrator1'
    int_T Integrator1_IWORK_o;         // '<S111>/Integrator1'
    int_T Integrator1_IWORK_b;         // '<S207>/Integrator1'
    int_T Integrator1_IWORK_c;         // '<S303>/Integrator1'
    int_T Integrator3_IWORK;           // '<S353>/Integrator3'
    int_T Integrator3_IWORK_i;         // '<S257>/Integrator3'
    int_T Integrator3_IWORK_c;         // '<S161>/Integrator3'
    int_T Integrator3_IWORK_l;         // '<S65>/Integrator3'
    int8_T If1_ActiveSubsystem;        // '<S285>/If1'
    int8_T If1_ActiveSubsystem_p;      // '<S189>/If1'
    int8_T If1_ActiveSubsystem_g;      // '<S93>/If1'
    int8_T If_ActiveSubsystem;         // '<S8>/If'
    int8_T If_ActiveSubsystem_l;       // '<S104>/If'
    int8_T If_ActiveSubsystem_a;       // '<S200>/If'
    int8_T If1_ActiveSubsystem_l;      // '<S381>/If1'
    int8_T If_ActiveSubsystem_k;       // '<S296>/If'
    int8_T If_ActiveSubsystem_i;       // '<S373>/If'
    int8_T If_ActiveSubsystem_j;       // '<S374>/If'
    int8_T If_ActiveSubsystem_d;       // '<S366>/If'
    int8_T If_ActiveSubsystem_h;       // '<S277>/If'
    int8_T If_ActiveSubsystem_m;       // '<S278>/If'
    int8_T If_ActiveSubsystem_c;       // '<S270>/If'
    int8_T If_ActiveSubsystem_e;       // '<S181>/If'
    int8_T If_ActiveSubsystem_g;       // '<S182>/If'
    int8_T If_ActiveSubsystem_dx;      // '<S174>/If'
    int8_T If_ActiveSubsystem_m5;      // '<S85>/If'
    int8_T If_ActiveSubsystem_ey;      // '<S86>/If'
    int8_T If_ActiveSubsystem_h0;      // '<S78>/If'
    DW_IfActionSubsystem_AHV_M_ci_T IfActionSubsystem_eo;// '<S296>/If Action Subsystem' 
    DW_Chart_AHV_Model_T sf_Chart_c;   // '<S354>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_j; // '<S298>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_o;       // '<S298>/Subsystem3'
    DW_IfActionSubsystem_AHV_M_ci_T IfActionSubsystem_e;// '<S200>/If Action Subsystem' 
    DW_Chart_AHV_Model_T sf_Chart_g;   // '<S258>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_p; // '<S202>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_h;       // '<S202>/Subsystem3'
    DW_IfActionSubsystem_AHV_M_ci_T IfActionSubsystem_c;// '<S104>/If Action Subsystem' 
    DW_Chart_AHV_Model_T sf_Chart_e;   // '<S162>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0_h; // '<S106>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3_d;       // '<S106>/Subsystem3'
    DW_Chart3_AHV_Model_T sf_Chart3;   // '<S94>/Chart3'
    DW_Chart1_AHV_Model_T sf_Chart2;   // '<S94>/Chart2'
    DW_Chart1_AHV_Model_T sf_Chart1;   // '<S94>/Chart1'
    DW_Chart_AHV_Model_T sf_Chart;     // '<S66>/Chart'
    DW_Wave_loads_fun_T WaveloadsU0;   // '<S10>/Wave loads (U=0)'
    DW_Wave_init_T Subsystem3;         // '<S10>/Subsystem3'
  } DW_AHV_Model_T;

  // Continuous states (default storage)
  typedef struct {
    real_T Integrator1_CSTATE[6];      // '<S15>/Integrator1'
    real_T Integrator1_CSTATE_n[6];    // '<S111>/Integrator1'
    real_T Integrator1_CSTATE_h[6];    // '<S207>/Integrator1'
    real_T Integrator1_CSTATE_hn[6];   // '<S303>/Integrator1'
    real_T Integrator_CSTATE[6];       // '<S303>/Integrator'
    real_T Integrator3_CSTATE[3];      // '<S353>/Integrator3'
    real_T Integrator4_CSTATE[3];      // '<S353>/Integrator4'
    real_T Dp11_CSTATE[5];             // '<S310>/Dp(1,1)'
    real_T Dp13_CSTATE[5];             // '<S310>/Dp(1,3)'
    real_T Dp15_CSTATE[5];             // '<S310>/Dp(1,5)'
    real_T Dp22_CSTATE[5];             // '<S310>/Dp(2,2)'
    real_T Dp24_CSTATE[5];             // '<S310>/Dp(2,4)'
    real_T Dp26_CSTATE[5];             // '<S310>/Dp(2,6)'
    real_T Dp31_CSTATE[5];             // '<S310>/Dp(3,1)'
    real_T Dp33_CSTATE[5];             // '<S310>/Dp(3,3)'
    real_T Dp35_CSTATE[5];             // '<S310>/Dp(3,5)'
    real_T Dp42_CSTATE[5];             // '<S310>/Dp(4,2)'
    real_T Integrator_CSTATE_d[5];     // '<S322>/Integrator'
    real_T Dp46_CSTATE[5];             // '<S310>/Dp(4,6)'
    real_T Dp51_CSTATE[5];             // '<S310>/Dp(5,1)'
    real_T Dp53_CSTATE[5];             // '<S310>/Dp(5,3)'
    real_T Dp55_CSTATE[5];             // '<S310>/Dp(5,5)'
    real_T Dp62_CSTATE[5];             // '<S310>/Dp(6,2)'
    real_T Dp64_CSTATE[5];             // '<S310>/Dp(6,4)'
    real_T Dp66_CSTATE[5];             // '<S310>/Dp(6,6)'
    real_T Integrator_CSTATE_h[6];     // '<S207>/Integrator'
    real_T Integrator3_CSTATE_c[3];    // '<S257>/Integrator3'
    real_T Integrator4_CSTATE_f[3];    // '<S257>/Integrator4'
    real_T Dp11_CSTATE_l[5];           // '<S214>/Dp(1,1)'
    real_T Dp13_CSTATE_c[5];           // '<S214>/Dp(1,3)'
    real_T Dp15_CSTATE_p[5];           // '<S214>/Dp(1,5)'
    real_T Dp22_CSTATE_j[5];           // '<S214>/Dp(2,2)'
    real_T Dp24_CSTATE_b[5];           // '<S214>/Dp(2,4)'
    real_T Dp26_CSTATE_n[5];           // '<S214>/Dp(2,6)'
    real_T Dp31_CSTATE_p[5];           // '<S214>/Dp(3,1)'
    real_T Dp33_CSTATE_m[5];           // '<S214>/Dp(3,3)'
    real_T Dp35_CSTATE_f[5];           // '<S214>/Dp(3,5)'
    real_T Dp42_CSTATE_b[5];           // '<S214>/Dp(4,2)'
    real_T Integrator_CSTATE_c[5];     // '<S226>/Integrator'
    real_T Dp46_CSTATE_o[5];           // '<S214>/Dp(4,6)'
    real_T Dp51_CSTATE_f[5];           // '<S214>/Dp(5,1)'
    real_T Dp53_CSTATE_i[5];           // '<S214>/Dp(5,3)'
    real_T Dp55_CSTATE_k[5];           // '<S214>/Dp(5,5)'
    real_T Dp62_CSTATE_g[5];           // '<S214>/Dp(6,2)'
    real_T Dp64_CSTATE_p[5];           // '<S214>/Dp(6,4)'
    real_T Dp66_CSTATE_f[5];           // '<S214>/Dp(6,6)'
    real_T Integrator_CSTATE_o[6];     // '<S111>/Integrator'
    real_T Integrator3_CSTATE_l[3];    // '<S161>/Integrator3'
    real_T Integrator4_CSTATE_p[3];    // '<S161>/Integrator4'
    real_T Dp11_CSTATE_o[5];           // '<S118>/Dp(1,1)'
    real_T Dp13_CSTATE_cn[5];          // '<S118>/Dp(1,3)'
    real_T Dp15_CSTATE_g[5];           // '<S118>/Dp(1,5)'
    real_T Dp22_CSTATE_jf[5];          // '<S118>/Dp(2,2)'
    real_T Dp24_CSTATE_j[5];           // '<S118>/Dp(2,4)'
    real_T Dp26_CSTATE_d[5];           // '<S118>/Dp(2,6)'
    real_T Dp31_CSTATE_j[5];           // '<S118>/Dp(3,1)'
    real_T Dp33_CSTATE_g[5];           // '<S118>/Dp(3,3)'
    real_T Dp35_CSTATE_e[5];           // '<S118>/Dp(3,5)'
    real_T Dp42_CSTATE_g[5];           // '<S118>/Dp(4,2)'
    real_T Integrator_CSTATE_k[5];     // '<S130>/Integrator'
    real_T Dp46_CSTATE_e[5];           // '<S118>/Dp(4,6)'
    real_T Dp51_CSTATE_l[5];           // '<S118>/Dp(5,1)'
    real_T Dp53_CSTATE_k[5];           // '<S118>/Dp(5,3)'
    real_T Dp55_CSTATE_i[5];           // '<S118>/Dp(5,5)'
    real_T Dp62_CSTATE_h[5];           // '<S118>/Dp(6,2)'
    real_T Dp64_CSTATE_h[5];           // '<S118>/Dp(6,4)'
    real_T Dp66_CSTATE_n[5];           // '<S118>/Dp(6,6)'
    real_T Integrator_CSTATE_e[6];     // '<S15>/Integrator'
    real_T Integrator3_CSTATE_n[3];    // '<S65>/Integrator3'
    real_T Integrator4_CSTATE_c[3];    // '<S65>/Integrator4'
    real_T Dp11_CSTATE_d[5];           // '<S22>/Dp(1,1)'
    real_T Dp13_CSTATE_d[5];           // '<S22>/Dp(1,3)'
    real_T Dp15_CSTATE_m[5];           // '<S22>/Dp(1,5)'
    real_T Dp22_CSTATE_h[5];           // '<S22>/Dp(2,2)'
    real_T Dp24_CSTATE_c[5];           // '<S22>/Dp(2,4)'
    real_T Dp26_CSTATE_h[5];           // '<S22>/Dp(2,6)'
    real_T Dp31_CSTATE_b[5];           // '<S22>/Dp(3,1)'
    real_T Dp33_CSTATE_o[5];           // '<S22>/Dp(3,3)'
    real_T Dp35_CSTATE_h[5];           // '<S22>/Dp(3,5)'
    real_T Dp42_CSTATE_e[5];           // '<S22>/Dp(4,2)'
    real_T Integrator_CSTATE_p[5];     // '<S34>/Integrator'
    real_T Dp46_CSTATE_h[5];           // '<S22>/Dp(4,6)'
    real_T Dp51_CSTATE_g[5];           // '<S22>/Dp(5,1)'
    real_T Dp53_CSTATE_m[5];           // '<S22>/Dp(5,3)'
    real_T Dp55_CSTATE_n[5];           // '<S22>/Dp(5,5)'
    real_T Dp62_CSTATE_a[5];           // '<S22>/Dp(6,2)'
    real_T Dp64_CSTATE_n[5];           // '<S22>/Dp(6,4)'
    real_T Dp66_CSTATE_fj[5];          // '<S22>/Dp(6,6)'
    real_T Integrator6_CSTATE[3];      // '<S65>/Integrator6'
    real_T Integrator1_CSTATE_a[3];    // '<S65>/Integrator1'
    real_T Integrator2_CSTATE[3];      // '<S65>/Integrator2'
    real_T Integrator6_CSTATE_p[3];    // '<S161>/Integrator6'
    real_T Integrator1_CSTATE_h0[3];   // '<S161>/Integrator1'
    real_T Integrator2_CSTATE_j[3];    // '<S161>/Integrator2'
    real_T Integrator6_CSTATE_i[3];    // '<S257>/Integrator6'
    real_T Integrator1_CSTATE_hl[3];   // '<S257>/Integrator1'
    real_T Integrator2_CSTATE_f[3];    // '<S257>/Integrator2'
    real_T Integrator6_CSTATE_ij[3];   // '<S353>/Integrator6'
    real_T Integrator1_CSTATE_nd[3];   // '<S353>/Integrator1'
    real_T Integrator2_CSTATE_p[3];    // '<S353>/Integrator2'
    real_T TransferFcn_CSTATE;         // '<S18>/Transfer Fcn'
    real_T TransferFcn_CSTATE_c;       // '<S114>/Transfer Fcn'
    real_T TransferFcn_CSTATE_d;       // '<S210>/Transfer Fcn'
    real_T TransferFcn_CSTATE_p;       // '<S306>/Transfer Fcn'
  } X_AHV_Model_T;

  // State derivatives (default storage)
  typedef struct {
    real_T Integrator1_CSTATE[6];      // '<S15>/Integrator1'
    real_T Integrator1_CSTATE_n[6];    // '<S111>/Integrator1'
    real_T Integrator1_CSTATE_h[6];    // '<S207>/Integrator1'
    real_T Integrator1_CSTATE_hn[6];   // '<S303>/Integrator1'
    real_T Integrator_CSTATE[6];       // '<S303>/Integrator'
    real_T Integrator3_CSTATE[3];      // '<S353>/Integrator3'
    real_T Integrator4_CSTATE[3];      // '<S353>/Integrator4'
    real_T Dp11_CSTATE[5];             // '<S310>/Dp(1,1)'
    real_T Dp13_CSTATE[5];             // '<S310>/Dp(1,3)'
    real_T Dp15_CSTATE[5];             // '<S310>/Dp(1,5)'
    real_T Dp22_CSTATE[5];             // '<S310>/Dp(2,2)'
    real_T Dp24_CSTATE[5];             // '<S310>/Dp(2,4)'
    real_T Dp26_CSTATE[5];             // '<S310>/Dp(2,6)'
    real_T Dp31_CSTATE[5];             // '<S310>/Dp(3,1)'
    real_T Dp33_CSTATE[5];             // '<S310>/Dp(3,3)'
    real_T Dp35_CSTATE[5];             // '<S310>/Dp(3,5)'
    real_T Dp42_CSTATE[5];             // '<S310>/Dp(4,2)'
    real_T Integrator_CSTATE_d[5];     // '<S322>/Integrator'
    real_T Dp46_CSTATE[5];             // '<S310>/Dp(4,6)'
    real_T Dp51_CSTATE[5];             // '<S310>/Dp(5,1)'
    real_T Dp53_CSTATE[5];             // '<S310>/Dp(5,3)'
    real_T Dp55_CSTATE[5];             // '<S310>/Dp(5,5)'
    real_T Dp62_CSTATE[5];             // '<S310>/Dp(6,2)'
    real_T Dp64_CSTATE[5];             // '<S310>/Dp(6,4)'
    real_T Dp66_CSTATE[5];             // '<S310>/Dp(6,6)'
    real_T Integrator_CSTATE_h[6];     // '<S207>/Integrator'
    real_T Integrator3_CSTATE_c[3];    // '<S257>/Integrator3'
    real_T Integrator4_CSTATE_f[3];    // '<S257>/Integrator4'
    real_T Dp11_CSTATE_l[5];           // '<S214>/Dp(1,1)'
    real_T Dp13_CSTATE_c[5];           // '<S214>/Dp(1,3)'
    real_T Dp15_CSTATE_p[5];           // '<S214>/Dp(1,5)'
    real_T Dp22_CSTATE_j[5];           // '<S214>/Dp(2,2)'
    real_T Dp24_CSTATE_b[5];           // '<S214>/Dp(2,4)'
    real_T Dp26_CSTATE_n[5];           // '<S214>/Dp(2,6)'
    real_T Dp31_CSTATE_p[5];           // '<S214>/Dp(3,1)'
    real_T Dp33_CSTATE_m[5];           // '<S214>/Dp(3,3)'
    real_T Dp35_CSTATE_f[5];           // '<S214>/Dp(3,5)'
    real_T Dp42_CSTATE_b[5];           // '<S214>/Dp(4,2)'
    real_T Integrator_CSTATE_c[5];     // '<S226>/Integrator'
    real_T Dp46_CSTATE_o[5];           // '<S214>/Dp(4,6)'
    real_T Dp51_CSTATE_f[5];           // '<S214>/Dp(5,1)'
    real_T Dp53_CSTATE_i[5];           // '<S214>/Dp(5,3)'
    real_T Dp55_CSTATE_k[5];           // '<S214>/Dp(5,5)'
    real_T Dp62_CSTATE_g[5];           // '<S214>/Dp(6,2)'
    real_T Dp64_CSTATE_p[5];           // '<S214>/Dp(6,4)'
    real_T Dp66_CSTATE_f[5];           // '<S214>/Dp(6,6)'
    real_T Integrator_CSTATE_o[6];     // '<S111>/Integrator'
    real_T Integrator3_CSTATE_l[3];    // '<S161>/Integrator3'
    real_T Integrator4_CSTATE_p[3];    // '<S161>/Integrator4'
    real_T Dp11_CSTATE_o[5];           // '<S118>/Dp(1,1)'
    real_T Dp13_CSTATE_cn[5];          // '<S118>/Dp(1,3)'
    real_T Dp15_CSTATE_g[5];           // '<S118>/Dp(1,5)'
    real_T Dp22_CSTATE_jf[5];          // '<S118>/Dp(2,2)'
    real_T Dp24_CSTATE_j[5];           // '<S118>/Dp(2,4)'
    real_T Dp26_CSTATE_d[5];           // '<S118>/Dp(2,6)'
    real_T Dp31_CSTATE_j[5];           // '<S118>/Dp(3,1)'
    real_T Dp33_CSTATE_g[5];           // '<S118>/Dp(3,3)'
    real_T Dp35_CSTATE_e[5];           // '<S118>/Dp(3,5)'
    real_T Dp42_CSTATE_g[5];           // '<S118>/Dp(4,2)'
    real_T Integrator_CSTATE_k[5];     // '<S130>/Integrator'
    real_T Dp46_CSTATE_e[5];           // '<S118>/Dp(4,6)'
    real_T Dp51_CSTATE_l[5];           // '<S118>/Dp(5,1)'
    real_T Dp53_CSTATE_k[5];           // '<S118>/Dp(5,3)'
    real_T Dp55_CSTATE_i[5];           // '<S118>/Dp(5,5)'
    real_T Dp62_CSTATE_h[5];           // '<S118>/Dp(6,2)'
    real_T Dp64_CSTATE_h[5];           // '<S118>/Dp(6,4)'
    real_T Dp66_CSTATE_n[5];           // '<S118>/Dp(6,6)'
    real_T Integrator_CSTATE_e[6];     // '<S15>/Integrator'
    real_T Integrator3_CSTATE_n[3];    // '<S65>/Integrator3'
    real_T Integrator4_CSTATE_c[3];    // '<S65>/Integrator4'
    real_T Dp11_CSTATE_d[5];           // '<S22>/Dp(1,1)'
    real_T Dp13_CSTATE_d[5];           // '<S22>/Dp(1,3)'
    real_T Dp15_CSTATE_m[5];           // '<S22>/Dp(1,5)'
    real_T Dp22_CSTATE_h[5];           // '<S22>/Dp(2,2)'
    real_T Dp24_CSTATE_c[5];           // '<S22>/Dp(2,4)'
    real_T Dp26_CSTATE_h[5];           // '<S22>/Dp(2,6)'
    real_T Dp31_CSTATE_b[5];           // '<S22>/Dp(3,1)'
    real_T Dp33_CSTATE_o[5];           // '<S22>/Dp(3,3)'
    real_T Dp35_CSTATE_h[5];           // '<S22>/Dp(3,5)'
    real_T Dp42_CSTATE_e[5];           // '<S22>/Dp(4,2)'
    real_T Integrator_CSTATE_p[5];     // '<S34>/Integrator'
    real_T Dp46_CSTATE_h[5];           // '<S22>/Dp(4,6)'
    real_T Dp51_CSTATE_g[5];           // '<S22>/Dp(5,1)'
    real_T Dp53_CSTATE_m[5];           // '<S22>/Dp(5,3)'
    real_T Dp55_CSTATE_n[5];           // '<S22>/Dp(5,5)'
    real_T Dp62_CSTATE_a[5];           // '<S22>/Dp(6,2)'
    real_T Dp64_CSTATE_n[5];           // '<S22>/Dp(6,4)'
    real_T Dp66_CSTATE_fj[5];          // '<S22>/Dp(6,6)'
    real_T Integrator6_CSTATE[3];      // '<S65>/Integrator6'
    real_T Integrator1_CSTATE_a[3];    // '<S65>/Integrator1'
    real_T Integrator2_CSTATE[3];      // '<S65>/Integrator2'
    real_T Integrator6_CSTATE_p[3];    // '<S161>/Integrator6'
    real_T Integrator1_CSTATE_h0[3];   // '<S161>/Integrator1'
    real_T Integrator2_CSTATE_j[3];    // '<S161>/Integrator2'
    real_T Integrator6_CSTATE_i[3];    // '<S257>/Integrator6'
    real_T Integrator1_CSTATE_hl[3];   // '<S257>/Integrator1'
    real_T Integrator2_CSTATE_f[3];    // '<S257>/Integrator2'
    real_T Integrator6_CSTATE_ij[3];   // '<S353>/Integrator6'
    real_T Integrator1_CSTATE_nd[3];   // '<S353>/Integrator1'
    real_T Integrator2_CSTATE_p[3];    // '<S353>/Integrator2'
    real_T TransferFcn_CSTATE;         // '<S18>/Transfer Fcn'
    real_T TransferFcn_CSTATE_c;       // '<S114>/Transfer Fcn'
    real_T TransferFcn_CSTATE_d;       // '<S210>/Transfer Fcn'
    real_T TransferFcn_CSTATE_p;       // '<S306>/Transfer Fcn'
  } XDot_AHV_Model_T;

  // State disabled
  typedef struct {
    boolean_T Integrator1_CSTATE[6];   // '<S15>/Integrator1'
    boolean_T Integrator1_CSTATE_n[6]; // '<S111>/Integrator1'
    boolean_T Integrator1_CSTATE_h[6]; // '<S207>/Integrator1'
    boolean_T Integrator1_CSTATE_hn[6];// '<S303>/Integrator1'
    boolean_T Integrator_CSTATE[6];    // '<S303>/Integrator'
    boolean_T Integrator3_CSTATE[3];   // '<S353>/Integrator3'
    boolean_T Integrator4_CSTATE[3];   // '<S353>/Integrator4'
    boolean_T Dp11_CSTATE[5];          // '<S310>/Dp(1,1)'
    boolean_T Dp13_CSTATE[5];          // '<S310>/Dp(1,3)'
    boolean_T Dp15_CSTATE[5];          // '<S310>/Dp(1,5)'
    boolean_T Dp22_CSTATE[5];          // '<S310>/Dp(2,2)'
    boolean_T Dp24_CSTATE[5];          // '<S310>/Dp(2,4)'
    boolean_T Dp26_CSTATE[5];          // '<S310>/Dp(2,6)'
    boolean_T Dp31_CSTATE[5];          // '<S310>/Dp(3,1)'
    boolean_T Dp33_CSTATE[5];          // '<S310>/Dp(3,3)'
    boolean_T Dp35_CSTATE[5];          // '<S310>/Dp(3,5)'
    boolean_T Dp42_CSTATE[5];          // '<S310>/Dp(4,2)'
    boolean_T Integrator_CSTATE_d[5];  // '<S322>/Integrator'
    boolean_T Dp46_CSTATE[5];          // '<S310>/Dp(4,6)'
    boolean_T Dp51_CSTATE[5];          // '<S310>/Dp(5,1)'
    boolean_T Dp53_CSTATE[5];          // '<S310>/Dp(5,3)'
    boolean_T Dp55_CSTATE[5];          // '<S310>/Dp(5,5)'
    boolean_T Dp62_CSTATE[5];          // '<S310>/Dp(6,2)'
    boolean_T Dp64_CSTATE[5];          // '<S310>/Dp(6,4)'
    boolean_T Dp66_CSTATE[5];          // '<S310>/Dp(6,6)'
    boolean_T Integrator_CSTATE_h[6];  // '<S207>/Integrator'
    boolean_T Integrator3_CSTATE_c[3]; // '<S257>/Integrator3'
    boolean_T Integrator4_CSTATE_f[3]; // '<S257>/Integrator4'
    boolean_T Dp11_CSTATE_l[5];        // '<S214>/Dp(1,1)'
    boolean_T Dp13_CSTATE_c[5];        // '<S214>/Dp(1,3)'
    boolean_T Dp15_CSTATE_p[5];        // '<S214>/Dp(1,5)'
    boolean_T Dp22_CSTATE_j[5];        // '<S214>/Dp(2,2)'
    boolean_T Dp24_CSTATE_b[5];        // '<S214>/Dp(2,4)'
    boolean_T Dp26_CSTATE_n[5];        // '<S214>/Dp(2,6)'
    boolean_T Dp31_CSTATE_p[5];        // '<S214>/Dp(3,1)'
    boolean_T Dp33_CSTATE_m[5];        // '<S214>/Dp(3,3)'
    boolean_T Dp35_CSTATE_f[5];        // '<S214>/Dp(3,5)'
    boolean_T Dp42_CSTATE_b[5];        // '<S214>/Dp(4,2)'
    boolean_T Integrator_CSTATE_c[5];  // '<S226>/Integrator'
    boolean_T Dp46_CSTATE_o[5];        // '<S214>/Dp(4,6)'
    boolean_T Dp51_CSTATE_f[5];        // '<S214>/Dp(5,1)'
    boolean_T Dp53_CSTATE_i[5];        // '<S214>/Dp(5,3)'
    boolean_T Dp55_CSTATE_k[5];        // '<S214>/Dp(5,5)'
    boolean_T Dp62_CSTATE_g[5];        // '<S214>/Dp(6,2)'
    boolean_T Dp64_CSTATE_p[5];        // '<S214>/Dp(6,4)'
    boolean_T Dp66_CSTATE_f[5];        // '<S214>/Dp(6,6)'
    boolean_T Integrator_CSTATE_o[6];  // '<S111>/Integrator'
    boolean_T Integrator3_CSTATE_l[3]; // '<S161>/Integrator3'
    boolean_T Integrator4_CSTATE_p[3]; // '<S161>/Integrator4'
    boolean_T Dp11_CSTATE_o[5];        // '<S118>/Dp(1,1)'
    boolean_T Dp13_CSTATE_cn[5];       // '<S118>/Dp(1,3)'
    boolean_T Dp15_CSTATE_g[5];        // '<S118>/Dp(1,5)'
    boolean_T Dp22_CSTATE_jf[5];       // '<S118>/Dp(2,2)'
    boolean_T Dp24_CSTATE_j[5];        // '<S118>/Dp(2,4)'
    boolean_T Dp26_CSTATE_d[5];        // '<S118>/Dp(2,6)'
    boolean_T Dp31_CSTATE_j[5];        // '<S118>/Dp(3,1)'
    boolean_T Dp33_CSTATE_g[5];        // '<S118>/Dp(3,3)'
    boolean_T Dp35_CSTATE_e[5];        // '<S118>/Dp(3,5)'
    boolean_T Dp42_CSTATE_g[5];        // '<S118>/Dp(4,2)'
    boolean_T Integrator_CSTATE_k[5];  // '<S130>/Integrator'
    boolean_T Dp46_CSTATE_e[5];        // '<S118>/Dp(4,6)'
    boolean_T Dp51_CSTATE_l[5];        // '<S118>/Dp(5,1)'
    boolean_T Dp53_CSTATE_k[5];        // '<S118>/Dp(5,3)'
    boolean_T Dp55_CSTATE_i[5];        // '<S118>/Dp(5,5)'
    boolean_T Dp62_CSTATE_h[5];        // '<S118>/Dp(6,2)'
    boolean_T Dp64_CSTATE_h[5];        // '<S118>/Dp(6,4)'
    boolean_T Dp66_CSTATE_n[5];        // '<S118>/Dp(6,6)'
    boolean_T Integrator_CSTATE_e[6];  // '<S15>/Integrator'
    boolean_T Integrator3_CSTATE_n[3]; // '<S65>/Integrator3'
    boolean_T Integrator4_CSTATE_c[3]; // '<S65>/Integrator4'
    boolean_T Dp11_CSTATE_d[5];        // '<S22>/Dp(1,1)'
    boolean_T Dp13_CSTATE_d[5];        // '<S22>/Dp(1,3)'
    boolean_T Dp15_CSTATE_m[5];        // '<S22>/Dp(1,5)'
    boolean_T Dp22_CSTATE_h[5];        // '<S22>/Dp(2,2)'
    boolean_T Dp24_CSTATE_c[5];        // '<S22>/Dp(2,4)'
    boolean_T Dp26_CSTATE_h[5];        // '<S22>/Dp(2,6)'
    boolean_T Dp31_CSTATE_b[5];        // '<S22>/Dp(3,1)'
    boolean_T Dp33_CSTATE_o[5];        // '<S22>/Dp(3,3)'
    boolean_T Dp35_CSTATE_h[5];        // '<S22>/Dp(3,5)'
    boolean_T Dp42_CSTATE_e[5];        // '<S22>/Dp(4,2)'
    boolean_T Integrator_CSTATE_p[5];  // '<S34>/Integrator'
    boolean_T Dp46_CSTATE_h[5];        // '<S22>/Dp(4,6)'
    boolean_T Dp51_CSTATE_g[5];        // '<S22>/Dp(5,1)'
    boolean_T Dp53_CSTATE_m[5];        // '<S22>/Dp(5,3)'
    boolean_T Dp55_CSTATE_n[5];        // '<S22>/Dp(5,5)'
    boolean_T Dp62_CSTATE_a[5];        // '<S22>/Dp(6,2)'
    boolean_T Dp64_CSTATE_n[5];        // '<S22>/Dp(6,4)'
    boolean_T Dp66_CSTATE_fj[5];       // '<S22>/Dp(6,6)'
    boolean_T Integrator6_CSTATE[3];   // '<S65>/Integrator6'
    boolean_T Integrator1_CSTATE_a[3]; // '<S65>/Integrator1'
    boolean_T Integrator2_CSTATE[3];   // '<S65>/Integrator2'
    boolean_T Integrator6_CSTATE_p[3]; // '<S161>/Integrator6'
    boolean_T Integrator1_CSTATE_h0[3];// '<S161>/Integrator1'
    boolean_T Integrator2_CSTATE_j[3]; // '<S161>/Integrator2'
    boolean_T Integrator6_CSTATE_i[3]; // '<S257>/Integrator6'
    boolean_T Integrator1_CSTATE_hl[3];// '<S257>/Integrator1'
    boolean_T Integrator2_CSTATE_f[3]; // '<S257>/Integrator2'
    boolean_T Integrator6_CSTATE_ij[3];// '<S353>/Integrator6'
    boolean_T Integrator1_CSTATE_nd[3];// '<S353>/Integrator1'
    boolean_T Integrator2_CSTATE_p[3]; // '<S353>/Integrator2'
    boolean_T TransferFcn_CSTATE;      // '<S18>/Transfer Fcn'
    boolean_T TransferFcn_CSTATE_c;    // '<S114>/Transfer Fcn'
    boolean_T TransferFcn_CSTATE_d;    // '<S210>/Transfer Fcn'
    boolean_T TransferFcn_CSTATE_p;    // '<S306>/Transfer Fcn'
  } XDis_AHV_Model_T;

  // Invariant block signals for system '<S17>/Wave_loads_for_heading1'
  typedef const struct tag_ConstB_Wave_loads_for_hea_T {
    real_T number_of_iterations;       // '<S37>/Width '
    real_T U;                     // '<S44>/Direct Lookup Table (n-D) Velocity'
    uint32_T Prelookup2;               // '<S44>/Prelookup2'
  } ConstB_Wave_loads_for_heading_T;

  // Invariant block signals for system '<S10>/Wave loads (U=0)'
  typedef const struct tag_ConstB_Wave_loads_fun_T {
    ConstB_Wave_loads_for_heading_T Wave_loads_for_heading2;// '<S17>/Wave_loads_for_heading2' 
    ConstB_Wave_loads_for_heading_T Wave_loads_for_heading1;// '<S17>/Wave_loads_for_heading1' 
  } ConstB_Wave_loads_fun_T;

  // Invariant block signals (default storage)
  typedef const struct tag_ConstB_AHV_Model_T {
    real_T N1;                         // '<S21>/Sum'
    real_T dx;                         // '<S21>/dx'
    real_T Gain2;                      // '<S21>/Gain2'
    real_T Product5;                   // '<S21>/Product5'
    real_T Sum2[36];                   // '<S15>/Sum2'
    real_T MathFunction[5];            // '<S34>/Math Function'
    real_T N1_f;                       // '<S117>/Sum'
    real_T dx_o;                       // '<S117>/dx'
    real_T Gain2_o;                    // '<S117>/Gain2'
    real_T Product5_h;                 // '<S117>/Product5'
    real_T Sum2_p[36];                 // '<S111>/Sum2'
    real_T MathFunction_l[5];          // '<S130>/Math Function'
    real_T N1_l;                       // '<S213>/Sum'
    real_T dx_g;                       // '<S213>/dx'
    real_T Gain2_h;                    // '<S213>/Gain2'
    real_T Product5_f;                 // '<S213>/Product5'
    real_T Sum2_g[36];                 // '<S207>/Sum2'
    real_T MathFunction_n[5];          // '<S226>/Math Function'
    real_T N1_lr;                      // '<S309>/Sum'
    real_T dx_f;                       // '<S309>/dx'
    real_T Gain2_hl;                   // '<S309>/Gain2'
    real_T Product5_hj;                // '<S309>/Product5'
    real_T Sum2_a[36];                 // '<S303>/Sum2'
    real_T MathFunction_e[5];          // '<S322>/Math Function'
    ConstB_Wave_loads_fun_T WaveloadsU0_j;// '<S298>/Wave loads (U=0)'
    ConstB_Wave_loads_fun_T WaveloadsU0_p;// '<S202>/Wave loads (U=0)'
    ConstB_Wave_loads_fun_T WaveloadsU0_h;// '<S106>/Wave loads (U=0)'
    ConstB_Wave_loads_fun_T WaveloadsU0;// '<S10>/Wave loads (U=0)'
  } ConstB_AHV_Model_T;

  // Constant parameters (default storage)
  typedef struct {
    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S78>/If Action Subsystem1'
    //    '<S174>/If Action Subsystem1'
    //    '<S270>/If Action Subsystem1'
    //    '<S366>/If Action Subsystem1'

    real_T pooled1[9];

    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S78>/If Action Subsystem'
    //    '<S174>/If Action Subsystem'
    //    '<S270>/If Action Subsystem'
    //    '<S366>/If Action Subsystem'

    real_T pooled2[9];

    // Pooled Parameter (Expression: phase_vector)
    //  Referenced by:
    //    '<S16>/Constant'
    //    '<S112>/Constant'
    //    '<S208>/Constant'
    //    '<S304>/Constant'

    real_T pooled9[200];

    // Pooled Parameter (Expression: psi_vector)
    //  Referenced by:
    //    '<S16>/Constant1'
    //    '<S112>/Constant1'
    //    '<S208>/Constant1'
    //    '<S304>/Constant1'

    real_T pooled10[200];

    // Pooled Parameter (Expression: vessel_dp.headings)
    //  Referenced by:
    //    '<S44>/Prelookup1'
    //    '<S44>/Direct Lookup Table (n-D) Beta'
    //    '<S51>/Prelookup1'
    //    '<S51>/Direct Lookup Table (n-D) Beta'
    //    '<S140>/Prelookup1'
    //    '<S140>/Direct Lookup Table (n-D) Beta'
    //    '<S147>/Prelookup1'
    //    '<S147>/Direct Lookup Table (n-D) Beta'
    //    '<S236>/Prelookup1'
    //    '<S236>/Direct Lookup Table (n-D) Beta'
    //    '<S243>/Prelookup1'
    //    '<S243>/Direct Lookup Table (n-D) Beta'
    //    '<S332>/Prelookup1'
    //    '<S332>/Direct Lookup Table (n-D) Beta'
    //    '<S339>/Prelookup1'
    //    '<S339>/Direct Lookup Table (n-D) Beta'
    //    '<S45>/Lookup Table (n-D) Amp1'
    //    '<S45>/Lookup Table (n-D) Amp2'
    //    '<S45>/Lookup Table (n-D) Amp3'
    //    '<S45>/Lookup Table (n-D) Amp4'
    //    '<S45>/Lookup Table (n-D) Amp5'
    //    '<S45>/Lookup Table (n-D) Amp6'
    //    '<S45>/Lookup Table (n-D) Phase'
    //    '<S45>/Lookup Table (n-D) Phase1'
    //    '<S45>/Lookup Table (n-D) Phase2'
    //    '<S45>/Lookup Table (n-D) Phase4'
    //    '<S45>/Lookup Table (n-D) Phase5'
    //    '<S45>/Lookup Table (n-D) Phase6'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S52>/Lookup Table (n-D) Amp1'
    //    '<S52>/Lookup Table (n-D) Amp2'
    //    '<S52>/Lookup Table (n-D) Amp3'
    //    '<S52>/Lookup Table (n-D) Amp4'
    //    '<S52>/Lookup Table (n-D) Amp5'
    //    '<S52>/Lookup Table (n-D) Amp6'
    //    '<S52>/Lookup Table (n-D) Phase'
    //    '<S52>/Lookup Table (n-D) Phase1'
    //    '<S52>/Lookup Table (n-D) Phase2'
    //    '<S52>/Lookup Table (n-D) Phase4'
    //    '<S52>/Lookup Table (n-D) Phase5'
    //    '<S52>/Lookup Table (n-D) Phase6'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S141>/Lookup Table (n-D) Amp1'
    //    '<S141>/Lookup Table (n-D) Amp2'
    //    '<S141>/Lookup Table (n-D) Amp3'
    //    '<S141>/Lookup Table (n-D) Amp4'
    //    '<S141>/Lookup Table (n-D) Amp5'
    //    '<S141>/Lookup Table (n-D) Amp6'
    //    '<S141>/Lookup Table (n-D) Phase'
    //    '<S141>/Lookup Table (n-D) Phase1'
    //    '<S141>/Lookup Table (n-D) Phase2'
    //    '<S141>/Lookup Table (n-D) Phase4'
    //    '<S141>/Lookup Table (n-D) Phase5'
    //    '<S141>/Lookup Table (n-D) Phase6'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S148>/Lookup Table (n-D) Amp1'
    //    '<S148>/Lookup Table (n-D) Amp2'
    //    '<S148>/Lookup Table (n-D) Amp3'
    //    '<S148>/Lookup Table (n-D) Amp4'
    //    '<S148>/Lookup Table (n-D) Amp5'
    //    '<S148>/Lookup Table (n-D) Amp6'
    //    '<S148>/Lookup Table (n-D) Phase'
    //    '<S148>/Lookup Table (n-D) Phase1'
    //    '<S148>/Lookup Table (n-D) Phase2'
    //    '<S148>/Lookup Table (n-D) Phase4'
    //    '<S148>/Lookup Table (n-D) Phase5'
    //    '<S148>/Lookup Table (n-D) Phase6'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S237>/Lookup Table (n-D) Amp1'
    //    '<S237>/Lookup Table (n-D) Amp2'
    //    '<S237>/Lookup Table (n-D) Amp3'
    //    '<S237>/Lookup Table (n-D) Amp4'
    //    '<S237>/Lookup Table (n-D) Amp5'
    //    '<S237>/Lookup Table (n-D) Amp6'
    //    '<S237>/Lookup Table (n-D) Phase'
    //    '<S237>/Lookup Table (n-D) Phase1'
    //    '<S237>/Lookup Table (n-D) Phase2'
    //    '<S237>/Lookup Table (n-D) Phase4'
    //    '<S237>/Lookup Table (n-D) Phase5'
    //    '<S237>/Lookup Table (n-D) Phase6'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S244>/Lookup Table (n-D) Amp1'
    //    '<S244>/Lookup Table (n-D) Amp2'
    //    '<S244>/Lookup Table (n-D) Amp3'
    //    '<S244>/Lookup Table (n-D) Amp4'
    //    '<S244>/Lookup Table (n-D) Amp5'
    //    '<S244>/Lookup Table (n-D) Amp6'
    //    '<S244>/Lookup Table (n-D) Phase'
    //    '<S244>/Lookup Table (n-D) Phase1'
    //    '<S244>/Lookup Table (n-D) Phase2'
    //    '<S244>/Lookup Table (n-D) Phase4'
    //    '<S244>/Lookup Table (n-D) Phase5'
    //    '<S244>/Lookup Table (n-D) Phase6'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S333>/Lookup Table (n-D) Amp1'
    //    '<S333>/Lookup Table (n-D) Amp2'
    //    '<S333>/Lookup Table (n-D) Amp3'
    //    '<S333>/Lookup Table (n-D) Amp4'
    //    '<S333>/Lookup Table (n-D) Amp5'
    //    '<S333>/Lookup Table (n-D) Amp6'
    //    '<S333>/Lookup Table (n-D) Phase'
    //    '<S333>/Lookup Table (n-D) Phase1'
    //    '<S333>/Lookup Table (n-D) Phase2'
    //    '<S333>/Lookup Table (n-D) Phase4'
    //    '<S333>/Lookup Table (n-D) Phase5'
    //    '<S333>/Lookup Table (n-D) Phase6'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S340>/Lookup Table (n-D) Amp1'
    //    '<S340>/Lookup Table (n-D) Amp2'
    //    '<S340>/Lookup Table (n-D) Amp3'
    //    '<S340>/Lookup Table (n-D) Amp4'
    //    '<S340>/Lookup Table (n-D) Amp5'
    //    '<S340>/Lookup Table (n-D) Amp6'
    //    '<S340>/Lookup Table (n-D) Phase'
    //    '<S340>/Lookup Table (n-D) Phase1'
    //    '<S340>/Lookup Table (n-D) Phase2'
    //    '<S340>/Lookup Table (n-D) Phase4'
    //    '<S340>/Lookup Table (n-D) Phase5'
    //    '<S340>/Lookup Table (n-D) Phase6'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled16[36];

    // Pooled Parameter (Expression: WD_Amp{1})
    //  Referenced by:
    //    '<S46>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 1'

    real_T pooled17[2808];

    // Pooled Parameter (Expression: WD_Freq)
    //  Referenced by:
    //    '<S46>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled18[39];

    // Pooled Parameter (Expression: RAO_Vel)
    //  Referenced by:
    //    '<S44>/Prelookup2'
    //    '<S44>/Direct Lookup Table (n-D) Velocity'
    //    '<S51>/Prelookup2'
    //    '<S51>/Direct Lookup Table (n-D) Velocity'
    //    '<S140>/Prelookup2'
    //    '<S140>/Direct Lookup Table (n-D) Velocity'
    //    '<S147>/Prelookup2'
    //    '<S147>/Direct Lookup Table (n-D) Velocity'
    //    '<S236>/Prelookup2'
    //    '<S236>/Direct Lookup Table (n-D) Velocity'
    //    '<S243>/Prelookup2'
    //    '<S243>/Direct Lookup Table (n-D) Velocity'
    //    '<S332>/Prelookup2'
    //    '<S332>/Direct Lookup Table (n-D) Velocity'
    //    '<S339>/Prelookup2'
    //    '<S339>/Direct Lookup Table (n-D) Velocity'
    //    '<S45>/Lookup Table (n-D) Amp1'
    //    '<S45>/Lookup Table (n-D) Amp2'
    //    '<S45>/Lookup Table (n-D) Amp3'
    //    '<S45>/Lookup Table (n-D) Amp4'
    //    '<S45>/Lookup Table (n-D) Amp5'
    //    '<S45>/Lookup Table (n-D) Amp6'
    //    '<S45>/Lookup Table (n-D) Phase'
    //    '<S45>/Lookup Table (n-D) Phase1'
    //    '<S45>/Lookup Table (n-D) Phase2'
    //    '<S45>/Lookup Table (n-D) Phase4'
    //    '<S45>/Lookup Table (n-D) Phase5'
    //    '<S45>/Lookup Table (n-D) Phase6'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S52>/Lookup Table (n-D) Amp1'
    //    '<S52>/Lookup Table (n-D) Amp2'
    //    '<S52>/Lookup Table (n-D) Amp3'
    //    '<S52>/Lookup Table (n-D) Amp4'
    //    '<S52>/Lookup Table (n-D) Amp5'
    //    '<S52>/Lookup Table (n-D) Amp6'
    //    '<S52>/Lookup Table (n-D) Phase'
    //    '<S52>/Lookup Table (n-D) Phase1'
    //    '<S52>/Lookup Table (n-D) Phase2'
    //    '<S52>/Lookup Table (n-D) Phase4'
    //    '<S52>/Lookup Table (n-D) Phase5'
    //    '<S52>/Lookup Table (n-D) Phase6'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S141>/Lookup Table (n-D) Amp1'
    //    '<S141>/Lookup Table (n-D) Amp2'
    //    '<S141>/Lookup Table (n-D) Amp3'
    //    '<S141>/Lookup Table (n-D) Amp4'
    //    '<S141>/Lookup Table (n-D) Amp5'
    //    '<S141>/Lookup Table (n-D) Amp6'
    //    '<S141>/Lookup Table (n-D) Phase'
    //    '<S141>/Lookup Table (n-D) Phase1'
    //    '<S141>/Lookup Table (n-D) Phase2'
    //    '<S141>/Lookup Table (n-D) Phase4'
    //    '<S141>/Lookup Table (n-D) Phase5'
    //    '<S141>/Lookup Table (n-D) Phase6'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S148>/Lookup Table (n-D) Amp1'
    //    '<S148>/Lookup Table (n-D) Amp2'
    //    '<S148>/Lookup Table (n-D) Amp3'
    //    '<S148>/Lookup Table (n-D) Amp4'
    //    '<S148>/Lookup Table (n-D) Amp5'
    //    '<S148>/Lookup Table (n-D) Amp6'
    //    '<S148>/Lookup Table (n-D) Phase'
    //    '<S148>/Lookup Table (n-D) Phase1'
    //    '<S148>/Lookup Table (n-D) Phase2'
    //    '<S148>/Lookup Table (n-D) Phase4'
    //    '<S148>/Lookup Table (n-D) Phase5'
    //    '<S148>/Lookup Table (n-D) Phase6'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S237>/Lookup Table (n-D) Amp1'
    //    '<S237>/Lookup Table (n-D) Amp2'
    //    '<S237>/Lookup Table (n-D) Amp3'
    //    '<S237>/Lookup Table (n-D) Amp4'
    //    '<S237>/Lookup Table (n-D) Amp5'
    //    '<S237>/Lookup Table (n-D) Amp6'
    //    '<S237>/Lookup Table (n-D) Phase'
    //    '<S237>/Lookup Table (n-D) Phase1'
    //    '<S237>/Lookup Table (n-D) Phase2'
    //    '<S237>/Lookup Table (n-D) Phase4'
    //    '<S237>/Lookup Table (n-D) Phase5'
    //    '<S237>/Lookup Table (n-D) Phase6'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S244>/Lookup Table (n-D) Amp1'
    //    '<S244>/Lookup Table (n-D) Amp2'
    //    '<S244>/Lookup Table (n-D) Amp3'
    //    '<S244>/Lookup Table (n-D) Amp4'
    //    '<S244>/Lookup Table (n-D) Amp5'
    //    '<S244>/Lookup Table (n-D) Amp6'
    //    '<S244>/Lookup Table (n-D) Phase'
    //    '<S244>/Lookup Table (n-D) Phase1'
    //    '<S244>/Lookup Table (n-D) Phase2'
    //    '<S244>/Lookup Table (n-D) Phase4'
    //    '<S244>/Lookup Table (n-D) Phase5'
    //    '<S244>/Lookup Table (n-D) Phase6'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S333>/Lookup Table (n-D) Amp1'
    //    '<S333>/Lookup Table (n-D) Amp2'
    //    '<S333>/Lookup Table (n-D) Amp3'
    //    '<S333>/Lookup Table (n-D) Amp4'
    //    '<S333>/Lookup Table (n-D) Amp5'
    //    '<S333>/Lookup Table (n-D) Amp6'
    //    '<S333>/Lookup Table (n-D) Phase'
    //    '<S333>/Lookup Table (n-D) Phase1'
    //    '<S333>/Lookup Table (n-D) Phase2'
    //    '<S333>/Lookup Table (n-D) Phase4'
    //    '<S333>/Lookup Table (n-D) Phase5'
    //    '<S333>/Lookup Table (n-D) Phase6'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S340>/Lookup Table (n-D) Amp1'
    //    '<S340>/Lookup Table (n-D) Amp2'
    //    '<S340>/Lookup Table (n-D) Amp3'
    //    '<S340>/Lookup Table (n-D) Amp4'
    //    '<S340>/Lookup Table (n-D) Amp5'
    //    '<S340>/Lookup Table (n-D) Amp6'
    //    '<S340>/Lookup Table (n-D) Phase'
    //    '<S340>/Lookup Table (n-D) Phase1'
    //    '<S340>/Lookup Table (n-D) Phase2'
    //    '<S340>/Lookup Table (n-D) Phase4'
    //    '<S340>/Lookup Table (n-D) Phase5'
    //    '<S340>/Lookup Table (n-D) Phase6'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled19[2];

    // Pooled Parameter (Expression: WD_Amp{2})
    //  Referenced by:
    //    '<S46>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 2'

    real_T pooled20[2808];

    // Pooled Parameter (Expression: WD_Amp{3})
    //  Referenced by:
    //    '<S46>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 3'

    real_T pooled21[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{1})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp1'
    //    '<S52>/Lookup Table (n-D) Amp1'
    //    '<S141>/Lookup Table (n-D) Amp1'
    //    '<S148>/Lookup Table (n-D) Amp1'
    //    '<S237>/Lookup Table (n-D) Amp1'
    //    '<S244>/Lookup Table (n-D) Amp1'
    //    '<S333>/Lookup Table (n-D) Amp1'
    //    '<S340>/Lookup Table (n-D) Amp1'

    real_T pooled22[2808];

    // Pooled Parameter (Expression: ForceRAO_Freq)
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp1'
    //    '<S45>/Lookup Table (n-D) Amp2'
    //    '<S45>/Lookup Table (n-D) Amp3'
    //    '<S45>/Lookup Table (n-D) Amp4'
    //    '<S45>/Lookup Table (n-D) Amp5'
    //    '<S45>/Lookup Table (n-D) Amp6'
    //    '<S45>/Lookup Table (n-D) Phase'
    //    '<S45>/Lookup Table (n-D) Phase1'
    //    '<S45>/Lookup Table (n-D) Phase2'
    //    '<S45>/Lookup Table (n-D) Phase4'
    //    '<S45>/Lookup Table (n-D) Phase5'
    //    '<S45>/Lookup Table (n-D) Phase6'
    //    '<S52>/Lookup Table (n-D) Amp1'
    //    '<S52>/Lookup Table (n-D) Amp2'
    //    '<S52>/Lookup Table (n-D) Amp3'
    //    '<S52>/Lookup Table (n-D) Amp4'
    //    '<S52>/Lookup Table (n-D) Amp5'
    //    '<S52>/Lookup Table (n-D) Amp6'
    //    '<S52>/Lookup Table (n-D) Phase'
    //    '<S52>/Lookup Table (n-D) Phase1'
    //    '<S52>/Lookup Table (n-D) Phase2'
    //    '<S52>/Lookup Table (n-D) Phase4'
    //    '<S52>/Lookup Table (n-D) Phase5'
    //    '<S52>/Lookup Table (n-D) Phase6'
    //    '<S141>/Lookup Table (n-D) Amp1'
    //    '<S141>/Lookup Table (n-D) Amp2'
    //    '<S141>/Lookup Table (n-D) Amp3'
    //    '<S141>/Lookup Table (n-D) Amp4'
    //    '<S141>/Lookup Table (n-D) Amp5'
    //    '<S141>/Lookup Table (n-D) Amp6'
    //    '<S141>/Lookup Table (n-D) Phase'
    //    '<S141>/Lookup Table (n-D) Phase1'
    //    '<S141>/Lookup Table (n-D) Phase2'
    //    '<S141>/Lookup Table (n-D) Phase4'
    //    '<S141>/Lookup Table (n-D) Phase5'
    //    '<S141>/Lookup Table (n-D) Phase6'
    //    '<S148>/Lookup Table (n-D) Amp1'
    //    '<S148>/Lookup Table (n-D) Amp2'
    //    '<S148>/Lookup Table (n-D) Amp3'
    //    '<S148>/Lookup Table (n-D) Amp4'
    //    '<S148>/Lookup Table (n-D) Amp5'
    //    '<S148>/Lookup Table (n-D) Amp6'
    //    '<S148>/Lookup Table (n-D) Phase'
    //    '<S148>/Lookup Table (n-D) Phase1'
    //    '<S148>/Lookup Table (n-D) Phase2'
    //    '<S148>/Lookup Table (n-D) Phase4'
    //    '<S148>/Lookup Table (n-D) Phase5'
    //    '<S148>/Lookup Table (n-D) Phase6'
    //    '<S237>/Lookup Table (n-D) Amp1'
    //    '<S237>/Lookup Table (n-D) Amp2'
    //    '<S237>/Lookup Table (n-D) Amp3'
    //    '<S237>/Lookup Table (n-D) Amp4'
    //    '<S237>/Lookup Table (n-D) Amp5'
    //    '<S237>/Lookup Table (n-D) Amp6'
    //    '<S237>/Lookup Table (n-D) Phase'
    //    '<S237>/Lookup Table (n-D) Phase1'
    //    '<S237>/Lookup Table (n-D) Phase2'
    //    '<S237>/Lookup Table (n-D) Phase4'
    //    '<S237>/Lookup Table (n-D) Phase5'
    //    '<S237>/Lookup Table (n-D) Phase6'
    //    '<S244>/Lookup Table (n-D) Amp1'
    //    '<S244>/Lookup Table (n-D) Amp2'
    //    '<S244>/Lookup Table (n-D) Amp3'
    //    '<S244>/Lookup Table (n-D) Amp4'
    //    '<S244>/Lookup Table (n-D) Amp5'
    //    '<S244>/Lookup Table (n-D) Amp6'
    //    '<S244>/Lookup Table (n-D) Phase'
    //    '<S244>/Lookup Table (n-D) Phase1'
    //    '<S244>/Lookup Table (n-D) Phase2'
    //    '<S244>/Lookup Table (n-D) Phase4'
    //    '<S244>/Lookup Table (n-D) Phase5'
    //    '<S244>/Lookup Table (n-D) Phase6'
    //    '<S333>/Lookup Table (n-D) Amp1'
    //    '<S333>/Lookup Table (n-D) Amp2'
    //    '<S333>/Lookup Table (n-D) Amp3'
    //    '<S333>/Lookup Table (n-D) Amp4'
    //    '<S333>/Lookup Table (n-D) Amp5'
    //    '<S333>/Lookup Table (n-D) Amp6'
    //    '<S333>/Lookup Table (n-D) Phase'
    //    '<S333>/Lookup Table (n-D) Phase1'
    //    '<S333>/Lookup Table (n-D) Phase2'
    //    '<S333>/Lookup Table (n-D) Phase4'
    //    '<S333>/Lookup Table (n-D) Phase5'
    //    '<S333>/Lookup Table (n-D) Phase6'
    //    '<S340>/Lookup Table (n-D) Amp1'
    //    '<S340>/Lookup Table (n-D) Amp2'
    //    '<S340>/Lookup Table (n-D) Amp3'
    //    '<S340>/Lookup Table (n-D) Amp4'
    //    '<S340>/Lookup Table (n-D) Amp5'
    //    '<S340>/Lookup Table (n-D) Amp6'
    //    '<S340>/Lookup Table (n-D) Phase'
    //    '<S340>/Lookup Table (n-D) Phase1'
    //    '<S340>/Lookup Table (n-D) Phase2'
    //    '<S340>/Lookup Table (n-D) Phase4'
    //    '<S340>/Lookup Table (n-D) Phase5'
    //    '<S340>/Lookup Table (n-D) Phase6'

    real_T pooled23[39];

    // Pooled Parameter (Expression: ForceRAO_Amp{2})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp2'
    //    '<S52>/Lookup Table (n-D) Amp2'
    //    '<S141>/Lookup Table (n-D) Amp2'
    //    '<S148>/Lookup Table (n-D) Amp2'
    //    '<S237>/Lookup Table (n-D) Amp2'
    //    '<S244>/Lookup Table (n-D) Amp2'
    //    '<S333>/Lookup Table (n-D) Amp2'
    //    '<S340>/Lookup Table (n-D) Amp2'

    real_T pooled24[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{3})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp3'
    //    '<S52>/Lookup Table (n-D) Amp3'
    //    '<S141>/Lookup Table (n-D) Amp3'
    //    '<S148>/Lookup Table (n-D) Amp3'
    //    '<S237>/Lookup Table (n-D) Amp3'
    //    '<S244>/Lookup Table (n-D) Amp3'
    //    '<S333>/Lookup Table (n-D) Amp3'
    //    '<S340>/Lookup Table (n-D) Amp3'

    real_T pooled25[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{4})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp4'
    //    '<S52>/Lookup Table (n-D) Amp4'
    //    '<S141>/Lookup Table (n-D) Amp4'
    //    '<S148>/Lookup Table (n-D) Amp4'
    //    '<S237>/Lookup Table (n-D) Amp4'
    //    '<S244>/Lookup Table (n-D) Amp4'
    //    '<S333>/Lookup Table (n-D) Amp4'
    //    '<S340>/Lookup Table (n-D) Amp4'

    real_T pooled26[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{5})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp5'
    //    '<S52>/Lookup Table (n-D) Amp5'
    //    '<S141>/Lookup Table (n-D) Amp5'
    //    '<S148>/Lookup Table (n-D) Amp5'
    //    '<S237>/Lookup Table (n-D) Amp5'
    //    '<S244>/Lookup Table (n-D) Amp5'
    //    '<S333>/Lookup Table (n-D) Amp5'
    //    '<S340>/Lookup Table (n-D) Amp5'

    real_T pooled27[2808];

    // Pooled Parameter (Expression: ForceRAO_Amp{6})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp6'
    //    '<S52>/Lookup Table (n-D) Amp6'
    //    '<S141>/Lookup Table (n-D) Amp6'
    //    '<S148>/Lookup Table (n-D) Amp6'
    //    '<S237>/Lookup Table (n-D) Amp6'
    //    '<S244>/Lookup Table (n-D) Amp6'
    //    '<S333>/Lookup Table (n-D) Amp6'
    //    '<S340>/Lookup Table (n-D) Amp6'

    real_T pooled28[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{1})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase1'
    //    '<S52>/Lookup Table (n-D) Phase1'
    //    '<S141>/Lookup Table (n-D) Phase1'
    //    '<S148>/Lookup Table (n-D) Phase1'
    //    '<S237>/Lookup Table (n-D) Phase1'
    //    '<S244>/Lookup Table (n-D) Phase1'
    //    '<S333>/Lookup Table (n-D) Phase1'
    //    '<S340>/Lookup Table (n-D) Phase1'

    real_T pooled29[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{2})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase2'
    //    '<S52>/Lookup Table (n-D) Phase2'
    //    '<S141>/Lookup Table (n-D) Phase2'
    //    '<S148>/Lookup Table (n-D) Phase2'
    //    '<S237>/Lookup Table (n-D) Phase2'
    //    '<S244>/Lookup Table (n-D) Phase2'
    //    '<S333>/Lookup Table (n-D) Phase2'
    //    '<S340>/Lookup Table (n-D) Phase2'

    real_T pooled30[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{3})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase'
    //    '<S52>/Lookup Table (n-D) Phase'
    //    '<S141>/Lookup Table (n-D) Phase'
    //    '<S148>/Lookup Table (n-D) Phase'
    //    '<S237>/Lookup Table (n-D) Phase'
    //    '<S244>/Lookup Table (n-D) Phase'
    //    '<S333>/Lookup Table (n-D) Phase'
    //    '<S340>/Lookup Table (n-D) Phase'

    real_T pooled31[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{4})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase4'
    //    '<S52>/Lookup Table (n-D) Phase4'
    //    '<S141>/Lookup Table (n-D) Phase4'
    //    '<S148>/Lookup Table (n-D) Phase4'
    //    '<S237>/Lookup Table (n-D) Phase4'
    //    '<S244>/Lookup Table (n-D) Phase4'
    //    '<S333>/Lookup Table (n-D) Phase4'
    //    '<S340>/Lookup Table (n-D) Phase4'

    real_T pooled32[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{5})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase5'
    //    '<S52>/Lookup Table (n-D) Phase5'
    //    '<S141>/Lookup Table (n-D) Phase5'
    //    '<S148>/Lookup Table (n-D) Phase5'
    //    '<S237>/Lookup Table (n-D) Phase5'
    //    '<S244>/Lookup Table (n-D) Phase5'
    //    '<S333>/Lookup Table (n-D) Phase5'
    //    '<S340>/Lookup Table (n-D) Phase5'

    real_T pooled33[2808];

    // Pooled Parameter (Expression: ForceRAO_Phase{6})
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Phase6'
    //    '<S52>/Lookup Table (n-D) Phase6'
    //    '<S141>/Lookup Table (n-D) Phase6'
    //    '<S148>/Lookup Table (n-D) Phase6'
    //    '<S237>/Lookup Table (n-D) Phase6'
    //    '<S244>/Lookup Table (n-D) Phase6'
    //    '<S333>/Lookup Table (n-D) Phase6'
    //    '<S340>/Lookup Table (n-D) Phase6'

    real_T pooled34[2808];

    // Pooled Parameter (Expression: ABC.G)
    //  Referenced by:
    //    '<S15>/Spring stiffness'
    //    '<S111>/Spring stiffness'
    //    '<S207>/Spring stiffness'
    //    '<S303>/Spring stiffness'

    real_T pooled55[36];

    // Pooled Parameter (Expression: ABC.Binf)
    //  Referenced by:
    //    '<S15>/damping'
    //    '<S111>/damping'
    //    '<S207>/damping'
    //    '<S303>/damping'

    real_T pooled58[36];

    // Pooled Parameter (Expression: vesselABC_dp.A44(:,:,1))
    //  Referenced by:
    //    '<S34>/A44'
    //    '<S130>/A44'
    //    '<S226>/A44'
    //    '<S322>/A44'

    real_T pooled59[25];

    // Pooled Parameter (Expression: vesselABC_dp.B44(:,:,1))
    //  Referenced by:
    //    '<S34>/B44'
    //    '<S130>/B44'
    //    '<S226>/B44'
    //    '<S322>/B44'

    real_T pooled60[5];

    // Pooled Parameter (Expression: Kd)
    //  Referenced by:
    //    '<S67>/Kd'
    //    '<S163>/Kd'
    //    '<S259>/Kd'
    //    '<S355>/Kd'

    real_T pooled70[9];

    // Pooled Parameter (Expression: [0 0 1 0 0 0])
    //  Referenced by:
    //    '<S10>/tau_cableGain1'
    //    '<S106>/tau_cableGain1'
    //    '<S202>/tau_cableGain1'
    //    '<S298>/tau_cableGain1'

    real_T pooled77[6];

    // Pooled Parameter (Expression: K4)
    //  Referenced by:
    //    '<S65>/K4'
    //    '<S161>/K4'
    //    '<S257>/K4'
    //    '<S353>/K4'

    real_T pooled114[9];

    // Pooled Parameter (Expression: D)
    //  Referenced by:
    //    '<S65>/Gain6'
    //    '<S161>/Gain6'
    //    '<S257>/Gain6'
    //    '<S353>/Gain6'

    real_T pooled115[9];

    // Pooled Parameter (Expression: inv(M))
    //  Referenced by:
    //    '<S65>/Gain3'
    //    '<S161>/Gain3'
    //    '<S257>/Gain3'
    //    '<S353>/Gain3'

    real_T pooled116[9];

    // Pooled Parameter (Expression: -2*(eye(3)-lambda)*diag([w_c(1,1)/w_o(1,1) w_c(2,2)/w_o(2,2) w_c(3,3)/w_o(3,3)]))
    //  Referenced by:
    //    '<S65>/K11'
    //    '<S161>/K11'
    //    '<S257>/K11'
    //    '<S353>/K11'

    real_T pooled117[9];

    // Pooled Parameter (Expression: 2*w_o*(eye(3)-lambda))
    //  Referenced by:
    //    '<S65>/K12'
    //    '<S161>/K12'
    //    '<S257>/K12'
    //    '<S353>/K12'

    real_T pooled118[9];

    // Pooled Parameter (Expression: 2*lambda*w_o)
    //  Referenced by:
    //    '<S65>/Gain1'
    //    '<S161>/Gain1'
    //    '<S257>/Gain1'
    //    '<S353>/Gain1'

    real_T pooled119[9];

    // Pooled Parameter (Expression: w_o*w_o)
    //  Referenced by:
    //    '<S65>/Gain2'
    //    '<S161>/Gain2'
    //    '<S257>/Gain2'
    //    '<S353>/Gain2'

    real_T pooled120[9];

    // Pooled Parameter (Expression: w_c)
    //  Referenced by:
    //    '<S65>/K2'
    //    '<S161>/K2'
    //    '<S257>/K2'
    //    '<S353>/K2'

    real_T pooled121[9];

    // Pooled Parameter (Expression: K3)
    //  Referenced by:
    //    '<S65>/K3'
    //    '<S161>/K3'
    //    '<S257>/K3'
    //    '<S353>/K3'

    real_T pooled122[9];

    // Pooled Parameter (Expression: diag([1/T_b(1,1) 1/T_b(2,2) 1/T_b(3,3)]))
    //  Referenced by:
    //    '<S65>/inv(T_b)'
    //    '<S161>/inv(T_b)'
    //    '<S257>/inv(T_b)'
    //    '<S353>/inv(T_b)'

    real_T pooled123[9];

    // Pooled Parameter (Expression: )
    //  Referenced by:
    //    '<S45>/Lookup Table (n-D) Amp1'
    //    '<S45>/Lookup Table (n-D) Amp2'
    //    '<S45>/Lookup Table (n-D) Amp3'
    //    '<S45>/Lookup Table (n-D) Amp4'
    //    '<S45>/Lookup Table (n-D) Amp5'
    //    '<S45>/Lookup Table (n-D) Amp6'
    //    '<S45>/Lookup Table (n-D) Phase'
    //    '<S45>/Lookup Table (n-D) Phase1'
    //    '<S45>/Lookup Table (n-D) Phase2'
    //    '<S45>/Lookup Table (n-D) Phase4'
    //    '<S45>/Lookup Table (n-D) Phase5'
    //    '<S45>/Lookup Table (n-D) Phase6'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S46>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S52>/Lookup Table (n-D) Amp1'
    //    '<S52>/Lookup Table (n-D) Amp2'
    //    '<S52>/Lookup Table (n-D) Amp3'
    //    '<S52>/Lookup Table (n-D) Amp4'
    //    '<S52>/Lookup Table (n-D) Amp5'
    //    '<S52>/Lookup Table (n-D) Amp6'
    //    '<S52>/Lookup Table (n-D) Phase'
    //    '<S52>/Lookup Table (n-D) Phase1'
    //    '<S52>/Lookup Table (n-D) Phase2'
    //    '<S52>/Lookup Table (n-D) Phase4'
    //    '<S52>/Lookup Table (n-D) Phase5'
    //    '<S52>/Lookup Table (n-D) Phase6'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S53>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S141>/Lookup Table (n-D) Amp1'
    //    '<S141>/Lookup Table (n-D) Amp2'
    //    '<S141>/Lookup Table (n-D) Amp3'
    //    '<S141>/Lookup Table (n-D) Amp4'
    //    '<S141>/Lookup Table (n-D) Amp5'
    //    '<S141>/Lookup Table (n-D) Amp6'
    //    '<S141>/Lookup Table (n-D) Phase'
    //    '<S141>/Lookup Table (n-D) Phase1'
    //    '<S141>/Lookup Table (n-D) Phase2'
    //    '<S141>/Lookup Table (n-D) Phase4'
    //    '<S141>/Lookup Table (n-D) Phase5'
    //    '<S141>/Lookup Table (n-D) Phase6'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S142>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S148>/Lookup Table (n-D) Amp1'
    //    '<S148>/Lookup Table (n-D) Amp2'
    //    '<S148>/Lookup Table (n-D) Amp3'
    //    '<S148>/Lookup Table (n-D) Amp4'
    //    '<S148>/Lookup Table (n-D) Amp5'
    //    '<S148>/Lookup Table (n-D) Amp6'
    //    '<S148>/Lookup Table (n-D) Phase'
    //    '<S148>/Lookup Table (n-D) Phase1'
    //    '<S148>/Lookup Table (n-D) Phase2'
    //    '<S148>/Lookup Table (n-D) Phase4'
    //    '<S148>/Lookup Table (n-D) Phase5'
    //    '<S148>/Lookup Table (n-D) Phase6'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S149>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S237>/Lookup Table (n-D) Amp1'
    //    '<S237>/Lookup Table (n-D) Amp2'
    //    '<S237>/Lookup Table (n-D) Amp3'
    //    '<S237>/Lookup Table (n-D) Amp4'
    //    '<S237>/Lookup Table (n-D) Amp5'
    //    '<S237>/Lookup Table (n-D) Amp6'
    //    '<S237>/Lookup Table (n-D) Phase'
    //    '<S237>/Lookup Table (n-D) Phase1'
    //    '<S237>/Lookup Table (n-D) Phase2'
    //    '<S237>/Lookup Table (n-D) Phase4'
    //    '<S237>/Lookup Table (n-D) Phase5'
    //    '<S237>/Lookup Table (n-D) Phase6'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S238>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S244>/Lookup Table (n-D) Amp1'
    //    '<S244>/Lookup Table (n-D) Amp2'
    //    '<S244>/Lookup Table (n-D) Amp3'
    //    '<S244>/Lookup Table (n-D) Amp4'
    //    '<S244>/Lookup Table (n-D) Amp5'
    //    '<S244>/Lookup Table (n-D) Amp6'
    //    '<S244>/Lookup Table (n-D) Phase'
    //    '<S244>/Lookup Table (n-D) Phase1'
    //    '<S244>/Lookup Table (n-D) Phase2'
    //    '<S244>/Lookup Table (n-D) Phase4'
    //    '<S244>/Lookup Table (n-D) Phase5'
    //    '<S244>/Lookup Table (n-D) Phase6'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S245>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S333>/Lookup Table (n-D) Amp1'
    //    '<S333>/Lookup Table (n-D) Amp2'
    //    '<S333>/Lookup Table (n-D) Amp3'
    //    '<S333>/Lookup Table (n-D) Amp4'
    //    '<S333>/Lookup Table (n-D) Amp5'
    //    '<S333>/Lookup Table (n-D) Amp6'
    //    '<S333>/Lookup Table (n-D) Phase'
    //    '<S333>/Lookup Table (n-D) Phase1'
    //    '<S333>/Lookup Table (n-D) Phase2'
    //    '<S333>/Lookup Table (n-D) Phase4'
    //    '<S333>/Lookup Table (n-D) Phase5'
    //    '<S333>/Lookup Table (n-D) Phase6'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S334>/Lookup Table (n-D)  Wavedrift 3'
    //    '<S340>/Lookup Table (n-D) Amp1'
    //    '<S340>/Lookup Table (n-D) Amp2'
    //    '<S340>/Lookup Table (n-D) Amp3'
    //    '<S340>/Lookup Table (n-D) Amp4'
    //    '<S340>/Lookup Table (n-D) Amp5'
    //    '<S340>/Lookup Table (n-D) Amp6'
    //    '<S340>/Lookup Table (n-D) Phase'
    //    '<S340>/Lookup Table (n-D) Phase1'
    //    '<S340>/Lookup Table (n-D) Phase2'
    //    '<S340>/Lookup Table (n-D) Phase4'
    //    '<S340>/Lookup Table (n-D) Phase5'
    //    '<S340>/Lookup Table (n-D) Phase6'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 1'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 2'
    //    '<S341>/Lookup Table (n-D)  Wavedrift 3'

    uint32_T pooled126[3];
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
    real_T odeY[472];
    real_T odeF[4][472];
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

  // private member function(s) for subsystem '<S21>/Cross-flow drag trapezoidal integration'
  void Crossflowdragtrapezoidali_Reset(real_T *memory2_PreviousInput, real_T
    *memory1_PreviousInput);
  void Crossflowdragtrapezoidalintegra(real_T rtu_N, real_T rtu_dx, real_T
    rtu_v_r, real_T rtu_r, real_T *rty_sum1, real_T *rty_sum2, real_T rtp_Lpp);

  // private member function(s) for subsystem '<S10>/Subsystem3'
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
  void AHV_Model_power_pn(const real_T a_data[], const int32_T a_size[2], real_T
    b, real_T y_data[], int32_T y_size[2]);
  void AHV_Model_power_p(real_T a, const real_T b_data[], const int32_T b_size[2],
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

  // private member function(s) for subsystem '<S17>/Wave_loads_for_heading1'
  void AHV_Mod_Wave_loads_for_heading1(const real_T rtu_eta[6], real_T rtu_psi,
    const real_T rtu_psi_wave[900], const real_T rtu_wavenum[900], const real_T
    rtu_Omega[900], const real_T rtu_Phase[900], const real_T rtu_Zeta_a[900],
    real_T rty_tau_WF[6], real_T rty_tau_WD[6], B_Wave_loads_for_heading1_AHV_T *
    localB, const ConstB_Wave_loads_for_heading_T *localC,
    DW_Wave_loads_for_heading1_AH_T *localDW);

  // private member function(s) for subsystem '<S10>/Wave loads (U=0)'
  void Wave_loads_fun(const real_T rtu_eta[6], const real_T rtu_waves[900],
                      const real_T rtu_waves_m[900], const real_T rtu_waves_p
                      [900], const real_T rtu_waves_pp[900], const real_T
                      rtu_waves_pz[900], real_T rty_tau_WF[6], real_T
                      rty_tau_WD[6], B_Wave_loads_fun_T *localB, const
                      ConstB_Wave_loads_fun_T *localC, DW_Wave_loads_fun_T
                      *localDW);

  // private member function(s) for subsystem '<S66>/Chart'
  void AHV_Model_Chart(boolean_T rtu_hold, real_T rtu_x_ref, real_T rtu_y_ref,
                       real_T rtu_x_hold, real_T rtu_y_hold, real_T rtu_yaw_ref,
                       real_T rtu_yaw_hold, real_T *rty_x_ref_rel, real_T
                       *rty_y_ref_rel, real_T *rty_yaw_ref_rel,
                       DW_Chart_AHV_Model_T *localDW);

  // private member function(s) for subsystem '<S78>/If Action Subsystem'
  void AHV_Model_IfActionSubsystem(real_T rtu_In1, real_T rtu_In1_b, real_T
    rtu_In1_l, real_T rty_Kp[3], const real_T rtp_Kp_A[9]);

  // private member function(s) for subsystem '<S78>/If Action Subsystem1'
  void AHV_Model_IfActionSubsystem1(real_T rtu_In1, real_T rtu_In1_i, real_T
    rtu_In1_k, real_T rty_Kp[3], const real_T rtp_Kp[9]);

  // private member function(s) for subsystem '<S85>/If Action Subsystem'
  void AHV_Model_IfActionSubsystem_f(real_T rtu_Increment, real_T rtu_main_T,
    real_T *rty_Out1);

  // private member function(s) for subsystem '<S94>/Chart1'
  void AHV_Model_Chart1_Init(real_T rtu_Rudder_angle, real_T rtu_eta_yaw_deg,
    real_T *rty_heading_deg, DW_Chart1_AHV_Model_T *localDW);
  void AHV_Model_Chart1(real_T rtu_Rudder_angle, real_T rtu_eta_yaw_deg, real_T *
                        rty_heading_deg, DW_Chart1_AHV_Model_T *localDW);

  // private member function(s) for subsystem '<S94>/Chart3'
  void AHV_Model_Chart3_Init(real_T rtu_eta_yaw_deg, real_T rtu_Rudder_angle,
    real_T *rty_heading_deg, DW_Chart3_AHV_Model_T *localDW);
  void AHV_Model_Chart3(real_T rtu_eta_yaw_deg, real_T rtu_Rudder_angle, real_T *
                        rty_heading_deg, DW_Chart3_AHV_Model_T *localDW);

  // private member function(s) for subsystem '<S8>/If Action Subsystem1'
  void AHV_Model_IfActionSubsystem1_k(real_T rtu_heading_mode, real_T
    rtu_Head_angel, real_T rtu_target_surge, real_T rtu_main_add, real_T
    rtu_fu_add, real_T rtu_target_sway, real_T *rty_Target_yaw, real_T
    *rty_Target_x, real_T *rty_surge_add, real_T *rty_sway_add, real_T
    *rty_Target_y, real_T *rty_heading_mode_dp, real_T *rty_SurgeForce, real_T
    *rty_SwayForce, real_T *rty_YawMoment);

  // private member function(s) for subsystem '<S93>/If Action Subsystem2'
  void AHV_Model_IfActionSubsystem2(boolean_T *rty_KeepPosition);

  // private member function(s) for subsystem '<S104>/If Action Subsystem'
  void AHV_Mo_IfActionSubsystem_j_Init(const real_T rtu_eta[6],
    B_IfActionSubsystem_AHV_Mod_d_T *localB, DW_IfActionSubsystem_AHV_M_ci_T
    *localDW);
  void AHV_Mo_IfActionSubsystem_Update(B_IfActionSubsystem_AHV_Mod_d_T *localB,
    DW_IfActionSubsystem_AHV_M_ci_T *localDW);
  void AHV_Model_IfActionSubsystem_c(const real_T rtu_eta[6], real_T
    rtu_Rudder_angle, real_T rtu_Thruster_percentage1, real_T rtu_heading_mode,
    boolean_T rtu_HoldPosition, real_T *rty_Target_yaw_headdeg, real_T
    *rty_Target_x_head, real_T *rty_surge_add_head, real_T *rty_sway_add_head,
    real_T *rty_Target_y_head, real_T *rty_heading_mode_head, real_T
    *rty_SurgeForce, real_T *rty_SwayForce, real_T *rty_YawMoment,
    B_IfActionSubsystem_AHV_Mod_d_T *localB, DW_IfActionSubsystem_AHV_M_ci_T
    *localDW);

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
//  Block '<S10>/Display' : Unused code path elimination
//  Block '<S10>/Display1' : Unused code path elimination
//  Block '<S10>/Display2' : Unused code path elimination
//  Block '<S10>/Display3' : Unused code path elimination
//  Block '<S10>/Display4' : Unused code path elimination
//  Block '<S10>/Display5' : Unused code path elimination
//  Block '<S56>/Gain' : Unused code path elimination
//  Block '<S57>/Gain' : Unused code path elimination
//  Block '<S58>/Gain' : Unused code path elimination
//  Block '<S59>/Gain' : Unused code path elimination
//  Block '<S60>/Gain' : Unused code path elimination
//  Block '<S61>/Gain' : Unused code path elimination
//  Block '<S19>/Gain' : Unused code path elimination
//  Block '<S10>/wave direction1' : Unused code path elimination
//  Block '<S94>/Display1' : Unused code path elimination
//  Block '<S94>/Display8' : Unused code path elimination
//  Block '<S106>/Display' : Unused code path elimination
//  Block '<S106>/Display1' : Unused code path elimination
//  Block '<S106>/Display2' : Unused code path elimination
//  Block '<S106>/Display3' : Unused code path elimination
//  Block '<S106>/Display4' : Unused code path elimination
//  Block '<S106>/Display5' : Unused code path elimination
//  Block '<S152>/Gain' : Unused code path elimination
//  Block '<S153>/Gain' : Unused code path elimination
//  Block '<S154>/Gain' : Unused code path elimination
//  Block '<S155>/Gain' : Unused code path elimination
//  Block '<S156>/Gain' : Unused code path elimination
//  Block '<S157>/Gain' : Unused code path elimination
//  Block '<S115>/Gain' : Unused code path elimination
//  Block '<S106>/wave direction1' : Unused code path elimination
//  Block '<S190>/Display1' : Unused code path elimination
//  Block '<S190>/Display8' : Unused code path elimination
//  Block '<S202>/Display' : Unused code path elimination
//  Block '<S202>/Display1' : Unused code path elimination
//  Block '<S202>/Display2' : Unused code path elimination
//  Block '<S202>/Display3' : Unused code path elimination
//  Block '<S202>/Display4' : Unused code path elimination
//  Block '<S202>/Display5' : Unused code path elimination
//  Block '<S248>/Gain' : Unused code path elimination
//  Block '<S249>/Gain' : Unused code path elimination
//  Block '<S250>/Gain' : Unused code path elimination
//  Block '<S251>/Gain' : Unused code path elimination
//  Block '<S252>/Gain' : Unused code path elimination
//  Block '<S253>/Gain' : Unused code path elimination
//  Block '<S211>/Gain' : Unused code path elimination
//  Block '<S202>/wave direction1' : Unused code path elimination
//  Block '<S286>/Display1' : Unused code path elimination
//  Block '<S286>/Display8' : Unused code path elimination
//  Block '<S298>/Display' : Unused code path elimination
//  Block '<S298>/Display1' : Unused code path elimination
//  Block '<S298>/Display2' : Unused code path elimination
//  Block '<S298>/Display3' : Unused code path elimination
//  Block '<S298>/Display4' : Unused code path elimination
//  Block '<S298>/Display5' : Unused code path elimination
//  Block '<S344>/Gain' : Unused code path elimination
//  Block '<S345>/Gain' : Unused code path elimination
//  Block '<S346>/Gain' : Unused code path elimination
//  Block '<S347>/Gain' : Unused code path elimination
//  Block '<S348>/Gain' : Unused code path elimination
//  Block '<S349>/Gain' : Unused code path elimination
//  Block '<S307>/Gain' : Unused code path elimination
//  Block '<S298>/wave direction1' : Unused code path elimination
//  Block '<S382>/Display1' : Unused code path elimination
//  Block '<S382>/Display8' : Unused code path elimination
//  Block '<S23>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S10>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S63>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S65>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S68>/Gain' : Eliminated nontunable gain of 1
//  Block '<S68>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S68>/Gain2' : Eliminated nontunable gain of 1
//  Block '<S68>/Gain3' : Eliminated nontunable gain of 1
//  Block '<S68>/Gain4' : Eliminated nontunable gain of 1
//  Block '<S68>/Gain5' : Eliminated nontunable gain of 1
//  Block '<S9>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S119>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S106>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S159>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S161>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S164>/Gain' : Eliminated nontunable gain of 1
//  Block '<S164>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S164>/Gain2' : Eliminated nontunable gain of 1
//  Block '<S164>/Gain3' : Eliminated nontunable gain of 1
//  Block '<S164>/Gain4' : Eliminated nontunable gain of 1
//  Block '<S164>/Gain5' : Eliminated nontunable gain of 1
//  Block '<S105>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S215>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S202>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S255>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S257>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S260>/Gain' : Eliminated nontunable gain of 1
//  Block '<S260>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S260>/Gain2' : Eliminated nontunable gain of 1
//  Block '<S260>/Gain3' : Eliminated nontunable gain of 1
//  Block '<S260>/Gain4' : Eliminated nontunable gain of 1
//  Block '<S260>/Gain5' : Eliminated nontunable gain of 1
//  Block '<S201>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S311>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S298>/Wave Drift on//off' : Eliminated nontunable gain of 1
//  Block '<S351>/Reshape 9x1->3x3' : Reshape block reduction
//  Block '<S353>/Signal Conversion' : Eliminate redundant signal conversion block
//  Block '<S356>/Gain' : Eliminated nontunable gain of 1
//  Block '<S356>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S356>/Gain2' : Eliminated nontunable gain of 1
//  Block '<S356>/Gain3' : Eliminated nontunable gain of 1
//  Block '<S356>/Gain4' : Eliminated nontunable gain of 1
//  Block '<S356>/Gain5' : Eliminated nontunable gain of 1
//  Block '<S297>/Gain1' : Eliminated nontunable gain of 1
//  Block '<S68>/Constant' : Unused code path elimination
//  Block '<S68>/Sway Failure' : Unused code path elimination
//  Block '<S68>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S164>/Constant' : Unused code path elimination
//  Block '<S164>/Sway Failure' : Unused code path elimination
//  Block '<S164>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S260>/Constant' : Unused code path elimination
//  Block '<S260>/Sway Failure' : Unused code path elimination
//  Block '<S260>/Yaw Moment Failure' : Unused code path elimination
//  Block '<S356>/Constant' : Unused code path elimination
//  Block '<S356>/Sway Failure' : Unused code path elimination
//  Block '<S356>/Yaw Moment Failure' : Unused code path elimination


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
//  '<S8>'   : 'AHV_Model/AH_Model1/Subsystem1'
//  '<S9>'   : 'AHV_Model/AH_Model1/sway_add'
//  '<S10>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1'
//  '<S11>'  : 'AHV_Model/AH_Model1/AHV_Model/Subsystem'
//  '<S12>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S13>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S14>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S15>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S16>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S17>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S18>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S19>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S20>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S21>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S22>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S23>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S24>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S25>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S26>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S27>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S28>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S29>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S30>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S31>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S32>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S33>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S34>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S35>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S36>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S37>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S38>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S39>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S40>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S41>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S42>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S43>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S44>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S45>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S46>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S47>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S48>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S49>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S50>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S51>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S52>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S53>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S54>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S55>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S56>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S57>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S58>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S59>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S60>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S61>'  : 'AHV_Model/AH_Model1/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S62>'  : 'AHV_Model/AH_Model1/AHV_Model/Subsystem/Cross product'
//  '<S63>'  : 'AHV_Model/AH_Model1/AHV_Model/Subsystem/Rbn_gnc'
//  '<S64>'  : 'AHV_Model/AH_Model1/DP_Thrust/Calculate Cable Angles'
//  '<S65>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter'
//  '<S66>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference1'
//  '<S67>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem'
//  '<S68>'  : 'AHV_Model/AH_Model1/DP_Thrust/Thrust Limitations'
//  '<S69>'  : 'AHV_Model/AH_Model1/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S70>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S71>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S72>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S73>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S74>'  : 'AHV_Model/AH_Model1/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S75>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference1/Chart'
//  '<S76>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference1/Degrees to Radians'
//  '<S77>'  : 'AHV_Model/AH_Model1/DP_Thrust/Position and Heading Reference1/[-inf inf] to [-pi pi]'
//  '<S78>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp2'
//  '<S79>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3'
//  '<S80>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S81>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S82>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S83>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp2/If Action Subsystem'
//  '<S84>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Kp2/If Action Subsystem1'
//  '<S85>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem'
//  '<S86>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem1'
//  '<S87>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem'
//  '<S88>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem1'
//  '<S89>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem'
//  '<S90>'  : 'AHV_Model/AH_Model1/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem1'
//  '<S91>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem'
//  '<S92>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem1'
//  '<S93>'  : 'AHV_Model/AH_Model1/Subsystem1/Subsystem'
//  '<S94>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem/Thruster_percentage'
//  '<S95>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem/Thruster_percentage/Chart1'
//  '<S96>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem/Thruster_percentage/Chart2'
//  '<S97>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem/Thruster_percentage/Chart3'
//  '<S98>'  : 'AHV_Model/AH_Model1/Subsystem1/If Action Subsystem/Thruster_percentage/deg-360'
//  '<S99>'  : 'AHV_Model/AH_Model1/Subsystem1/Subsystem/If Action Subsystem2'
//  '<S100>' : 'AHV_Model/AH_Model1/Subsystem1/Subsystem/If Action Subsystem3'
//  '<S101>' : 'AHV_Model/AH_Model2/AHV_Model'
//  '<S102>' : 'AHV_Model/AH_Model2/DP_Thrust'
//  '<S103>' : 'AHV_Model/AH_Model2/Manual Control'
//  '<S104>' : 'AHV_Model/AH_Model2/Subsystem1'
//  '<S105>' : 'AHV_Model/AH_Model2/sway_add'
//  '<S106>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1'
//  '<S107>' : 'AHV_Model/AH_Model2/AHV_Model/Subsystem'
//  '<S108>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S109>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S110>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S111>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S112>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S113>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S114>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S115>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S116>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S117>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S118>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S119>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S120>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S121>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S122>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S123>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S124>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S125>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S126>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S127>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S128>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S129>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S130>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S131>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S132>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S133>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S134>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S135>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S136>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S137>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S138>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S139>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S140>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S141>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S142>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S143>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S144>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S145>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S146>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S147>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S148>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S149>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S150>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S151>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S152>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S153>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S154>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S155>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S156>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S157>' : 'AHV_Model/AH_Model2/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S158>' : 'AHV_Model/AH_Model2/AHV_Model/Subsystem/Cross product'
//  '<S159>' : 'AHV_Model/AH_Model2/AHV_Model/Subsystem/Rbn_gnc'
//  '<S160>' : 'AHV_Model/AH_Model2/DP_Thrust/Calculate Cable Angles'
//  '<S161>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter'
//  '<S162>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference1'
//  '<S163>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem'
//  '<S164>' : 'AHV_Model/AH_Model2/DP_Thrust/Thrust Limitations'
//  '<S165>' : 'AHV_Model/AH_Model2/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S166>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S167>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S168>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S169>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S170>' : 'AHV_Model/AH_Model2/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S171>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference1/Chart'
//  '<S172>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference1/Degrees to Radians'
//  '<S173>' : 'AHV_Model/AH_Model2/DP_Thrust/Position and Heading Reference1/[-inf inf] to [-pi pi]'
//  '<S174>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp2'
//  '<S175>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3'
//  '<S176>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S177>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S178>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S179>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp2/If Action Subsystem'
//  '<S180>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Kp2/If Action Subsystem1'
//  '<S181>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem'
//  '<S182>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem1'
//  '<S183>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem'
//  '<S184>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem1'
//  '<S185>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem'
//  '<S186>' : 'AHV_Model/AH_Model2/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem1'
//  '<S187>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem'
//  '<S188>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem1'
//  '<S189>' : 'AHV_Model/AH_Model2/Subsystem1/Subsystem'
//  '<S190>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem/Thruster_percentage'
//  '<S191>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem/Thruster_percentage/Chart1'
//  '<S192>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem/Thruster_percentage/Chart2'
//  '<S193>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem/Thruster_percentage/Chart3'
//  '<S194>' : 'AHV_Model/AH_Model2/Subsystem1/If Action Subsystem/Thruster_percentage/deg-360'
//  '<S195>' : 'AHV_Model/AH_Model2/Subsystem1/Subsystem/If Action Subsystem2'
//  '<S196>' : 'AHV_Model/AH_Model2/Subsystem1/Subsystem/If Action Subsystem3'
//  '<S197>' : 'AHV_Model/AH_Model3/AHV_Model'
//  '<S198>' : 'AHV_Model/AH_Model3/DP_Thrust'
//  '<S199>' : 'AHV_Model/AH_Model3/Manual Control'
//  '<S200>' : 'AHV_Model/AH_Model3/Subsystem1'
//  '<S201>' : 'AHV_Model/AH_Model3/sway_add'
//  '<S202>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1'
//  '<S203>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem'
//  '<S204>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S205>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S206>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S207>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S208>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S209>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S210>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S211>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S212>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S213>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S214>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S215>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S216>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S217>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S218>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S219>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S220>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S221>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S222>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S223>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S224>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S225>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S226>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S227>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S228>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S229>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S230>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S231>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S232>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S233>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S234>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S235>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S236>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S237>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S238>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S239>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S240>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S241>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S242>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S243>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S244>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S245>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S246>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S247>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S248>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S249>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S250>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S251>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S252>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S253>' : 'AHV_Model/AH_Model3/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S254>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem/Cross product'
//  '<S255>' : 'AHV_Model/AH_Model3/AHV_Model/Subsystem/Rbn_gnc'
//  '<S256>' : 'AHV_Model/AH_Model3/DP_Thrust/Calculate Cable Angles'
//  '<S257>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter'
//  '<S258>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference1'
//  '<S259>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem'
//  '<S260>' : 'AHV_Model/AH_Model3/DP_Thrust/Thrust Limitations'
//  '<S261>' : 'AHV_Model/AH_Model3/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S262>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S263>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S264>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S265>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S266>' : 'AHV_Model/AH_Model3/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S267>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference1/Chart'
//  '<S268>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference1/Degrees to Radians'
//  '<S269>' : 'AHV_Model/AH_Model3/DP_Thrust/Position and Heading Reference1/[-inf inf] to [-pi pi]'
//  '<S270>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp2'
//  '<S271>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3'
//  '<S272>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S273>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S274>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S275>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp2/If Action Subsystem'
//  '<S276>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Kp2/If Action Subsystem1'
//  '<S277>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem'
//  '<S278>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem1'
//  '<S279>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem'
//  '<S280>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem1'
//  '<S281>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem'
//  '<S282>' : 'AHV_Model/AH_Model3/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem1'
//  '<S283>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem'
//  '<S284>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem1'
//  '<S285>' : 'AHV_Model/AH_Model3/Subsystem1/Subsystem'
//  '<S286>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem/Thruster_percentage'
//  '<S287>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem/Thruster_percentage/Chart1'
//  '<S288>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem/Thruster_percentage/Chart2'
//  '<S289>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem/Thruster_percentage/Chart3'
//  '<S290>' : 'AHV_Model/AH_Model3/Subsystem1/If Action Subsystem/Thruster_percentage/deg-360'
//  '<S291>' : 'AHV_Model/AH_Model3/Subsystem1/Subsystem/If Action Subsystem2'
//  '<S292>' : 'AHV_Model/AH_Model3/Subsystem1/Subsystem/If Action Subsystem3'
//  '<S293>' : 'AHV_Model/AH_Model4/AHV_Model'
//  '<S294>' : 'AHV_Model/AH_Model4/DP_Thrust'
//  '<S295>' : 'AHV_Model/AH_Model4/Manual Control'
//  '<S296>' : 'AHV_Model/AH_Model4/Subsystem1'
//  '<S297>' : 'AHV_Model/AH_Model4/sway_add'
//  '<S298>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1'
//  '<S299>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem'
//  '<S300>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/3 to 6 DOF'
//  '<S301>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Model Info1'
//  '<S302>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/North-East current'
//  '<S303>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem'
//  '<S304>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem3'
//  '<S305>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)'
//  '<S306>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots'
//  '<S307>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/rad2deg1'
//  '<S308>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation'
//  '<S309>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1'
//  '<S310>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping'
//  '<S311>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to rotation matrix1'
//  '<S312>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/6 DOF transformation/Euler angles to attitude transformation matrix'
//  '<S313>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration'
//  '<S314>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 0.5'
//  '<S315>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/Cross-flow drag and surge resistance1/Cross-flow drag trapezoidal integration/multiply with 1'
//  '<S316>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/heave'
//  '<S317>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/pitch'
//  '<S318>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/roll'
//  '<S319>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/surge'
//  '<S320>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/sway'
//  '<S321>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/yaw'
//  '<S322>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem/h-frame  potential//viscous  damping/zero speed viiscous  roll damping'
//  '<S323>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Subsystem3/Wave'
//  '<S324>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Compare To Constant '
//  '<S325>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1'
//  '<S326>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2'
//  '<S327>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/deg to -180//180'
//  '<S328>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi1'
//  '<S329>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/psi2'
//  '<S330>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Position phase'
//  '<S331>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Subsystem'
//  '<S332>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop'
//  '<S333>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/Force RAO tables'
//  '<S334>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/WD RAO tables'
//  '<S335>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S336>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading1/Wave component loop/x to y'
//  '<S337>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Position phase'
//  '<S338>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Subsystem'
//  '<S339>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop'
//  '<S340>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/Force RAO tables'
//  '<S341>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/WD RAO tables'
//  '<S342>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/[-pi pi] to [0 2pi]'
//  '<S343>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/Wave loads (U=0)/Wave_loads_for_heading2/Wave component loop/x to y'
//  '<S344>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees1'
//  '<S345>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees2'
//  '<S346>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees3'
//  '<S347>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees4'
//  '<S348>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees5'
//  '<S349>' : 'AHV_Model/AH_Model4/AHV_Model/Modified MSS Vessel Model1/plots/Radians to Degrees6'
//  '<S350>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem/Cross product'
//  '<S351>' : 'AHV_Model/AH_Model4/AHV_Model/Subsystem/Rbn_gnc'
//  '<S352>' : 'AHV_Model/AH_Model4/DP_Thrust/Calculate Cable Angles'
//  '<S353>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter'
//  '<S354>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference1'
//  '<S355>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem'
//  '<S356>' : 'AHV_Model/AH_Model4/DP_Thrust/Thrust Limitations'
//  '<S357>' : 'AHV_Model/AH_Model4/DP_Thrust/Calculate Cable Angles/Radians to Degrees1'
//  '<S358>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Rotation matrix in yaw 1'
//  '<S359>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw'
//  '<S360>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/Transposed rotation  matrix in yaw1'
//  '<S361>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]'
//  '<S362>' : 'AHV_Model/AH_Model4/DP_Thrust/Passive DP wave filter/[-inf inf] to [-pi pi]1'
//  '<S363>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference1/Chart'
//  '<S364>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference1/Degrees to Radians'
//  '<S365>' : 'AHV_Model/AH_Model4/DP_Thrust/Position and Heading Reference1/[-inf inf] to [-pi pi]'
//  '<S366>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp2'
//  '<S367>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3'
//  '<S368>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Transposed rotation  matrix in yaw1'
//  '<S369>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]'
//  '<S370>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/[-inf inf] to [-pi pi]1'
//  '<S371>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp2/If Action Subsystem'
//  '<S372>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Kp2/If Action Subsystem1'
//  '<S373>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem'
//  '<S374>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem1'
//  '<S375>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem'
//  '<S376>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem/If Action Subsystem1'
//  '<S377>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem'
//  '<S378>' : 'AHV_Model/AH_Model4/DP_Thrust/Subsystem/Subsystem3/Subsystem1/If Action Subsystem1'
//  '<S379>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem'
//  '<S380>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem1'
//  '<S381>' : 'AHV_Model/AH_Model4/Subsystem1/Subsystem'
//  '<S382>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem/Thruster_percentage'
//  '<S383>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem/Thruster_percentage/Chart1'
//  '<S384>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem/Thruster_percentage/Chart2'
//  '<S385>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem/Thruster_percentage/Chart3'
//  '<S386>' : 'AHV_Model/AH_Model4/Subsystem1/If Action Subsystem/Thruster_percentage/deg-360'
//  '<S387>' : 'AHV_Model/AH_Model4/Subsystem1/Subsystem/If Action Subsystem2'
//  '<S388>' : 'AHV_Model/AH_Model4/Subsystem1/Subsystem/If Action Subsystem3'

#endif                                 // RTW_HEADER_AHV_Model_h_

//
// File trailer for generated code.
//
// [EOF]
//
