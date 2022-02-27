//
// File: Wave_loads.cpp
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
#include "Wave_loads.h"

// Include model header file for global data
#include "AHV_Model.h"
#include "AHV_Model_private.h"

//
// Outputs for atomic system:
//    '<S15>/Wave_loads_for_heading1'
//    '<S15>/Wave_loads_for_heading2'
//    '<S99>/Wave_loads_for_heading1'
//    '<S99>/Wave_loads_for_heading2'
//    '<S183>/Wave_loads_for_heading1'
//    '<S183>/Wave_loads_for_heading2'
//    '<S267>/Wave_loads_for_heading1'
//    '<S267>/Wave_loads_for_heading2'
//
void AH_Model_v1ModelClass::AHV_Mod_Wave_loads_for_heading1(const real_T
  rtu_eta[6], real_T rtu_psi, const real_T rtu_psi_wave[900], const real_T
  rtu_wavenum[900], const real_T rtu_Omega[900], const real_T rtu_Phase[900],
  const real_T rtu_Zeta_a[900], real_T rty_tau_WF[6], real_T rty_tau_WD[6],
  B_Wave_loads_for_heading1_AHV_T *localB, const ConstB_Wave_loads_for_heading_T
  *localC, DW_Wave_loads_for_heading1_AH_T *localDW)
{
  real_T frac;
  uint32_T bpIndices[3];
  real_T fractions[3];
  uint32_T bpIndices_0[3];
  real_T fractions_0[3];
  uint32_T bpIndices_1[3];
  real_T fractions_1[3];
  uint32_T bpIndices_2[3];
  real_T fractions_2[3];
  uint32_T bpIndices_3[3];
  real_T fractions_3[3];
  uint32_T bpIndices_4[3];
  real_T fractions_4[3];
  uint32_T bpIndices_5[3];
  real_T fractions_5[3];
  uint32_T bpIndices_6[3];
  real_T fractions_6[3];
  uint32_T bpIndices_7[3];
  real_T fractions_7[3];
  uint32_T bpIndices_8[3];
  real_T fractions_8[3];
  uint32_T bpIndices_9[3];
  real_T fractions_9[3];
  uint32_T bpIndices_a[3];
  real_T fractions_a[3];
  uint32_T bpIndices_b[3];
  real_T fractions_b[3];
  uint32_T bpIndices_c[3];
  real_T fractions_c[3];
  uint32_T bpIndices_d[3];
  real_T fractions_d[3];
  int32_T tmp;
  int32_T i_1_N;
  real_T rtb_LookupTablenDAmp3;
  real_T rtb_phase;
  real_T rtb_LookupTablenDPhase5;
  real_T rtb_LookupTablenDPhase;
  real_T rtb_LookupTablenDPhase2;
  real_T rtb_LookupTablenDPhase1;
  real_T Memoryampzeta_PreviousInput[6];
  real_T Memory2_PreviousInput[6];
  int32_T i;
  real_T rtb_LookupTablenDPhase1_0[6];
  uint32_T u0;
  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    for (i = 0; i < 900; i++) {
      // Trigonometry: '<S40>/sin'
      localB->sin_h[i] = std::sin(rtu_psi_wave[i]);

      // Trigonometry: '<S40>/cos'
      localB->cos_d[i] = std::cos(rtu_psi_wave[i]);

      // Product: '<S35>/Product1' incorporates:
      //   Delay: '<S41>/Delay1'

      localB->clock_Omega[i] = rtu_Omega[i] * localDW->Delay1_DSTATE;
    }
  }

  for (i = 0; i < 900; i++) {
    // Sum: '<S35>/Sum' incorporates:
    //   Product: '<S40>/Product2'
    //   Product: '<S40>/Product6'
    //   Product: '<S40>/Product7'
    //   Sum: '<S40>/Sum5'

    localB->Phase_tot[i] = (localB->clock_Omega[i] - (localB->sin_h[i] *
      rtu_eta[1] + localB->cos_d[i] * rtu_eta[0]) * rtu_wavenum[i]) +
      rtu_Phase[i];

    // Sum: '<S35>/Sum1'
    localB->psi_r[i] = rtu_psi_wave[i] - rtu_psi;
  }

  if (rtmIsMajorTimeStep((&AHV_Model_M))) {
    // Outputs for Iterator SubSystem: '<S35>/Wave component loop' incorporates:
    //   ForIterator: '<S42>/For Iterator(i=1:N)'

    for (i = 0; i < 6; i++) {
      // InitializeConditions for Memory: '<S42>/Memory(amp*zeta)'
      Memoryampzeta_PreviousInput[i] = 0.0;

      // InitializeConditions for Memory: '<S42>/Memory2'
      Memory2_PreviousInput[i] = 0.0;
    }

    if (localC->number_of_iterations < 2.147483648E+9) {
      if (localC->number_of_iterations >= -2.147483648E+9) {
        tmp = static_cast<int32_T>(localC->number_of_iterations);
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

    // Sum: '<S42>/Sum5' incorporates:
    //   Memory: '<S42>/Memory2'
    //   Selector: '<S42>/Omega_Selector'
    //   Selector: '<S42>/Phase_tot_Selector'
    //   Selector: '<S42>/Psi_r_Selector'
    //   Selector: '<S42>/Zeta_a_Selector'

    for (i_1_N = 0; i_1_N < tmp; i_1_N++) {
      // Saturate: '<S46>/x_Saturation'
      if (localB->psi_r[i_1_N] > 1.0E+10) {
        rtb_LookupTablenDAmp3 = 1.0E+10;
      } else if (localB->psi_r[i_1_N] < -1.0E+10) {
        rtb_LookupTablenDAmp3 = -1.0E+10;
      } else {
        rtb_LookupTablenDAmp3 = localB->psi_r[i_1_N];
      }

      // End of Saturate: '<S46>/x_Saturation'

      // Signum: '<S46>/x_Sign'
      if (rtb_LookupTablenDAmp3 < 0.0) {
        rtb_phase = -1.0;
      } else if (rtb_LookupTablenDAmp3 > 0.0) {
        rtb_phase = 1.0;
      } else if (rtb_LookupTablenDAmp3 == 0.0) {
        rtb_phase = 0.0;
      } else {
        rtb_phase = (rtNaN);
      }

      // End of Signum: '<S46>/x_Sign'

      // Gain: '<S46>/pi'
      rtb_phase *= 3.1415926535897931;

      // Sum: '<S46>/Sum1'
      rtb_LookupTablenDAmp3 += rtb_phase;

      // Math: '<S46>/Math Function' incorporates:
      //   Constant: '<S46>/Constant'

      rtb_LookupTablenDAmp3 = rt_remd_snf(rtb_LookupTablenDAmp3,
        6.2831853071795862);

      // Sum: '<S46>/Sum'
      rtb_LookupTablenDAmp3 -= rtb_phase;

      // LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta' incorporates:
      //   Constant: '<S45>/Constant'
      //   Constant: '<S45>/Constant1'
      //   PreLookup: '<S42>/Prelookup1'
      //   Signum: '<S45>/x_Sign'
      //   Sum: '<S45>/Sum'
      //   Switch: '<S45>/Switch'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      if (rtb_LookupTablenDAmp3 < 0.0) {
        rtb_phase = -1.0;
      } else if (rtb_LookupTablenDAmp3 > 0.0) {
        rtb_phase = 1.0;
      } else if (rtb_LookupTablenDAmp3 == 0.0) {
        rtb_phase = 0.0;
      } else {
        rtb_phase = (rtNaN);
      }

      if (rtb_phase >= 0.0) {
        rtb_phase = 0.0;
      } else {
        rtb_phase = 6.2831853071795862;
      }

      u0 = plook_u32d_evencka(rtb_phase + rtb_LookupTablenDAmp3, 0.0,
        0.17453292519943295, 35U);
      if (u0 >= 35U) {
        u0 = 35U;
      }

      // Lookup_n-D: '<S44>/Lookup Table (n-D)  Wavedrift 1' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled16,
        38U, &frac, &localDW->m_bpIndex[0U]);
      fractions[0U] = frac;
      bpIndices[1U] = plook_lincp(AHV_Model_ConstP.pooled14[static_cast<int32_T>
        (u0)], AHV_Model_ConstP.pooled14, 35U, &frac, &localDW->m_bpIndex[1U]);
      fractions[1U] = frac;
      bpIndices[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex[2U]);
      fractions[2U] = frac;
      rtb_LookupTablenDAmp3 = intrp3d_l_pw(bpIndices, fractions,
        AHV_Model_ConstP.pooled15, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S44>/Lookup Table (n-D)  Wavedrift 2' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_0[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled16,
        38U, &frac, &localDW->m_bpIndex_n[0U]);
      fractions_0[0U] = frac;
      bpIndices_0[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_n[1U]);
      fractions_0[1U] = frac;
      bpIndices_0[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_n[2U]);
      fractions_0[2U] = frac;
      rtb_LookupTablenDPhase5 = intrp3d_l_pw(bpIndices_0, fractions_0,
        AHV_Model_ConstP.pooled18, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S44>/Lookup Table (n-D)  Wavedrift 3' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_1[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled16,
        38U, &frac, &localDW->m_bpIndex_l[0U]);
      fractions_1[0U] = frac;
      bpIndices_1[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_l[1U]);
      fractions_1[1U] = frac;
      bpIndices_1[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_l[2U]);
      fractions_1[2U] = frac;
      rtb_phase = intrp3d_l_pw(bpIndices_1, fractions_1,
        AHV_Model_ConstP.pooled19, AHV_Model_ConstP.pooled119);

      // Sum: '<S42>/Sum1' incorporates:
      //   Memory: '<S42>/Memory(amp*zeta)'
      //   Product: '<S42>/amp2*zeta'

      rty_tau_WD[0] = rtb_LookupTablenDAmp3 * rtu_Zeta_a[i_1_N] *
        rtu_Zeta_a[i_1_N] + Memoryampzeta_PreviousInput[0];
      rty_tau_WD[1] = rtb_LookupTablenDPhase5 * rtu_Zeta_a[i_1_N] *
        rtu_Zeta_a[i_1_N] + Memoryampzeta_PreviousInput[1];

      // Product: '<S42>/amp2*zeta' incorporates:
      //   Constant: '<S44>/Constant'

      rtb_LookupTablenDAmp3 = 0.0 * rtu_Zeta_a[i_1_N] * rtu_Zeta_a[i_1_N];

      // Sum: '<S42>/Sum1' incorporates:
      //   Memory: '<S42>/Memory(amp*zeta)'
      //   Product: '<S42>/amp2*zeta'

      rty_tau_WD[2] = rtb_LookupTablenDAmp3 + Memoryampzeta_PreviousInput[2];
      rty_tau_WD[3] = rtb_LookupTablenDAmp3 + Memoryampzeta_PreviousInput[3];
      rty_tau_WD[4] = rtb_LookupTablenDAmp3 + Memoryampzeta_PreviousInput[4];
      rty_tau_WD[5] = rtb_phase * rtu_Zeta_a[i_1_N] * rtu_Zeta_a[i_1_N] +
        Memoryampzeta_PreviousInput[5];

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp1' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_2[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_h[0U]);
      fractions_2[0U] = frac;
      bpIndices_2[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_h[1U]);
      fractions_2[1U] = frac;
      bpIndices_2[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_h[2U]);
      fractions_2[2U] = frac;
      rtb_phase = intrp3d_l_pw(bpIndices_2, fractions_2,
        AHV_Model_ConstP.pooled20, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp2' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_3[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_b[0U]);
      fractions_3[0U] = frac;
      bpIndices_3[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_b[1U]);
      fractions_3[1U] = frac;
      bpIndices_3[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_b[2U]);
      fractions_3[2U] = frac;
      rtb_LookupTablenDPhase5 = intrp3d_l_pw(bpIndices_3, fractions_3,
        AHV_Model_ConstP.pooled22, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp3' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_4[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_hm[0U]);
      fractions_4[0U] = frac;
      bpIndices_4[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_hm[1U]);
      fractions_4[1U] = frac;
      bpIndices_4[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_hm[2U]);
      fractions_4[2U] = frac;
      rtb_LookupTablenDAmp3 = intrp3d_l_pw(bpIndices_4, fractions_4,
        AHV_Model_ConstP.pooled23, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp4' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_5[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_c[0U]);
      fractions_5[0U] = frac;
      bpIndices_5[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_c[1U]);
      fractions_5[1U] = frac;
      bpIndices_5[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_c[2U]);
      fractions_5[2U] = frac;
      rtb_LookupTablenDPhase = intrp3d_l_pw(bpIndices_5, fractions_5,
        AHV_Model_ConstP.pooled24, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp5' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_6[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_ni[0U]);
      fractions_6[0U] = frac;
      bpIndices_6[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_ni[1U]);
      fractions_6[1U] = frac;
      bpIndices_6[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_ni[2U]);
      fractions_6[2U] = frac;
      rtb_LookupTablenDPhase2 = intrp3d_l_pw(bpIndices_6, fractions_6,
        AHV_Model_ConstP.pooled25, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Amp6' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_7[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_nf[0U]);
      fractions_7[0U] = frac;
      bpIndices_7[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_nf[1U]);
      fractions_7[1U] = frac;
      bpIndices_7[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_nf[2U]);
      fractions_7[2U] = frac;
      rtb_LookupTablenDPhase1 = intrp3d_l_pw(bpIndices_7, fractions_7,
        AHV_Model_ConstP.pooled26, AHV_Model_ConstP.pooled119);

      // Product: '<S42>/phase*amp1'
      rty_tau_WF[0] = rtb_phase * rtu_Zeta_a[i_1_N];
      rty_tau_WF[1] = rtb_LookupTablenDPhase5 * rtu_Zeta_a[i_1_N];
      rty_tau_WF[2] = rtb_LookupTablenDAmp3 * rtu_Zeta_a[i_1_N];
      rty_tau_WF[3] = rtb_LookupTablenDPhase * rtu_Zeta_a[i_1_N];
      rty_tau_WF[4] = rtb_LookupTablenDPhase2 * rtu_Zeta_a[i_1_N];
      rty_tau_WF[5] = rtb_LookupTablenDPhase1 * rtu_Zeta_a[i_1_N];

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase1' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_8[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_nq[0U]);
      fractions_8[0U] = frac;
      bpIndices_8[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_nq[1U]);
      fractions_8[1U] = frac;
      bpIndices_8[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_nq[2U]);
      fractions_8[2U] = frac;
      rtb_LookupTablenDPhase1 = intrp3d_l_pw(bpIndices_8, fractions_8,
        AHV_Model_ConstP.pooled27, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase2' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_9[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_cl[0U]);
      fractions_9[0U] = frac;
      bpIndices_9[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_cl[1U]);
      fractions_9[1U] = frac;
      bpIndices_9[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_cl[2U]);
      fractions_9[2U] = frac;
      rtb_LookupTablenDPhase2 = intrp3d_l_pw(bpIndices_9, fractions_9,
        AHV_Model_ConstP.pooled28, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_a[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_j[0U]);
      fractions_a[0U] = frac;
      bpIndices_a[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_j[1U]);
      fractions_a[1U] = frac;
      bpIndices_a[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_j[2U]);
      fractions_a[2U] = frac;
      rtb_LookupTablenDPhase = intrp3d_l_pw(bpIndices_a, fractions_a,
        AHV_Model_ConstP.pooled29, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase4' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_b[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_a[0U]);
      fractions_b[0U] = frac;
      bpIndices_b[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_a[1U]);
      fractions_b[1U] = frac;
      bpIndices_b[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_a[2U]);
      fractions_b[2U] = frac;
      rtb_phase = intrp3d_l_pw(bpIndices_b, fractions_b,
        AHV_Model_ConstP.pooled30, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase5' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_c[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_d[0U]);
      fractions_c[0U] = frac;
      bpIndices_c[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_d[1U]);
      fractions_c[1U] = frac;
      bpIndices_c[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_d[2U]);
      fractions_c[2U] = frac;
      rtb_LookupTablenDPhase5 = intrp3d_l_pw(bpIndices_c, fractions_c,
        AHV_Model_ConstP.pooled31, AHV_Model_ConstP.pooled119);

      // Lookup_n-D: '<S43>/Lookup Table (n-D) Phase6' incorporates:
      //   LookupNDDirect: '<S42>/Direct Lookup Table (n-D) Beta'
      //
      //  About '<S42>/Direct Lookup Table (n-D) Beta':
      //   1-dimensional Direct Look-Up returning a Scalar,
      //
      //      Remove protection against out-of-range input in generated code: 'off'

      bpIndices_d[0U] = plook_lincp(rtu_Omega[i_1_N], AHV_Model_ConstP.pooled21,
        38U, &frac, &localDW->m_bpIndex_f[0U]);
      fractions_d[0U] = frac;
      bpIndices_d[1U] = plook_lincp(AHV_Model_ConstP.pooled14
        [static_cast<int32_T>(u0)], AHV_Model_ConstP.pooled14, 35U, &frac,
        &localDW->m_bpIndex_f[1U]);
      fractions_d[1U] = frac;
      bpIndices_d[2U] = plook_lincp(localC->U, AHV_Model_ConstP.pooled17, 1U,
        &frac, &localDW->m_bpIndex_f[2U]);
      fractions_d[2U] = frac;
      rtb_LookupTablenDAmp3 = intrp3d_l_pw(bpIndices_d, fractions_d,
        AHV_Model_ConstP.pooled32, AHV_Model_ConstP.pooled119);
      rtb_LookupTablenDPhase1_0[0] = rtb_LookupTablenDPhase1 + localB->
        Phase_tot[i_1_N];
      rtb_LookupTablenDPhase1_0[1] = rtb_LookupTablenDPhase2 + localB->
        Phase_tot[i_1_N];
      rtb_LookupTablenDPhase1_0[2] = rtb_LookupTablenDPhase + localB->
        Phase_tot[i_1_N];
      rtb_LookupTablenDPhase1_0[3] = rtb_phase + localB->Phase_tot[i_1_N];
      rtb_LookupTablenDPhase1_0[4] = rtb_LookupTablenDPhase5 + localB->
        Phase_tot[i_1_N];
      rtb_LookupTablenDPhase1_0[5] = rtb_LookupTablenDAmp3 + localB->
        Phase_tot[i_1_N];
      for (i = 0; i < 6; i++) {
        // Sum: '<S42>/Sum4' incorporates:
        //   Memory: '<S42>/Memory2'
        //   Product: '<S42>/Product1'
        //   Trigonometry: '<S42>/                '

        rty_tau_WF[i] = std::sin(rtb_LookupTablenDPhase1_0[i]) * rty_tau_WF[i] +
          Memory2_PreviousInput[i];

        // Update for Memory: '<S42>/Memory(amp*zeta)'
        Memoryampzeta_PreviousInput[i] = rty_tau_WD[i];
        Memory2_PreviousInput[i] = rty_tau_WF[i];
      }
    }

    // End of Sum: '<S42>/Sum5'
    // End of Outputs for SubSystem: '<S35>/Wave component loop'

    // Sum: '<S41>/Sum2' incorporates:
    //   Constant: '<S41>/Constant9'
    //   Delay: '<S41>/Delay1'

    localDW->Delay1_DSTATE += 0.05;
  }
}

//
// Outputs for atomic system:
//    '<S8>/Wave loads (U=0)'
//    '<S92>/Wave loads (U=0)'
//    '<S176>/Wave loads (U=0)'
//    '<S260>/Wave loads (U=0)'
//
void AH_Model_v1ModelClass::Wave_loads_fun(const real_T rtu_eta[6], const real_T
  rtu_waves[900], const real_T rtu_waves_p[900], const real_T rtu_waves_c[900],
  const real_T rtu_waves_m[900], const real_T rtu_waves_mv[900], real_T
  rty_tau_WF[6], real_T rty_tau_WD[6], B_Wave_loads_fun_T *localB, const
  ConstB_Wave_loads_fun_T *localC, DW_Wave_loads_fun_T *localDW)
{
  uint32_T rtb_k;
  real_T rtb_f;
  uint8_T rtb_Compare;
  real_T rtb_Product_g;
  int32_T i;

  // Gain: '<S15>/R2D'
  rtb_Product_g = 57.295779513082323 * rtu_eta[5];

  // Saturate: '<S37>/Saturation'
  if (rtb_Product_g > 1.0E+10) {
    rtb_Product_g = 1.0E+10;
  } else {
    if (rtb_Product_g < -1.0E+10) {
      rtb_Product_g = -1.0E+10;
    }
  }

  // End of Saturate: '<S37>/Saturation'

  // Signum: '<S37>/Sign'
  if (rtb_Product_g < 0.0) {
    rtb_f = -1.0;
  } else if (rtb_Product_g > 0.0) {
    rtb_f = 1.0;
  } else if (rtb_Product_g == 0.0) {
    rtb_f = 0.0;
  } else {
    rtb_f = (rtNaN);
  }

  // End of Signum: '<S37>/Sign'

  // Gain: '<S37>/Gain'
  rtb_f *= 180.0;

  // Sum: '<S37>/Sum1'
  rtb_Product_g += rtb_f;

  // Math: '<S37>/Math Function' incorporates:
  //   Constant: '<S37>/Constant'

  rtb_Product_g = rt_remd_snf(rtb_Product_g, 360.0);

  // Sum: '<S37>/Sum'
  rtb_Product_g -= rtb_f;

  // Abs: '<S15>/Abs_eta'
  rtb_f = std::abs(rtb_Product_g);

  // PreLookup: '<S15>/Prelookup [0,180]'
  rtb_k = plook_evenc(rtb_f, 0.0, 10.0, 18U, &rtb_f);

  // Signum: '<S15>/Sign(eta)'
  if (rtb_Product_g < 0.0) {
    rtb_Product_g = -1.0;
  } else if (rtb_Product_g > 0.0) {
    rtb_Product_g = 1.0;
  } else if (rtb_Product_g == 0.0) {
    rtb_Product_g = 0.0;
  } else {
    rtb_Product_g = (rtNaN);
  }

  // End of Signum: '<S15>/Sign(eta)'

  // RelationalOperator: '<S34>/Compare' incorporates:
  //   Constant: '<S34>/Constant'

  rtb_Compare = (rtb_Product_g < 0.0);

  // Outputs for Atomic SubSystem: '<S15>/Wave_loads_for_heading1'
  // Gain: '<S15>/D2R' incorporates:
  //   Constant: '<S38>/const'
  //   DataTypeConversion: '<S15>/Conversion'
  //   Gain: '<S38>/Gain1'
  //   Gain: '<S38>/Gain2'
  //   Product: '<S38>/Product'
  //   Sum: '<S38>/Sum7'
  //   Sum: '<S38>/Sum8'

  AHV_Mod_Wave_loads_for_heading1(rtu_eta, 0.017453292519943295 * (10.0 *
    static_cast<real_T>(rtb_k) + 20.0 * (18.0 - static_cast<real_T>(rtb_k)) *
    static_cast<real_T>(rtb_Compare)), rtu_waves, rtu_waves_p, rtu_waves_c,
    rtu_waves_m, rtu_waves_mv, localB->tau_WF_j, localB->tau_WD_k,
    &localB->Wave_loads_for_heading1, &localC->Wave_loads_for_heading1,
    &localDW->Wave_loads_for_heading1);

  // End of Outputs for SubSystem: '<S15>/Wave_loads_for_heading1'

  // Outputs for Atomic SubSystem: '<S15>/Wave_loads_for_heading2'
  // Gain: '<S15>/D2R ' incorporates:
  //   Constant: '<S15>/one'
  //   Constant: '<S39>/const'
  //   DataTypeConversion: '<S15>/Conversion'
  //   Gain: '<S39>/Gain1'
  //   Gain: '<S39>/Gain2'
  //   Product: '<S39>/Product'
  //   Sum: '<S15>/Sum2'
  //   Sum: '<S39>/Sum7'
  //   Sum: '<S39>/Sum8'

  AHV_Mod_Wave_loads_for_heading1(rtu_eta, 0.017453292519943295 * (10.0 * (
    static_cast<real_T>(rtb_k) + 1.0) + 20.0 * (18.0 - (static_cast<real_T>
    (rtb_k) + 1.0)) * static_cast<real_T>(rtb_Compare)), rtu_waves, rtu_waves_p,
    rtu_waves_c, rtu_waves_m, rtu_waves_mv, localB->tau_WF, localB->tau_WD,
    &localB->Wave_loads_for_heading2, &localC->Wave_loads_for_heading2,
    &localDW->Wave_loads_for_heading2);

  // End of Outputs for SubSystem: '<S15>/Wave_loads_for_heading2'
  for (i = 0; i < 6; i++) {
    // Sum: '<S15>/Sum5' incorporates:
    //   Constant: '<S15>/one '
    //   Product: '<S15>/(1-f)*tau_WD'
    //   Product: '<S15>/f*tau_WD'
    //   Sum: '<S15>/Sum4'

    rty_tau_WD[i] = (1.0 - rtb_f) * localB->tau_WD_k[i] + rtb_f * localB->
      tau_WD[i];

    // Sum: '<S15>/Sum1' incorporates:
    //   Constant: '<S15>/one '
    //   Product: '<S15>/(1-f)*tau_WF'
    //   Product: '<S15>/f*tau_WF'
    //   Sum: '<S15>/Sum4'

    rty_tau_WF[i] = (1.0 - rtb_f) * localB->tau_WF_j[i] + rtb_f * localB->
      tau_WF[i];
  }
}

//
// File trailer for generated code.
//
// [EOF]
//
