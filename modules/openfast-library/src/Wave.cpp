//
// File: Wave.cpp
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
#include "Wave.h"

// Include model header file for global data
#include "AHV_Model.h"
#include "AHV_Model_private.h"

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else {
    if (u1 < 0.0) {
      u1_0 = std::ceil(u1);
    } else {
      u1_0 = std::floor(u1);
    }

    if ((u1 != 0.0) && (u1 != u1_0)) {
      u1_0 = std::abs(u0 / u1);
      if (!(std::abs(u1_0 - std::floor(u1_0 + 0.5)) > DBL_EPSILON * u1_0)) {
        y = 0.0 * u0;
      } else {
        y = std::fmod(u0, u1);
      }
    } else {
      y = std::fmod(u0, u1);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S16>/Wave'
real_T AH_Model_v1ModelClass::AHV_Model_rad2pipi2(real_T angle)
{
  real_T y;
  real_T y_tmp;
  real_T u;
  if (angle < 0.0) {
    y_tmp = -1.0;
  } else if (angle > 0.0) {
    y_tmp = 1.0;
  } else if (angle == 0.0) {
    y_tmp = 0.0;
  } else {
    y_tmp = (rtNaN);
  }

  u = std::abs(rt_remd_snf(angle + 3.1415926535897931, 6.2831853071795862) /
               6.2831853071795862);
  if (u > 0.0) {
    u = 1.0;
  } else if (u == 0.0) {
    u = 0.0;
  } else {
    u = (rtNaN);
  }

  u = (u - 1.0) * 2.0 + y_tmp;
  if (u < 0.0) {
    u = -1.0;
  } else if (u > 0.0) {
    u = 1.0;
  } else if (u == 0.0) {
    u = 0.0;
  } else {
    u = (rtNaN);
  }

  y = rt_remd_snf(y_tmp * 3.1415926535897931 + angle, 6.2831853071795862) - u *
    3.1415926535897931;
  return y;
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void AH_Model_v1ModelClass::AHV_Model_emxInit_real_T(emxArray_real_T_AHV_Model_T
  **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T_AHV_Model_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T_AHV_Model_T *)std::malloc(sizeof
    (emxArray_real_T_AHV_Model_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)std::malloc(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void AH_Model_v1ModelClass::AHV_Mo_emxEnsureCapacity_real_T
  (emxArray_real_T_AHV_Model_T *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = std::calloc(static_cast<uint32_T>(i), sizeof(real_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(real_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (real_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

// Function for MATLAB Function: '<S16>/Wave'
real_T AH_Model_v1ModelClass::AHV_Model_sum(const real_T x_data[], const int32_T
  x_size[2])
{
  real_T y;
  int32_T k;
  if (x_size[1] == 0) {
    y = 0.0;
  } else {
    y = x_data[0];
    for (k = 2; k <= x_size[1]; k++) {
      y += x_data[k - 1];
    }
  }

  return y;
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_gamma(real_T *x)
{
  real_T n;
  real_T r;
  real_T xnum;
  real_T xden;
  real_T sum;
  int32_T b_i;
  static const real_T gam[23] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
    40320.0, 362880.0, 3.6288E+6, 3.99168E+7, 4.790016E+8, 6.2270208E+9,
    8.71782912E+10, 1.307674368E+12, 2.0922789888E+13, 3.55687428096E+14,
    6.402373705728E+15, 1.21645100408832E+17, 2.43290200817664E+18,
    5.109094217170944E+19, 1.1240007277776077E+21 };

  static const real_T c[7] = { -0.001910444077728, 0.00084171387781295,
    -0.00059523799130430121, 0.0007936507935003503, -0.0027777777777776816,
    0.083333333333333329, 0.0057083835261 };

  static const real_T p[8] = { -1.716185138865495, 24.76565080557592,
    -379.80425647094563, 629.3311553128184, 866.96620279041326,
    -31451.272968848367, -36144.413418691176, 66456.143820240541 };

  static const real_T q[8] = { -30.840230011973897, 315.35062697960416,
    -1015.1563674902192, -3107.7716715723109, 22538.11842098015,
    4755.8462775278813, -134659.95986496931, -115132.25967555349 };

  if ((*x >= 1.0) && (*x <= 23.0) && (*x == std::floor(*x))) {
    *x = gam[static_cast<int32_T>(*x) - 1];
  } else if ((*x < 1.0) && (*x == std::floor(*x))) {
    *x = (rtInf);
  } else if (rtIsNaN(*x)) {
    *x = (rtNaN);
  } else if (rtIsInf(*x)) {
    *x = (rtInf);
  } else {
    n = 0.0;
    if (*x < 12.0) {
      sum = *x;
      if (*x < 1.0) {
        r = *x;
        (*x)++;
      } else {
        n = std::floor(*x) - 1.0;
        *x -= n;
        r = *x - 1.0;
      }

      xnum = 0.0 * r;
      xden = 1.0;
      for (b_i = 0; b_i < 8; b_i++) {
        xnum = (xnum + p[b_i]) * r;
        xden = xden * r + q[b_i];
      }

      r = xnum / xden + 1.0;
      if (sum < *x) {
        r /= sum;
      } else {
        if (sum > *x) {
          for (b_i = 0; b_i < static_cast<int32_T>(n); b_i++) {
            r *= *x;
            (*x)++;
          }
        }
      }
    } else {
      n = *x * *x;
      sum = 0.0057083835261;
      for (b_i = 0; b_i < 6; b_i++) {
        sum = sum / n + c[b_i];
      }

      sum = (sum / *x - *x) + 0.91893853320467278;
      sum += (*x - 0.5) * std::log(*x);
      r = std::exp(sum);
    }

    *x = r;
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_power(const real_T a_data[], const int32_T
  a_size[2], real_T y_data[], int32_T y_size[2])
{
  real_T z1_data[127];
  int32_T loop_ub;
  y_size[1] = static_cast<uint8_T>(a_size[1]);
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&z1_data[0], &y_data[0], (loop_ub + 1) * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub < y_size[1]; loop_ub++) {
    z1_data[loop_ub] = rt_powd_snf(a_data[loop_ub], 2.0);
  }

  y_size[0] = 1;
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&y_data[0], &z1_data[0], (loop_ub + 1) * sizeof(real_T));
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_exp(real_T x_data[], const int32_T x_size
  [2])
{
  int32_T k;
  for (k = 0; k < x_size[1]; k++) {
    x_data[k] = std::exp(x_data[k]);
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_power_pn(const real_T a_data[], const
  int32_T a_size[2], real_T b, real_T y_data[], int32_T y_size[2])
{
  real_T z1_data[127];
  int32_T loop_ub;
  y_size[1] = static_cast<uint8_T>(a_size[1]);
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&z1_data[0], &y_data[0], (loop_ub + 1) * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub < y_size[1]; loop_ub++) {
    z1_data[loop_ub] = rt_powd_snf(a_data[loop_ub], b);
  }

  y_size[0] = 1;
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&y_data[0], &z1_data[0], (loop_ub + 1) * sizeof(real_T));
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_power_p(real_T a, const real_T b_data[],
  const int32_T b_size[2], real_T y_data[], int32_T y_size[2])
{
  real_T z1_data[127];
  int32_T loop_ub;
  y_size[1] = static_cast<uint8_T>(b_size[1]);
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&z1_data[0], &y_data[0], (loop_ub + 1) * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub < y_size[1]; loop_ub++) {
    z1_data[loop_ub] = rt_powd_snf(a, b_data[loop_ub]);
  }

  y_size[0] = 1;
  loop_ub = y_size[1] - 1;
  if (0 <= loop_ub) {
    std::memcpy(&y_data[0], &z1_data[0], (loop_ub + 1) * sizeof(real_T));
  }
}

void AH_Model_v1ModelClass::AHV_Model_emxFree_real_T(emxArray_real_T_AHV_Model_T
  **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T_AHV_Model_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_real_T_AHV_Model_T *)NULL;
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_torset_spec(real_T Hs, real_T wo, const
  real_T omg_data[], const int32_T omg_size[2], real_T S_data[], int32_T *S_size)
{
  real_T Tp;
  real_T tf;
  real_T hsw;
  real_T hss;
  real_T tps;
  real_T mw;
  real_T g_argw;
  real_T g_args;
  real_T g0w;
  real_T agammaw;
  real_T epsu;
  real_T rps;
  real_T f_data[127];
  real_T fnw_data[127];
  int32_T in;
  real_T ftest1_data[127];
  real_T b_S_data[127];
  boolean_T x_data[127];
  uint8_T ii_data[127];
  int32_T idx;
  emxArray_real_T_AHV_Model_T *y_tmp;
  emxArray_real_T_AHV_Model_T *b_y_tmp;
  int32_T loop_ub;
  real_T tmp_data[127];
  real_T tmp_data_0[127];
  real_T tmp_data_1[127];
  int32_T f_size[2];
  int32_T fnw_size[2];
  int32_T ftest1_size[2];
  int32_T b_S_size[2];
  int32_T b_S_size_0[2];
  int32_T y_tmp_size[2];
  int32_T tmp_size[2];
  int32_T y_tmp_size_0[2];
  int32_T b_y_tmp_size[2];
  int32_T b_y_tmp_size_0[2];
  int32_T x_size_idx_1;
  uint8_T k_idx_1;
  boolean_T exitg1;
  Tp = 6.2831853071795862 / wo;
  if (Hs > 0.0) {
    tf = 6.6 * rt_powd_snf(Hs, 0.33333333333333331);
    if (Tp < tf) {
      rps = (tf - Tp) / (tf - 2.0 * std::sqrt(Hs)) / 0.5;
      rps = std::exp(-(rps * rps)) * 0.30000000000000004 + 0.7;
      hsw = rps * Hs;
      hss = std::sqrt(1.0 - rps * rps) * Hs;
      tps = tf + 2.0;
      g_argw = (3.5 * std::exp(-Hs) + 1.0) * 35.0 * rt_powd_snf
        (0.64048779889700158 * hsw / (Tp * Tp), 0.857);
      tf = 1.0;
      rps = 0.5 * std::sqrt(Hs) + 3.2;
      mw = 4.0;
      epsu = (rps - 1.0) / 4.0;
      g_args = epsu;
      AHV_Model_gamma(&g_args);
      g_args = 1.0 / (0.25 * g_args / rt_powd_snf(rps / 4.0, epsu));
      g0w = epsu;
      AHV_Model_gamma(&g0w);
      g0w = 1.0 / (0.25 * g0w / rt_powd_snf(rps / 4.0, epsu));
      agammaw = (4.1 / rt_powd_snf(rps - -2.3514615654177975,
                  0.70561261474570092) * rt_powd_snf(std::log(g_argw),
                  0.5926790422164091 / rt_powd_snf(rps, 0.43870198653690839) +
                  0.86533460692519271) + 1.0) / g_argw;
      epsu = 1.0;
    } else {
      epsu = (Tp - tf) / (25.0 - tf);
      rps = epsu / 0.3;
      rps = std::exp(-(rps * rps)) * 0.4 + 0.6;
      hss = rps * Hs;
      hsw = std::sqrt(1.0 - rps * rps) * Hs;
      tps = Tp;
      rps = 0.5 * std::sqrt(Hs) + 3.2;
      Tp = std::exp(-Hs / 3.0);
      mw = (1.0 - 0.7 * Tp) * 4.0;
      g_argw = (rps - 1.0) / mw;
      g_args = (rps - 1.0) / 4.0;
      g0w = g_args;
      AHV_Model_gamma(&g0w);
      g_args = 1.0 / (0.25 * g0w / rt_powd_snf(rps / 4.0, g_args));
      g0w = g_argw;
      AHV_Model_gamma(&g0w);
      g0w = 1.0 / (1.0 / mw * g0w / rt_powd_snf(rps / mw, g_argw));
      Tp = rt_powd_snf(hsw * hsw * g0w / ((1.0 - Tp) * 0.08 * 16.0 * rt_powd_snf
        (0.4, rps)), 1.0 / (rps - 1.0));
      g_argw = 1.0;
      tf = (3.5 * std::exp(-Hs) + 1.0) * 35.0 * rt_powd_snf(0.64048779889700158 *
        Hs / (tf * tf), 0.857) * (6.0 * epsu + 1.0);
      agammaw = 1.0;
      epsu = (4.1 / rt_powd_snf(rps - -2.3514615654177975, 0.70561261474570092) *
              rt_powd_snf(std::log(tf), 0.5926790422164091 / rt_powd_snf(rps,
                0.43870198653690839) + 0.86533460692519271) + 1.0) / tf;
    }

    hss = hss * hss * tps / 16.0;
    hsw = hsw * hsw * Tp / 16.0;
    loop_ub = omg_size[0] * omg_size[1] - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      f_data[idx] = omg_data[idx] / 6.2831853071795862;
    }

    AHV_Model_emxInit_real_T(&y_tmp, 2);
    idx = y_tmp->size[0] * y_tmp->size[1];
    y_tmp->size[0] = 1;
    y_tmp->size[1] = omg_size[1];
    AHV_Mo_emxEnsureCapacity_real_T(y_tmp, idx);
    loop_ub = omg_size[1] - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      y_tmp->data[idx] = f_data[idx] * Tp;
    }

    x_size_idx_1 = y_tmp->size[1];
    loop_ub = y_tmp->size[0] * y_tmp->size[1] - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      x_data[idx] = (y_tmp->data[idx] < 1.0);
    }

    idx = 0;
    in = 0;
    exitg1 = false;
    while ((!exitg1) && (in <= x_size_idx_1 - 1)) {
      if (x_data[in]) {
        idx++;
        ii_data[idx - 1] = static_cast<uint8_T>(in + 1);
        if (idx >= x_size_idx_1) {
          exitg1 = true;
        } else {
          in++;
        }
      } else {
        in++;
      }
    }

    if (x_size_idx_1 == 1) {
      if (idx == 0) {
        x_size_idx_1 = 0;
      }
    } else if (1 > idx) {
      x_size_idx_1 = 0;
    } else {
      x_size_idx_1 = idx;
    }

    loop_ub = x_size_idx_1 - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      b_S_data[idx] = ii_data[idx];
    }

    if (x_size_idx_1 <= 2) {
      if (x_size_idx_1 == 1) {
        in = static_cast<int32_T>(b_S_data[0]);
      } else {
        loop_ub = x_size_idx_1 - 1;
        for (idx = 0; idx <= loop_ub; idx++) {
          ii_data[idx] = static_cast<uint8_T>(b_S_data[idx]);
        }

        if (ii_data[0] < ii_data[1]) {
          in = ii_data[1];
        } else {
          in = ii_data[0];
        }
      }
    } else {
      loop_ub = x_size_idx_1 - 1;
      for (idx = 0; idx <= loop_ub; idx++) {
        ii_data[idx] = static_cast<uint8_T>(b_S_data[idx]);
      }

      in = ii_data[0];
      for (idx = 1; idx < x_size_idx_1; idx++) {
        if (in < ii_data[idx]) {
          in = ii_data[idx];
        }
      }
    }

    ftest1_size[0] = 1;
    ftest1_size[1] = omg_size[1];
    loop_ub = omg_size[1] - 1;
    if (0 <= loop_ub) {
      std::memset(&ftest1_data[0], 0, (loop_ub + 1) * sizeof(real_T));
    }

    y_tmp_size[0] = 1;
    y_tmp_size[1] = in;
    for (idx = 0; idx < in; idx++) {
      b_S_data[idx] = y_tmp->data[idx] - 1.0;
    }

    AHV_Model_power(b_S_data, y_tmp_size, tmp_data, tmp_size);
    b_S_size[0] = 1;
    b_S_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = -(tmp_data[idx] / 0.0098000000000000014);
    }

    AHV_Model_exp(b_S_data, b_S_size);
    loop_ub = b_S_size[1];
    if (0 <= loop_ub - 1) {
      std::memcpy(&ftest1_data[0], &b_S_data[0], loop_ub * sizeof(real_T));
    }

    if (in + 1 > omg_size[1]) {
      in = 0;
      idx = 0;
      x_size_idx_1 = 0;
    } else {
      idx = omg_size[1];
      x_size_idx_1 = in;
    }

    y_tmp_size_0[0] = 1;
    loop_ub = idx - in;
    y_tmp_size_0[1] = loop_ub;
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = y_tmp->data[in + idx] - 1.0;
    }

    AHV_Model_power(b_S_data, y_tmp_size_0, tmp_data, tmp_size);
    b_S_size[0] = 1;
    b_S_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = -(tmp_data[idx] / 0.0162);
    }

    AHV_Model_exp(b_S_data, b_S_size);
    loop_ub = b_S_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      ftest1_data[x_size_idx_1 + idx] = b_S_data[idx];
    }

    AHV_Model_emxInit_real_T(&b_y_tmp, 2);
    Tp = -(rps / mw);
    g0w *= agammaw;
    idx = b_y_tmp->size[0] * b_y_tmp->size[1];
    b_y_tmp->size[0] = 1;
    b_y_tmp->size[1] = omg_size[1];
    AHV_Mo_emxEnsureCapacity_real_T(b_y_tmp, idx);
    loop_ub = omg_size[1] - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      b_y_tmp->data[idx] = f_data[idx] * tps;
    }

    x_size_idx_1 = b_y_tmp->size[1];
    loop_ub = b_y_tmp->size[0] * b_y_tmp->size[1] - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      x_data[idx] = (b_y_tmp->data[idx] < 1.0);
    }

    idx = 0;
    in = 0;
    exitg1 = false;
    while ((!exitg1) && (in <= x_size_idx_1 - 1)) {
      if (x_data[in]) {
        idx++;
        ii_data[idx - 1] = static_cast<uint8_T>(in + 1);
        if (idx >= x_size_idx_1) {
          exitg1 = true;
        } else {
          in++;
        }
      } else {
        in++;
      }
    }

    if (x_size_idx_1 == 1) {
      if (idx == 0) {
        x_size_idx_1 = 0;
      }
    } else if (1 > idx) {
      x_size_idx_1 = 0;
    } else {
      x_size_idx_1 = idx;
    }

    loop_ub = x_size_idx_1 - 1;
    for (idx = 0; idx <= loop_ub; idx++) {
      b_S_data[idx] = ii_data[idx];
    }

    if (x_size_idx_1 <= 2) {
      if (x_size_idx_1 == 1) {
        in = static_cast<int32_T>(b_S_data[0]);
      } else {
        loop_ub = x_size_idx_1 - 1;
        for (idx = 0; idx <= loop_ub; idx++) {
          ii_data[idx] = static_cast<uint8_T>(b_S_data[idx]);
        }

        if (ii_data[0] < ii_data[1]) {
          in = ii_data[1];
        } else {
          in = ii_data[0];
        }
      }
    } else {
      loop_ub = x_size_idx_1 - 1;
      for (idx = 0; idx <= loop_ub; idx++) {
        ii_data[idx] = static_cast<uint8_T>(b_S_data[idx]);
      }

      in = ii_data[0];
      for (idx = 1; idx < x_size_idx_1; idx++) {
        if (in < ii_data[idx]) {
          in = ii_data[idx];
        }
      }
    }

    f_size[0] = 1;
    f_size[1] = omg_size[1];
    loop_ub = omg_size[1] - 1;
    if (0 <= loop_ub) {
      std::memset(&f_data[0], 0, (loop_ub + 1) * sizeof(real_T));
    }

    b_y_tmp_size[0] = 1;
    b_y_tmp_size[1] = in;
    for (idx = 0; idx < in; idx++) {
      b_S_data[idx] = b_y_tmp->data[idx] - 1.0;
    }

    AHV_Model_power(b_S_data, b_y_tmp_size, tmp_data, tmp_size);
    b_S_size[0] = 1;
    b_S_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = -(tmp_data[idx] / 0.0098000000000000014);
    }

    AHV_Model_exp(b_S_data, b_S_size);
    loop_ub = b_S_size[1];
    if (0 <= loop_ub - 1) {
      std::memcpy(&f_data[0], &b_S_data[0], loop_ub * sizeof(real_T));
    }

    if (in + 1 > omg_size[1]) {
      in = 0;
      idx = 0;
      x_size_idx_1 = 0;
    } else {
      idx = omg_size[1];
      x_size_idx_1 = in;
    }

    b_y_tmp_size_0[0] = 1;
    loop_ub = idx - in;
    b_y_tmp_size_0[1] = loop_ub;
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = b_y_tmp->data[in + idx] - 1.0;
    }

    AHV_Model_power(b_S_data, b_y_tmp_size_0, tmp_data, tmp_size);
    b_S_size[0] = 1;
    b_S_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = -(tmp_data[idx] / 0.0162);
    }

    AHV_Model_exp(b_S_data, b_S_size);
    loop_ub = b_S_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      f_data[x_size_idx_1 + idx] = b_S_data[idx];
    }

    tps = -(rps / 4.0);
    g_args *= epsu;
    AHV_Model_power_pn(y_tmp->data, y_tmp->size, -mw, tmp_data, tmp_size);
    b_S_size[0] = 1;
    b_S_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_S_data[idx] = Tp * tmp_data[idx];
    }

    AHV_Model_exp(b_S_data, b_S_size);
    AHV_Model_power_pn(b_y_tmp->data, b_y_tmp->size, -4.0, tmp_data, tmp_size);
    fnw_size[0] = 1;
    fnw_size[1] = tmp_size[1];
    loop_ub = tmp_size[0] * tmp_size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      fnw_data[idx] = tps * tmp_data[idx];
    }

    AHV_Model_exp(fnw_data, fnw_size);
    AHV_Model_power_pn(y_tmp->data, y_tmp->size, -rps, tmp_data, tmp_size);
    AHV_Model_power_p(g_argw, ftest1_data, ftest1_size, tmp_data_0, y_tmp_size);
    AHV_Model_power_pn(b_y_tmp->data, b_y_tmp->size, -rps, ftest1_data,
                       y_tmp_size);
    AHV_Model_power_p(tf, f_data, f_size, tmp_data_1, y_tmp_size);
    loop_ub = tmp_size[0] * tmp_size[1] - 1;
    b_S_size[1] = tmp_size[1];
    AHV_Model_emxFree_real_T(&b_y_tmp);
    AHV_Model_emxFree_real_T(&y_tmp);
    for (idx = 0; idx <= loop_ub; idx++) {
      b_S_data[idx] = tmp_data[idx] * b_S_data[idx] * g0w * tmp_data_0[idx] *
        hsw / 6.2831853071795862 + ftest1_data[idx] * fnw_data[idx] * g_args *
        tmp_data_1[idx] * hss / 6.2831853071795862;
    }
  } else {
    k_idx_1 = static_cast<uint8_T>(omg_size[1]);
    b_S_size[1] = k_idx_1;
    loop_ub = k_idx_1 - 1;
    if (0 <= loop_ub) {
      std::memset(&b_S_data[0], 0, (loop_ub + 1) * sizeof(real_T));
    }
  }

  b_S_size_0[0] = 1;
  b_S_size_0[1] = b_S_size[1];
  loop_ub = b_S_size[1];
  if (0 <= loop_ub - 1) {
    std::memset(&f_data[0], 0, loop_ub * sizeof(real_T));
  }

  if (AHV_Model_sum(f_data, b_S_size_0) != 0.0) {
    k_idx_1 = static_cast<uint8_T>(omg_size[1]);
    b_S_size[1] = k_idx_1;
    loop_ub = k_idx_1 - 1;
    if (0 <= loop_ub) {
      std::memset(&b_S_data[0], 0, (loop_ub + 1) * sizeof(real_T));
    }
  }

  *S_size = b_S_size[1];
  loop_ub = b_S_size[1];
  if (0 <= loop_ub - 1) {
    std::memcpy(&S_data[0], &b_S_data[0], loop_ub * sizeof(real_T));
  }
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_wavespec2(real_T SpecType, const real_T
  Par_data[], const real_T W_data[], const int32_T *W_size,
  emxArray_real_T_AHV_Model_T *S)
{
  real_T A;
  real_T B;
  real_T b_gamma;
  int32_T b_k;
  real_T b_S_data[127];
  real_T c_a;
  int32_T loop_ub;
  real_T W_data_0[127];
  int32_T W_size_0[2];
  real_T W_data_1;
  b_k = S->size[0] * S->size[1];
  S->size[0] = 1;
  S->size[1] = *W_size;
  AHV_Mo_emxEnsureCapacity_real_T(S, b_k);
  loop_ub = *W_size - 1;
  for (b_k = 0; b_k <= loop_ub; b_k++) {
    S->data[b_k] = 0.0;
  }

  switch (static_cast<int32_T>(SpecType)) {
   case 1:
    for (b_k = 0; b_k < *W_size; b_k++) {
      S->data[b_k] = Par_data[0] * rt_powd_snf(W_data[b_k], -5.0) * std::exp
        (-Par_data[1] / rt_powd_snf(W_data[b_k], 4.0));
    }
    break;

   case 2:
    B = rt_powd_snf(9.81 / Par_data[0], 4.0) * 0.74;
    for (b_k = 0; b_k < *W_size; b_k++) {
      S->data[b_k] = 0.77951241 * rt_powd_snf(W_data[b_k], -5.0) * std::exp(-B /
        rt_powd_snf(W_data[b_k], 4.0));
    }
    break;

   case 3:
    B = rt_powd_snf(Par_data[1], 4.0);
    A = Par_data[0] * Par_data[0] * 487.0 / B;
    B = 1949.0 / B;
    for (b_k = 0; b_k < *W_size; b_k++) {
      S->data[b_k] = A * rt_powd_snf(W_data[b_k], -5.0) * std::exp(-B /
        rt_powd_snf(W_data[b_k], 4.0));
    }
    break;

   case 4:
    B = rt_powd_snf(Par_data[1], 4.0);
    A = Par_data[0] * Par_data[0] * 173.0 / B;
    B = 691.0 / B;
    for (b_k = 0; b_k < *W_size; b_k++) {
      S->data[b_k] = A * rt_powd_snf(W_data[b_k], -5.0) * std::exp(-B /
        rt_powd_snf(W_data[b_k], 4.0));
    }
    break;

   case 5:
    B = rt_powd_snf(Par_data[1], 4.0);
    A = Par_data[0] * Par_data[0] * 123.0 / B;
    B = 495.0 / B;
    for (b_k = 0; b_k < *W_size; b_k++) {
      S->data[b_k] = A * rt_powd_snf(W_data[b_k], -5.0) * std::exp(-B /
        rt_powd_snf(W_data[b_k], 4.0));
    }
    break;

   case 6:
    B = 9.81 * Par_data[1] / (Par_data[0] * Par_data[0]);
    A = 9.81 / Par_data[0] * 3.5 * rt_powd_snf(B, -0.33) * 6.2831853071795862;
    B = 0.076 * rt_powd_snf(B, -0.22);
    for (b_k = 0; b_k < *W_size; b_k++) {
      b_gamma = W_data[b_k] - A;
      if (W_data[b_k] < A) {
        W_data_1 = 0.07;
      } else {
        W_data_1 = 0.09;
      }

      W_data_1 *= A;
      S->data[b_k] = B * 96.236100000000008 * rt_powd_snf(W_data[b_k], -5.0) *
        std::exp(rt_powd_snf(A / W_data[b_k], 4.0) * -1.25) * rt_powd_snf(3.3,
        std::exp(-(b_gamma * b_gamma) / (W_data_1 * W_data_1 * 2.0)));
    }
    break;

   case 7:
    A = Par_data[1];
    b_gamma = Par_data[2];
    B = Par_data[0] * Par_data[0] * 0.2 * rt_powd_snf(Par_data[1], 4.0) /
      96.236100000000008;
    if ((Par_data[2] < 1.0) || (Par_data[2] > 7.0)) {
      W_data_1 = 6.2831853071795862 / (Par_data[1] * std::sqrt(Par_data[0]));
      if (W_data_1 <= 3.6) {
        b_gamma = 5.0;
      } else if (W_data_1 <= 5.0) {
        b_gamma = std::exp(5.75 - 1.15 * W_data_1);
      } else {
        b_gamma = 1.0;
      }
    }

    for (b_k = 0; b_k < *W_size; b_k++) {
      c_a = W_data[b_k] - A;
      if (W_data[b_k] < A) {
        W_data_1 = 0.07;
      } else {
        W_data_1 = 0.09;
      }

      W_data_1 *= A;
      S->data[b_k] = B * 96.236100000000008 * rt_powd_snf(W_data[b_k], -5.0) *
        std::exp(rt_powd_snf(A / W_data[b_k], 4.0) * -1.25) * (1.0 - 0.287 * std::
        log(b_gamma)) * rt_powd_snf(b_gamma, std::exp(-(c_a * c_a) / (W_data_1 *
        W_data_1 * 2.0)));
    }
    break;

   case 8:
    W_size_0[0] = 1;
    W_size_0[1] = *W_size;
    if (0 <= *W_size - 1) {
      std::memcpy(&W_data_0[0], &W_data[0], *W_size * sizeof(real_T));
    }

    AHV_Model_torset_spec(Par_data[0], Par_data[1], W_data_0, W_size_0, b_S_data,
                          &loop_ub);
    b_k = S->size[0] * S->size[1];
    S->size[0] = loop_ub;
    S->size[1] = 1;
    AHV_Mo_emxEnsureCapacity_real_T(S, b_k);
    loop_ub--;
    for (b_k = 0; b_k <= loop_ub; b_k++) {
      S->data[b_k] = b_S_data[b_k];
    }
    break;
  }
}

// Function for MATLAB Function: '<S16>/Wave'
real_T AH_Model_v1ModelClass::AHV_Model_eml_rand_mt19937ar(uint32_T state[625])
{
  real_T r;
  uint32_T u[2];
  uint32_T mti;
  uint32_T y;
  int32_T kk;
  int32_T k;
  boolean_T b_isvalid;
  int32_T exitg1;
  boolean_T exitg2;

  // ========================= COPYRIGHT NOTICE ============================
  //  This is a uniform (0,1) pseudorandom number generator based on:
  //
  //  A C-program for MT19937, with initialization improved 2002/1/26.
  //  Coded by Takuji Nishimura and Makoto Matsumoto.
  //
  //  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
  //  All rights reserved.
  //
  //  Redistribution and use in source and binary forms, with or without
  //  modification, are permitted provided that the following conditions
  //  are met:
  //
  //    1. Redistributions of source code must retain the above copyright
  //       notice, this list of conditions and the following disclaimer.
  //
  //    2. Redistributions in binary form must reproduce the above copyright
  //       notice, this list of conditions and the following disclaimer
  //       in the documentation and/or other materials provided with the
  //       distribution.
  //
  //    3. The names of its contributors may not be used to endorse or
  //       promote products derived from this software without specific
  //       prior written permission.
  //
  //  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  //  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  //  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  //  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
  //  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  //  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  //  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  //  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  //  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  //  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  //  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  //
  // =============================   END   =================================
  do {
    exitg1 = 0;
    for (k = 0; k < 2; k++) {
      mti = state[624] + 1U;
      if (mti >= 625U) {
        for (kk = 0; kk < 227; kk++) {
          y = (state[kk + 1] & 2147483647U) | (state[kk] & 2147483648U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }

          state[kk] = state[kk + 397] ^ y;
        }

        for (kk = 0; kk < 396; kk++) {
          y = (state[kk + 227] & 2147483648U) | (state[kk + 228] & 2147483647U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }

          state[kk + 227] = state[kk] ^ y;
        }

        y = (state[623] & 2147483648U) | (state[0] & 2147483647U);
        if ((y & 1U) == 0U) {
          y >>= 1U;
        } else {
          y = y >> 1U ^ 2567483615U;
        }

        state[623] = state[396] ^ y;
        mti = 1U;
      }

      y = state[static_cast<int32_T>(mti) - 1];
      state[624] = mti;
      y ^= y >> 11U;
      y ^= y << 7U & 2636928640U;
      y ^= y << 15U & 4022730752U;
      u[k] = y >> 18U ^ y;
    }

    r = (static_cast<real_T>(u[0] >> 5U) * 6.7108864E+7 + static_cast<real_T>(u
          [1] >> 6U)) * 1.1102230246251565E-16;
    if (r == 0.0) {
      b_isvalid = ((state[624] >= 1U) && (state[624] < 625U));
      if (b_isvalid) {
        b_isvalid = false;
        k = 1;
        exitg2 = false;
        while ((!exitg2) && (k < 625)) {
          if (state[k - 1] == 0U) {
            k++;
          } else {
            b_isvalid = true;
            exitg2 = true;
          }
        }
      }

      if (!b_isvalid) {
        mti = 5489U;
        state[0] = 5489U;
        for (k = 0; k < 623; k++) {
          mti = ((mti >> 30U ^ mti) * 1812433253U + k) + 1U;
          state[k + 1] = mti;
        }

        state[624] = 624U;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return r;
}

// Function for MATLAB Function: '<S16>/Wave'
void AH_Model_v1ModelClass::AHV_Model_rand(real_T varargin_2,
  emxArray_real_T_AHV_Model_T *r, DW_Wave_init_T *localDW)
{
  int32_T c;
  int32_T k;
  c = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = static_cast<int32_T>(varargin_2);
  AHV_Mo_emxEnsureCapacity_real_T(r, c);
  c = r->size[1];
  for (k = 0; k < c; k++) {
    r->data[k] = AHV_Model_eml_rand_mt19937ar(localDW->state);
  }
}

// Function for MATLAB Function: '<S16>/Wave'
real_T AH_Model_v1ModelClass::AHV_Model_fzero(real_T
  FunFcn_tunableEnvironment_f1, const real_T FunFcn_tunableEnvironment_f2[900],
  real_T FunFcn_tunableEnvironment_f3, const real_T x[2])
{
  real_T b;
  real_T a;
  real_T fa;
  real_T fb;
  real_T fc;
  real_T c;
  real_T e;
  real_T d;
  real_T m;
  real_T toler;
  real_T s;
  real_T r;
  int32_T fa_tmp;
  real_T fa_tmp_0;
  boolean_T exitg1;
  a = x[0];
  b = 1.0E+10;
  fa_tmp = static_cast<int32_T>(FunFcn_tunableEnvironment_f3) - 1;
  fa_tmp_0 = FunFcn_tunableEnvironment_f2[fa_tmp] *
    FunFcn_tunableEnvironment_f2[fa_tmp] / 9.81;
  fa = std::tanh(x[0] * FunFcn_tunableEnvironment_f1) * x[0] - fa_tmp_0;
  fb = std::tanh(1.0E+10 * FunFcn_tunableEnvironment_f1) * 1.0E+10 - fa_tmp_0;
  if (fa == 0.0) {
    b = x[0];
  } else {
    if (!(fb == 0.0)) {
      fc = fb;
      c = 1.0E+10;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }

        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }

        m = (c - b) * 0.5;
        toler = std::abs(b);
        if (!(toler > 1.0)) {
          toler = 1.0;
        }

        toler *= 4.4408920985006262E-16;
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            s = fb / fa;
            if (a == c) {
              a = 2.0 * m * s;
              fa = 1.0 - s;
            } else {
              fa /= fc;
              r = fb / fc;
              a = (2.0 * m * fa * (fa - r) - (b - a) * (r - 1.0)) * s;
              fa = (fa - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (a > 0.0) {
              fa = -fa;
            } else {
              a = -a;
            }

            if ((2.0 * a < 3.0 * m * fa - std::abs(toler * fa)) && (a < std::abs
                 (0.5 * e * fa))) {
              e = d;
              d = a / fa;
            } else {
              d = m;
              e = m;
            }
          }

          a = b;
          fa = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }

          fb = std::tanh(b * FunFcn_tunableEnvironment_f1) * b - fa_tmp_0;
        }
      }
    }
  }

  return b;
}

//
// System initialize for atomic system:
//    '<S10>/Subsystem3'
//    '<S106>/Subsystem3'
//    '<S202>/Subsystem3'
//    '<S298>/Subsystem3'
//
void AH_Model_v1ModelClass::Wave_init_Init(DW_Wave_init_T *localDW)
{
  uint32_T r;
  int32_T mti;

  // SystemInitialize for MATLAB Function: '<S16>/Wave'
  std::memset(&localDW->state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  localDW->state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = ((r >> 30U ^ r) * 1812433253U + mti) + 1U;
    localDW->state[mti + 1] = r;
  }

  localDW->state[624] = 624U;

  // End of SystemInitialize for MATLAB Function: '<S16>/Wave'
}

//
// Output and update for atomic system:
//    '<S10>/Subsystem3'
//    '<S106>/Subsystem3'
//    '<S202>/Subsystem3'
//    '<S298>/Subsystem3'
//
void AH_Model_v1ModelClass::Wave_init(real_T rtu_spectrum_type, real_T rtu_hs,
  real_T rtu_omega_peak, real_T rtu_psi_mean, real_T rtu_gamma, real_T
  rtu_spread, real_T rtu_depth, real_T rtu_nfreq, real_T rtu_ndir, real_T
  rtu_energylim, real_T rtu_freq_cutoff, real_T rtu_dir_cutoff, real_T
  rtu_rand_freq, real_T rtu_rand_dir, real_T rty_Zeta_a[900], real_T rty_Omega
  [900], real_T rty_Phase[900], real_T rty_Wavenum[900], real_T rty_Psi[900],
  real_T *rty_wave_direction, B_Wave_init_T *localB, DW_Wave_init_T *localDW,
  real_T rtp_nfreq, real_T rtp_ndir)
{
  real_T hs_0;
  real_T spread_0;
  real_T depth_0;
  real_T ndir_0;
  real_T energylim_0;
  real_T psi_mean_0;
  int32_T SpecType;
  real_T Par_data[3];
  real_T delta_omega;
  real_T Omega_vec_data[127];
  real_T delta_psi;
  real_T psi_start;
  real_T psi_max;
  emxArray_real_T_AHV_Model_T *psi_loop;
  emxArray_real_T_AHV_Model_T *S_vec;
  emxArray_real_T_AHV_Model_T *rand_omega_vector;
  emxArray_real_T_AHV_Model_T *phase_vector;
  emxArray_real_T_AHV_Model_T *rand_psi_vector;
  real_T r;
  real_T k;
  emxArray_real_T_AHV_Model_T *b_a;
  int32_T nm1d2;
  int32_T e_n;
  int32_T b_k;
  real_T ndbl;
  real_T apnd;
  real_T cdiff;
  real_T Omega_vec_data_0[127];
  real_T Omega[2];
  int32_T Omega_vec_size_idx_1;
  int8_T tmp;
  real_T u0;
  real_T u1;
  static const real_T b[170] = { 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
    40320.0, 362880.0, 3.6288E+6, 3.99168E+7, 4.790016E+8, 6.2270208E+9,
    8.71782912E+10, 1.307674368E+12, 2.0922789888E+13, 3.55687428096E+14,
    6.402373705728E+15, 1.21645100408832E+17, 2.43290200817664E+18,
    5.109094217170944E+19, 1.1240007277776077E+21, 2.5852016738884978E+22,
    6.2044840173323941E+23, 1.5511210043330986E+25, 4.0329146112660565E+26,
    1.0888869450418352E+28, 3.0488834461171384E+29, 8.8417619937397008E+30,
    2.6525285981219103E+32, 8.2228386541779224E+33, 2.6313083693369352E+35,
    8.6833176188118859E+36, 2.9523279903960412E+38, 1.0333147966386144E+40,
    3.7199332678990118E+41, 1.3763753091226343E+43, 5.23022617466601E+44,
    2.0397882081197442E+46, 8.1591528324789768E+47, 3.3452526613163803E+49,
    1.4050061177528798E+51, 6.0415263063373834E+52, 2.6582715747884485E+54,
    1.1962222086548019E+56, 5.5026221598120885E+57, 2.5862324151116818E+59,
    1.2413915592536073E+61, 6.0828186403426752E+62, 3.0414093201713376E+64,
    1.5511187532873822E+66, 8.0658175170943877E+67, 4.2748832840600255E+69,
    2.3084369733924138E+71, 1.2696403353658276E+73, 7.1099858780486348E+74,
    4.0526919504877221E+76, 2.3505613312828789E+78, 1.3868311854568986E+80,
    8.3209871127413916E+81, 5.0758021387722484E+83, 3.1469973260387939E+85,
    1.98260831540444E+87, 1.2688693218588417E+89, 8.2476505920824715E+90,
    5.4434493907744307E+92, 3.6471110918188683E+94, 2.4800355424368305E+96,
    1.711224524281413E+98, 1.197857166996989E+100, 8.5047858856786218E+101,
    6.1234458376886077E+103, 4.4701154615126834E+105, 3.3078854415193856E+107,
    2.4809140811395391E+109, 1.8854947016660498E+111, 1.4518309202828584E+113,
    1.1324281178206295E+115, 8.9461821307829729E+116, 7.1569457046263779E+118,
    5.7971260207473655E+120, 4.75364333701284E+122, 3.9455239697206569E+124,
    3.314240134565352E+126, 2.8171041143805494E+128, 2.4227095383672724E+130,
    2.1077572983795269E+132, 1.8548264225739836E+134, 1.6507955160908452E+136,
    1.4857159644817607E+138, 1.3520015276784023E+140, 1.24384140546413E+142,
    1.1567725070816409E+144, 1.0873661566567424E+146, 1.0329978488239052E+148,
    9.916779348709491E+149, 9.6192759682482062E+151, 9.426890448883242E+153,
    9.33262154439441E+155, 9.33262154439441E+157, 9.4259477598383536E+159,
    9.6144667150351211E+161, 9.9029007164861754E+163, 1.0299016745145622E+166,
    1.0813967582402903E+168, 1.1462805637347078E+170, 1.2265202031961373E+172,
    1.3246418194518284E+174, 1.4438595832024928E+176, 1.5882455415227421E+178,
    1.7629525510902437E+180, 1.9745068572210728E+182, 2.2311927486598123E+184,
    2.5435597334721862E+186, 2.9250936934930141E+188, 3.3931086844518965E+190,
    3.969937160808719E+192, 4.6845258497542883E+194, 5.5745857612076033E+196,
    6.6895029134491239E+198, 8.09429852527344E+200, 9.8750442008335976E+202,
    1.2146304367025325E+205, 1.5061417415111404E+207, 1.8826771768889254E+209,
    2.3721732428800459E+211, 3.0126600184576582E+213, 3.8562048236258025E+215,
    4.9745042224772855E+217, 6.4668554892204716E+219, 8.4715806908788174E+221,
    1.1182486511960039E+224, 1.4872707060906852E+226, 1.9929427461615181E+228,
    2.6904727073180495E+230, 3.6590428819525472E+232, 5.01288874827499E+234,
    6.9177864726194859E+236, 9.6157231969410859E+238, 1.346201247571752E+241,
    1.89814375907617E+243, 2.6953641378881614E+245, 3.8543707171800706E+247,
    5.5502938327393013E+249, 8.0479260574719866E+251, 1.17499720439091E+254,
    1.7272458904546376E+256, 2.5563239178728637E+258, 3.8089226376305671E+260,
    5.7133839564458505E+262, 8.6272097742332346E+264, 1.3113358856834518E+267,
    2.0063439050956811E+269, 3.0897696138473489E+271, 4.7891429014633912E+273,
    7.47106292628289E+275, 1.1729568794264138E+278, 1.8532718694937338E+280,
    2.9467022724950369E+282, 4.714723635992059E+284, 7.5907050539472148E+286,
    1.2296942187394488E+289, 2.0044015765453015E+291, 3.2872185855342945E+293,
    5.423910666131586E+295, 9.0036917057784329E+297, 1.5036165148649983E+300,
    2.5260757449731969E+302, 4.2690680090047027E+304, 7.257415615307994E+306 };

  boolean_T guard1 = false;

  // MATLAB Function: '<S16>/Wave' incorporates:
  //   Constant: '<S16>/Constant'
  //   Constant: '<S16>/Constant1'
  //   Constant: '<S16>/Constant3'
  //   Constant: '<S16>/Constant4'

  *rty_wave_direction = rtu_psi_mean;
  hs_0 = rtu_hs;
  r = rtu_omega_peak;
  depth_0 = rtu_depth;
  k = rtu_nfreq;
  ndir_0 = rtu_ndir;
  energylim_0 = rtu_energylim;
  delta_omega = rtu_freq_cutoff;
  apnd = rtu_dir_cutoff;
  psi_mean_0 = *rty_wave_direction;
  std::memset(&localB->Zeta_a[0], 0, 900U * sizeof(real_T));
  std::memset(&localB->Omega[0], 0, 900U * sizeof(real_T));
  std::memset(&localB->Phase[0], 0, 900U * sizeof(real_T));
  std::memset(&localB->Wavenum[0], 0, 900U * sizeof(real_T));
  std::memset(&localB->Psi[0], 0, 900U * sizeof(real_T));
  SpecType = 0;
  Par_data[0] = 0.0;
  if (rtu_spectrum_type < 4.0) {
    if (rtu_hs < 0.0) {
      hs_0 = 0.0;
    }

    psi_mean_0 = AHV_Model_rad2pipi2(*rty_wave_direction);
    spread_0 = rt_roundd_snf(rtu_spread);
    if (spread_0 < 1.0) {
      spread_0 = 1.0;
    } else {
      if (spread_0 > 5.0) {
        spread_0 = 5.0;
      }
    }

    if (rtu_depth < 0.001) {
      depth_0 = 0.001;
    }

    if (rtu_nfreq <= 0.0) {
      k = 1.0;
    }

    if (rtu_ndir <= 0.0) {
      ndir_0 = 1.0;
    }

    if (rtu_freq_cutoff < 2.0) {
      delta_omega = 2.0;
    }

    if (rtu_dir_cutoff < 0.0) {
      apnd = 0.0;
    } else {
      if (rtu_dir_cutoff > 1.1780972450961724) {
        apnd = 1.1780972450961724;
      }
    }

    if (rtu_omega_peak <= 0.0) {
      r = 6.2831853071795862 / (2.68 * rt_powd_snf(hs_0, 0.54) + 4.883);
    }

    if (rtu_energylim < 0.0) {
      energylim_0 = 0.0;
    } else {
      if (rtu_energylim > 1.0) {
        energylim_0 = rt_roundd_snf(rtu_energylim);
        delta_psi = k * ndir_0;
        if (energylim_0 > delta_psi) {
          energylim_0 = delta_psi;
        }
      }
    }

    delta_omega = delta_omega * r / k;
    delta_psi = rt_roundd_snf(k);
    if (delta_psi < 128.0) {
      if (delta_psi >= -128.0) {
        tmp = static_cast<int8_T>(delta_psi);
        Omega_vec_size_idx_1 = static_cast<int8_T>(delta_psi);
      } else {
        tmp = MIN_int8_T;
        Omega_vec_size_idx_1 = -128;
      }
    } else {
      tmp = MAX_int8_T;
      Omega_vec_size_idx_1 = 127;
    }

    for (nm1d2 = 0; nm1d2 < static_cast<int32_T>(k); nm1d2++) {
      Omega_vec_data[nm1d2] = (static_cast<real_T>(nm1d2) + 1.0) * delta_omega;
    }

    delta_psi = (3.1415926535897931 - 2.0 * apnd) / ndir_0;
    psi_start = ((psi_mean_0 - 1.5707963267948966) + delta_psi / 2.0) + apnd;
    psi_max = ((psi_mean_0 + 1.5707963267948966) - delta_psi / 2.0) - apnd;
    AHV_Model_emxInit_real_T(&psi_loop, 2);
    if (ndir_0 > 1.0) {
      if (rtIsNaN(psi_start)) {
        nm1d2 = psi_loop->size[0] * psi_loop->size[1];
        psi_loop->size[0] = 1;
        psi_loop->size[1] = 1;
        AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
        psi_loop->data[0] = (rtNaN);
      } else if (rtIsNaN(delta_psi)) {
        nm1d2 = psi_loop->size[0] * psi_loop->size[1];
        psi_loop->size[0] = 1;
        psi_loop->size[1] = 1;
        AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
        psi_loop->data[0] = (rtNaN);
      } else if (rtIsNaN(psi_max)) {
        nm1d2 = psi_loop->size[0] * psi_loop->size[1];
        psi_loop->size[0] = 1;
        psi_loop->size[1] = 1;
        AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
        psi_loop->data[0] = (rtNaN);
      } else if ((delta_psi == 0.0) || ((psi_max < psi_start) && (delta_psi >
                   0.0))) {
        psi_loop->size[0] = 1;
        psi_loop->size[1] = 0;
      } else {
        guard1 = false;
        if (rtIsInf(psi_start) || rtIsInf(psi_max)) {
          if (rtIsInf(delta_psi)) {
            nm1d2 = psi_loop->size[0] * psi_loop->size[1];
            psi_loop->size[0] = 1;
            psi_loop->size[1] = 1;
            AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
            psi_loop->data[0] = (rtNaN);
          } else if (psi_start == psi_max) {
            nm1d2 = psi_loop->size[0] * psi_loop->size[1];
            psi_loop->size[0] = 1;
            psi_loop->size[1] = 1;
            AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
            psi_loop->data[0] = (rtNaN);
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          if (rtIsInf(delta_psi)) {
            nm1d2 = psi_loop->size[0] * psi_loop->size[1];
            psi_loop->size[0] = 1;
            psi_loop->size[1] = 1;
            AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
            psi_loop->data[0] = psi_start;
          } else if ((std::floor(psi_start) == psi_start) && (std::floor
                      (delta_psi) == delta_psi)) {
            nm1d2 = psi_loop->size[0] * psi_loop->size[1];
            psi_loop->size[0] = 1;
            e_n = static_cast<int32_T>(std::floor((psi_max - psi_start) /
              delta_psi));
            psi_loop->size[1] = e_n + 1;
            AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
            for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
              psi_loop->data[nm1d2] = delta_psi * static_cast<real_T>(nm1d2) +
                psi_start;
            }
          } else {
            ndbl = std::floor((psi_max - psi_start) / delta_psi + 0.5);
            apnd = ndbl * delta_psi + psi_start;
            if (delta_psi > 0.0) {
              cdiff = apnd - psi_max;
            } else {
              cdiff = psi_max - apnd;
            }

            u0 = std::abs(psi_start);
            u1 = std::abs(psi_max);
            if ((u0 > u1) || rtIsNaN(u1)) {
              u1 = u0;
            }

            if (std::abs(cdiff) < 4.4408920985006262E-16 * u1) {
              ndbl++;
              apnd = psi_max;
            } else if (cdiff > 0.0) {
              apnd = (ndbl - 1.0) * delta_psi + psi_start;
            } else {
              ndbl++;
            }

            if (ndbl >= 0.0) {
              e_n = static_cast<int32_T>(ndbl) - 1;
            } else {
              e_n = -1;
            }

            nm1d2 = psi_loop->size[0] * psi_loop->size[1];
            psi_loop->size[0] = 1;
            psi_loop->size[1] = e_n + 1;
            AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
            if (e_n + 1 > 0) {
              psi_loop->data[0] = psi_start;
              if (e_n + 1 > 1) {
                psi_loop->data[e_n] = apnd;
                nm1d2 = e_n / 2;
                for (b_k = 1; b_k - 1 <= nm1d2 - 2; b_k++) {
                  psi_max = static_cast<real_T>(b_k) * delta_psi;
                  psi_loop->data[b_k] = psi_start + psi_max;
                  psi_loop->data[e_n - b_k] = apnd - psi_max;
                }

                if (nm1d2 << 1 == e_n) {
                  psi_loop->data[nm1d2] = (psi_start + apnd) / 2.0;
                } else {
                  psi_max = static_cast<real_T>(nm1d2) * delta_psi;
                  psi_loop->data[nm1d2] = psi_start + psi_max;
                  psi_loop->data[nm1d2 + 1] = apnd - psi_max;
                }
              }
            }
          }
        }
      }
    } else {
      nm1d2 = psi_loop->size[0] * psi_loop->size[1];
      psi_loop->size[0] = 1;
      psi_loop->size[1] = 1;
      AHV_Mo_emxEnsureCapacity_real_T(psi_loop, nm1d2);
      psi_loop->data[0] = psi_mean_0;
    }

    if (hs_0 != 0.0) {
      psi_start = 4.0 * hs_0;
      psi_start = psi_start * psi_start / (k * ndir_0);
    } else {
      psi_start = 1.0E-5;
    }

    if (rtu_spectrum_type == 1.0) {
      SpecType = 3;
      Par_data[0] = hs_0;
      Par_data[1] = 6.2831853071795862 / r;
    } else if (rtu_spectrum_type == 2.0) {
      SpecType = 7;
      Par_data[0] = hs_0;
      Par_data[1] = r;
      Par_data[2] = rtu_gamma;
    } else {
      if (rtu_spectrum_type == 3.0) {
        SpecType = 8;
        Par_data[0] = hs_0;
        Par_data[1] = r;
      }
    }

    nm1d2 = tmp;
    if (0 <= Omega_vec_size_idx_1 - 1) {
      std::memcpy(&Omega_vec_data_0[0], &Omega_vec_data[0], Omega_vec_size_idx_1
                  * sizeof(real_T));
    }

    AHV_Model_emxInit_real_T(&S_vec, 2);
    AHV_Model_emxInit_real_T(&rand_omega_vector, 2);
    AHV_Model_wavespec2(static_cast<real_T>(SpecType), Par_data,
                        Omega_vec_data_0, &nm1d2, S_vec);
    localB->Zeta_a[0] = -1.0;
    AHV_Model_rand(k * ndir_0, rand_omega_vector, localDW);
    nm1d2 = rand_omega_vector->size[0] * rand_omega_vector->size[1];
    SpecType = rand_omega_vector->size[0] * rand_omega_vector->size[1];
    rand_omega_vector->size[0] = 1;
    AHV_Mo_emxEnsureCapacity_real_T(rand_omega_vector, SpecType);
    e_n = nm1d2 - 1;
    for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
      rand_omega_vector->data[nm1d2] -= 0.5;
    }

    nm1d2 = rand_omega_vector->size[0] * rand_omega_vector->size[1];
    SpecType = rand_omega_vector->size[0] * rand_omega_vector->size[1];
    rand_omega_vector->size[0] = 1;
    AHV_Mo_emxEnsureCapacity_real_T(rand_omega_vector, SpecType);
    e_n = nm1d2 - 1;
    for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
      rand_omega_vector->data[nm1d2] *= delta_omega;
    }

    AHV_Model_emxInit_real_T(&phase_vector, 2);
    AHV_Model_emxInit_real_T(&rand_psi_vector, 2);
    if ((k != rtp_nfreq) || (ndir_0 != rtp_ndir)) {
      AHV_Model_emxInit_real_T(&b_a, 2);
      AHV_Model_rand(k * ndir_0, b_a, localDW);
      nm1d2 = phase_vector->size[0] * phase_vector->size[1];
      phase_vector->size[0] = 1;
      phase_vector->size[1] = b_a->size[1];
      AHV_Mo_emxEnsureCapacity_real_T(phase_vector, nm1d2);
      e_n = b_a->size[0] * b_a->size[1] - 1;
      for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
        phase_vector->data[nm1d2] = b_a->data[nm1d2] * 2.0 * 3.1415926535897931;
      }

      AHV_Model_rand(k * ndir_0, b_a, localDW);
      nm1d2 = b_a->size[0] * b_a->size[1];
      SpecType = b_a->size[0] * b_a->size[1];
      b_a->size[0] = 1;
      AHV_Mo_emxEnsureCapacity_real_T(b_a, SpecType);
      e_n = nm1d2 - 1;
      for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
        b_a->data[nm1d2] -= 0.5;
      }

      nm1d2 = rand_psi_vector->size[0] * rand_psi_vector->size[1];
      rand_psi_vector->size[0] = 1;
      rand_psi_vector->size[1] = b_a->size[1];
      AHV_Mo_emxEnsureCapacity_real_T(rand_psi_vector, nm1d2);
      e_n = b_a->size[0] * b_a->size[1] - 1;
      for (nm1d2 = 0; nm1d2 <= e_n; nm1d2++) {
        rand_psi_vector->data[nm1d2] = delta_psi * b_a->data[nm1d2];
      }

      AHV_Model_emxFree_real_T(&b_a);
    } else {
      nm1d2 = phase_vector->size[0] * phase_vector->size[1];
      phase_vector->size[0] = 200;
      phase_vector->size[1] = 1;
      AHV_Mo_emxEnsureCapacity_real_T(phase_vector, nm1d2);
      for (nm1d2 = 0; nm1d2 < 200; nm1d2++) {
        phase_vector->data[nm1d2] = AHV_Model_ConstP.pooled9[nm1d2];
      }

      nm1d2 = rand_psi_vector->size[0] * rand_psi_vector->size[1];
      rand_psi_vector->size[0] = 200;
      rand_psi_vector->size[1] = 1;
      AHV_Mo_emxEnsureCapacity_real_T(rand_psi_vector, nm1d2);
      for (nm1d2 = 0; nm1d2 < 200; nm1d2++) {
        rand_psi_vector->data[nm1d2] = AHV_Model_ConstP.pooled10[nm1d2];
      }
    }

    r = 0.0;
    k = 0.0;
    apnd = 0.0;
    for (nm1d2 = 0; nm1d2 < tmp; nm1d2++) {
      apnd++;
      psi_max = S_vec->data[static_cast<int32_T>(apnd) - 1];
      for (SpecType = 0; SpecType < psi_loop->size[1]; SpecType++) {
        r++;
        if (ndir_0 > 1.0) {
          ndbl = 2.0 * spread_0 - 1.0;
          if (spread_0 != spread_0) {
            cdiff = (rtNaN);
          } else if (rtIsNaN(spread_0)) {
            cdiff = (rtNaN);
          } else {
            cdiff = b[static_cast<int32_T>(spread_0) - 1];
          }

          if (spread_0 - 1.0 != spread_0 - 1.0) {
            u0 = (rtNaN);
          } else if (rtIsNaN(spread_0 - 1.0)) {
            u0 = (rtNaN);
          } else if (spread_0 - 1.0 < 1.0) {
            u0 = 1.0;
          } else {
            u0 = b[static_cast<int32_T>(spread_0 - 1.0) - 1];
          }

          if (ndbl != ndbl) {
            u1 = (rtNaN);
          } else if (rtIsNaN(ndbl)) {
            u1 = (rtNaN);
          } else {
            u1 = b[static_cast<int32_T>(ndbl) - 1];
          }

          ndbl = rt_powd_snf(2.0, ndbl) * cdiff * u0 / (3.1415926535897931 * u1)
            * rt_powd_snf(std::cos(psi_loop->data[SpecType] - psi_mean_0), 2.0 *
                          spread_0);
        } else {
          ndbl = 1.0 / delta_psi;
        }

        if ((psi_max * ndbl * delta_omega * delta_psi / psi_start >= energylim_0)
            || (energylim_0 > 1.0)) {
          k++;
          Omega_vec_size_idx_1 = static_cast<int32_T>(k) - 1;
          localB->Zeta_a[Omega_vec_size_idx_1] = std::sqrt(2.0 * psi_max * ndbl *
            delta_omega * delta_psi);
          if (rtu_rand_freq == 1.0) {
            localB->Omega[Omega_vec_size_idx_1] = rand_omega_vector->data[
              static_cast<int32_T>(r) - 1] + Omega_vec_data[nm1d2];
          } else {
            localB->Omega[Omega_vec_size_idx_1] = Omega_vec_data[nm1d2];
          }

          if (depth_0 == (rtInf)) {
            localB->Wavenum[Omega_vec_size_idx_1] = localB->
              Omega[Omega_vec_size_idx_1] * localB->Omega[Omega_vec_size_idx_1] /
              9.81;
          } else {
            Omega[0] = localB->Omega[Omega_vec_size_idx_1] * localB->
              Omega[Omega_vec_size_idx_1] / 9.81;
            Omega[1] = 1.0E+10;
            localB->Wavenum[Omega_vec_size_idx_1] = AHV_Model_fzero(depth_0,
              localB->Omega, k, Omega);
          }

          if ((ndir_0 > 1.0) && (rtu_rand_dir == 1.0)) {
            ndbl = psi_loop->data[SpecType];
            cdiff = psi_loop->data[SpecType];
            u0 = std::abs(rt_remd_snf(psi_loop->data[SpecType] +
              3.1415926535897931, 6.2831853071795862) / 6.2831853071795862);
            if (cdiff < 0.0) {
              cdiff = -1.0;
            } else if (cdiff > 0.0) {
              cdiff = 1.0;
            } else if (cdiff == 0.0) {
              cdiff = 0.0;
            } else {
              cdiff = (rtNaN);
            }

            if (u0 > 0.0) {
              u0 = 1.0;
            } else if (u0 == 0.0) {
              u0 = 0.0;
            } else {
              u0 = (rtNaN);
            }

            cdiff += (u0 - 1.0) * 2.0;
            if (ndbl < 0.0) {
              ndbl = -1.0;
            } else if (ndbl > 0.0) {
              ndbl = 1.0;
            } else if (ndbl == 0.0) {
              ndbl = 0.0;
            } else {
              ndbl = (rtNaN);
            }

            if (cdiff < 0.0) {
              cdiff = -1.0;
            } else if (cdiff > 0.0) {
              cdiff = 1.0;
            } else if (cdiff == 0.0) {
              cdiff = 0.0;
            } else {
              cdiff = (rtNaN);
            }

            localB->Psi[Omega_vec_size_idx_1] = (rt_remd_snf(ndbl *
              3.1415926535897931 + psi_loop->data[SpecType], 6.2831853071795862)
              - cdiff * 3.1415926535897931) + rand_psi_vector->data
              [static_cast<int32_T>(r) - 1];
          } else {
            ndbl = psi_loop->data[SpecType];
            cdiff = psi_loop->data[SpecType];
            u0 = std::abs(rt_remd_snf(psi_loop->data[SpecType] +
              3.1415926535897931, 6.2831853071795862) / 6.2831853071795862);
            if (cdiff < 0.0) {
              cdiff = -1.0;
            } else if (cdiff > 0.0) {
              cdiff = 1.0;
            } else if (cdiff == 0.0) {
              cdiff = 0.0;
            } else {
              cdiff = (rtNaN);
            }

            if (u0 > 0.0) {
              u0 = 1.0;
            } else if (u0 == 0.0) {
              u0 = 0.0;
            } else {
              u0 = (rtNaN);
            }

            cdiff += (u0 - 1.0) * 2.0;
            if (ndbl < 0.0) {
              ndbl = -1.0;
            } else if (ndbl > 0.0) {
              ndbl = 1.0;
            } else if (ndbl == 0.0) {
              ndbl = 0.0;
            } else {
              ndbl = (rtNaN);
            }

            if (cdiff < 0.0) {
              cdiff = -1.0;
            } else if (cdiff > 0.0) {
              cdiff = 1.0;
            } else if (cdiff == 0.0) {
              cdiff = 0.0;
            } else {
              cdiff = (rtNaN);
            }

            localB->Psi[Omega_vec_size_idx_1] = rt_remd_snf(ndbl *
              3.1415926535897931 + psi_loop->data[SpecType], 6.2831853071795862)
              - cdiff * 3.1415926535897931;
          }

          localB->Phase[Omega_vec_size_idx_1] = phase_vector->data
            [static_cast<int32_T>(r) - 1];
        }
      }
    }

    AHV_Model_emxFree_real_T(&rand_psi_vector);
    AHV_Model_emxFree_real_T(&phase_vector);
    AHV_Model_emxFree_real_T(&rand_omega_vector);
    AHV_Model_emxFree_real_T(&S_vec);
    AHV_Model_emxFree_real_T(&psi_loop);
    if (localB->Zeta_a[0] < 0.0) {
      if (hs_0 > 0.0) {
        std::memset(&localB->Zeta_a[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Omega[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Phase[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Wavenum[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Psi[0], 0, 900U * sizeof(real_T));
      } else {
        std::memset(&localB->Zeta_a[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Omega[0], 0, 900U * sizeof(real_T));
        localB->Omega[0] = 1.0;
        std::memset(&localB->Phase[0], 0, 900U * sizeof(real_T));
        std::memset(&localB->Wavenum[0], 0, 900U * sizeof(real_T));
        localB->Wavenum[0] = 0.1019367991845056;
        std::memset(&localB->Psi[0], 0, 900U * sizeof(real_T));
        localB->Psi[0] = psi_mean_0;
      }
    }
  }

  *rty_wave_direction = psi_mean_0;
  std::memcpy(&rty_Psi[0], &localB->Psi[0], 900U * sizeof(real_T));
  std::memcpy(&rty_Wavenum[0], &localB->Wavenum[0], 900U * sizeof(real_T));
  std::memcpy(&rty_Phase[0], &localB->Phase[0], 900U * sizeof(real_T));
  std::memcpy(&rty_Omega[0], &localB->Omega[0], 900U * sizeof(real_T));
  std::memcpy(&rty_Zeta_a[0], &localB->Zeta_a[0], 900U * sizeof(real_T));

  // End of MATLAB Function: '<S16>/Wave'
}

//
// File trailer for generated code.
//
// [EOF]
//
