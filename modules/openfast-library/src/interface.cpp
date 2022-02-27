
#include <stdio.h>

#include "AHV_Model.h"  // Model's header file

static AH_Model_v1ModelClass AHV_Model_Obj;  // Instance of model class

double turn[6] = {1, 1, -1, 1, 1, 1};

double V1_Fairlead_pos[3] = {-60, 0, 0};  //导缆孔位置

double V2_Fairlead_pos[3] = {-60, 0, 0};  //导缆孔位置

double V3_Fairlead_pos[3] = {-60, 0, 0};  //导缆孔位置

double V4_Fairlead_pos[3] = {-60, 0, 0};  //导缆孔位置
extern "C" {
void ship_init_fun_(
    __float128* dtime, double* V1_surge_init, double* V1_sway_init,
    double* V1_heave_init, double* V1_roll_init, double* V1_pitch_init,
    double* V1_yaw_init

    ,
    double* V2_surge_init, double* V2_sway_init, double* V2_heave_init,
    double* V2_roll_init, double* V2_pitch_init, double* V2_yaw_init

    ,
    double* V3_surge_init, double* V3_sway_init, double* V3_heave_init,
    double* V3_roll_init, double* V3_pitch_init, double* V3_yaw_init

    ,
    double* V4_surge_init, double* V4_sway_init, double* V4_heave_init,
    double* V4_roll_init, double* V4_pitch_init, double* V4_yaw_init)

{
  //海况参数
  printf("\tInitialization of tugboat mathematical model。\n");
  printf("\ttime: %.36Qg\n", *dtime);

  spectrum_type = 1;

  omega_peak = -1;

  gamma_value = 0;

  spread = 2;

  depth = 287;

  nfreq = 20;

  ndir = 10;

  energylim = 0.005;

  freq_cutoff = 3;

  dir_cutoff = 0;

  rand_freq = 0;

  rand_dir = 0;

  AHV_Model_Obj.initialize(*dtime);

  Vessel_init1[0] = *V1_surge_init;

  Vessel_init1[1] = *V1_sway_init;

  Vessel_init1[2] = *V1_heave_init;

  Vessel_init1[3] = *V1_roll_init;

  Vessel_init1[4] = *V1_pitch_init;

  Vessel_init1[5] = *V1_yaw_init;

  Vessel_init2[0] = *V2_surge_init;

  Vessel_init2[1] = *V2_sway_init;

  Vessel_init2[2] = *V2_heave_init;

  Vessel_init2[3] = *V2_roll_init;

  Vessel_init2[4] = *V2_pitch_init;

  Vessel_init2[5] = *V2_yaw_init;

  Vessel_init3[0] = *V3_surge_init;

  Vessel_init3[1] = *V3_sway_init;

  Vessel_init3[2] = *V3_heave_init;

  Vessel_init3[3] = *V3_roll_init;

  Vessel_init3[4] = *V3_pitch_init;

  Vessel_init3[5] = *V3_yaw_init;

  Vessel_init4[0] = *V4_surge_init;

  Vessel_init4[1] = *V4_sway_init;

  Vessel_init4[2] = *V4_heave_init;

  Vessel_init4[3] = *V4_roll_init;

  Vessel_init4[4] = *V4_pitch_init;

  Vessel_init4[5] = *V4_yaw_init;
}

void shipsetfun1_(double* MeaningWaveHeight, int* heading_contral_mode11,
                  double* WaveDirection, double* CurrentDirection,
                  double* CurrentSpeed, int* V1_IfTugKeepInitPos,

                  double* V1_TugTargetPosX, double* V1_TugTargetPosY,
                  double* heading_angle1, double V1_Flines[])

{
  Hold_Position1 = *V1_IfTugKeepInitPos;

  hs = *MeaningWaveHeight;

  psi_mean = *WaveDirection * 3.14 / 180;

  Current_direction = *CurrentDirection;

  Current_speed = *CurrentSpeed;

  heading_mode1 = *heading_contral_mode11;  // 1自动，0手动

  heading_angle_ref1 = *heading_angle1;

  Vessel_X_Ref1 = *V1_TugTargetPosX;

  Vessel_Y_Ref1 = *V1_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead1[i] = V1_Fairlead_pos[i];

    tau_cable1[i] = V1_Flines[i] * turn[i];
  }
}

void shipsetfun2_(int* V2_IfTugKeepInitPos, int* heading_contral_mode12,
                  double* V2_TugTargetPosX, double* V2_TugTargetPosY,
                  double* heading_angle2, double V2_Flines[])

{
  Hold_Position2 = *V2_IfTugKeepInitPos;

  heading_mode2 = *heading_contral_mode12;

  heading_angle_ref2 = *heading_angle2;

  Vessel_X_Ref2 = *V2_TugTargetPosX;

  Vessel_Y_Ref2 = *V2_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead2[i] = V2_Fairlead_pos[i];

    tau_cable2[i] = V2_Flines[i] * turn[i];
  }
}

void shipsetfun3_(int* V3_IfTugKeepInitPos, int* heading_contral_mode13,
                  double* V3_TugTargetPosX, double* V3_TugTargetPosY,
                  double* heading_angle3, double V3_Flines[])

{
  Hold_Position3 = *V3_IfTugKeepInitPos;

  heading_mode3 = *heading_contral_mode13;

  heading_angle_ref3 = *heading_angle3;

  Vessel_X_Ref3 = *V3_TugTargetPosX;

  Vessel_Y_Ref3 = *V3_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead3[i] = V3_Fairlead_pos[i];

    tau_cable3[i] = V3_Flines[i] * turn[i];
  }
}

void shipsetfun4_(int* V4_IfTugKeepInitPos, int* heading_contral_mode14,
                  double* V4_TugTargetPosX, double* V4_TugTargetPosY,
                  double* heading_angle4, double V4_Flines[])

{
  Hold_Position4 = *V4_IfTugKeepInitPos;

  heading_mode4 = *heading_contral_mode14;

  heading_angle_ref4 = *heading_angle4;

  Vessel_X_Ref4 = *V4_TugTargetPosX;

  Vessel_Y_Ref4 = *V4_TugTargetPosY;

  // printf("\tshipsetfunc target position x: %f\n", Vessel_X_Ref4);

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead4[i] = V4_Fairlead_pos[i];

    tau_cable4[i] = V4_Flines[i] * turn[i];
  }

  AHV_Model_Obj.step();
}

void shipsportfun1_(double* V1_Surge, double* V1_Sway, double* V1_Heave,
                    double* V1_Roll, double* V1_Pitch, double* V1_Yaw,
                    double V1_Freedom[], double V1_velocity[])

{
  *V1_Surge = eta_AHV1[0];

  *V1_Sway = eta_AHV1[1];

  *V1_Heave = -eta_AHV1[2];

  *V1_Roll = eta_AHV1[3];

  *V1_Pitch = eta_AHV1[4];

  *V1_Yaw = eta_AHV1[5];

  for (int j = 0; j < 6; j++) {
    V1_Freedom[j] = eta_AHV1[j] * turn[j];

    if (j > 2)

    {
      V1_velocity[j] = nu1[j] * turn[j] * 3.14 / 180;

    }

    else

    {
      V1_velocity[j] = nu1[j] * turn[j];
    }
  }
}

void shipsportfun2_(double* V2_Surge, double* V2_Sway, double* V2_Heave,
                    double* V2_Roll, double* V2_Pitch, double* V2_Yaw,
                    double V2_Freedom[], double V2_velocity[])

{
  *V2_Surge = eta_AHV2[0];

  *V2_Sway = eta_AHV2[1];

  *V2_Heave = -eta_AHV2[2];

  *V2_Roll = eta_AHV2[3];

  *V2_Pitch = eta_AHV2[4];

  *V2_Yaw = eta_AHV2[5];

  for (int j = 0; j < 6; j++)

  {
    V2_Freedom[j] = eta_AHV2[j] * turn[j];

    if (j > 2)

    {
      V2_velocity[j] = nu2[j] * turn[j] * 3.14 / 180;

    }

    else

    {
      V2_velocity[j] = nu2[j] * turn[j];
    }
  }
}

void shipsportfun3_(double* V3_Surge, double* V3_Sway, double* V3_Heave,
                    double* V3_Roll, double* V3_Pitch, double* V3_Yaw,
                    double V3_Freedom[], double V3_velocity[])

{
  *V3_Surge = eta_AHV3[0];

  *V3_Sway = eta_AHV3[1];

  *V3_Heave = -eta_AHV3[2];

  *V3_Roll = eta_AHV3[3];

  *V3_Pitch = eta_AHV3[4];

  *V3_Yaw = eta_AHV3[5];

  for (int j = 0; j < 6; j++)

  {
    V3_Freedom[j] = eta_AHV3[j] * turn[j];

    if (j > 2)

    {
      V3_velocity[j] = nu3[j] * turn[j] * 3.14 / 180;

    }

    else

    {
      V3_velocity[j] = nu3[j] * turn[j];
    }
  }

  // printf("\t Ship3 Roll: %f\n", eta_AHV3[1]);
}

void shipsportfun4_(double* V4_Surge, double* V4_Sway, double* V4_Heave,
                    double* V4_Roll, double* V4_Pitch, double* V4_Yaw,
                    double V4_Freedom[], double V4_velocity[])

{
  *V4_Surge = eta_AHV4[0];

  *V4_Sway = eta_AHV4[1];
  *V4_Heave = -eta_AHV4[2];
  *V4_Roll = eta_AHV4[3];
  *V4_Pitch = eta_AHV4[4];
  *V4_Yaw = eta_AHV4[5];

  for (int j = 0; j < 6; j++) {
    V4_Freedom[j] = eta_AHV4[j] * turn[j];

    if (j > 2) {
      V4_velocity[j] = nu4[j] * turn[j] * 3.14 / 180;
    } else {
      V4_velocity[j] = nu4[j] * turn[j];
    }
  }
}

void shipterminatefun_()

{
  AHV_Model_Obj.terminate();
}
}