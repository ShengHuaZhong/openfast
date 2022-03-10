
#include <stdio.h>

#include "AHV_Model.h"  // Model's header file
extern "C" {
#include "log.h"
}
static AH_Model_v1ModelClass AHV_Model_Obj;  // Instance of model class

double turn[6] = {1, 1, -1, 1, 1, 1};

double V1_Fairlead_pos[3] = {-60, 0, 0};  //���¿�λ��

double V2_Fairlead_pos[3] = {-60, 0, 0};  //���¿�λ��

double V3_Fairlead_pos[3] = {-60, 0, 0};  //���¿�λ��

double V4_Fairlead_pos[3] = {-60, 0, 0};  //���¿�λ��

// RuddeAngle1[-35 35]deg
// ThrusterPercentage1[-0.01 0.01]
// TragetSpeed1-m/s
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

void shipsetfun1_(double* MeaningWaveHeight, double* WaveDirection,
                  double* CurrentDirection, double* CurrentSpeed,
                  int* DrvingMode1, int* HeadingMode11,
                  int* V1_IfTugKeepInitPos, double* RuddeAngle1,
                  double* ThrusterPercentage1, double* TragetSpeed1,
                  double* V1_TugTargetPosX, double* V1_TugTargetPosY,
                  double* heading_angle1, double V1_Flines[])

{
  log_trace(
      "Call Function(shipsetfun1_) (MeaningWaveHeight %f) (WaveDirection %f) "
      "(CurrentDirection %f) (CurrentSpeed %f) (RuddeAngle1 %f) "
      "(Thruste %f)",
      *MeaningWaveHeight, *WaveDirection, *CurrentDirection, *CurrentSpeed,
      *RuddeAngle1, *ThrusterPercentage1);
  Drving_Mode1 = *DrvingMode1;

  heading_mode1 = *HeadingMode11;

  Rudde_angle1 = *RuddeAngle1;

  Thruster_percentage1 = *ThrusterPercentage1;

  traget_speed1 = *TragetSpeed1;

  hs = *MeaningWaveHeight;

  psi_mean = *WaveDirection * 3.14 / 180;

  Current_direction = *CurrentDirection;

  Current_speed = *CurrentSpeed;

  heading_angle_ref1 = *heading_angle1;

  Vessel_X_Ref1 = *V1_TugTargetPosX;

  Vessel_Y_Ref1 = *V1_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead1[i] = V1_Fairlead_pos[i];

    tau_cable1[i] = V1_Flines[i];
  }
}

void shipsetfun2_(int* DrvingMode2, int* HeadingMode12,
                  int* V2_IfTugKeepInitPos, double* RuddeAngle2,
                  double* ThrusterPercentage2, double* TragetSpeed2,
                  double* V2_TugTargetPosX, double* V2_TugTargetPosY,
                  double* heading_angle2, double V2_Flines[])

{
  log_trace(
      "Call funcion (shipsetfun2_) (DrvingMode %d) (HeadingMode %d) "
      "(KeepInitPos %d) (Rudde %f) (Thruster %f) (TargetSpeed %f) "
      "(TugTargetPosX %f) (TugTargetPosY %f ) (HeadAngle %f)",
      *DrvingMode2, *HeadingMode12, *V2_IfTugKeepInitPos, *RuddeAngle2,
      *ThrusterPercentage2, *TragetSpeed2, *V2_TugTargetPosX, *V2_TugTargetPosY,
      *heading_angle2);
  Drving_Mode2 = *DrvingMode2;

  heading_mode2 = *HeadingMode12;

  Rudde_angle2 = *RuddeAngle2;

  Thruster_percentage2 = *ThrusterPercentage2;

  traget_speed2 = *TragetSpeed2;

  Hold_Position2 = *V2_IfTugKeepInitPos;

  heading_angle_ref2 = *heading_angle2;

  Vessel_X_Ref2 = *V2_TugTargetPosX;

  Vessel_Y_Ref2 = *V2_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead2[i] = V2_Fairlead_pos[i];

    tau_cable2[i] = V2_Flines[i] * turn[i];
  }
}

void shipsetfun3_(int* DrvingMode3, int* HeadingMode13,
                  int* V3_IfTugKeepInitPos, double* RuddeAngle3,
                  double* ThrusterPercentage3, double* TragetSpeed3,
                  double* V3_TugTargetPosX, double* V3_TugTargetPosY,
                  double* heading_angle3, double V3_Flines[])

{
  log_trace(
      "Call funcion (shipsetfun3_) (DrvingMode %d) (HeadingMode %d) "
      "(KeepInitPos %d) (Rudde %f) (Thruster %f) (TargetSpeed %f) "
      "(TugTargetPosX %f) (TugTargetPosY %f ) (HeadAngle %f)",
      *DrvingMode3, *HeadingMode13, *V3_IfTugKeepInitPos, *RuddeAngle3,
      *ThrusterPercentage3, *TragetSpeed3, *V3_TugTargetPosX, *V3_TugTargetPosY,
      *heading_angle3);
  Drving_Mode3 = *DrvingMode3;

  heading_mode3 = *HeadingMode13;

  Rudde_angle3 = *RuddeAngle3;

  Thruster_percentage3 = *ThrusterPercentage3;

  traget_speed3 = *TragetSpeed3;

  Hold_Position3 = *V3_IfTugKeepInitPos;

  heading_angle_ref3 = *heading_angle3;

  Vessel_X_Ref3 = *V3_TugTargetPosX;

  Vessel_Y_Ref3 = *V3_TugTargetPosY;

  for (int i = 0; i < 3; i++)

  {
    ahv_fairlead3[i] = V3_Fairlead_pos[i];

    tau_cable3[i] = V3_Flines[i] * turn[i];
  }
}

void shipsetfun4_(int* DrvingMode4, int* HeadingMode14,
                  int* V4_IfTugKeepInitPos, double* RuddeAngle4,
                  double* ThrusterPercentage4, double* TragetSpeed4,
                  double* V4_TugTargetPosX, double* V4_TugTargetPosY,
                  double* heading_angle4, double V4_Flines[])

{
  log_trace(
      "Call funcion (shipsetfun4_) (DrvingMode %d) (HeadingMode %d) "
      "(KeepInitPos %d) (Rudde %f) (Thruster %f) (TargetSpeed %f) "
      "(TugTargetPosX %f) (TugTargetPosY %f) (HeadAngle %f)",
      *DrvingMode4, *HeadingMode14, *V4_IfTugKeepInitPos, *RuddeAngle4,
      *ThrusterPercentage4, *TragetSpeed4, *V4_TugTargetPosX, *V4_TugTargetPosY,
      *heading_angle4);
  Drving_Mode4 = *DrvingMode4;

  heading_mode4 = *HeadingMode14;

  Rudde_angle4 = *RuddeAngle4;

  Thruster_percentage4 = *ThrusterPercentage4;

  traget_speed4 = *TragetSpeed4;

  Hold_Position4 = *V4_IfTugKeepInitPos;

  heading_angle_ref4 = *heading_angle4;

  Vessel_X_Ref4 = *V4_TugTargetPosX;

  Vessel_Y_Ref4 = *V4_TugTargetPosY;

  for (int i = 0; i < 3; i++) {
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