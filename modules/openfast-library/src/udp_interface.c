#include <netinet/in.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>

#include "cjson/cJSON.h"

#define SHIP_PID_PACKET (1)

#define WINCH_CONTROL_ORDER_PACKET (2)

#define SEA_ENV_PACKET (3)

typedef struct {
  int keep_pos_;
  int kepp_head_;
  double target_x_;
  double target_y_;
  double target_head_;
} ship_pid;

typedef struct {
  int winch_id_;
  int order_;
  double winch_speed_;
} winch_control_order;

typedef struct {
  double wave_height_;
  double wave_direction_;
  double current_direction_;
  double current_speed_;
} sea_env;

static int socket_fd = 0;
static struct sockaddr_in server_addr;
static struct sockaddr_in remote_addr;

static ship_pid ship1;
static ship_pid ship2;
static ship_pid ship3;
static ship_pid ship4;
static winch_control_order winch1;
static winch_control_order winch2;
static sea_env env;

int init_socket_();

void* recv_data_(void*);

void update_ship_control_(ship_pid* ship1_fortran, ship_pid* ship2_fortran,
                          ship_pid* ship3_fortran, ship_pid* ship4_fortran);

void update_sea_env_(double* wave_height, double* wave_direction,
                     double* current_direction, double* current_speed);

int send_data_(long double* time, double* ptfm_surge, double* ptfm_sway,
               double* ptfm_heave, double* ptfm_roll, double* ptfm_pitch,
               double* ptfm_yaw);

void init_ship(ship_pid*);

void debug_ship_pid_(ship_pid* ship) {
  if (ship != NULL) {
    printf("\t\tkeep_pos     : %d\n", ship->keep_pos_);
    printf("\t\tkepp_head    : %d\n", ship->kepp_head_);
    printf("\t\ttarget_x_    : %f\n", ship->target_x_);
    printf("\t\ttarget_y_    : %f\n", ship->target_y_);
    printf("\t\ttarget_head  : %f\n", ship->target_head_);
    printf("\n");
  }
}

int init_socket_() {
  init_ship(&ship1);
  init_ship(&ship2);
  init_ship(&ship3);
  init_ship(&ship4);
  printf("socket !\n");
  int ret = 0;
  socket_fd = socket(AF_INET, SOCK_DGRAM, 0);

  memset(&server_addr, 0, sizeof(struct sockaddr_in));
  server_addr.sin_family = AF_INET;
  server_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  server_addr.sin_port = htons(9000);

  memset(&remote_addr, 0, sizeof(struct sockaddr_in));
  remote_addr.sin_family = AF_INET;
  remote_addr.sin_addr.s_addr = inet_addr("8.140.105.185");
  remote_addr.sin_port = htons(7023);

  ret = bind(socket_fd, (struct sockaddr*)&server_addr,
             sizeof(struct sockaddr_in));
  if (ret < 0) {
    printf("socket bind fail!\n");
    return -1;
  }

  // create recv thread
  pthread_t ptid = 0;
  int error = pthread_create(&ptid, NULL, recv_data_, NULL);
  if (error != 0) {
    printf("Can't create thread\n");
  }

  pthread_detach(ptid);
  return 1;
}

int sned_data_(long double* time, double* ptfm_surge, double* ptfm_sway,
               double* ptfm_heave, double* ptfm_roll, double* ptfm_pitch,
               double* ptfm_yaw) {
  char buf[4096];
  memset(buf, 0, 4096);
  sprintf(buf, "surge:%f,sway:%f,heave:%f,roll:%f,pitch:%f,yaw:%f,",
          *ptfm_surge, *ptfm_sway, *ptfm_heave, *ptfm_roll, *ptfm_pitch,
          *ptfm_yaw);
  if (socket_fd > 0) {
    int buf_len = strnlen(buf, 4096);
    socklen_t socklen = sizeof(server_addr);
    sendto(socket_fd, buf, buf_len, 0, (struct sockaddr*)&remote_addr, socklen);
  }
}

int send_tug_data_(int* tug_id, double* tug_surge, double* tug_sway,
                   double* tug_heave, double* tug_roll, double* tug_pitch,
                   double* tug_yaw) {
  char buf[4096];
  memset(buf, 0, 4096);
  double roll_deg = (*tug_roll) * 57.324;
  double pitch_deg = *tug_pitch * 57.324;
  double yaw_deg = *tug_yaw * 57.324;
  sprintf(buf,
          "tug%d_surge:%f,tug%d_sway:%f,tug%d_heave:%f,tug%d_roll:%f"
          ",tug%d_pitch:%f,tug%d_yaw:%f,",
          *tug_id, *tug_surge, *tug_id, *tug_sway, *tug_id, *tug_heave, *tug_id,
          roll_deg, *tug_id, pitch_deg, *tug_id, yaw_deg);
  if (socket_fd > 0) {
    int buf_len = strnlen(buf, 4096) + 1;
    socklen_t socklen = sizeof(server_addr);
    sendto(socket_fd, buf, buf_len, 0, (struct sockaddr*)&remote_addr, socklen);
  }
}

void* recv_data_(void* args) {
  char buf[4096];
  struct sockaddr_in client_addres;
  socklen_t socklen = 0;
  ship_pid* ships[4];
  ships[0] = &ship1;
  ships[1] = &ship2;
  ships[2] = &ship3;
  ships[3] = &ship4;

  winch_control_order* winchs[2];
  winchs[0] = &winch1;
  winchs[1] = &winch2;

  if (socket_fd > 0)
    while (1) {
      memset(buf, 0, 4096);
      int recv_size = recvfrom(socket_fd, buf, 4096, 0,
                               (struct sockaddr*)&client_addres, &socklen);
      if (recv_size > 2) {
        cJSON* packet_json = cJSON_Parse(buf);
        if (packet_json == NULL) {
          printf("JSON parse fail\n");
          continue;
        }

        cJSON* type_json = cJSON_GetObjectItem(packet_json, "type");
        // printf("\ttype  :%d\n", type_json->valueint);
        if (type_json != NULL) {
          const int packet_type = type_json->valueint;
          if (packet_type == SHIP_PID_PACKET) {
            cJSON* ship_id = cJSON_GetObjectItem(packet_json, "ship_id");
            if (ship_id == NULL) {
              printf("Don't have field (ship_id)\n");
              continue;
            }
            cJSON* keep_pos_json = cJSON_GetObjectItem(packet_json, "keep_pos");
            if (keep_pos_json == NULL) {
              printf("Don't have field (keep_pos)\n");
              continue;
            }
            cJSON* keep_head_json =
                cJSON_GetObjectItem(packet_json, "keep_head");
            if (keep_head_json == NULL) {
              printf("Don't have field (keep_head)\n");
              continue;
            }
            cJSON* target_x_json = cJSON_GetObjectItem(packet_json, "target_x");
            if (target_x_json == NULL) {
              printf("Don't have field (target_x)\n");
              continue;
            }
            cJSON* target_y_json = cJSON_GetObjectItem(packet_json, "target_y");
            if (target_y_json == NULL) {
              printf("Don't have field (target_y)\n");
              continue;
            }
            cJSON* target_head_json =
                cJSON_GetObjectItem(packet_json, "target_head");
            if (target_head_json == NULL) {
              printf("Don't have field (target_head)\n");
              continue;
            }

            if (ship_id->valueint > 0 && ship_id->valueint < 5) {
              const int ship_index = ship_id->valueint - 1;
              ships[ship_index]->keep_pos_ = keep_pos_json->valueint;
              ships[ship_index]->kepp_head_ = keep_head_json->valueint;
              ships[ship_index]->target_x_ = target_x_json->valuedouble;
              ships[ship_index]->target_y_ = target_y_json->valuedouble;
              ships[ship_index]->target_head_ = target_head_json->valuedouble;
              printf("\t\tship_id      : %d\n", ship_id->valueint);
              debug_ship_pid_(ships[ship_index]);
            }
          } else if (packet_type == WINCH_CONTROL_ORDER_PACKET) {
            cJSON* winch_id = cJSON_GetObjectItem(packet_json, "winch_id");
            if (winch_id == NULL) {
              printf("Don't have field (winch_id)\n");
              continue;
            }
            cJSON* winch_speed_json =
                cJSON_GetObjectItem(packet_json, "winch_speed");
            if (winch_speed_json == NULL) {
              printf("Don't have field (winch_speed)\n");
              continue;
            }
            cJSON* winch_order_json =
                cJSON_GetObjectItem(packet_json, "winch_order");
            if (winch_order_json == NULL) {
              printf("Don't have field (winch_order)\n");
              continue;
            }
            if (winch_id->valueint >= 0 && winch_id->valueint < 2) {
              winchs[winch_id->valueint]->order_ = winch_order_json->valueint;
              winchs[winch_id->valueint]->winch_speed_ =
                  winch_speed_json->valuedouble;
            }
          } else if (packet_type == SEA_ENV_PACKET) {
            cJSON* wave_height_json =
                cJSON_GetObjectItem(packet_json, "wave_height");
            if (wave_height_json == NULL) {
              printf("Don't have field (wave_height)\n");
              continue;
            }
            cJSON* wave_direction_json =
                cJSON_GetObjectItem(packet_json, "wave_direction");
            if (wave_direction_json == NULL) {
              printf("Don't have field (wave_direction)\n");
              continue;
            }
            cJSON* current_direction_json =
                cJSON_GetObjectItem(packet_json, "current_direction");
            if (current_direction_json == NULL) {
              printf("Don't have field (current_direction)\n");
              continue;
            }
            cJSON* current_speed_json =
                cJSON_GetObjectItem(packet_json, "current_speed");
            if (current_speed_json == NULL) {
              printf("Don't have field (current_speed)\n");
              continue;
            }

            env.current_direction_ = current_direction_json->valuedouble;
            env.current_speed_ = current_speed_json->valuedouble;
            env.wave_direction_ = wave_direction_json->valuedouble;
            env.wave_height_ = wave_height_json->valuedouble;
          }
        } else {
          printf("No type field exists\n");
        }
        cJSON_Delete(packet_json);
      }
    }
}

void update_winch_order(int* id, double* winch_speed) {
  if (id == NULL || winch_speed == NULL) {
    return;
  }
  if (*id == 5) {
    *winch_speed = winch1.winch_speed_;
  }
  if (*id == 6) {
    *winch_speed = winch2.winch_speed_;
  }
}

void update_ship_control_(ship_pid* ship1_fortran, ship_pid* ship2_fortran,
                          ship_pid* ship3_fortran, ship_pid* ship4_fortran) {
  if (ship1_fortran != NULL) {
    ship1_fortran->keep_pos_ = ship1.keep_pos_;
    ship1_fortran->kepp_head_ = ship1.kepp_head_;
    ship1_fortran->target_x_ = ship1.target_x_;
    ship1_fortran->target_x_ = ship1.target_x_;
    ship1_fortran->target_head_ = ship1.target_head_;
  }

  if (ship2_fortran != NULL) {
    ship2_fortran->keep_pos_ = ship2.keep_pos_;
    ship2_fortran->kepp_head_ = ship2.kepp_head_;
    ship2_fortran->target_x_ = ship2.target_x_;
    ship2_fortran->target_x_ = ship2.target_x_;
    ship2_fortran->target_head_ = ship2.target_head_;
  }

  if (ship3_fortran != NULL) {
    ship3_fortran->keep_pos_ = ship3.keep_pos_;
    ship3_fortran->kepp_head_ = ship3.kepp_head_;
    ship3_fortran->target_x_ = ship3.target_x_;
    ship3_fortran->target_x_ = ship3.target_x_;
    ship3_fortran->target_head_ = ship3.target_head_;
  }

  if (ship4_fortran != NULL) {
    ship4_fortran->keep_pos_ = ship4.keep_pos_;
    ship4_fortran->kepp_head_ = ship4.kepp_head_;
    ship4_fortran->target_x_ = ship4.target_x_;
    ship4_fortran->target_x_ = ship4.target_x_;
    ship4_fortran->target_head_ = ship4.target_head_;
  }
}

void init_ship(ship_pid* ship) {
  ship->keep_pos_ = 1;
  ship->kepp_head_ = 1;
}

void update_sea_env_(double* wave_height, double* wave_direction,
                     double* current_direction, double* current_speed) {
  if (wave_height != NULL) *wave_height = env.wave_height_;
  if (wave_direction != NULL) *wave_direction = env.wave_direction_;
  if (current_direction != NULL) *current_direction = env.current_direction_;
  if (current_speed != NULL) *current_speed = env.current_speed_;
}