#include <netinet/in.h>
#include <pthread.h>
#include <quadmath.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>

#include "cjson/cJSON.h"
#include "proj.h"
#define SHIP_PID_PACKET (1)

#define WINCH_CONTROL_ORDER_PACKET (2)

#define SEA_ENV_PACKET (3)

#define LINE_FORCE_PACKET (4)

#define LINE_LEN_PACKET (5)

#define DISPLAY_PACKET (5)

typedef struct {
  int keep_pos_;
  int kepp_head_;
  double target_x_;
  double target_y_;
  double target_head_;
} ship_pid;

typedef struct {
  double sway_;
  double surge_;
  double heave_;
  double pitch_;
  double roll_;
  double yaw_;
} six_freedom;

typedef struct {
  double fairlead1_tension_;
  double fairlead2_tension_;
  double fairlead3_tension_;
  double fairlead4_tension_;
  double fairlead9_tension_;
  double fairlead10_tension_;
  __float128 chain_aht_len_;
  __float128 chain_bw_len_;
} display_pack;

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

typedef struct {
  double base_point_x_;
  double base_point_y_;
  struct sockaddr_in addr_;
  six_freedom* ship_;
} location_monitor_pack;

static int socket_fd = 0;
static struct sockaddr_in server_addr;
static struct sockaddr_in remote_addr;
static struct sockaddr_in display_remote_addr;
static location_monitor_pack location_monitor_pack_buf[5];

static ship_pid ship1;
static ship_pid ship2;
static ship_pid ship3;
static ship_pid ship4;
static six_freedom ship1_freedom;
static six_freedom ship2_freedom;
static six_freedom ship3_freedom;
static six_freedom ship4_freedom;
static six_freedom platform_freedom;
static display_pack pack;
static winch_control_order winch1;
static winch_control_order winch2;
static sea_env env;
static PJ_CONTEXT* C;
static PJ* P;
static PJ* norm;
static PJ_COORD a, b;
void read_profile_and_init(const char* filename);

int init_socket_();

void* recv_data_(void*);

void update_ship_control_(ship_pid* ship1_fortran, ship_pid* ship2_fortran,
                          ship_pid* ship3_fortran, ship_pid* ship4_fortran);

void update_sea_env_(double* wave_height, double* wave_direction,
                     double* current_direction, double* current_speed);

void send_line_force_(int* id, double* force);

int send_data_(long double* time, double* ptfm_surge, double* ptfm_sway,
               double* ptfm_heave, double* ptfm_roll, double* ptfm_pitch,
               double* ptfm_yaw);
void send_dispaly_pack(double* Fairten1, double* Fairten2, double* Fairten3,
                       double* Fairten4, double* BWFairten9,
                       double* BWFairten10, __float128* AHTChainLength,
                       __float128* BWChainLength);
void init_ship(ship_pid*);
void send_gpgga_();

void send_gpgga_() {
  time_t times = time(NULL);
  struct tm* utcTime = gmtime(&times);
  char time_format[30];
  sprintf(time_format, "%02d%02d%02d.000", utcTime->tm_hour + 8,
          utcTime->tm_min, utcTime->tm_sec);
  const int send_buf_len = 200;
  char send_buf[send_buf_len];
  socklen_t socklen = sizeof(server_addr);
  for (int i = 0; i < 5; ++i) {
    memset(send_buf, 0, send_buf_len);
    location_monitor_pack* pack_buf_ptr = &location_monitor_pack_buf[i];
    if (location_monitor_pack_buf[i].ship_ != NULL) {
      a = proj_coord(pack_buf_ptr->base_point_x_ + pack_buf_ptr->ship_->surge_,
                     pack_buf_ptr->base_point_y_ + pack_buf_ptr->ship_->sway_,
                     0, 0);
      b = proj_trans(P, PJ_INV, a);
    }
    int send_len = snprintf(send_buf, send_buf_len,
                            "$GPGGA,%s,%08.4f,N,%9.4f,E,1,04,24.4,19.7,M,,,,"
                            "0000*1F\r\n",
                            time_format, b.lp.phi * 100, b.lp.lam * 100);
    if (send_len > 0) {
      sendto(socket_fd, send_buf, send_len, 0,
             (struct sockaddr*)&pack_buf_ptr->addr_, socklen);
    }

    send_len = snprintf(send_buf, send_buf_len, "$GPHDT,%f*1F\r\n",
                        pack_buf_ptr->ship_->roll_);
    if (send_len > 0) {
      sendto(socket_fd, send_buf, send_len, 0,
             (struct sockaddr*)&pack_buf_ptr->addr_, socklen);
    }
  }
}

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

void debug_winch_order(winch_control_order* order) {
  if (order != NULL) {
    printf("\t\twinch id      : %d\n", order->winch_id_);
    printf("\t\twinch speed   : %f\n", order->winch_speed_);
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
  ret = bind(socket_fd, (struct sockaddr*)&server_addr,
             sizeof(struct sockaddr_in));
  read_profile_and_init("profile.json");
  if (ret < 0) {
    printf("socket bind fail!\n");
    return -1;
  }

  // init proj
  C = proj_context_create();

  P = proj_create_crs_to_crs(C, "EPSG:4326", "EPSG:32650", /* or EPSG:32632 */
                             NULL);

  norm = proj_normalize_for_visualization(C, P);
  if (0 == norm) {
    fprintf(stderr, "Failed to normalize transformation object.\n");
    return 1;
  }
  proj_destroy(P);
  P = norm;

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
  platform_freedom.surge_ = *ptfm_surge;
  platform_freedom.sway_ = *ptfm_sway;
  platform_freedom.heave_ = *ptfm_heave;
  platform_freedom.roll_ = *ptfm_roll;
  platform_freedom.pitch_ = *ptfm_pitch;
  platform_freedom.yaw_ = *ptfm_yaw;
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

  // save tug data to array
  six_freedom* ship_array[4];
  ship_array[0] = &ship1_freedom;
  ship_array[1] = &ship2_freedom;
  ship_array[2] = &ship3_freedom;
  ship_array[3] = &ship4_freedom;
  if ((*tug_id) >= 1 && (*tug_id) <= 4) {
    six_freedom* update_ship_ptr = ship_array[(*tug_id) - 1];
    if (update_ship_ptr != NULL) {
      update_ship_ptr->heave_ = *tug_heave;
      update_ship_ptr->surge_ = *tug_surge;
      update_ship_ptr->sway_ = *tug_sway;
      update_ship_ptr->roll_ = roll_deg;
      update_ship_ptr->pitch_ = pitch_deg;
      update_ship_ptr->yaw_ = yaw_deg;
    }
  }
}

void send_line_force_(int* id, double* force) {
  cJSON* root = cJSON_CreateObject();
  if (NULL != root) {
    cJSON_AddItemToObject(root, "id", cJSON_CreateNumber(*id));
    cJSON_AddItemToObject(root, "x_force", cJSON_CreateNumber(force[0]));
    cJSON_AddItemToObject(root, "y_force", cJSON_CreateNumber(force[1]));
    cJSON_AddItemToObject(root, "z_force", cJSON_CreateNumber(force[2]));
    cJSON_AddItemToObject(root, "type", cJSON_CreateNumber(LINE_FORCE_PACKET));
    const char* json_str = cJSON_Print(root);
    if (json_str != NULL) {
      int buf_len = strnlen(json_str, 4096) + 1;
      socklen_t socklen = sizeof(server_addr);
      sendto(socket_fd, json_str, buf_len, 0,
             (struct sockaddr*)&display_remote_addr, socklen);
    }
  }
  cJSON_Delete(root);
}

void update_dispaly_pack_(double* Fairten1, double* Fairten2, double* Fairten3,
                          double* Fairten4, double* BWFairten9,
                          double* BWFairten10, __float128* AHTChainLength,
                          __float128* BWChainLength) {
  pack.chain_aht_len_ = *AHTChainLength;
  pack.chain_bw_len_ = *BWChainLength;
  pack.fairlead1_tension_ = *Fairten1;
  pack.fairlead2_tension_ = *Fairten2;
  pack.fairlead3_tension_ = *Fairten3;
  pack.fairlead4_tension_ = *Fairten4;
  pack.fairlead9_tension_ = *BWFairten9;
  pack.fairlead10_tension_ = *BWFairten10;
}

void send_line_length_(int* id, double* len) {
  cJSON* root = cJSON_CreateObject();
  if (NULL != root) {
    cJSON_AddItemToObject(root, "id", cJSON_CreateNumber(*id));
    cJSON_AddItemToObject(root, "len", cJSON_CreateNumber(*len));
    cJSON_AddItemToObject(root, "type", cJSON_CreateNumber(LINE_LEN_PACKET));
    const char* json_str = cJSON_Print(root);
    if (json_str != NULL) {
      int buf_len = strnlen(json_str, 4096) + 1;
      socklen_t socklen = sizeof(server_addr);
      sendto(socket_fd, json_str, buf_len, 0,
             (struct sockaddr*)&display_remote_addr, socklen);
    }
  }
  cJSON_Delete(root);
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
              winchs[winch_id->valueint]->order_ =
                  winch_order_json->valuedouble;
              winchs[winch_id->valueint]->winch_speed_ =
                  winch_speed_json->valuedouble;
              // debug_winch_order(winchs[winch_id->valueint]);
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

void update_winch_order_(int* id, double* winch_speed) {
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
// TODO

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

void send_line_node_pos_(int* id, __float128* r, int* n) {
  const int conver_buf_len = 40;
  if (id != NULL && r != NULL && n != NULL) {
    int malloc_size = sizeof(char) * conver_buf_len * 3 * (*n);
    char* send_buf = (char*)malloc(malloc_size);
    memset(send_buf, 0, malloc_size);
    char buf[conver_buf_len];
    memset(buf, 0, conver_buf_len);
    sprintf(buf, "%d", *id);
    strcat(send_buf, buf);
    strcat(send_buf, ":");
    for (int i = 0; i < (*n); ++i) {
      for (int j = 0; j < 3; ++j) {
        memset(buf, 0, conver_buf_len);
        double val = r[3 * i + j];
        sprintf(buf, "%f", val);
        strcat(send_buf, buf);
        strcat(send_buf, ",");
      }
    }
    int buf_len = strnlen(send_buf, malloc_size) + 1;
    socklen_t socklen = sizeof(remote_addr);
    sendto(socket_fd, send_buf, buf_len, 0, (struct sockaddr*)&remote_addr,
           socklen);
    free(send_buf);
  }
}

void send_display_data_(int* time) {
  six_freedom* ship_array[5];
  ship_array[0] = &ship1_freedom;
  ship_array[1] = &ship2_freedom;
  ship_array[2] = &ship3_freedom;
  ship_array[3] = &ship4_freedom;
  ship_array[4] = &platform_freedom;
  cJSON* root = cJSON_CreateObject();
  for (int i = 0; i < 5; ++i) {
    six_freedom* ship_six_freedom_ptr = ship_array[i];
    cJSON* ship_root = cJSON_CreateObject();
    if (ship_root != NULL) {
      cJSON_AddItemToObject(ship_root, "id", cJSON_CreateNumber(i));
      cJSON_AddItemToObject(ship_root, "surge",
                            cJSON_CreateNumber(ship_six_freedom_ptr->surge_));
      cJSON_AddItemToObject(ship_root, "sway",
                            cJSON_CreateNumber(ship_six_freedom_ptr->sway_));
      cJSON_AddItemToObject(ship_root, "heave",
                            cJSON_CreateNumber(ship_six_freedom_ptr->heave_));
      cJSON_AddItemToObject(ship_root, "roll",
                            cJSON_CreateNumber(ship_six_freedom_ptr->roll_));
      cJSON_AddItemToObject(ship_root, "pitch",
                            cJSON_CreateNumber(ship_six_freedom_ptr->pitch_));
      cJSON_AddItemToObject(ship_root, "yaw",
                            cJSON_CreateNumber(ship_six_freedom_ptr->yaw_));
      cJSON_AddItemToObject(ship_root, "type",
                            cJSON_CreateNumber(SHIP_PID_PACKET));
      if (root != NULL) {
        char buf[4];
        sprintf(buf, "%d", i);
        cJSON_AddItemToObject(root, buf, ship_root);
      }
    }
  }
  if (root != NULL) {
    cJSON_AddItemToObject(root, "Fairten1",
                          cJSON_CreateNumber(pack.fairlead1_tension_));
    cJSON_AddItemToObject(root, "Fairten2",
                          cJSON_CreateNumber(pack.fairlead2_tension_));
    cJSON_AddItemToObject(root, "Fairten3",
                          cJSON_CreateNumber(pack.fairlead3_tension_));
    cJSON_AddItemToObject(root, "Fairten4",
                          cJSON_CreateNumber(pack.fairlead4_tension_));
    cJSON_AddItemToObject(root, "BWFairten9",
                          cJSON_CreateNumber(pack.fairlead9_tension_));
    cJSON_AddItemToObject(root, "BWFairten10",
                          cJSON_CreateNumber(pack.fairlead10_tension_));
    cJSON_AddItemToObject(root, "AHTChainLength",
                          cJSON_CreateNumber(pack.chain_aht_len_));
    cJSON_AddItemToObject(root, "BWChainLength",
                          cJSON_CreateNumber(pack.chain_bw_len_));
    cJSON_AddItemToObject(root, "type", cJSON_CreateNumber(DISPLAY_PACKET));
    cJSON_AddItemToObject(root, "time", cJSON_CreateNumber(*time));
  }

  if (root != NULL) {
    const char* json_str = cJSON_Print(root);
    if (json_str != NULL) {
      int buf_len = strnlen(json_str, 4096) + 1;
      socklen_t socklen = sizeof(server_addr);
      int rt = sendto(socket_fd, json_str, buf_len, 0,
                      (struct sockaddr*)&display_remote_addr, socklen);
    }
  }
}
void update_sea_env_(double* wave_height, double* wave_direction,
                     double* current_direction, double* current_speed) {
  if (wave_height != NULL) *wave_height = env.wave_height_;
  if (wave_direction != NULL) *wave_direction = env.wave_direction_;
  if (current_direction != NULL) *current_direction = env.current_direction_;
  if (current_speed != NULL) *current_speed = env.current_speed_;
}
void read_profile_and_init(const char* filename) {
  FILE* profile_file = fopen(filename, "r");
  if (profile_file != NULL) {
    //获取文件长度
    fseek(profile_file, 0, SEEK_END);            //定位到文件末
    int profile_file_len = ftell(profile_file);  //文件长度
    printf("profile file len:%d\n ", profile_file_len);
    char* profile_file_buf = (void*)malloc(profile_file_len + 1);
    memset(profile_file_buf, 0, profile_file_len + 1);
    fseek(profile_file, 0, SEEK_SET);
    fread(profile_file_buf, 1, profile_file_len, profile_file);
    fclose(profile_file);

    cJSON* profile_json_root = cJSON_Parse(profile_file_buf);
    if (profile_json_root != NULL) {
      six_freedom* ship_array[5];
      ship_array[0] = &ship1_freedom;
      ship_array[1] = &ship2_freedom;
      ship_array[2] = &ship3_freedom;
      ship_array[3] = &ship4_freedom;
      ship_array[4] = &platform_freedom;
      for (int i = 0; i < 5; ++i) {
        memset(&location_monitor_pack_buf[i], 0, sizeof(location_monitor_pack));

        const int field_name_len = 20;
        char buf[field_name_len];
        memset(buf, 0, field_name_len);
        snprintf(buf, field_name_len, "location_monitor_%d", i);
        cJSON* location_monitor_json =
            cJSON_GetObjectItem(profile_json_root, buf);
        if (location_monitor_json != NULL) {
          cJSON* network_addres_json =
              cJSON_GetObjectItem(location_monitor_json, "network_addres");
          if (network_addres_json == NULL) {
            printf(
                "\tMissing field from reading profile file: "
                "(network_addres),using 127.0.0.1\n");
            location_monitor_pack_buf[i].addr_.sin_addr.s_addr =
                inet_addr("127.0.0.1");
          } else {
            location_monitor_pack_buf[i].addr_.sin_addr.s_addr =
                inet_addr(network_addres_json->valuestring);
          }
          cJSON* network_port_json =
              cJSON_GetObjectItem(location_monitor_json, "network_port");
          if (network_port_json == NULL) {
            printf(
                "\tMissing field from reading profile file: "
                "(network_port), using 9000\n");
            location_monitor_pack_buf[i].addr_.sin_port = htons(9000);
          } else {
            location_monitor_pack_buf[i].addr_.sin_port =
                htons(network_port_json->valueint);
          }

          cJSON* ship_base_x_json =
              cJSON_GetObjectItem(location_monitor_json, "ship_base_x");
          if (ship_base_x_json == NULL) {
            printf(
                "\tMissing field from reading profile file: "
                "(ship_base_x), using 0\n");
          } else {
            location_monitor_pack_buf[i].base_point_x_ =
                ship_base_x_json->valuedouble;
          }
          cJSON* ship_base_y_json =
              cJSON_GetObjectItem(location_monitor_json, "ship_base_y");
          if (ship_base_y_json == NULL) {
            printf(
                "\tMissing field from reading profile file: "
                "(ship_base_y), using 0\n");
          } else {
            location_monitor_pack_buf[i].base_point_y_ =
                ship_base_y_json->valuedouble;
          }
          location_monitor_pack_buf[i].ship_ = &ship_array[i];
        }
      }
      //
      memset(&display_remote_addr, 0, sizeof(struct sockaddr_in));
      display_remote_addr.sin_family = AF_INET;
      cJSON* data_display_network_address =
          cJSON_GetObjectItem(profile_json_root, "data_display_network");
      if (data_display_network_address == NULL) {
        printf(
            "\tMissing field from reading profile file: "
            "(data_display_network), using 127.0.0.1\n");
        display_remote_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
      } else {
        display_remote_addr.sin_addr.s_addr =
            inet_addr(data_display_network_address->valuestring);
      }
      cJSON* data_display_network_port =
          cJSON_GetObjectItem(profile_json_root, "data_display_port");
      if (data_display_network_port == NULL) {
        printf(
            "\tMissing field from reading profile file: "
            "(data_display_port), using 9000\n");
        display_remote_addr.sin_port = htons(9000);
      } else {
        display_remote_addr.sin_port =
            htons(data_display_network_port->valueint);
      }

      //
      memset(&remote_addr, 0, sizeof(struct sockaddr_in));
      remote_addr.sin_family = AF_INET;
      cJSON* remote_network_address =
          cJSON_GetObjectItem(profile_json_root, "visual_network");
      if (remote_network_address == NULL) {
        printf(
            "\tMissing field from reading profile file: "
            "(visual_network), using 127.0.0.1\n");
        remote_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
      } else {
        remote_addr.sin_addr.s_addr =
            inet_addr(remote_network_address->valuestring);
      }
      cJSON* remote_network_port =
          cJSON_GetObjectItem(profile_json_root, "visual_port");
      if (remote_network_port == NULL) {
        printf(
            "\tMissing field from reading profile file: "
            "(visual_port), using 9000\n");
        remote_addr.sin_port = htons(9000);
      } else {
        remote_addr.sin_port = htons(remote_network_port->valueint);
        printf("\t\t\t%d\n", remote_network_port->valueint);
      }
    } else {
      printf(
          "\t The content of the profile file does not conform to JSON format "
          "\n");
      exit(-1);
    }
  } else {
    printf("\tThe profile cannot exist\n");
    exit(-1);
  }
}
