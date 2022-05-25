#include <fcntl.h>
#include <math.h>
#include <netinet/in.h>
#include <pthread.h>
#include <quadmath.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/epoll.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>

#include "cjson/cJSON.h"
#include "light_matrix.h"
#include "log.h"
#include "proj.h"
#define SHIP_PID_PACKET (1)

#define WINCH_CONTROL_ORDER_PACKET (2)

#define SEA_ENV_PACKET (3)

#define LINE_FORCE_PACKET (4)

#define TUG_DRIVER_CONSOLE (6)

#define DIRECTOR_INIT_PACK (7)

#define LINE_LEN_PACKET (5)

#define DISPLAY_PACKET (5)

typedef struct {
  int driver_mode_;
  int keep_pos_;
  int kepp_head_;
  double target_x_;
  double target_y_;
  double target_head_;
  double target_velocity_;
  double rudder_;
  double thrust_;
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
static struct sockaddr_in director_addr;
static location_monitor_pack location_monitor_pack_buf[5];

static int tcp_socket_fd = 0;

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
static int epoll_fd;
static sea_env env;
static PJ_CONTEXT* C;
static PJ* P;
static PJ* norm;
static PJ_COORD a, b;
const double ARC_TO_DEG = 57.29577951308238;
const double DEG_TO_ARC = 0.0174532925199433;

void read_profile_and_init(const char* filename);

int init_socket_();

void* recv_data_(void*);

static void update_tug_driver(ship_pid*, cJSON* root);

void update_ship_control_(ship_pid* ship1_fortran, ship_pid* ship2_fortran,
                          ship_pid* ship3_fortran, ship_pid* ship4_fortran);

void update_sea_env_(double* wave_height, double* wave_direction,
                     double* current_direction, double* current_speed);

void send_line_force_(int* id, double* force);
double eulerAnglesToRotationMatrix(double roll, double pitch, double yaw);
int send_data_(long double* time, double* ptfm_surge, double* ptfm_sway,
               double* ptfm_heave, double* ptfm_roll, double* ptfm_pitch,
               double* ptfm_yaw);
void send_dispaly_pack(double* Fairten1, double* Fairten2, double* Fairten3,
                       double* Fairten4, double* BWFairten9,
                       double* BWFairten10, __float128* AHTChainLength,
                       __float128* BWChainLength);
void init_ship(ship_pid*);
void send_gpgga_();

static void parse_director_json(cJSON*);

static void cJSON_GetDouble(double*, cJSON*, const char*);

static void cJSON_GetInt(int* value, cJSON* item, const char* error);

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
      // Openfast has no way to rotate, so rotate 135 degrees here
      double x = pack_buf_ptr->ship_->surge_;
      double y = pack_buf_ptr->ship_->sway_;

      double tmp_x, tmp_y;

      tmp_x = x * (-0.7071) + y * (-0.7071);
      tmp_y = x * 0.7071 + y * (-0.7071);

      x = pack_buf_ptr->base_point_x_ + tmp_x;
      y = pack_buf_ptr->base_point_y_ + tmp_y;
      // printf("id:%d , x:%f, y:%f, surge:%f, sway:%f\n", i, x, y,
      //        pack_buf_ptr->ship_->surge_, pack_buf_ptr->ship_->sway_);
      a = proj_coord(x, y, 0, 0);
      b = proj_trans(P, PJ_INV, a);
    }

    int latitude_deg = floor(b.lp.phi);
    double latitude_minute = (b.lp.phi - latitude_deg) * 60;
    int longitude_deg = floor(b.lp.lam);
    double longitude_minute = (b.lp.lam - longitude_deg) * 60;
    int send_len =
        snprintf(send_buf, send_buf_len,
                 "$GPGGA,%s,%02d%6.4f,N,%03d%6.4f,E,1,04,24.4,19.7,M,,,,"
                 "0000*1F\r\n",
                 time_format, latitude_deg, latitude_minute, longitude_deg,
                 longitude_minute);
    if (send_len > 0) {
      sendto(socket_fd, send_buf, send_len, 0,
             (struct sockaddr*)&pack_buf_ptr->addr_, socklen);
    }
    double yaw_ang = 0;
    if (i < 4) {
      // ship
      yaw_ang = eulerAnglesToRotationMatrix(pack_buf_ptr->ship_->roll_,
                                            pack_buf_ptr->ship_->pitch_,
                                            pack_buf_ptr->ship_->yaw_) *
                ARC_TO_DEG;
    } else {
      // platfor
      yaw_ang =
          eulerAnglesToRotationMatrix(
              pack_buf_ptr->ship_->roll_ * DEG_TO_ARC,
              pack_buf_ptr->ship_->pitch_ * DEG_TO_ARC,
              (pack_buf_ptr->ship_->yaw_ * DEG_TO_ARC) + 225 * DEG_TO_ARC) *
          ARC_TO_DEG;
    }
    if (yaw_ang > 0) {
      yaw_ang = 360 - yaw_ang;
    } else {
      yaw_ang = -yaw_ang;
    }

    send_len = snprintf(send_buf, send_buf_len, "$GPHDT,%f*1F\r\n", yaw_ang);
    if (send_len > 0) {
      sendto(socket_fd, send_buf, send_len, 0,
             (struct sockaddr*)&pack_buf_ptr->addr_, socklen);
    }
  }
}

void debug_ship_pid_(ship_pid* ship) {
  if (ship != NULL) {
    log_trace("\t\tkeep_pos     : %d", ship->keep_pos_);
    log_trace("\t\tkepp_head    : %d", ship->kepp_head_);
    log_trace("\t\ttarget_x_    : %f", ship->target_x_);
    log_trace("\t\ttarget_y_    : %f", ship->target_y_);
    log_trace("\t\ttarget_head  : %f", ship->target_head_);
  }
}

void debug_winch_order(winch_control_order* order) {
  if (order != NULL) {
    log_trace("\t\twinch id      : %d", order->winch_id_);
    log_trace("\t\twinch speed   : %f", order->winch_speed_);
  }
}

static void update_ship_pid_from_udp(ship_pid*, cJSON*);

int init_socket_() {
  log_set_level(LOG_TRACE);
  init_ship(&ship1);
  init_ship(&ship2);
  init_ship(&ship3);
  init_ship(&ship4);
  printf("socket !\n");

  epoll_fd = epoll_create(1024);

  int ret = 0;

  {
    socket_fd = socket(AF_INET, SOCK_DGRAM, 0);
    memset(&server_addr, 0, sizeof(struct sockaddr_in));
    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    server_addr.sin_port = htons(9000);
    ret = bind(socket_fd, (struct sockaddr*)&server_addr,
               sizeof(struct sockaddr_in));

    if (setsockopt(socket_fd, SOL_SOCKET, SO_REUSEADDR, &(int){1},
                   sizeof(int)) < 0) {
      log_error("setsockopt(SO_REUSEADDR) failed");
    }
    struct ip_mreqn group;

    inet_pton(AF_INET, "224.0.1.0", &group.imr_multiaddr);
    group.imr_address.s_addr = htonl(INADDR_ANY);
    inet_pton(AF_INET, "10.192.102.198", &group.imr_address);
    int BROADCAST_ENABLE = 1;

    int r = setsockopt(socket_fd, IPPROTO_IP, IP_ADD_MEMBERSHIP, &group,
                       sizeof(group));
    log_trace("setsockopt return %d", r);

    struct epoll_event event;
    event.events = EPOLLIN;
    event.data.fd = socket_fd;
    int ret = epoll_ctl(epoll_fd, EPOLL_CTL_ADD, socket_fd, &event);
    if (ret != 0) {
      log_error("epoll_ctl return value:%d", ret);
    }
  }

  read_profile_and_init("profile.json");
  if (ret < 0) {
    log_fatal("socket bind fail!");
    return -1;
  }

  // Initialize the connection with TCP
  {
    tcp_socket_fd = socket(AF_INET, SOCK_STREAM, 0);
    memset(&server_addr, 0, sizeof(struct sockaddr_in));
    int flags = fcntl(tcp_socket_fd, F_GETFL, 0);
    fcntl(tcp_socket_fd, F_SETFL, flags | O_NONBLOCK);
    socklen_t sock_len = sizeof(struct sockaddr_in);
    int ret =
        connect(tcp_socket_fd, (struct sockaddr*)&director_addr, sock_len);

    struct epoll_event event;
    event.events = EPOLLIN;
    event.data.fd = tcp_socket_fd;
    ret = epoll_ctl(epoll_fd, EPOLL_CTL_ADD, tcp_socket_fd, &event);
    if (ret != 0) {
      log_error("epoll_ctl return value:%d", ret);
    }
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
      update_ship_ptr->roll_ = *tug_roll;
      update_ship_ptr->pitch_ = *tug_pitch;
      update_ship_ptr->yaw_ = *tug_yaw;
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
  const int pack_buf_len = 4096;
  char buf[pack_buf_len];
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

  const int max_event_num = 1024;
  struct epoll_event events[max_event_num];
  memset(events, 0, sizeof(struct epoll_event) * max_event_num);
  int16_t pack_len = 0;
  int tcp_flag = 0;
  if (socket_fd > 0)
    while (1) {
      int recv_size = 0;

      int fds = epoll_wait(epoll_fd, events, max_event_num, -1);
      memset(buf, 0, 4096);
      for (int i = 0; i < fds; ++i) {
        int fd = events[i].data.fd;
        if (events[i].events & EPOLLIN == EPOLLIN) {
          if (fd == tcp_socket_fd) {
            if (tcp_flag == 0) {
              int r = read(fd, &pack_len, sizeof(pack_len));
              pack_len = ntohs(pack_len);
              log_debug("TCP data received %d, packet length %u", r, pack_len);
              if (r == 0) {
                log_warn(
                    "Disconnect the TCP connection with the director station");
                struct epoll_event event;
                epoll_ctl(epoll_fd, EPOLL_CTL_DEL, tcp_socket_fd, &event);
              }
              tcp_flag = 1;
            } else if (tcp_flag == 1) {
              if (pack_len < pack_buf_len) {
                recv_size = read(fd, buf, pack_len);
                log_debug("TCP data received %d, packet: %s", recv_size, buf);
                tcp_flag = 0;
                pack_len = 0;
              }
            }

          } else if (fd == socket_fd) {
            recv_size = recvfrom(socket_fd, buf, 4096, 0,
                                 (struct sockaddr*)&client_addres, &socklen);
          }
        }
        if (recv_size > 2) {
          cJSON* packet_json = cJSON_Parse(buf);
          if (packet_json == NULL) {
            log_warn("JSON parse fail");
            continue;
          }

          cJSON* type_json = cJSON_GetObjectItem(packet_json, "type");
          if (type_json == NULL) {
            log_error("Don't have field (type)");
            continue;
          }
          log_trace("json packet type:%d", type_json->valueint);
          if (type_json != NULL) {
            const int packet_type = type_json->valueint;

            if (packet_type == SHIP_PID_PACKET) {
              cJSON* ship_id = cJSON_GetObjectItem(packet_json, "ship_id");
              if (ship_id == NULL) {
                log_warn("Don't have field (ship_id)");
                return;
              }
              const int ship_index = ship_id->valueint - 1;
              if (ship_index > 0 && ship_index < 4) {
                update_ship_pid_from_udp(ships[ship_index], packet_json);
              } else {
                log_warn("Don't have field (ship_id %d)", ship_index);
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
              parse_director_json(packet_json);
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
            } else if (packet_type == TUG_DRIVER_CONSOLE) {
              log_trace("recive TUG_DRIVER_CONSOLE");
              cJSON* console_id_item = cJSON_GetObjectItem(packet_json, "ID");
              if (console_id_item != NULL) {
                int console_id = console_id_item->valueint;
                if (console_id == 1) {
                  // 360 degree ball screen
                  log_trace("recive TUG 1");
                  update_tug_driver(&ship4, packet_json);
                } else if (console_id == 2) {
                  // 102 room
                  log_trace("recive TUG 2");
                  update_tug_driver(&ship1, packet_json);
                } else if (console_id == 3) {
                  update_tug_driver(&ship3, packet_json);
                } else if (console_id == 4) {
                  update_tug_driver(&ship2, packet_json);
                }
              }
            } else if (packet_type == DIRECTOR_INIT_PACK) {
              parse_director_json(packet_json);
            }
          } else {
            printf("No type field exists\n");
          }
          cJSON_Delete(packet_json);
        }
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
    ship1_fortran->driver_mode_ = ship1.driver_mode_;
    ship1_fortran->rudder_ = ship1.rudder_;
    ship1_fortran->target_velocity_ = ship1.target_velocity_;
    ship1_fortran->thrust_ = ship1.thrust_;
  }

  if (ship2_fortran != NULL) {
    ship2_fortran->keep_pos_ = ship2.keep_pos_;
    ship2_fortran->kepp_head_ = ship2.kepp_head_;
    ship2_fortran->target_x_ = ship2.target_x_;
    ship2_fortran->target_x_ = ship2.target_x_;
    ship2_fortran->target_head_ = ship2.target_head_;
    ship2_fortran->driver_mode_ = ship2.driver_mode_;
    ship2_fortran->rudder_ = ship2.rudder_;
    ship2_fortran->target_velocity_ = ship2.target_velocity_;
    ship2_fortran->thrust_ = ship2.thrust_;
  }

  if (ship3_fortran != NULL) {
    ship3_fortran->keep_pos_ = ship3.keep_pos_;
    ship3_fortran->kepp_head_ = ship3.kepp_head_;
    ship3_fortran->target_x_ = ship3.target_x_;
    ship3_fortran->target_x_ = ship3.target_x_;
    ship3_fortran->target_head_ = ship3.target_head_;
    ship3_fortran->driver_mode_ = ship3.driver_mode_;
    ship3_fortran->rudder_ = ship3.rudder_;
    ship3_fortran->target_velocity_ = ship3.target_velocity_;
    ship3_fortran->thrust_ = ship3.thrust_;
  }

  if (ship4_fortran != NULL) {
    ship4_fortran->keep_pos_ = ship4.keep_pos_;
    ship4_fortran->kepp_head_ = ship4.kepp_head_;
    ship4_fortran->target_x_ = ship4.target_x_;
    ship4_fortran->target_x_ = ship4.target_x_;
    ship4_fortran->target_head_ = ship4.target_head_;
    ship4_fortran->driver_mode_ = ship4.driver_mode_;
    ship4_fortran->rudder_ = ship4.rudder_;
    ship4_fortran->target_velocity_ = ship4.target_velocity_;
    ship4_fortran->thrust_ = ship4.thrust_;
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
            log_trace(
                "Missing field from reading profile file: "
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
            log_trace(
                "Missing field from reading profile file: "
                "(network_port), using 9000\n");
            location_monitor_pack_buf[i].addr_.sin_port = htons(9000);
          } else {
            location_monitor_pack_buf[i].addr_.sin_port =
                htons(network_port_json->valueint);
          }

          cJSON* ship_base_x_json =
              cJSON_GetObjectItem(location_monitor_json, "ship_base_x");
          if (ship_base_x_json == NULL) {
            log_warn(
                "Missing field from reading profile file: "
                "(ship_base_x), using 0\n");
          } else {
            location_monitor_pack_buf[i].base_point_x_ =
                ship_base_x_json->valuedouble;
          }
          cJSON* ship_base_y_json =
              cJSON_GetObjectItem(location_monitor_json, "ship_base_y");
          if (ship_base_y_json == NULL) {
            log_warn(
                "Missing field from reading profile file: "
                "(ship_base_y), using 0\n");
          } else {
            location_monitor_pack_buf[i].base_point_y_ =
                ship_base_y_json->valuedouble;
          }
          location_monitor_pack_buf[i].ship_ = ship_array[i];
        }
      }
      //
      memset(&display_remote_addr, 0, sizeof(struct sockaddr_in));
      display_remote_addr.sin_family = AF_INET;
      cJSON* data_display_network_address =
          cJSON_GetObjectItem(profile_json_root, "data_display_network");
      if (data_display_network_address == NULL) {
        log_warn(
            "Missing field from reading profile file: "
            "(data_display_network), using 127.0.0.1");
        display_remote_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
      } else {
        display_remote_addr.sin_addr.s_addr =
            inet_addr(data_display_network_address->valuestring);
      }
      cJSON* data_display_network_port =
          cJSON_GetObjectItem(profile_json_root, "data_display_port");
      if (data_display_network_port == NULL) {
        log_warn(
            "Missing field from reading profile file: "
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
        log_warn(
            "Missing field from reading profile file: "
            "(visual_network), using 127.0.0.1");
        remote_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
      } else {
        remote_addr.sin_addr.s_addr =
            inet_addr(remote_network_address->valuestring);
      }
      cJSON* remote_network_port =
          cJSON_GetObjectItem(profile_json_root, "visual_port");
      if (remote_network_port == NULL) {
        log_warn(
            "Missing field from reading profile file: "
            "(visual_port), using 9000");
        remote_addr.sin_port = htons(9000);
      } else {
        remote_addr.sin_port = htons(remote_network_port->valueint);
      }

      //
      memset(&director_addr, 0, sizeof(director_addr));
      director_addr.sin_family = AF_INET;
      cJSON* director_network_address =
          cJSON_GetObjectItem(profile_json_root, "director_network");
      if (director_network_address == NULL) {
        log_warn(
            "Missing field from reading profile file: "
            "(director_network), using 127.0.0.1");
        director_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
      } else {
        director_addr.sin_addr.s_addr =
            inet_addr(director_network_address->valuestring);
      }
      cJSON* director_port =
          cJSON_GetObjectItem(profile_json_root, "director_port");
      if (director_port == NULL) {
        log_warn(
            "Missing field from reading profile file: "
            "(director_port), using 9000");
        director_addr.sin_port = htons(9000);
      } else {
        director_addr.sin_port = htons(director_port->valueint);
      }

      cJSON* logger_json =
          cJSON_GetObjectItem(profile_json_root, "logger_level");
      if (logger_json == NULL)
        log_set_level(LOG_INFO);
      else
        log_set_level(logger_json->valueint);

    } else {
      log_error(
          "The content of the profile file does not conform to JSON format ");
      exit(-1);
    }
  } else {
    log_debug("The profile cannot exist");
    exit(-1);
  }
}

static void update_tug_driver(ship_pid* ship, cJSON* root) {
  cJSON* velocity_item = cJSON_GetObjectItem(root, "ZT_L_Y");
  cJSON* rudder_item = cJSON_GetObjectItem(root, "ZT_L_W");
  cJSON_GetDouble(&ship->thrust_, velocity_item, "Missing field (ZT_L_W)");
  cJSON_GetDouble(&ship->rudder_, rudder_item, "Missing field (ZT_L_Y)");
  ship->thrust_ = ship->thrust_ / 100;
  log_trace("update_tug_driver:(thrust %f) (ruder %f)", ship->thrust_,
            ship->rudder_, ship->target_velocity_, ship->driver_mode_);
}

static void cJSON_GetDouble(double* value, cJSON* item, const char* error) {
  if (item == NULL) {
    log_warn("cJSON_GetDouble (value null): %s", error);
    return;
  }
  if (value == NULL) {
    log_warn("cJSON_GetDouble (item null): %s", error);
    return;
  }
  if (error == NULL) {
    log_warn("cJSON_GetDouble (error null)");
    return;
  }
  *value = item->valuedouble;
}

static void cJSON_GetInt(int* value, cJSON* item, const char* error) {
  if (item == NULL) {
    log_warn("cJSON_GetInt (value null): %s", error);
    return;
  }
  if (value == NULL) {
    log_warn("cJSON_GetInt (item null): %s", error);
    return;
  }
  if (error == NULL) {
    log_warn("cJSON_GetDouble (error null)");
    return;
  }
  *value = item->valueint;
}

static void parse_director_json(cJSON* root) {
  log_trace("Call function (parse_director_json)");
  cJSON* current_direction_json = cJSON_GetObjectItem(root, "CurrentDirection");
  cJSON* current_velocity_json = cJSON_GetObjectItem(root, "CurrentSpeed");
  cJSON* wave_direction_json = cJSON_GetObjectItem(root, "WaveDirection");
  cJSON* wave_height_json = cJSON_GetObjectItem(root, "WaveHeight");

  cJSON_GetDouble(&env.current_speed_, current_velocity_json,
                  "Miss field (CurrentSpeed)");
  cJSON_GetDouble(&env.current_direction_, current_direction_json,
                  "Miss field (CurrentDirection)");
  cJSON_GetDouble(&env.wave_direction_, wave_direction_json,
                  "Miss field (WaveDirection)");
  cJSON_GetDouble(&env.wave_height_, wave_height_json,
                  "Miss field (WaveHeight)");
  log_trace("Finash call cuntion (parse_director_json)");
}

static void update_ship_pid_from_udp(ship_pid* ship, cJSON* packet_json) {
  if (ship == NULL) {
    log_error("update_ship_pid_from_udp (ship null)");
  }
  if (packet_json == NULL) {
    log_error("update_ship_pid_from_udp (packet_json null)");
  }
  cJSON* keep_pos_json = cJSON_GetObjectItem(packet_json, "keep_pos");
  cJSON* keep_head_json = cJSON_GetObjectItem(packet_json, "keep_head");
  cJSON* target_x_json = cJSON_GetObjectItem(packet_json, "target_x");
  cJSON* target_y_json = cJSON_GetObjectItem(packet_json, "target_y");
  cJSON* target_head_json = cJSON_GetObjectItem(packet_json, "target_head");
  cJSON* target_velocity = cJSON_GetObjectItem(packet_json, "target_velocity");
  cJSON* driver_mode_item = cJSON_GetObjectItem(packet_json, "driver_mode");
  cJSON_GetInt(&ship->keep_pos_, keep_pos_json, "Missing field (keep_pos)");
  cJSON_GetInt(&ship->kepp_head_, keep_pos_json, "Missing field (keep_head)");
  cJSON_GetDouble(&ship->target_x_, target_x_json, "Missing field (target_x)");
  cJSON_GetDouble(&ship->target_y_, target_y_json, "Missing field (target_y)");
  cJSON_GetDouble(&ship->target_head_, target_head_json,
                  "Missing field (target_head)");

  cJSON_GetDouble(&ship->target_velocity_, target_velocity,
                  "Missing field (target_velocity)");
  cJSON_GetInt(&ship->driver_mode_, driver_mode_item,
               "Missing field (driver_mode)");

  log_trace(
      "update_ship_pid_from_udp (keep_pos %d) (keep_head %d) (target_x %f) "
      "(target_y %f) (target_head %f) (target_velocity %f) (driver_mode %d)",
      ship->keep_pos_, ship->kepp_head_, ship->target_x_, ship->target_y_,
      ship->target_head_, ship->target_velocity_, ship->driver_mode_);
  debug_ship_pid_(ship);
}

double eulerAnglesToRotationMatrix(double roll, double pitch, double yaw) {
  log_trace(
      "Call function (eulerAnglesToRotationMatrix)(roll %f)(pitch %f)(yaw angle"
      " %f)(yaw %f)",
      roll, pitch, yaw * ARC_TO_DEG, yaw);
  Mat r_x;
  Mat r_y;
  Mat r_z;
  Mat r_x_det;
  Mat r_z_det;
  Mat r;

  float val_x[9] = {1, 0, 0, 0, cos(roll), -sin(roll), 0, sin(roll), cos(roll)};

  float val_y[9] = {cos(pitch), 0,           sin(pitch), 0,         1,
                    0,          -sin(pitch), 0,          cos(pitch)};

  float val_z[9] = {cos(yaw), -sin(yaw), 0, sin(yaw), cos(yaw), 0, 0, 0, 1};

  float val_x_det[9] = {1,
                        0,
                        0,
                        0,
                        cos(-180 * DEG_TO_ARC),
                        -sin(-180 * DEG_TO_ARC),
                        0,
                        sin(-180 * DEG_TO_ARC),
                        cos(-180 * DEG_TO_ARC)};

  float val_z_det[9] = {cos(-90 * DEG_TO_ARC),
                        -sin(-90 * DEG_TO_ARC),
                        0,
                        sin(-90 * DEG_TO_ARC),
                        cos(-90 * DEG_TO_ARC),
                        0,
                        0,
                        0,
                        1};

  MatCreate(&r_x, 3, 3);
  MatSetVal(&r_x, val_x);

  MatCreate(&r_y, 3, 3);
  MatSetVal(&r_y, val_y);

  MatCreate(&r_z, 3, 3);
  MatSetVal(&r_z, val_z);

  MatCreate(&r_x_det, 3, 3);
  MatSetVal(&r_x_det, val_x_det);

  MatCreate(&r_z_det, 3, 3);
  MatSetVal(&r_z_det, val_z_det);

  MatCreate(&r, 3, 3);
  Mat tmp_value;
  MatCreate(&tmp_value, 3, 3);

  MatMul(&r_z_det, &r_x_det, &r);
  MatMul(&r, &r_z, &tmp_value);
  MatMul(&tmp_value, &r_y, &r);
  MatMul(&r, &r_x, &tmp_value);
  MatTrans(&tmp_value, &r);
  double sy = sqrt(r.element[0][0] * r.element[0][0] +
                   r.element[1][0] * r.element[1][0]);

  int i = sy < 1e-6;
  double x = 0;
  double y = 0;
  double z = 0;

  if (!i) {
    x = atan2(r.element[2][1], r.element[2][2]);
    y = atan2(-r.element[2][0], sy);
    z = atan2(r.element[1][0], r.element[0][0]);
  } else {
    x = atan2(-r.element[1][2], r.element[1][1]);
    y = atan2(-r.element[2][0], sy);
    z = 0;
  }
  log_trace("(roll %f)(ptich %f)(yaw angle %f)", x, y, z * ARC_TO_DEG);
  return z;
}