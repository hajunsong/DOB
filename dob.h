#ifndef DOB_H
#define DOB_H

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
using namespace std;

typedef unsigned int uint;
inline void tilde(double *a, double *b) {
    *(b++) = 0;	*(b++) = -a[2];	*(b++) = a[1];
    *(b++) = a[2];	*(b++) = 0;	*(b++) = -a[0];
    *(b++) = -a[1];	*(b++) = a[0];	*(b++) = 0;
}

inline void ang2mat(double ang_z1, double ang_x, double ang_z2, double *mat, bool deg_flag = true){
    double z1, x, z2;
    if (deg_flag){
        z1 = ang_z1*M_PI/180.0;
        x = ang_x*M_PI/180.0;
        z2 = ang_z2*M_PI/180.0;
    }
    else{
        z1 = ang_z1;
        x = ang_x;
        z2 = ang_z2;
    }

    double Rz1[9] = {cos(z1), -sin(z1), 0, sin(z1), cos(z1), 0, 0, 0, 1};
    double Rx[9] = {1, 0, 0, 0, cos(x), -sin(x), 0, sin(x), cos(x)};
    double Rz2[9] = {cos(z2), -sin(z2), 0, sin(z2), cos(z2), 0, 0, 0, 1};
    double Rz1Rx[9] = {0,};
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                Rz1Rx[i*3+j] += Rz1[i*3+k]*Rx[k*3+j];
            }
        }
    }

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            mat[i*3+j] = 0;
            for(int k = 0; k < 3; k++){
                mat[i*3+j] += Rz1Rx[i*3+k]*Rz2[k*3+j];
            }
        }
    }
}

class Body
{
public:
    Body();
    ~Body();
    // base body information
    double A0[9], C01[9], s01p[3], J0p[9], r0[3], s01[3], A0_C01[9];
    // body initial data
    double qi, qi_dot, mi;
    double ri[3], ri_dot[3], wi[3], rhoip[3], sijp[3], Jip[9], Cii[9], Cij[9];
    // Orientation
    double Aijpp[9], Ai[9], Hi[3], u_vec[3];
    // Position
    double sij[3], rhoi[3], ric[3], rit[9];
    // Velocity State
    double Bi[6], Yih[6];
    // Cartesian velocity state
    double Ti[36], wit[9], Yib[6], ric_dot[3];
    // Mass & Force
    double Jic[9], rict[9], rict_dot[9], Mih[36], Fic[3], Tic[3], Qih[6], Qih_g[6], Qih_c[6];
    // Velocity Coupling
    double rit_dot[9], dHi[3], Di[6];
    // System EQM
    double Ki[36], Li[6], Li_g[6], Li_c[6];
    // Acceleration
    double qi_ddot;
    // Control
    double des_vel, err_vel, err_vel_accum, err_vel_prev, T_control, T_control_vel, Kp_vel, Ki_vel, Kd_vel;
    // Residual
    double Tg, Tc, alpha;
    double r_hat, K, p, Ta, yp, y, tor;
    double NThreshold, PThreshold;
};

inline Body::Body(){
    u_vec[0] = 0;
    u_vec[1] = 0;
    u_vec[2] = 1;

    qi = 0;
    qi_dot = 0;
    qi_ddot = 0;
    mi = 0;

    Tg = 0;
    Tc = 0;

    des_vel = 0;
    err_vel = 0;
    err_vel_accum = 0;
    err_vel_prev = 0;
    T_control = 0;
    T_control_vel = 0;
    Kp_vel = 0;
    Ki_vel = 0;
    Kd_vel = 0;

    r_hat = 0;
    K = 0;
    p = 0;
    Ta = 0;
    yp = 0;

    NThreshold = -999;
    PThreshold = 999;
}

inline Body::~Body(){}

class DOB
{
public:
    DOB(uint num_body, uint payload);
    ~DOB();
    void run(double *q, double *q_dot, double *tau, int *collision_flag);
    double total_time, average_time;
private:
    uint num_body, payload;
    Body* body;

    double start_time, end_time, h, g, t_current;

    double *Y, *Yp, *Y_old, *Yp_old;

    void analysis();
        void kinematics_analysis();
        void dynamic_analysis();
        void residual_analysis();

    FILE *fp;
    char file_name[256];
};

#endif // DOB_H
