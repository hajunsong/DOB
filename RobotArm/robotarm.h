#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>

#include "FileIO/fileio.h"
#include "Numerical/numerical.h"

using namespace std;

class RobotArm
{
public:
    RobotArm(uint numbody, uint DOF);
    ~RobotArm();
#ifdef FILEIO_H_
    void run_kinematics();
    void run_inverse_kinematics();
    void run_dynamics();
    void run_dynamics_with_DOB();
    void run_dynamics_with_DOB(double *qi, double *qi_dot, double *Ta);
#endif
    void run_kinematics(double *q, double *pose);
    void run_inverse_kinematics(double* cur_joint, double* des_pose, double* res_joint, double* res_pose);

private:
    inline void tilde(double *a, double *b) {
        *(b++) = 0;	*(b++) = -a[2];	*(b++) = a[1];
        *(b++) = a[2];	*(b++) = 0;	*(b++) = -a[0];
        *(b++) = -a[1];	*(b++) = a[0];	*(b++) = 0;
    }

    class Body
    {
    public:
        Body();
        Body(double psi, double theta, double phi, double sijp_x, double sijp_y, double sijp_z);
        ~Body();
        // body initial data
        double qi, qi_dot, mi, qi_init;
        double ri[3], ri_dot[3], wi[3], rhoip[3], sijp[3], Jip[9], Cii[9], Cij[9], Ai_Cij[9], Cij_Aijpp[9], Ai_Cii[9];
        // Orientation
        double Aijpp[9], Ai[9], Hi[3], u_vec[3];
        // Position
        double sij[3], rhoi[3], ric[3], rit[9];
        // End point
        double Ce[9], sep[3], se[3], re[3], Ae[9], roll, pitch, yaw;
        // Jacobian
        double Jvi[3], Jwi[3], re_qi[3], Ae_qi[9], r6_qi[3], A6_qi[9], Aijpp_qi[9];
        double Ae_qi_31, Ae_qi_32, Ae_qi_33, Ae_qi_21, Ae_qi_11, roll_qi, pitch_qi, yaw_qi;
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
        // DOB residual
        double y, yp, Ta, r_hat, K, p_linear, p_rotate, p, y_old, yp_old, alpha;

        static void ang2mat(double ang_z1, double ang_x, double ang_z2, double* mat, bool deg_flag = true);
    };

    double DH[6*4];

    uint num_body, dof;
    double *PH, *PH_pos, *PH_ori, *delta_q, *J, *JD;

    // system variable
    double start_time, end_time, h, t_current;
    double g;

    // EQM Mass & Force
    double M, Q_g, Q_c, Q;

    // file
    char file_name[256];
    FILE *fp;

    Body *body;
    Numerical *numeric;

    double lamda;

    double *Y, *Yp;

    void kinematics();
    void inverse_kinematics(double pos_d[3], double ori_d[3]=nullptr);
        void jacobian();
    void dynamics();
    void systemEQM();
    void residual();

    void save_data();
};

