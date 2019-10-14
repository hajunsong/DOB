#include "robotarm.h"

RobotArm::Body::Body(){
    u_vec[0] = 0;
    u_vec[1] = 0;
    u_vec[2] = 1;
}

RobotArm::Body::~Body(){}

void RobotArm::Body::ang2mat(double ang_z1, double ang_x, double ang_z2, double *mat, bool deg_flag)
{
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

RobotArm::RobotArm(uint numbody, uint DOF) {
    num_body = numbody;
    dof = DOF;

    PH = new double[dof];
    PH_pos = new double[3 * num_body];
    PH_ori = new double[3 * num_body];
    delta_q = new double[dof];
    J = new double[num_body * dof];
    JD = new double[dof * num_body];

    body = new Body[num_body+1];

    Y = new double[dof*2];
    Yp = new double[dof*2];

    lamda = 0.0001;

    // read data
    start_time = 0;
    h = 0.001;
    g = -9.80665;
    end_time = 2;

    // body 0 variable
    body[0].Ai[0] = 1; body[0].Ai[1] = 0; body[0].Ai[2] = 0;
    body[0].Ai[3] = 0; body[0].Ai[4] = 1; body[0].Ai[5] = 0;
    body[0].Ai[6] = 0; body[0].Ai[7] = 0; body[0].Ai[8] = 1;

    memset(body[0].ri, 0, sizeof(double)*3);
    memset(body[0].ri_dot, 0, sizeof(double)*3);

    Body::ang2mat(M_PI, M_PI_2, M_PI_2, body[0].Cij, false);
    body[0].sijp[0] = 1.37633e-6; body[0].sijp[1] = 0.0546755; body[0].sijp[2] = 0.0116372;

    memset(body[0].Yih, 0, sizeof(double)*6);
    memset(body[0].wit, 0, sizeof(double)*3);

    // body 1 variable
    body[1].qi = 0;
    body[1].qi_dot = 0;

    body[1].u_vec[0] = 0;   body[1].u_vec[1] = 0;   body[1].u_vec[2] = 1;

    body[1].mi = 1.15828469407282;

    body[1].Jip[0] =  3.33759019862649e-004; body[1].Jip[1] = -2.09252503527173e-021; body[1].Jip[2] = -2.15421743813125e-008;
    body[1].Jip[3] = -2.09252503527173e-021; body[1].Jip[4] =  1.75972275017953e-002; body[1].Jip[5] = -3.70356905869084e-020;
    body[1].Jip[6] = -2.15421743813125e-008; body[1].Jip[7] = -3.70356905869084e-020; body[1].Jip[8] =  1.72798722147092e-002;
    Body::ang2mat(M_PI_2, M_PI_2, 0, body[1].Cii, false);
    body[1].rhoip[0] = 1.13687e-16; body[1].rhoip[1] = 4.70024e-15; body[1].rhoip[2] = 0.0114834;

    Body::ang2mat(M_PI_2, M_PI_2, 0, body[1].Cij, false);
    body[1].sijp[0] = 0; body[1].sijp[1] = -0.2; body[1].sijp[2] = 0.0147;

    // DOB residual
    body[1].r_hat = 0;
    body[1].y = 0;
    body[1].yp = 0;
    body[1].K = 150;
    body[1].p_linear = 0;
    body[1].p_rotate = 0;
    body[1].p = 0;
    body[1].y_old = 0;
    body[1].yp_old = 0;

    // body 2 variable
    body[2].qi = 0;
    body[2].qi_dot = 0;

    body[2].mi = 0.489451904456537;

    body[2].u_vec[0] = 0;   body[2].u_vec[1] = 0;   body[2].u_vec[2] = 0;

    body[2].Jip[0] = 2.34784629597537e-004;  body[2].Jip[1] = 0;                     body[2].Jip[2] = 0;
    body[2].Jip[3] = 0;                      body[2].Jip[4] = 4.60887362792117e-004; body[2].Jip[5] = 2.38198731578976e-022;
    body[2].Jip[6] = 0;                      body[2].Jip[7] = 2.38198731578976e-022; body[2].Jip[8] = 2.35989837216437e-004;
    Body::ang2mat(0, 0, 0, body[2].Cii, false);
    body[2].rhoip[0] = -1.98952e-16; body[2].rhoip[1] = 0.00546345; body[2].rhoip[2] = 2.84217e-17;

    Body::ang2mat(0, 0, 0, body[2].Cij, false);
    body[2].sijp[0] = 0; body[2].sijp[1] = 0; body[2].sijp[2] = 0;

    numeric = new Numerical();
}

RobotArm::~RobotArm() {
    delete[] PH;
    delete[] PH_pos;
    delete[] PH_ori;
    delete[] delta_q;
    delete[] J;
    delete[] JD;
    delete[] Y;
    delete[] Yp;

    delete[] body;
    delete numeric;
}

#ifdef FILEIO_H_
void RobotArm::run_kinematics()
{
    sprintf(file_name, "../data/hj_kinematics_result.txt");
    fp = fopen(file_name, "w+");

    //    uint data_size = 2001;
    //    double *input = new double[data_size * col_size];
    vector<double> input;
    load_data("../data/kinematics_input_q.txt", &input, "\t");
    uint row = 8001;
    uint col = 13;

    for(uint indx = 0; indx < row; indx++) {
        for (uint i = 1; i <= 6; i++) {
            body[i].qi = input[indx*col+i];
        }

        kinematics();

        save_data();

        //        printf("Time : %.3f[s]\n", static_cast<double>(t_current));
        cout << t_current << endl;

        t_current += h;
    }

    input.clear();
    fclose(fp);
}
#endif

void RobotArm::run_kinematics(double *q, double *des_pose){
    for(int i = 0; i < 6; i++){
        body[i+1].qi = q[i];
    }

    kinematics();

    des_pose[0] = body[6].re[0];
    des_pose[1] = body[6].re[1];
    des_pose[2] = body[6].re[2];
    des_pose[3] = body[6].roll;
    des_pose[4] = body[6].pitch;
    des_pose[5] = body[6].yaw;
}

#ifdef FILEIO_H_
void RobotArm::run_inverse_kinematics() {
    sprintf(file_name, "../data/hj_inverse_kinematics_result.txt");
    fp = fopen(file_name, "w+");

    //    uint data_size = 4001;
    //    uint col_size = 7;

    //    uint sim_flag = 1;

    //    double *q_data = new double[data_size * col_size];
    //    char q_file_name[256];
    //sprintf(q_file_name, "../data/inverse_kinematics_output_q%d.txt", sim_flag);
    //    sprintf(q_file_name, "../data/inverse_kinematics_output_q_pick.txt");
    vector<double> q_data;
    load_data("../data/inverse_kinematics_output_q_pick.txt", &q_data, "\t");

    for (uint i = 1; i <= 6; i++) {
        body[i].qi = q_data[i];
    }
    //    delete[] q_data;
    q_data.clear();

    //    double *input = new double[data_size * col_size];
    //    char end_file_name[256];
    //sprintf_s(end_file_name, "../data/inverse_kinematics_input_end%d.txt", sim_flag);
    //    sprintf(end_file_name, "../data/inverse_kinematics_input_end_pick.txt");
    vector<double> input;
    load_data("../data/inverse_kinematics_input_end_pick.txt", &input, "\t");
    uint row = 8001;
    uint col = 14;

    double pos_d[3], ori_d[3];
    for (uint indx = 0; indx < row; indx++) {
        pos_d[0] = input[indx*col+1];
        pos_d[1] = input[indx*col+2];
        pos_d[2] = input[indx*col+3];
        ori_d[0] = input[indx*col+4];
        ori_d[1] = input[indx*col+5];
        ori_d[2] = input[indx*col+6];

        kinematics();

        inverse_kinematics(pos_d, ori_d);

        save_data();

        //        printf("Time : %.3f[s]\n", static_cast<double>(t_current));
        cout << t_current << endl;

        t_current += h;
    }

    input.clear();
    fclose(fp);
}
#endif

void RobotArm::run_inverse_kinematics(double* input_q, double* des_pose, double* cur_joint, double* cur_pose){
    double epsilon_pos = 0.05;
    double epsilon_ang = 1;

    for (uint i = 1; i <= num_body; i++) {
        body[i].qi = input_q[i - 1];
    }

    double pos_d[3], ori_d[3];

    for(uint i = 0; i < 5; i++){
        pos_d[0] = des_pose[0];
        pos_d[1] = des_pose[1];
        pos_d[2] = des_pose[2];
        ori_d[0] = des_pose[3];
        ori_d[1] = des_pose[4];
        ori_d[2] = des_pose[5];

        kinematics();

        inverse_kinematics(pos_d, ori_d);

        for(uint i = 1; i <= num_body; i++){
            cur_joint[i - 1] = body[i].qi;
        }

        kinematics();

        cur_pose[0] = body[num_body].re[0];
        cur_pose[1] = body[num_body].re[1];
        cur_pose[2] = body[num_body].re[2];
        cur_pose[3] = body[num_body].roll;
        cur_pose[4] = body[num_body].pitch;
        cur_pose[5] = body[num_body].yaw;

        double pos = sqrt(pow(des_pose[0] - cur_pose[0], 2) + pow(des_pose[1] - cur_pose[1], 2) + pow(des_pose[2] - cur_pose[2], 2));
        double ang_r = abs(des_pose[3] - cur_pose[3]);
        double ang_p = abs(des_pose[4] - cur_pose[4]);
        double ang_y = abs(des_pose[5] - cur_pose[5]);

        printf("[IK]pos : %f\n", pos);
        printf("[IK]ang_r : %f\t ang_p : %f\t ang_y : %f\n", ang_r, ang_p, ang_y);

        if (pos < epsilon_pos && ang_r < epsilon_ang && ang_p < epsilon_ang && ang_y < epsilon_ang){
            printf("[IK]iteration : %d\n", i);
            break;
        }
    }
}

void RobotArm::run_dynamics(){
    Y[0] = body[1].qi;
    Y[1] = body[1].qi_dot;

    t_current = start_time;
    numeric->absh3Initialize(h, dof*2);

    fp = fopen("/home/hajun/Desktop/DOB/analysis_result_c.txt", "w+");

    while(t_current <= end_time){
        body[1].qi = Y[0];
        body[1].qi_dot = Y[1];

        kinematics();
        dynamics();
        systemEQM();

        Yp[0] = body[1].qi_dot;
        Yp[1] = body[1].qi_ddot;

        save_data();

        t_current = numeric->absh3(Y, Yp, t_current);
        numeric->getY_next(Y);

        printf("%.7f\t%.7f\t%.7f\t%.7f\t%.7f\n",t_current, body[1].qi, body[2].ri[0], body[2].ri[1], body[2].ri[2]);
    }

    fclose(fp);
}

void RobotArm::run_dynamics_with_DOB(){

}

void RobotArm::run_dynamics_with_DOB(double *qi, double *qi_dot, double *Ta){
    body[1].qi = qi[0];
    body[1].qi_dot = qi_dot[0];
    body[1].Ta = Ta[0];

    fp = fopen("/home/hajun/Desktop/DOB/analysis_result_c.txt", "a+");

    kinematics();
    dynamics();
    systemEQM();
    residual();

    save_data();

    fclose(fp);
}

void RobotArm::kinematics()
{
    for (uint indx = 1; indx <= num_body; indx++) {
        // Orientation relationship
        double *Aijpp_ptr = body[indx].Aijpp;
        *(Aijpp_ptr++) = cos(body[indx].qi);	*(Aijpp_ptr++) = -sin(body[indx].qi);	*(Aijpp_ptr++) = 0;
        *(Aijpp_ptr++) = sin(body[indx].qi);	*(Aijpp_ptr++) = cos(body[indx].qi);	*(Aijpp_ptr++) = 0;
        *(Aijpp_ptr++) = 0;						*(Aijpp_ptr++) = 0;						*(Aijpp_ptr++) = 1;
        memset(body[indx].Hi, 0, sizeof(double) * 3);
        memset(body[indx].Ai, 0, sizeof(double) * 9);
        memset(body[indx].Ai_Cij, 0, sizeof(double) * 9);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                for (uint k = 0; k < 3; k++) {
                    body[indx].Ai_Cij[i * 3 + j] += body[indx - 1].Ai[i * 3 + k] * body[indx - 1].Cij[k * 3 + j];
                }
            }
        }
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Hi[i] += body[indx].Ai_Cij[i * 3 + j] * body[indx].u_vec[j];
                for (uint k = 0; k < 3; k++) {
                    body[indx].Ai[i * 3 + j] += body[indx].Ai_Cij[i * 3 + k] * body[indx].Aijpp[k * 3 + j];
                }
            }
        }
        // Position relationship
        memset(body[indx - 1].sij, 0, sizeof(double) * 3);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx - 1].sij[i] += body[indx - 1].Ai[i * 3 + j] * body[indx - 1].sijp[j];
            }
            body[indx].ri[i] = body[indx - 1].ri[i] + body[indx - 1].sij[i];
        }
    }

    // End point
    for (uint i = 0; i < 3; i++) {
        body[num_body].sij[i] = 0;
        for (uint j = 0; j < 3; j++) {
            body[num_body].sij[i] += body[num_body].Ai[i * 3 + j] * body[num_body].sijp[j];
        }
        body[num_body].re[i] = body[num_body].ri[i] + body[num_body].sij[i];
    }
    for (uint i = 0; i < 3; i++) {
        for (uint j = 0; j < 3; j++) {
            body[num_body].Ae[i * 3 + j] = 0;
            for (uint k = 0; k < 3; k++) {
                body[num_body].Ae[i * 3 + j] += body[num_body].Ai[i * 3 + k] * body[num_body].Cij[k * 3 + j];
            }
        }
    }
    body[num_body].roll = atan2(body[num_body].Ae[2 * 3 + 1], body[num_body].Ae[2 * 3 + 2]);
    body[num_body].pitch = atan2(-body[num_body].Ae[2 * 3 + 0], sqrt(pow(body[num_body].Ae[2 * 3 + 1], 2.0) + pow(body[num_body].Ae[2 * 3 + 2], 2.0)));
    body[num_body].yaw = atan2(body[num_body].Ae[1 * 3 + 0], body[num_body].Ae[0 * 3 + 0]);
}

void RobotArm::inverse_kinematics(double des_pos[3], double des_ang[3]) {
//    for (uint i = 0; i < 3; i++) {
//        PH_pos[i] = des_pos[i] - body[num_body].re[i];
//    }
//    PH_ori[0] = des_ang[0] - body[num_body].roll;
//    PH_ori[1] = des_ang[1] - body[num_body].pitch;
//    PH_ori[2] = des_ang[2] - body[num_body].yaw;

//    for (uint i = 0; i < 3; i++) {
//        PH[i] = PH_pos[i];
//        PH[i + 3] = PH_ori[i];
//    }

//    jacobian();

#if 0
    double *U, *s, *V;
    U = new double[dof * dof];
    s = new double[MIN(dof, num_body)];
    V = new double[num_body*num_body];

    numeric->svdcmp(J, dof, num_body, U, s, V);

    memset(JD, 0, sizeof(double) * num_body*dof);
    double *temp = new double[num_body*dof];
    double lamda = 1e-5;
    for (uint i = 0; i < dof; i++) {
        for (uint j = 0; j < num_body; j++) {
            for (uint k = 0; k < dof; k++) {
                temp[j * dof + k] = V[j * num_body + i] * U[k * num_body + i];
            }
        }
        for (uint j = 0; j < num_body; j++) {
            for (uint k = 0; k < dof; k++) {
                JD[j * dof + k] += (s[i] / (s[i]*s[i] +lamda*lamda))*(temp[j * dof + k]);
            }
        }
    }

    delete[] s;
    delete[] U;
    delete[] V;
    delete[] temp;


    memset(delta_q, 0, sizeof(double) * 6);
    for (uint i = 0; i < num_body; i++) {
        for (uint j = 0; j < num_body; j++) {
            delta_q[i] += JD[i * num_body + j] * PH[j];
        }
    }
#else

    int *indx = new int[6];
    double *fac = new double[6*6];
    double errmax = 0;
    int NRcount = 0;

    do{
        for (uint i = 0; i < 3; i++) {
            PH_pos[i] = des_pos[i] - body[num_body].re[i];
        }
        PH_ori[0] = des_ang[0] - body[num_body].roll;
        PH_ori[1] = des_ang[1] - body[num_body].pitch;
        PH_ori[2] = des_ang[2] - body[num_body].yaw;

        for (uint i = 0; i < 3; i++) {
            PH[i] = PH_pos[i];
            PH[i + 3] = PH_ori[i];
        }

        jacobian();

        numeric->ludcmp(J, 6, indx, 0.0, fac);
        memset(delta_q, 0, sizeof(double) * 6);
        numeric->lubksb(fac, 6, indx, PH, delta_q);

        for (uint i = 0; i < num_body; i++) {
            body[i + 1].qi += delta_q[i];
        }

        kinematics();

        for (uint i = 0; i < 3; i++) {
            PH_pos[i] = des_pos[i] - body[num_body].re[i];
        }
        PH_ori[0] = des_ang[0] - body[num_body].roll;
        PH_ori[1] = des_ang[1] - body[num_body].pitch;
        PH_ori[2] = des_ang[2] - body[num_body].yaw;

        for (uint i = 0; i < 3; i++) {
            PH[i] = PH_pos[i];
            PH[i + 3] = PH_ori[i];
        }

        errmax = PH[0];
        for(int i = 1; i < 6;i++){
            errmax = errmax > PH[i] ? errmax : PH[i];
        }
        printf("[IK]Err Max : %f\t : Iteration : %d\n", errmax, NRcount++);
    }while(errmax > 1e-3 && NRcount < 5);

    delete[] indx;
    delete[] fac;
#endif
}

void RobotArm::jacobian()
{
    double *Jv = new double[3 * num_body];
    double *Jw = new double[3 * num_body];

    //for (uint indx = 1; indx <= num_body; indx++) {
    //	for (uint i = 0; i < 3; i++) {
    //		body[indx].oi[i] = des_pos[i] - body[indx].ri[i];
    //	}
    //	tilde(body[indx].Hi, body[indx].zit);
    //	for (uint i = 0; i < 3; i++) {
    //		body[indx].Jvi[i] = 0;
    //		for (uint j = 0; j < 3; j++) {
    //			body[indx].Jvi[i] += body[indx].zit[i * 3 + j] * body[indx].oi[j];
    //		}
    //	}
    //	for (uint j = 0; j < 3; j++) {
    //		Jv[j * num_body + indx - 1] = body[indx].Jvi[j];
    //		Jw[j * num_body + indx - 1] = body[indx].Hi[j];
    //	}
    //}

    for (uint indx = 1; indx <= num_body; indx++) {
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx - 1].Cij_Aijpp[i * 3 + j] = 0;
                for (uint k = 0; k < 3; k++) {
                    body[indx - 1].Cij_Aijpp[i * 3 + j] += body[indx - 1].Cij[i * 3 + k] * body[indx].Aijpp[k * 3 + j];
                }
            }
        }
    }

    for (uint indx = 1; indx <= num_body; indx++) {
        double *Aijpp_qi_ptr = body[indx].Aijpp_qi;
        *(Aijpp_qi_ptr++) = -sin(body[indx].qi);	*(Aijpp_qi_ptr++) = -cos(body[indx].qi);	*(Aijpp_qi_ptr++) = 0;
        *(Aijpp_qi_ptr++) =  cos(body[indx].qi);	*(Aijpp_qi_ptr++) = -sin(body[indx].qi);	*(Aijpp_qi_ptr++) = 0;
        *(Aijpp_qi_ptr++) = 0;						*(Aijpp_qi_ptr++) = 0;						*(Aijpp_qi_ptr) = 0;
        memset(body[indx].A6_qi, 0, sizeof(double) * 9);
        memset(body[indx].Ae_qi, 0, sizeof(double) * 9);
        memset(body[indx].r6_qi, 0, sizeof(double) * 3);
        memset(body[indx].re_qi, 0, sizeof(double) * 3);
        double temp[9] = { 0, };
        for (uint indx2 = indx; indx2 <= num_body; indx2++) {
            if (indx2 == indx) {
                for (uint i = 0; i < 3; i++) {
                    for (uint j = 0; j < 3; j++) {
                        for (uint k = 0; k < 3; k++) {
                            body[indx].A6_qi[i * 3 + j] += body[indx2].Ai_Cij[i * 3 + k] * body[indx2].Aijpp_qi[k * 3 + j];
                        }
                    }
                }
            }
            else {
                for (uint i = 0; i < 3; i++) {
                    for (uint j = 0; j < 3; j++) {
                        temp[i * 3 + j] = 0;
                        for (uint k = 0; k < 3; k++) {
                            temp[i * 3 + j] += body[indx].A6_qi[i * 3 + k] * body[indx2 - 1].Cij_Aijpp[k * 3 + j];
                        }
                    }
                }
                memcpy(body[indx].A6_qi, temp, sizeof(double) * 9);
            }
            if (indx2 < num_body) {
                for (uint i = 0; i < 3; i++) {
                    for (uint j = 0; j < 3; j++) {
                        body[indx].r6_qi[i] += body[indx].A6_qi[i * 3 + j] * body[indx2].sijp[j];
                    }
                }
            }
        }
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                for (uint k = 0; k < 3; k++) {
                    body[indx].Ae_qi[i * 3 + j] += body[indx].A6_qi[i * 3 + k] * body[num_body].Cij[k * 3 + j];
                }
            }
        }
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].re_qi[i] += body[indx].A6_qi[i * 3 + j] * body[num_body].sijp[j];
            }
            if (indx < num_body) {
                body[indx].re_qi[i] += body[indx].r6_qi[i];
            }
        }
    }
    for (uint i = 0; i < 3; i++) {
        for (uint j = 0; j < num_body; j++) {
            Jv[i*num_body + j] = body[j + 1].re_qi[i];
        }
    }

    double Ae_31 = body[num_body].Ae[6];
    double Ae_32 = body[num_body].Ae[7];
    double Ae_33 = body[num_body].Ae[8];
    double Ae_21 = body[num_body].Ae[3];
    double Ae_11 = body[num_body].Ae[0];
    for (uint indx = 1; indx <= num_body; indx++) {
        body[indx].Ae_qi_31 = body[indx].Ae_qi[6];
        body[indx].Ae_qi_32 = body[indx].Ae_qi[7];
        body[indx].Ae_qi_33 = body[indx].Ae_qi[8];
        body[indx].Ae_qi_21 = body[indx].Ae_qi[3];
        body[indx].Ae_qi_11 = body[indx].Ae_qi[0];
    }

    double roll_q_temp1 = Ae_32 * Ae_32 + Ae_33 * Ae_33;
    double roll_q_temp2 = sqrt(roll_q_temp1);
    double roll_q_temp3 = Ae_33 + roll_q_temp2;
    double roll_q_temp4 = roll_q_temp2 * (roll_q_temp1 + Ae_33*roll_q_temp2);

    double pitch_q_temp1 = sqrt(Ae_32*Ae_32 + Ae_33*Ae_33);
    double pitch_q_temp2 = Ae_31 * Ae_31 + pitch_q_temp1 * pitch_q_temp1;
    double pitch_q_temp3 = sqrt(pitch_q_temp2);
    double pitch_q_temp4 = pitch_q_temp3 * (pitch_q_temp2 + pitch_q_temp1 * pitch_q_temp3);

    double yaw_q_temp1 = Ae_21 * Ae_21 + Ae_11 * Ae_11;
    double yaw_q_temp2 = sqrt(yaw_q_temp1);
    double yaw_q_temp3 = Ae_11 + yaw_q_temp2;
    double yaw_q_temp4 = yaw_q_temp2 * (yaw_q_temp1 + Ae_11*yaw_q_temp2);

    for (uint indx = 1; indx <= num_body; indx++) {
        body[indx].roll_qi = (roll_q_temp3*(body[indx].Ae_qi_32*Ae_33 - Ae_32*body[indx].Ae_qi_33)) / roll_q_temp4;
        body[indx].pitch_qi = -((pitch_q_temp3 + pitch_q_temp1)*(body[indx].Ae_qi_31*pitch_q_temp1 - Ae_31 * (Ae_32*body[indx].Ae_qi_32 + Ae_33 * body[indx].Ae_qi_33)/pitch_q_temp1))/ pitch_q_temp4;
        body[indx].yaw_qi = (yaw_q_temp3*(body[indx].Ae_qi_21*Ae_11 - Ae_21*body[indx].Ae_qi_11)) / yaw_q_temp4;
    }

    for (uint i = 0; i < num_body; i++) {
        Jw[0 * num_body + i] = body[i + 1].roll_qi;
        Jw[1 * num_body + i] = body[i + 1].pitch_qi;
        Jw[2 * num_body + i] = body[i + 1].yaw_qi;
    }

    memcpy(J, Jv, sizeof(double) * 3 * num_body);
    memcpy(J + 3 * num_body, Jw, sizeof(double) * 3 * num_body);

    delete[] Jv;
    delete[] Jw;
}

void RobotArm::dynamics(){
    for(uint indx = 1; indx <= num_body; indx++){
        memset(body[indx].rhoi, 0, sizeof(double) * 3);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].rhoi[i] += body[indx].Ai[i * 3 + j] * body[indx].rhoip[j];
            }
            body[indx].ric[i] = body[indx].ri[i] + body[indx].rhoi[i];
        }
        // Inertial matrix
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Ai_Cii[i*3+j] = 0;
                for (uint k = 0; k < 3; k++) {
                    body[indx].Ai_Cii[i * 3 + j] += body[indx].Ai[i * 3 + k] * body[indx].Cii[k * 3 + j];
                }
            }
        }
        memset(body[indx].Jic, 0, sizeof(double) * 9);
        double temp[9] = { 0, };
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                for (uint k = 0; k < 3; k++) {
                    temp[i * 3 + j] += body[indx].Ai_Cii[i * 3 + k] * body[indx].Jip[k * 3 + j];
                }
            }
        }
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Jic[i * 3 + j] = 0;
                for (uint k = 0; k < 3; k++) {
                    body[indx].Jic[i * 3 + j] += temp[i * 3 + k] * body[indx].Ai_Cii[j * 3 + k];
                }
            }
        }

        // Velocity State
        tilde(body[indx].ri, body[indx].rit);
        memset(body[indx].Bi, 0, sizeof(double) * 3);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Bi[i] += body[indx].rit[i * 3 + j] * body[indx].Hi[j];
            }
        }
        memcpy(body[indx].Bi + 3, body[indx].Hi, sizeof(double) * 3);

        for (uint i = 0; i < 6; i++) {
            body[indx].Yih[i] = body[indx - 1].Yih[i] + body[indx].Bi[i] * body[indx].qi_dot;
        }

        // Cartesian Velocity
        for (uint i = 0; i < 6; i++) {
            for (uint j = 0; j < 6; j++) {
                if (i < 3 && j >= 3) {
                    body[indx].Ti[i * 6 + j] = -body[indx].rit[i * 3 + (j - 3)];
                }
                else {
                    body[indx].Ti[i * 6 + j] = i == j ? 1 : 0;
                }
            }
        }
        memset(body[indx].Yib, 0, sizeof(double) * 6);
        for (uint i = 0; i < 6; i++) {
            for (uint j = 0; j < 6; j++) {
                body[indx].Yib[i] += body[indx].Ti[i * 6 + j] * body[indx].Yih[j];
            }
        }
        memcpy(body[indx].ri_dot, body[indx].Yib, sizeof(double) * 3);
        memcpy(body[indx].wi, body[indx].Yib + 3, sizeof(double) * 3);
        tilde(body[indx].wi, body[indx].wit);
        memset(body[indx].ric_dot, 0, sizeof(double) * 3);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].ric_dot[i] += body[indx].wit[i * 3 + j] * body[indx].rhoi[j];
            }
            body[indx].ric_dot[i] += body[indx].ri_dot[i];
        }
        // Mass & Force
        tilde(body[indx].ric, body[indx].rict);
        tilde(body[indx].ric_dot, body[indx].rict_dot);
        double mi_rict[9] = { 0, };
        for (uint i = 0; i < 9; i++) {
            mi_rict[i] = body[indx].mi*body[indx].rict[i];
        }
        memset(body[indx].Mih, 0, sizeof(double) * 36);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Mih[i * 6 + j] = i == j ? body[indx].mi : 0;
                body[indx].Mih[i * 6 + (j + 3)] = -mi_rict[i * 3 + j];
                body[indx].Mih[(i + 3) * 6 + j] = mi_rict[i * 3 + j];
                for (uint k = 0; k < 3; k++) {
                    body[indx].Mih[(i + 3) * 6 + (j + 3)] += -mi_rict[i * 3 + k] * body[indx].rict[k * 3 + j];
                }
                body[indx].Mih[(i + 3) * 6 + (j + 3)] += body[indx].Jic[i * 3 + j];
            }
        }
        memset(body[indx].Fic, 0, sizeof(double) * 3);
        body[indx].Fic[2] = body[indx].mi*g;

        double mi_rict_drict_wi[3] = { 0, }, wit_Jic_wi[3] = { 0, };
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                double temp1[3] = { 0, }, temp2[3] = { 0, };
                for (uint k = 0; k < 3; k++) {
                    for (uint m = 0; m < 3; m++) {
                        temp1[k] += body[indx].rict_dot[k * 3 + m] * body[indx].wi[m];
                        temp2[k] += body[indx].Jic[k * 3 + m] * body[indx].wi[m];
                    }
                }
                mi_rict_drict_wi[i] += mi_rict[i * 3 + j] * temp1[j];
                wit_Jic_wi[i] += body[indx].wit[i * 3 + j] * temp2[j];
            }
        }
        memset(body[indx].Qih_g, 0, sizeof(double) * 6);
        memset(body[indx].Qih_c, 0, sizeof(double) * 6);
        for (uint i = 0; i < 3; i++) {
            for (uint j = 0; j < 3; j++) {
                body[indx].Qih_g[i + 3] += body[indx].rict[i * 3 + j] * body[indx].Fic[j];
                body[indx].Qih_c[i] += body[indx].rict_dot[i * 3 + j] * body[indx].wi[j];
            }
            body[indx].Qih_g[i] = body[indx].Fic[i];
            body[indx].Qih_c[i] *= body[indx].mi;
            body[indx].Qih_c[i + 3] += (mi_rict_drict_wi[i] - wit_Jic_wi[i]);
        }
        // Velocity Coupling
        tilde(body[indx].ri_dot, body[indx].rit_dot);
        memset(body[indx].dHi, 0, sizeof(double) * 3);
        if (indx != 0) {
            for (uint i = 0; i < 3; i++) {
                for (uint j = 0; j < 3; j++) {
                    body[indx].dHi[i] += body[indx - 1].wit[i * 3 + j] * body[indx].Hi[j];
                }
            }
        }
        memset(body[indx].Di, 0, sizeof(double) * 3);
        for (uint i = 0; i < 3; i++) {
            double temp1 = 0, temp2 = 0;
            for (uint j = 0; j < 3; j++) {
                temp1 += body[indx].rit_dot[i * 3 + j] * body[indx].Hi[j];
                temp2 += body[indx].rit[i * 3 + j] * body[indx].dHi[j];
            }
            body[indx].Di[i] = (temp1 + temp2)*body[indx].qi_dot;
            body[indx].Di[i + 3] = body[indx].dHi[i] * body[indx].qi_dot;
        }
    }
}

void RobotArm::systemEQM(){
    memcpy(body[2].Ki, body[2].Mih, sizeof(double)*36);
    for(uint i = 0; i < 36; i++){
        body[1].Ki[i] = body[2].Ki[i] + body[1].Mih[i];
    }

    double K2_D2[6] = {0,};
    for(uint i = 0; i <6; i++){
        K2_D2[i] = 0;
        for(uint j = 0; j < 6; j++){
            K2_D2[i] += body[2].Ki[i*6+j]*body[2].Di[j];
        }
    }
    memcpy(body[2].Li_g, body[2].Qih_g, sizeof(double)*6);
    memcpy(body[2].Li_c, body[2].Qih_c, sizeof(double)*6);
    for(uint i = 0; i < 6; i++){
        body[1].Li_g[i] = body[1].Qih_g[i] + body[2].Li_g[i] - K2_D2[i];
        body[1].Li_c[i] = body[1].Qih_c[i] + body[2].Li_c[i] - K2_D2[i];
    }

    M = 0;
    double K1_B1[6] = {0,};
    for(uint i = 0; i < 6; i++){
        K1_B1[i] = 0;
        for(uint j = 0; j < 6; j++){
            K1_B1[i] += body[1].Ki[i*6+j]*body[1].Bi[j];
        }
    }
    for(uint i = 0; i < 6; i++){
        M += body[1].Bi[i]*K1_B1[i];
    }

    double K1_D1[6] = {0,};
    for(uint i = 0; i < 6; i++){
        K1_D1[i] = 0;
        for(uint j = 0; j < 6; j++){
            K1_D1[i] += body[1].Ki[i*6+j]*body[1].Di[j];
        }
    }

    Q_g = 0;
    Q_c = 0;
    for(uint i = 0; i < 6; i++){
        Q_g += body[1].Bi[i]*(body[1].Li_g[i] - K1_D1[i]);
        Q_c += body[1].Bi[i]*(body[1].Li_c[i] - K1_D1[i]);
    }

    Q = 0;
    Q = Q_g + Q_c;
    body[1].qi_ddot = Q/M;
}

void RobotArm::residual(){
    body[1].alpha = Q_g - Q_c;
    body[1].yp = body[1].Ta - body[1].alpha - body[1].r_hat;

    body[1].y = body[1].y_old + body[1].yp*h + 0.5*h*h*(body[1].yp - body[1].yp_old);

    for(uint i = 1; i <= num_body; i++){
        body[i].p_linear = 0;
        for (uint j = 0; j < 3; j++) {
            body[i].p_linear += pow(body[i].ric_dot[j], 2);
        }
        body[i].p_linear *= 0.5*body[i].mi;
        double temp[3] = { 0, };
        for (uint j = 0; j < 3; j++) {
            for (uint k = 0; k < 3; k++) {
                temp[k] += body[i].Jic[j * 3 + k] * body[i].wi[k];
            }
        }
        body[i].p_rotate = 0;
        for (uint j = 0; j < 3; j++) {
            body[i].p_rotate += body[i].wi[j] * temp[j];
        }
        body[i].p = body[i].p_linear + 0.5*body[i].p_rotate;
    }
    body[1].r_hat = body[1].K*(body[1].y - (body[1].p + body[2].p));

    body[1].y_old = body[1].y;
    body[1].yp_old = body[1].yp;

    t_current += h;
}

void RobotArm::save_data() {
//    fprintf(fp, "%.7f\t", t_current);
//    for (uint i = 1; i <= num_body; i++) {
//        fprintf(fp, "%.7f\t", body[i].qi);
//    }
//    kinematics();
//    fprintf(fp, "%.7f\t%.7f\t%.7f\t", body[num_body].re[0], body[num_body].re[1], body[num_body].re[2]);
//    fprintf(fp, "%.7f\t%.7f\t%.7f", body[num_body].roll, body[num_body].pitch, body[num_body].yaw);
//    fprintf(fp, "\n");

    fprintf(fp, "%.7f\t", t_current);
    fprintf(fp, "%.7f\t", body[1].qi);
    fprintf(fp, "%.7f\t%.7f\t%.7f\t", body[2].ri[0], body[2].ri[1], body[2].ri[2]);
    fprintf(fp, "%.7f\n", body[1].r_hat);
}

