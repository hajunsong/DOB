clc; clear all; close all;

global t_current start_time end_time h Y Yp intcount r_hat q1 dq1 dr1c dr2c w1 w2 m1 m2 J1c J2c

read_data;
r_hat = 0;
Y = 0;
K = 150;

define_Y_vector;

t_current = start_time;
intcount = 1;

% q_data = load('recurdyn_result_collision.txt');

indx = 1;

while(1)
%     q1 = q_data(indx, 3);
%     dq1 = q_data(indx, 7);
    Y2qdq;
    
    kinematics;
    dynamics;
%     residual;
    
    dqddq2Yp;
    
    save_data;
    
    [Y_next, t_next] = absh3(t_current, Y, Yp, h);
    
    if intcount == 0 || intcount == 5 || intcount >= 7
        indx = indx + 1;
    end
    
    t_current = t_next;
    Y = Y_next;
    disp(t_current);
    
%     p = 0.5*m1*(dr1c'*dr1c) + 0.5*m2*(dr2c'*dr2c) + 0.5*w1'*J1c*w1 + 0.5*w2'*J2c*w2;
%     r_hat = K*(Y(1) - p);
    
    if t_current > end_time; break; end
end

plotting2;