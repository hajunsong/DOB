function J = jacobian2

    global A1 C12 A2 C23 A3 C34 A4 C45 A5 C56 A6 C6e A0 C01
    global r1 r2 r3 r4 r5 r6 re r0 Ae
    global q1 q2 q3 q4 q5 q6 A12pp A23pp A34pp A45pp A56pp s12p s23p s34p s45p s56p s6ep 
    
%     z1 = A0*C01*[0;0;1];
%     z2 = A1*C12*[0;0;1];
%     z3 = A2*C23*[0;0;1];
%     z4 = A3*C34*[0;0;1];
%     z5 = A4*C45*[0;0;1];
%     z6 = A5*C56*[0;0;1];
%     
%     o1 = r_e - r1;
%     o2 = r_e - r2;
%     o3 = r_e - r3;
%     o4 = r_e - r4;
%     o5 = r_e - r5;
%     o6 = r_e - r6;
%     
%     z1t = tilde(z1);
%     z2t = tilde(z2);
%     z3t = tilde(z3);
%     z4t = tilde(z4);
%     z5t = tilde(z5);
%     z6t = tilde(z6);
    
%     Jv = [z1t*o1 z2t*o2 z3t*o3 z4t*o4 z5t*o5 z6t*o6];
%     Jw = [z1 z2 z3 z4 z5 z6];
    
%     J = [Jv;Jw];
    
    A01pp_q1 = [-sin(q1) -cos(q1) 0;cos(q1) -sin(q1) 0; 0 0 0];
    A12pp_q2 = [-sin(q2) -cos(q2) 0;cos(q2) -sin(q2) 0; 0 0 0];
    A23pp_q3 = [-sin(q3) -cos(q3) 0;cos(q3) -sin(q3) 0; 0 0 0];
    A34pp_q4 = [-sin(q4) -cos(q4) 0;cos(q4) -sin(q4) 0; 0 0 0];
    A45pp_q5 = [-sin(q5) -cos(q5) 0;cos(q5) -sin(q5) 0; 0 0 0];
    A56pp_q6 = [-sin(q6) -cos(q6) 0;cos(q6) -sin(q6) 0; 0 0 0];
    
    C23_A23pp = C23*A23pp; C34_A34pp = C34*A34pp; C45_A45pp = C45*A45pp; C56_A56pp = C56*A56pp;
    
    A1_q1 = A0*C01*A01pp_q1;
    A2_q1 = A1_q1*C12*A12pp; A2_q2 = A1*C12*A12pp_q2;
    A3_q1 = A2_q1*C23_A23pp; A3_q2 = A2_q2*C23_A23pp; A3_q3 = A2*C23*A23pp_q3;
    A4_q1 = A3_q1*C34_A34pp; A4_q2 = A3_q2*C34_A34pp; A4_q3 = A3_q3*C34_A34pp; A4_q4 = A3*C34*A34pp_q4;
    A5_q1 = A4_q1*C45_A45pp; A5_q2 = A4_q2*C45_A45pp; A5_q3 = A4_q3*C45_A45pp; A5_q4 = A4_q4*C45_A45pp; A5_q5 = A4*C45*A45pp_q5;
    A6_q1 = A5_q1*C56_A56pp; A6_q2 = A5_q2*C56_A56pp; A6_q3 = A5_q3*C56_A56pp; A6_q4 = A5_q4*C56_A56pp; A6_q5 = A5_q5*C56_A56pp; A6_q6 = A5*C56*A56pp_q6;
    
    s12_q1 = A1_q1*s12p;
    s23_q1 = A2_q1*s23p; s23_q2 = A2_q2*s23p;
    s34_q1 = A3_q1*s34p; s34_q2 = A3_q2*s34p; s34_q3 = A3_q3*s34p;
    s45_q1 = A4_q1*s45p; s45_q2 = A4_q2*s45p; s45_q3 = A4_q3*s45p; s45_q4 = A4_q4*s45p;
    s56_q1 = A5_q1*s56p; s56_q2 = A5_q2*s56p; s56_q3 = A5_q3*s56p; s56_q4 = A5_q4*s56p; s56_q5 = A5_q5*s56p;
    
    r6_q1 = s12_q1 + s23_q1 + s34_q1 + s45_q1 + s56_q1;
    r6_q2 = s23_q2 + s34_q2 + s45_q2 + s56_q2;
    r6_q3 = s34_q3 + s45_q3 + s56_q3;
    r6_q4 = s45_q4 + s56_q4;
    r6_q5 = s56_q5;
    
    re_q1 = r6_q1 + A6_q1*s6ep;
    re_q2 = r6_q2 + A6_q2*s6ep;
    re_q3 = r6_q3 + A6_q3*s6ep;
    re_q4 = r6_q4 + A6_q4*s6ep;
    re_q5 = r6_q5 + A6_q5*s6ep;
    re_q6 = A6_q6*s6ep;
    
    Jv = [re_q1 re_q2 re_q3 re_q4 re_q5 re_q6];
    
%     Ae = A6*C6e;
    Ae_q1 = A6_q1*C6e; Ae_q2 = A6_q2*C6e; Ae_q3 = A6_q3*C6e; Ae_q4 = A6_q4*C6e; Ae_q5 = A6_q5*C6e; Ae_q6 = A6_q6*C6e;
    Ae_31 = Ae(3,1); Ae_32 = Ae(3,2); Ae_33 = Ae(3,3); Ae_21 = Ae(2,1);	Ae_11 = Ae(1,1);
    Ae_q1_31 = Ae_q1(3,1);	Ae_q1_32 = Ae_q1(3,2);	Ae_q1_33 = Ae_q1(3,3);
	Ae_q2_31 = Ae_q2(3,1);	Ae_q2_32 = Ae_q2(3,2);	Ae_q2_33 = Ae_q2(3,3);
	Ae_q3_31 = Ae_q3(3,1);	Ae_q3_32 = Ae_q3(3,2);	Ae_q3_33 = Ae_q3(3,3);
	Ae_q4_31 = Ae_q4(3,1);	Ae_q4_32 = Ae_q4(3,2);	Ae_q4_33 = Ae_q4(3,3);
	Ae_q5_31 = Ae_q5(3,1);	Ae_q5_32 = Ae_q5(3,2);	Ae_q5_33 = Ae_q5(3,3);
	Ae_q6_31 = Ae_q6(3,1);	Ae_q6_32 = Ae_q6(3,2);	Ae_q6_33 = Ae_q6(3,3);
    Ae_q1_21 = Ae_q1(2,1);	Ae_q1_11 = Ae_q1(1,1);
	Ae_q2_21 = Ae_q2(2,1);	Ae_q2_11 = Ae_q2(1,1);
	Ae_q3_21 = Ae_q3(2,1);	Ae_q3_11 = Ae_q3(1,1);
	Ae_q4_21 = Ae_q4(2,1);	Ae_q4_11 = Ae_q4(1,1);
	Ae_q5_21 = Ae_q5(2,1);	Ae_q5_11 = Ae_q5(1,1);
	Ae_q6_21 = Ae_q6(2,1);	Ae_q6_11 = Ae_q6(1,1);
    
    % roll = atan2(Ae(3,2),Ae(3,3));
    roll_q_temp1 = Ae_32^2 + Ae_33^2;
    roll_q_temp2 = sqrt(roll_q_temp1);
    roll_q1 = ((roll_q_temp2 + Ae_33)*(Ae_q1_32*Ae_33 - Ae_32*Ae_q1_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));
    roll_q2 = ((roll_q_temp2 + Ae_33)*(Ae_q2_32*Ae_33 - Ae_32*Ae_q2_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));
    roll_q3 = ((roll_q_temp2 + Ae_33)*(Ae_q3_32*Ae_33 - Ae_32*Ae_q3_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));
    roll_q4 = ((roll_q_temp2 + Ae_33)*(Ae_q4_32*Ae_33 - Ae_32*Ae_q4_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));
    roll_q5 = ((roll_q_temp2 + Ae_33)*(Ae_q5_32*Ae_33 - Ae_32*Ae_q5_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));
    roll_q6 = ((roll_q_temp2 + Ae_33)*(Ae_q6_32*Ae_33 - Ae_32*Ae_q6_33))/(roll_q_temp2*(roll_q_temp1 + Ae_33*roll_q_temp2));

    % pitch = atan2(-Ae(3,1),sqrt(Ae(3,2)^2 + Ae(3,3)^2));
    pitch_q_temp1 = sqrt(Ae_32^2 + Ae_33^2);
    pitch_q_temp2 = Ae_31^2 + pitch_q_temp1^2;
    pitch_q_temp3 = sqrt(pitch_q_temp2);
    pitch_q1 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q1_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q1_32 + Ae_33*Ae_q1_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));
    pitch_q2 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q2_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q2_32 + Ae_33*Ae_q2_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));
    pitch_q3 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q3_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q3_32 + Ae_33*Ae_q3_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));
    pitch_q4 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q4_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q4_32 + Ae_33*Ae_q4_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));
    pitch_q5 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q5_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q5_32 + Ae_33*Ae_q5_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));
    pitch_q6 = -((pitch_q_temp3 + pitch_q_temp1)*(Ae_q6_31*pitch_q_temp1 - Ae_31*(Ae_32*Ae_q6_32 + Ae_33*Ae_q6_33)/pitch_q_temp1))/(pitch_q_temp3*(pitch_q_temp2 + pitch_q_temp1*pitch_q_temp3));

    % yaw = atan2(Ae(2,1),Ae(1,1));
    yaw_q_temp1 = Ae_21^2 + Ae_11^2;
    yaw_q_temp2 = sqrt(yaw_q_temp1);
    yaw_q1 = ((yaw_q_temp2 + Ae_11)*(Ae_q1_21*Ae_11 - Ae_21*Ae_q1_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    yaw_q2 = ((yaw_q_temp2 + Ae_11)*(Ae_q2_21*Ae_11 - Ae_21*Ae_q2_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    yaw_q3 = ((yaw_q_temp2 + Ae_11)*(Ae_q3_21*Ae_11 - Ae_21*Ae_q3_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    yaw_q4 = ((yaw_q_temp2 + Ae_11)*(Ae_q4_21*Ae_11 - Ae_21*Ae_q4_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    yaw_q5 = ((yaw_q_temp2 + Ae_11)*(Ae_q5_21*Ae_11 - Ae_21*Ae_q5_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    yaw_q6 = ((yaw_q_temp2 + Ae_11)*(Ae_q6_21*Ae_11 - Ae_21*Ae_q6_11))/(yaw_q_temp2*(yaw_q_temp1 + Ae_11*yaw_q_temp2));
    
	Jw = [roll_q1 roll_q2 roll_q3 roll_q4 roll_q5 roll_q6;
        pitch_q1 pitch_q2 pitch_q3 pitch_q4 pitch_q5 pitch_q6;
        yaw_q1 yaw_q2 yaw_q3 yaw_q4 yaw_q5 yaw_q6];
    
    J = [Jv;Jw];
    
end