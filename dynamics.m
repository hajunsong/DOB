function dynamics

    global A0 A1 A2 A3 A4 A5 A6 C11 C22 C33 C44 C55 C66 J1p J2p J3p J4p J5p J6p
    global rho1p rho2p rho3p rho4p rho5p rho6p r1 r2 r3 r4 r5 r6 g
    global C01 C12 C23 C34 C45 C56 w0 w1 w2 w3 w4 w5 w6 dq1 dq2 dq3 dq4 dq5 dq6
    global dr0 dr1 dr2 dr3 dr4 dr5 dr6 m1 m2 m3 m4 m5 m6 s01 s12 s23 s34 s45 s56
    global Q_g Q_c Ta dr1c dr2c J1c J2c
    
    global ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 r1t r2t r3t r4t r5t r6t
    global L1 L2 L3 L4 L5 L6 K1 K2 K3 K4 K5 K6 B1 B2 B3 B4 B5 B6 D1 D2 D3 D4 D5 D6

    J1c = A1*C11*J1p*(A1*C11)';
    J2c = A2*C22*J2p*(A2*C22)';
%     J3c = A3*C33*J3p*(A3*C33)';
%     J4c = A4*C44*J4p*(A4*C44)';
%     J5c = A5*C55*J5p*(A5*C55)';
%     J6c = A6*C66*J6p*(A6*C66)';
    
    rho1 = A1*rho1p;
    rho2 = A2*rho2p;
%     rho3 = A3*rho3p;
%     rho4 = A4*rho4p;
%     rho5 = A5*rho5p;
%     rho6 = A6*rho6p;
    
    r1c = r1 + rho1;
    r2c = r2 + rho2;
%     r3c = r3 + rho3;
%     r4c = r4 + rho4;
%     r5c = r5 + rho5;
%     r6c = r6 + rho6;
    
    H1 = A0*C01*[0;0;1];
    H2 = zeros(3,1);%A1*C12*[0;0;1];
%     H3 = A2*C23*[0;0;1];
%     H4 = A3*C34*[0;0;1];
%     H5 = A4*C45*[0;0;1];
%     H6 = A5*C56*[0;0;1];
    
    w1 = w0 + H1*dq1;
    w2 = w1;
%     w3 = w2 + H3*dq3;
%     w4 = w3 + H4*dq4;
%     w5 = w4 + H5*dq5;
%     w6 = w5 + H6*dq6;
    
    w0t = tilde(w0);
    w1t = tilde(w1);
    w2t = tilde(w2);
%     w3t = tilde(w3);
%     w4t = tilde(w4);
%     w5t = tilde(w5);
%     w6t = tilde(w6);
    
    dr1 = dr0 + w0t*s01;
    dr2 = dr1 + w1t*s12;
%     dr3 = dr2 + w2t*s23;
%     dr4 = dr3 + w3t*s34;
%     dr5 = dr4 + w4t*s45;
%     dr6 = dr5 + w5t*s56;
    
    r1t = tilde(r1);
    r2t = tilde(r2);
%     r3t = tilde(r3);
%     r4t = tilde(r4);
%     r5t = tilde(r5);
%     r6t = tilde(r6);
    
    B1 = [r1t*H1;H1];
    B2 = [r2t*H2;H2];
%     B3 = [r3t*H3;H3];
%     B4 = [r4t*H4;H4];
%     B5 = [r5t*H5;H5];
%     B6 = [r6t*H6;H6];
    
    dr1t = tilde(dr1);
    dr2t = tilde(dr2);
%     dr3t = tilde(dr3);
%     dr4t = tilde(dr4);
%     dr5t = tilde(dr5);
%     dr6t = tilde(dr6);
    
    dr1c = dr1 + w1t*rho1;
    dr2c = dr2 + w2t*rho2;
%     dr3c = dr3 + w3t*rho3;
%     dr4c = dr4 + w4t*rho4;
%     dr5c = dr5 + w5t*rho5;
%     dr6c = dr6 + w6t*rho6;

	dr1ct = tilde(dr1c);
    dr2ct = tilde(dr2c);
%     dr3ct = tilde(dr3c);
%     dr4ct = tilde(dr4c);
%     dr5ct = tilde(dr5c);
%     dr6ct = tilde(dr6c);
    
    dH1 = w0t*H1;
    dH2 = w1t*H2;
%     dH3 = w2t*H3;
%     dH4 = w3t*H4;
%     dH5 = w4t*H5;
%     dH6 = w5t*H6;
    
    D1 = [dr1t*H1 + r1t*dH1;dH1]*dq1;
    D2 = [dr2t*H2 + r2t*dH2;dH2]*dq2;
%     D3 = [dr3t*H3 + r3t*dH3;dH3]*dq3;
%     D4 = [dr4t*H4 + r4t*dH4;dH4]*dq4;
%     D5 = [dr5t*H5 + r5t*dH5;dH5]*dq5;
%     D6 = [dr6t*H6 + r6t*dH6;dH6]*dq6;
    
    r1ct = tilde(r1c);
    r2ct = tilde(r2c);
%     r3ct = tilde(r3c);
%     r4ct = tilde(r4c);
%     r5ct = tilde(r5c);
%     r6ct = tilde(r6c);
	
	M1h = [m1*eye(3) -m1*r1ct;m1*r1ct J1c - m1*r1ct*r1ct];
	M2h = [m2*eye(3) -m2*r2ct;m2*r2ct J2c - m2*r2ct*r2ct];
% 	M3h = [m3*eye(3) -m3*r3ct;m3*r3ct J3c - m3*r3ct*r3ct];
% 	M4h = [m4*eye(3) -m4*r4ct;m4*r4ct J4c - m4*r4ct*r4ct];
% 	M5h = [m5*eye(3) -m5*r5ct;m5*r5ct J5c - m5*r5ct*r5ct];
% 	M6h = [m6*eye(3) -m6*r6ct;m6*r6ct J6c - m6*r6ct*r6ct];
	
    F1c = [0;0;m1*g];
    F2c = [0;0;m2*g];
%     F3c = [0;m3*g;0];
%     F4c = [0;m4*g;0];
%     F5c = [0;m5*g;0];
%     F6c = [0;m6*g;0];

    T1c = [0;1;0];
    T2c = [0;0;0];
%     T3c = [0;0;0];
%     T4c = [0;0;0];
%     T5c = [0;0;0];
%     T6c = [0;0;0];

% 	Q1h = [F1c + m1*dr1ct*w1;T1c + r1ct*F1c + m1*r1ct*dr1ct*w1 - w1t*J1c*w1];
% 	Q2h = [F2c + m2*dr2ct*w2;T2c + r2ct*F2c + m2*r2ct*dr2ct*w2 - w2t*J2c*w2];
% 	Q3h = [F3c + m3*dr3ct*w3;T3c + r3ct*F3c + m3*r3ct*dr3ct*w3 - w3t*J3c*w3];
% 	Q4h = [F4c + m4*dr4ct*w4;T4c + r4ct*F4c + m4*r4ct*dr4ct*w4 - w4t*J4c*w4];
% 	Q5h = [F5c + m5*dr5ct*w5;T5c + r5ct*F5c + m5*r5ct*dr5ct*w5 - w5t*J5c*w5];
% 	Q6h = [F6c + m6*dr6ct*w6;T6c + r6ct*F6c + m6*r6ct*dr6ct*w6 - w6t*J6c*w6];

    Q1h_g = [F1c;r1ct*F1c];
    Q1h_c = [F2c;r2ct*F2c];
    Q2h_g = [m1*dr1ct*w1;m1*r1ct*dr1ct*w1 - w1t*J1c*w1];
    Q2h_c = [m2*dr2ct*w2;m2*r2ct*dr2ct*w2 - w2t*J2c*w2];
    
%     K6 = M6h;
%     K5 = K6 + M5h;
%     K4 = K5 + M4h;
%     K3 = K4 + M3h;
    K2 = M2h;
    K1 = K2 + M1h;

% 	L6 = Q6h;
% 	L5 = Q5h + L6 - K6*D6;
% 	L4 = Q4h + L5 - K5*D5;
% 	L3 = Q3h + L4 - K4*D4;
% 	L2 = Q2h;
% 	L1 = Q1h + L2 - K2*D2;

    L2_g = Q2h_g;
    L1_g = Q1h_g + L2_g - K2*D2;
    L2_c = Q2h_c;
    L1_c = Q1h_c + L2_c - K2*D2;
    
%     M = [B1'*K1*B1	B1'*K2*B2	B1'*K3*B3	B1'*K4*B4	B1'*K5*B5	B1'*K6*B6;
%          B2'*K2*B1	B2'*K2*B2	B2'*K3*B3	B2'*K4*B4	B2'*K5*B5	B2'*K6*B6;
%          B3'*K3*B1	B3'*K3*B2	B3'*K3*B3	B3'*K4*B4	B3'*K5*B5	B3'*K6*B6;
%          B4'*K4*B1	B4'*K4*B2	B4'*K4*B3	B4'*K4*B4	B4'*K5*B5	B4'*K6*B6;
%          B5'*K5*B1	B5'*K5*B2	B5'*K5*B3	B5'*K5*B4	B5'*K5*B5	B5'*K6*B6;
%          B6'*K6*B1	B6'*K6*B2	B6'*K6*B3	B6'*K6*B4	B6'*K6*B5	B6'*K6*B6];
%     Q = [B1'*(L1 - K1*(D1));
% 		 B2'*(L2 - K2*(D1 + D2));
% 		 B3'*(L3 - K3*(D1 + D2 + D3));
% 		 B4'*(L4 - K4*(D1 + D2 + D3 + D4));
% 		 B5'*(L5 - K5*(D1 + D2 + D3 + D4 + D5));
% 		 B6'*(L6 - K6*(D1 + D2 + D3 + D4 + D5 + D6))];

    M = B1'*K1*B1;
    
%     Q = B1'*(L1 - K1*(D1));
    
    Q_g = B1'*(L1_g - K1*D1);
    Q_c = B1'*(L1_c - K1*D1);
    Ta = 2;
    Q = Ta + Q_g + Q_c;
    
    q_ddot = M\Q;
    
    ddq1 = q_ddot(1,1);
%     ddq2 = q_ddot(2,1);
%     ddq3 = q_ddot(3,1);
%     ddq4 = q_ddot(4,1);
%     ddq5 = q_ddot(5,1);
%     ddq6 = q_ddot(6,1);

end