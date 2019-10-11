function kinematics

    global p0 q1 q2 q3 q4 q5 q6 C01 C12 C23 C34 C45 C56 r0 s01p s12p s23p s34p s45p s56p
    global A0 A1 A2 A3 A4 A5 A6 s01 s12 s23 s34 s45 s56 r1 r2 r3 r4 r5 r6 re Ae roll pitch yaw
    global A01pp A12pp A23pp A34pp A45pp A56pp s2ep C2e
    global s6ep C6e
    
    E0 = [-p0(2:4,1),tilde(p0(2:4,1))+p0(1)*eye(3)];
    G0 = [-p0(2:4,1),-tilde(p0(2:4,1))+p0(1)*eye(3)];
    A0 = E0*G0';
%     psi = pi; theta = pi/2; phi = 0;
%     A0 = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1]...
%         *[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)]...
%         *[cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];
        
    A01pp = [cos(q1) -sin(q1) 0;sin(q1) cos(q1) 0;0 0 1];
%     A12pp = [cos(q2) -sin(q2) 0;sin(q2) cos(q2) 0;0 0 1];
%     A23pp = [cos(q3) -sin(q3) 0;sin(q3) cos(q3) 0;0 0 1];
%     A34pp = [cos(q4) -sin(q4) 0;sin(q4) cos(q4) 0;0 0 1];
%     A45pp = [cos(q5) -sin(q5) 0;sin(q5) cos(q5) 0;0 0 1];
%     A56pp = [cos(q6) -sin(q6) 0;sin(q6) cos(q6) 0;0 0 1];
    
    A1 = A0*C01*A01pp;
    A2 = A1*C12;
%     A3 = A2*C23*A23pp;
%     A4 = A3*C34*A34pp;
%     A5 = A4*C45*A45pp;
%     A6 = A5*C56*A56pp;
    
    s01 = A0*s01p;
    s12 = A1*s12p;
%     s23 = A2*s23p;
%     s34 = A3*s34p;
%     s45 = A4*s45p;
%     s56 = A5*s56p;
    
    r1 = r0 + s01;
    r2 = r1 + s12;
%     r3 = r2 + s23;
%     r4 = r3 + s34;
%     r5 = r4 + s45;
%     r6 = r5 + s56;
    
%     re = r2 + A2*s2ep;
%     Ae = A2*C2e;
%     
%     roll = atan2(Ae(3,2),Ae(3,3));
%     pitch = atan2(-Ae(3,1),sqrt(Ae(3,2)^2 + Ae(3,3)^2));
%     yaw = atan2(Ae(2,1),Ae(1,1));

end