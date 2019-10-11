function inverse_kinematics(r_e, A_e)

    global re r6 A6 C6e s6ep PH_q q1 q2 q3 q4 q5 q6 NRcount roll pitch yaw Ae manipulability
    global manipulability_pos manipulability_ori
        
    PH_pos = r_e - re;
    PH_ori = A_e - [roll;pitch;yaw];
    PH = [PH_pos;PH_ori];
    
    errmax = max(abs(PH));
    
    NRcount = 0;
        
%         PH_q = jacobian;
        PH_q = jacobian2;


%         manipulability = sqrt(det(PH_q*PH_q'));
%         manipulability_pos = sqrt(det(PH_q(1:3,:)*PH_q(1:3,:)'));
%         manipulability_ori = sqrt(det(PH_q(4:6,:)*PH_q(4:6,:)'));
%         manipulability = det(PH_q);
    
    while errmax > 1e-5 && NRcount < 10
        NRcount = NRcount + 1;
        
        PH_q = jacobian;
%         PH_q = jacobian2(r_e);


%         manipulability = sqrt(det(PH_q*PH_q'));
%         manipulability_pos = sqrt(det(PH_q(1:3,:)*PH_q(1:3,:)'));
%         manipulability_ori = sqrt(det(PH_q(4:6,:)*PH_q(4:6,:)'));
        
        delta_q = inv(PH_q)*PH;
        
        q1 = q1 + delta_q(1);
        q2 = q2 + delta_q(2);
        q3 = q3 + delta_q(3);
        q4 = q4 + delta_q(4);
        q5 = q5 + delta_q(5);
        q6 = q6 + delta_q(6);
        
        kinematics;

        PH_pos = r_e - re;
        PH_ori = A_e - [roll;pitch;yaw];
        PH = [PH_pos;PH_ori];
        
        errmax = max(abs(PH));
    end

end