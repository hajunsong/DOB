function inverse_kinematics2(r_e, A_e)

    global re roll pitch yaw PH PH_q q1 q2 q3 q4 q5 q6 manipulability manipulability_pos manipulability_ori
    global q1_limit q2_limit q3_limit q4_limit q5_limit q6_limit

    PH_pos = r_e - re;
    PH_ori = A_e - [roll;pitch;yaw];
    PH = [PH_pos;PH_ori];
    
    PH_q = jacobian2;
    
    [U, S, V] = svd(PH_q);
    
    lamda = 1e-5;
    JD = zeros(6,6);
    for i = 1 : 6
        JD = JD + (S(i,i)/(S(i,i)^2 + lamda^2))*V(:,i)*U(:,i)';
    end
    
    delta_q = JD*PH;
    
    q1 = q1 + delta_q(1);
    q2 = q2 + delta_q(2);
    q3 = q3 + delta_q(3);
    q4 = q4 + delta_q(4);
    q5 = q5 + delta_q(5);
    q6 = q6 + delta_q(6);
    
%     if  q1*180/pi < q1_limit(1) || q1*180/pi > q1_limit(2)
%         q1 = q1 - delta_q(1);
%     end
%     if  q2*180/pi < q2_limit(1) || q2*180/pi > q2_limit(2)
%         q2 = q2 - delta_q(2);
%     end
%     if  q3*180/pi < q3_limit(1) || q3*180/pi > q3_limit(2)
%         q3 = q3 - delta_q(3);
%     end
%     if  q4*180/pi < q4_limit(1) || q4*180/pi > q4_limit(2)
%         q4 = q4 - delta_q(4);
%     end
%     if  q5*180/pi < q5_limit(1) || q5*180/pi > q5_limit(2)
%         q5 = q5 - delta_q(5);
%     end
%     if  q6*180/pi < q6_limit(1) || q6*180/pi > q6_limit(2)
%         q6 = q6 - delta_q(6);
%     end
    
    manipulability = sqrt(det(PH_q*PH_q'));
    manipulability_pos = sqrt(det(PH_q(1:3,:)*PH_q(1:3,:)'));
    manipulability_ori = sqrt(det(PH_q(4:6,:)*PH_q(4:6,:)'));
    manipulability = det(PH_q);

end