function inverse_kinematics_main

    global start_time t_current end_time h NRcount
    global q1 q2 q3 q4 q5 q6 s6ep sim_flag

%     for sim_flag = 1 : 5
    t_current = start_time;
    NRcount = 0;
        
    read_data;
    
    IK_data = load('ear3_motion_pick.txt');
    IK_input = [IK_data(:,2), IK_data(:,9:14)];
    IK_output = [IK_data(:,2), IK_data(:,3:8)];
    
    q1 = IK_output(1,2);
    q2 = IK_output(1,3);
    q3 = IK_output(1,4);
    q4 = IK_output(1,5);
    q5 = IK_output(1,6);
    q6 = IK_output(1,7);
    
    for indx = 1 : size(IK_input,1)
        
        kinematics;
        
        inverse_kinematics(IK_input(indx,2:4)', IK_input(indx,5:7)');
        save_data;

        fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
        t_current = t_current + h;

    end
    

plotting;
%     end
    
%     for indx = 1 : 750
%         q1 = q1 + h;
%         q5 = q5 - h;
%         
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     for indx = 1 : 250        
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     for indx = 1 : 750
%         q1 = q1 - h;
%         q5 = q5 + h;
%         
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     for indx = 1 : 500
%         if q6 <= pi
%             q6 = q6 + h*10;
%         end
%         
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     
%     for indx = 1 : 750
%         q1 = q1 + h;
%         q5 = q5 - h;
%         
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     for indx = 1 : 250        
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end
%     
%     for indx = 1 : 750
%         q1 = q1 - h;
%         q5 = q5 + h;
%         
%         kinematics;
%         
%         save_data;
%         
%         fprintf('%s %f %d\n','inverse kinematics analysis',t_current, NRcount);
%         t_current = t_current + h;
%     end

end