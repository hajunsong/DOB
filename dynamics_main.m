function dynamics_main

    global start_time t_current end_time h Y Yp intcount

    t_current = start_time;
    intcount = 1;
    
    read_data;
    
    define_Y_vector;
    Yp = zeros(12,1);
    
    while t_current <= end_time
    
%         [t_current, Y_next] = ode23s(@analysis, [t_current, t_current+h], Y);
%         Y = Y_next(end,:)';
%         t_current = t_current(end);
%         save_data;
%         fprintf('%s %f\n','dynamics analysis',t_current);
        
        Y2qdq;
    
        kinematics;

        dynamics;

        dqddq2Yp;
        
        [Y_next, t_next] = absh3(t_current, Y, Yp, h);
        Y = Y_next;
        
        if intcount == 2 || intcount == 6 || intcount >= 8
            save_data;
            fprintf('%s %f\n','dynamics analysis',t_current);
        end
        t_current = t_next;
        
    end

end

% function Y_next = analysis(tc, Yc)
%     global t_current Y
%     global Yp
%     
%     t_current = tc;
%     Y = Yc;
%     
%     Y2qdq;
%     
%     forward_kinematics;
% 
%     dynamics_analysis;
%     
%     dqddq2Yp;
%     
%     Y_next = Yp;
%     
%     if t_current == 0
%         save_data;
%     end
% end