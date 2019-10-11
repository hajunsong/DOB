function inverse_dynamics_main

    global start_time t_current end_time h Y Yp intcount

    t_current = start_time;
    intcount = 1;
    
    read_data;
    
    define_Y_vector;
    Yp = zeros(12,1);
    
    while t_current <= end_time
        
        Y2qdq;
    
        kinematics;

        dynamics;
        
        inverse_dynamics;

        dqddq2Yp;
        
        [Y_next, t_next] = absh3(t_current, Y, Yp, h);
        Y = Y_next;
        
        if intcount == 2 || intcount == 6 || intcount >= 8
            save_data;
            fprintf('%s %f\n','inverse dynamics analysis',t_current);
        end
        t_current = t_next;
        
    end

end