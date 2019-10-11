function inverse_kinematics_main2

    global start_time t_current end_time h NRcount

    index = 1;
    t_current = start_time;
    NRcount = 0;
    
    IK_input = load('inverse_kinematics/inverse_kinematics_input.txt');
    
    read_data;
    
    kinematics;
    
    while t_current <= end_time
        
        inverse_kinematics2(IK_input(index,2:4)', IK_input(index,5:7)');

        save_data;
        fprintf('%s %f %d\n','inverse kinematics analysis',t_current,NRcount);

        t_current = t_current + h;
        index = index + 1;

    end

end