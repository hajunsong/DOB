function kinematics_main

    global end_time t_current start_time NRcount h 
    global q1 q2 q3 q4 q5 q6 q1_init q2_init q3_init q4_init q5_init q6_init
    global r6 A6 C6e s6ep re


    t_current = start_time;
    NRcount = 0;
    
    read_data;
    
    data = load('ear3.txt');
    
    for indx = 1 : size(data,1)
        q1 = data(indx, 3);
        q2 = data(indx, 4);
        q3 = data(indx, 5);
        q4 = data(indx, 6);
        q5 = data(indx, 7);
        q6 = data(indx, 8);

        kinematics;

        save_data;
        fprintf('%s %f\n','kinematics analysis',t_current);

        t_current = t_current + h;

    end

end