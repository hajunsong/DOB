function save_data
    
    global t_current re NRcount flag q1 q2 q3 q4 q5 q6 roll pitch yaw
    global dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6
    global f10 f21 f32 f43 f54 f65 n10 n21 n32 n43 n54 n65 manipulability sim_flag
    global manipulability_pos manipulability_ori r2 r_hat
    
    if t_current == 0
        fp = fopen('analysis_result.txt','w+');
    else
        fp = fopen('analysis_result.txt','a+');
    end
    
    fprintf(fp, '%3.7f\t%7.7f\t%7.7f\t%7.7f\t%7.7f\t%7.7f\n',t_current, q1, r2',r_hat);
    fclose(fp);

end