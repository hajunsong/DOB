function read_body2

    global m2 q2 dq2 J2p s12p C12 rho2p C22 q2_init s23p C23 q2_limit

    q2 = 0;
    C23 = [0 1 0;
            -1 0 0;
            0 0 1];
    s23p = [0;-0.15175;0];
    q2_limit = [-90 60];

end