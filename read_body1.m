function read_body1

    global m1 q1 dq1 J1p s01p C01 rho1p C11 q1_init s12p C12 q1_limit

    q1 = 0;
    C12 = [1 0 0;
            0 0 1;
            0 -1 0];
    s12p = [0;0.0245;0.042];
    q1_limit = [-180,180];

end