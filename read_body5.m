function read_body5

    global m5 q5 dq5 J5p s45p C45 rho5p C55 q5_init s56p C56 q5_limit

    q5 = 0;
    s56p = [0.0825;0;-0.04475];
    C56 = [0 0 1;
            1 0 0
            0 1 0];
    q5_limit = [-45 100];

end