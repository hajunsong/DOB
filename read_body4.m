function read_body4

    global m4 q4 dq4 J4p s34p C34 rho4p C44 q4_init s45p C45 q4_limit

    q4 = 0;
    C45 = [0 0 1;
            1 0 0
            0 1 0];
    s45p = [0.0245;0.08675;-0.0245];
    q4_limit = [-180 45];

end