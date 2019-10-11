function inverse_dynamics

    global ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 r1t r2t r3t r4t r5t r6t
    global L1 L2 L3 L4 L5 L6 K1 K2 K3 K4 K5 K6 B1 B2 B3 B4 B5 B6 D1 D2 D3 D4 D5 D6
    
    global f10 f21 f32 f43 f54 f65 n10 n21 n32 n43 n54 n65
    
    dY1h = B1*ddq1 + D1;
    dY2h = dY1h + B2*ddq2 + D2;
    dY3h = dY2h + B3*ddq3 + D3;
    dY4h = dY3h + B4*ddq4 + D4;
    dY5h = dY4h + B5*ddq5 + D5;
    dY6h = dY5h + B6*ddq6 + D6;
    
    R10h = -L1 + K1*dY1h + K2*B2*ddq2 + K3*B3*ddq3 + K4*B4*ddq4 + K5*B5*ddq5 + K6*B6*ddq6;
    R21h = -L2 + K2*dY2h + K3*B3*ddq3 + K4*B4*ddq4 + K5*B5*ddq5 + K6*B6*ddq6;
    R32h = -L3 + K3*dY3h + K4*B4*ddq4 + K5*B5*ddq5 + K6*B6*ddq6;
    R43h = -L4 + K4*dY4h + K5*B5*ddq5 + K6*B6*ddq6;
    R54h = -L5 + K5*dY5h + K6*B6*ddq6;
    R65h = -L6 + K6*dY6h;
    
    T1 = [eye(3), zeros(3); -r1t, eye(3)];
    T2 = [eye(3), zeros(3); -r2t, eye(3)];
    T3 = [eye(3), zeros(3); -r3t, eye(3)];
    T4 = [eye(3), zeros(3); -r4t, eye(3)];
    T5 = [eye(3), zeros(3); -r5t, eye(3)];
    T6 = [eye(3), zeros(3); -r6t, eye(3)];
    
    R10b = T1*R10h;
    R21b = T2*R21h;
    R32b = T3*R32h;
    R43b = T4*R43h;
    R54b = T5*R54h;
    R65b = T6*R65h;
    
    f10 = R10b(1:3,1); n10 = R10b(4:6,1);
    f21 = R21b(1:3,1); n21 = R21b(4:6,1);
    f32 = R32b(1:3,1); n32 = R32b(4:6,1);
    f43 = R43b(1:3,1); n43 = R43b(4:6,1);
    f54 = R54b(1:3,1); n54 = R54b(4:6,1);
    f65 = R65b(1:3,1); n65 = R65b(4:6,1);

end