function residual

    global Q_g Q_c Ta r_hat
    global Yp

    alpha = Q_g - Q_c;
    Yp = Ta - alpha - r_hat;
    

end