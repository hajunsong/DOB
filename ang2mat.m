function mat = ang2mat(ang_z1, ang_x, ang_z2, flag)

    if flag == 1
        z1 = ang_z1*pi/180;
        x =  ang_x*pi/180;
        z2 = ang_z2*pi/180;
    elseif flag == 0
        z1 = ang_z1;
        x =  ang_x;
        z2 = ang_z2;
    end

    Rz1 = [cos(z1) -sin(z1) 0;sin(z1) cos(z1) 0;0 0 1];
    Rx = [1 0 0;0 cos(x) -sin(x);0 sin(x) cos(x)];
    Rz2 = [cos(z2) -sin(z2) 0;sin(z2) cos(z2) 0;0 0 1];

    mat = Rz1*Rx*Rz2;
    
end