function [QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx)

for i=1:4 % for M1,M2,M3,Cam
    
    dz=vx((i-1)*6+1);
    dx=vx((i-1)*6+2);
    dy=vx((i-1)*6+3);
    tx=vx((i-1)*6+4)/3600/180*pi; %convert from arcsec to radian
    ty=vx((i-1)*6+5)/3600/180*pi; %convert from arcsec to radian
    tz=vx((i-1)*6+6)/3600/180*pi; %convert from arcsec to radian
    %the following T, Rx, Ry, Rz matrices are for rotation operations in TCRS
    % the matrix elements have inversed sign, compared to those for coordinate system changes
    T=[1 0 0 dx;
        0 1 0 dy;
        0 0 1 dz;
        0 0 0 1];
    Rx=[1 0 0 0;
        0 1 -tx 0;
        0 tx 1 0;
        0 0 0 1];
    Ry=[1 0 ty 0;
        0 1 0 0;
        -ty 0 1 0;
        0 0 0 1];
    Rz=[1 -tz 0 0;
        tz 1 0 0;
        0 0 1 0;
        0 0 0 1];

    Q=T*Rx*Ry*Rz;
    if i==1
        QM1 = Q;
    elseif i==2
        QM2 = Q;
    elseif i==3
        QM3 = Q;
    elseif i==4
        QCam = Q;
    end
end

end
