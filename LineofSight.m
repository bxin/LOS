function [LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam)

% LINEOFSIGHT calculates the line of sight deviation from the overall 
% sensitivity matrix senM

[SxFinal, SyFinal] = make_LOS_senM();

%alpha (deg) is the athimuth angle of the telescope (by convention, counter-RH rule)
%beta (deg) is the elevation angle of the telescope
alphaRad=alpha/180*pi; %in rad now
betaRad = beta/180*pi; %in rad now

D_EL_AZ = 5425 * 1000; %in um
D_EL_M1 = 1895 * 1000; %in um
% D_M2_M1 = (6628-472) * 1000; %in um
% D_Cam_M1 = (5336-1938) * 1000; %in um
% D_M1_M3 = 233.8 * 1000; %in um

T_TCRS_ACRS = [cos(alphaRad) -sin(alphaRad) 0 0;
               sin(alphaRad) cos(alphaRad) 0 0;
               0             0       1 0;
               0             0       0 1];
T_ACRS_RCRS=[1 0 0 0; 0 1 0 0; 0 0 1 -D_EL_AZ; 0 0 0 1];         
T_RCRS_ACRS=[1 0 0 0; 0 1 0 0; 0 0 1 D_EL_AZ; 0 0 0 1];         

T_RCRS_ECRS=[1     0         0      0;
             0 cos(betaRad) sin(betaRad) 0;
             0 -sin(betaRad) cos(betaRad) 0;
             0   0         0        1];
         
T_ECRS_OCRS = [1 0 0 0; 0 1 0 0 ; 0 0 1 D_EL_M1; 0 0 0 1];        
T_OCRS_ECRS = [1 0 0 0; 0 1 0 0 ; 0 0 1 -D_EL_M1; 0 0 0 1];        

% T_OCRS_M2 = [1 0 0 0; 0 1 0 0; 0 0 1 -D_M2_M1; 0 0 0 1]; 
% T_M2_OCRS = [1 0 0 0; 0 1 0 0; 0 0 1 D_M2_M1; 0 0 0 1]; 
% 
% T_OCRS_M3 = [1 0 0 0; 0 1 0 0; 0 0 1  D_M1_M3; 0 0 0 1]; 
% T_M3_OCRS = [1 0 0 0; 0 1 0 0; 0 0 1  -D_M1_M3; 0 0 0 1]; 
% 
% T_OCRS_Cam = [1 0 0 0; 0 1 0 0; 0 0 1 -D_Cam_M1; 0 0 0 1]; 
% T_Cam_OCRS = [1 0 0 0; 0 1 0 0; 0 0 1 D_Cam_M1; 0 0 0 1]; 

T_reverse_z = [1 0 0 0; 0 -1 0 0 ; 0 0 -1 0; 0 0 0 1];

for i=1:4
    if i==1
        Q = QM1;
    elseif i==2
        Q = QM2;
    elseif i==3
        Q = QM3;
    elseif i==4
        Q = QCam;
    end
    QA=T_TCRS_ACRS*Q*transpose(T_TCRS_ACRS); %now in ACRS
    QR=T_ACRS_RCRS*QA*T_RCRS_ACRS; %now in RCRS
    QE=T_RCRS_ECRS*QR*transpose(T_RCRS_ECRS); %now in ECRS
    QO=T_ECRS_OCRS*QE*T_OCRS_ECRS; %now in OCRS
    if i==1
        QM1 = T_reverse_z*QO*transpose(T_reverse_z);
    elseif i==2
        %         QM2 = T_reverse_z*T_OCRS_M2*QO*T_M2_OCRS*T_reverse_z';
        QM2 = T_reverse_z*QO*transpose(T_reverse_z);
    elseif i==3
        %         QM3 = T_reverse_z*T_OCRS_M3*QO*T_M3_OCRS*T_reverse_z';
        QM3 = T_reverse_z*QO*transpose(T_reverse_z);
    elseif i==4
        %         QCam = T_reverse_z*T_OCRS_Cam*QO*T_Cam_OCRS*T_reverse_z';
        QCam = T_reverse_z*QO*(T_reverse_z);
    end
end
Qfinal = blkdiag(QM1,QM2,QM3,QCam);

LOSx = -sum(sum(Qfinal.*SxFinal)); %by convention, clockwise rotation around z is positive azimuth angle
LOSy = sum(sum(Qfinal.*SyFinal));

end




