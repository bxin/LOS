function [] = LineofSight_Test(iTest)

%iTest = 1, 2

%% test input #1
% this is an example 24x1 x-vector, dy=1um for all 4 optical elements
if iTest==1
    vx = [0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 ]';
    % vx = [0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 ]';
    % vx = [1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 ]';

%% test input #2
elseif iTest==2
    dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    dM1M3 = -233.8; %in mm
    dM1EL = -1895+5425;
    dM2EL = dM1EL + dM1M2;
    dM3EL = dM1EL + dM1M3;
    dCamEL = dM1EL + dM1Cam;
    
    
    angleAS = 1; %angle in arcsec
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian
    M1Motion = [-dM1EL*(1-cos(angleRad))*1000 0 -dM1EL*sin(angleRad)*1000 angleAS 0 0];
    M2Motion = [-dM2EL*(1-cos(angleRad))*1000 0 -dM2EL*sin(angleRad)*1000 angleAS 0 0];
    M3Motion = [-dM3EL*(1-cos(angleRad))*1000 0 -dM3EL*sin(angleRad)*1000 angleAS 0 0];
    CamMotion = [-dCamEL*(1-cos(angleRad))*1000 0 -dCamEL*sin(angleRad)*1000 angleAS 0 0];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
end

%% coordinate transforms

alpha=0; %for tests
beta=0; 

[QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);
[LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);

fprintf('LoS = (%5.2f, %5.2f) arcsec\n',LOSx, LOSy);

end




