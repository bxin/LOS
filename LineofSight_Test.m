function [] = LineofSight_Test(iTest)

%iTest = 1, 2   for two test examples

%% test input #1
% this is an example 24x1 x-vector, dy=1um for all 4 optical elements
if iTest==1
    angleAS = 0;
    d = 1;
    vx = [0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 ]';
    % vx = [d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 ]';
    % vx = [0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 ]';
    alpha=0; %for tests
    beta=0;

%% test input #2
% a system-wide tilt around x-axis of TCRS.
elseif iTest==2
    dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    dM1M3 = -233.8; %in mm
    dM1AZ = -1895+5425;
    dM2AZ = dM1AZ + dM1M2;
    dM3AZ = dM1AZ + dM1M3;
    dCamAZ = dM1AZ + dM1Cam;
    
    angleAS = 1; %angle in arcsec
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian
    M1Motion = [0 -dM1AZ*sin(angleRad)*1000 -dM1AZ*(1-cos(angleRad))*1000 angleAS 0 0];
    M2Motion = [0 -dM2AZ*sin(angleRad)*1000 -dM2AZ*(1-cos(angleRad))*1000 angleAS 0 0];
    M3Motion = [0 -dM3AZ*sin(angleRad)*1000 -dM3AZ*(1-cos(angleRad))*1000 angleAS 0 0];
    CamMotion = [0 -dCamAZ*sin(angleRad)*1000 -dCamAZ*(1-cos(angleRad))*1000 angleAS 0 0];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
    
    alpha=0; %for tests
    beta=0;

elseif iTest==3
    
    % a system-wide tilt around elevation-axis.
    
    dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    dM1M3 = -233.8; %in mm
    dM1EL = -1895;
    dM2EL = dM1EL + dM1M2;
    dM3EL = dM1EL + dM1M3;
    dCamEL = dM1EL + dM1Cam;
    
    angleAS = 1; %angle in arcsec
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian
    M1Motion = [0 -dM1EL*sin(angleRad)*1000 -dM1EL*(1-cos(angleRad))*1000 angleAS 0 0];
    M2Motion = [0 -dM2EL*sin(angleRad)*1000 -dM2EL*(1-cos(angleRad))*1000 angleAS 0 0];
    M3Motion = [0 -dM3EL*sin(angleRad)*1000 -dM3EL*(1-cos(angleRad))*1000 angleAS 0 0];
    CamMotion = [0 -dCamEL*sin(angleRad)*1000 -dCamEL*(1-cos(angleRad))*1000 angleAS 0 0];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
    
    alpha=0; %for tests
    beta=0;

elseif iTest ==4 % FEA output from Christoph
    load('christoph/undeformedTelescope45');
    alpha=asin(mirrorAxis(1))/pi*180;
    beta = acos(mirrorAxis(3))/pi*180;
    
    load('christoph/deformedTelescopeElevationOnly');
    angleAS = elevationRotAS;
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian    
    
    M1Motion = [];
end

%% coordinate transforms


[QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);
[LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);

fprintf('angle in arcsec: %3.1f\nLoS = (%5.2f, %5.2f) arcsec\n',angleAS, LOSx, LOSy);

end




