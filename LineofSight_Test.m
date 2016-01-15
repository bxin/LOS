function [] = LineofSight_Test(iTest)

%iTest = 1, 2   for two test examples

%% test input #1
% this is an example 24x1 x-vector, dy=1um for all 4 optical elements
if iTest==1
    alpha=0; %for tests
    beta=0;

    angleAS = 0;
    d = 1;
    vx = [0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 ]';
    % vx = [d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 ]';
    % vx = [0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 ]';

%% test input #2
% a system-wide tilt around x-axis of TCRS.
elseif iTest==2
    
    alpha=0; %for tests
    beta=0;

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

elseif iTest==3
    
    % a system-wide tilt around elevation-axis.
    
    alpha=0; %for tests
    beta=0;
    
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

elseif iTest==4
    
    % a system-wide tilt around elevation-axis.
    alpha=0; %for tests
    beta=45;

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
    e = [1 0 0; 0 cos(beta/180*pi) -sin(beta/180*pi); 0 sin(beta/180*pi) cos(beta/180*pi)];
    M1Motion = [(e*[0 -dM1EL*sin(angleRad)*1000 -dM1EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    M2Motion = [(e*[0 -dM2EL*sin(angleRad)*1000 -dM2EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    M3Motion = [(e*[0 -dM3EL*sin(angleRad)*1000 -dM3EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    CamMotion = [(e*[0 -dCamEL*sin(angleRad)*1000 -dCamEL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
    
elseif iTest ==5 % FEA output from Christoph

    load('FEAdata/deformedTelescopeElevationOnly');
    
    %Note: in the current FEA, the telescope pointing is a negative x-rotation from the zenith pointing.
    % so this is really -45 elevation angle, rather than +45 elevation angle.  
    alpha=asin(mirrorAxis(1))/pi*180;
    beta = -acos(mirrorAxis(3))/pi*180;
    
    angleAS = elevationRotAS;

    %original: dx dy dz in meter, Rx,Ry,Rz in Rad
    % new: dx dy dz in um, Rx, Ry, Rz in arcsec

    M1Motion = [m1m3Trans*1e6 m1m3Rot/pi*180*3600]; 
    M2Motion = [m2Trans*1e6 m2Rot/pi*180*3600];
    CamMotion = [cameraTrans*1e6 cameraRot/pi*180*3600];
    vxCG = [M1Motion M2Motion CamMotion]';

    vx = shift_CG2Vtx(vxCG);

    %% independently check FEA output
    dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    dM1M3 = -233.8; %in mm
    dM1EL = -1895;
    dM2EL = dM1EL + dM1M2;
    dM3EL = dM1EL + dM1M3;
    dCamEL = dM1EL + dM1Cam;
    
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian
    e = [1 0 0; 0 cos(beta/180*pi) -sin(beta/180*pi); 0 sin(beta/180*pi) cos(beta/180*pi)];
    M1MotionC = [(e*[0 -dM1EL*sin(angleRad)*1000 -dM1EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    M2MotionC = [(e*[0 -dM2EL*sin(angleRad)*1000 -dM2EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    M3MotionC = [(e*[0 -dM3EL*sin(angleRad)*1000 -dM3EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    CamMotionC = [(e*[0 -dCamEL*sin(angleRad)*1000 -dCamEL*(1-cos(angleRad))*1000]')' angleAS 0 0];

elseif iTest ==6 % FEA output from Christoph %same as iTest=5, but we use matrix form below
    load('FEAdata/deformedTelescopeElevationOnly');
    
    %Note: in the current FEA, the telescope pointing is a negative x-rotation from the zenith pointing.
    % so this is really -45 elevation angle, rather than +45 elevation angle.  
    alpha=asin(mirrorAxis(1))/pi*180;
    beta = -acos(mirrorAxis(3))/pi*180;
    
    angleAS = elevationRotAS;

    %original: dx dy dz in meter, Rx,Ry,Rz in Rad
    % new: dx dy dz in um, Rx, Ry, Rz in arcsec

    M1Motion = [m1m3Trans*1e6 m1m3Rot/pi*180*3600]; 
    M2Motion = [m2Trans*1e6 m2Rot/pi*180*3600];
    CamMotion = [cameraTrans*1e6 cameraRot/pi*180*3600];
    vxCG = [M1Motion M2Motion CamMotion]';

end

if iTest<=5

    [QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);
    [LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);
else
    LOSM = LOS_matrix(alpha,beta);
    LOS = LOSM * vxCG;
    LOSx = double(LOS(1));
    LOSy = double(LOS(2));
end

fprintf('angle in arcsec: %6.3f\nLoS = (%6.3f, %6.3f) arcsec\n',angleAS, LOSx, LOSy);

end




