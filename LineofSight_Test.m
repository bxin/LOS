function [] = LineofSight_Test(iTest)

% with telescope pointing at zenith
%iTest = 1: move all optics (as defined by Zemax, as opposed to FEA dummy mass) along y for 1um
%           tests along x and z can be easily done by uncommenting the commented lines
%      = 11: same as 1. but move along line that is 45 deg to x-axis, and in x-y plane
%iTest = 2: tilt all optics (as defined by Zemax, as opposed to FEA dummy mass) around 
%           x-axis of TCRS for 1 arcsec
%      = 12: same as 2. but tilt around y-axis
%      = 22: same as 2. but tilt around z-axis
%iTest = 3: tilt all optics (as defined by Zemax, as opposed to FEA dummy mass) around 
%           x-axis of ECRS for 1 arcsec

% now at 45 deg zenith angle
%iTest = 4: tilt all optics (as defined by Zemax, as opposed to FEA dummy mass) around 
%           x-axis of ECRS for 1 arcsec
%      = 14: same as 4. but rotate around z-axis.
%iTest = 5: tilt all optics (defined as FEA dummy mass) around 
%           x-axis of ECRS for 1 arcsec. The input is from FEA
%        15: same as 5. but rotate around z-axis.
%iTest = 6: Same as 5. But here we derive the 2x18 matrix, so that LOS=LOSM*vxCG;
%      = 16: Same as 15. But here we derive the 2x18 matrix, so that LOS=LOSM*vxCG;


if iTest==1
%% test input #1
% a system-wide translation in TCRS.
    alpha=0; %for tests
    beta=0;

    angleAS = 0;
    d = 1;
    vx = [0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 ]';
    % vx = [d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 ]';
    % vx = [0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 ]';
elseif iTest == 11
    %dy = dx = d
    alpha=0; %for tests
    beta=0;

    angleAS = 0;
    d = 1;
    vx = [d d 0 0 0 0 d d 0 0 0 0 d d 0 0 0 0 d d 0 0 0 0 ]'; 

elseif iTest==2
%% test input #2
% a system-wide tilt around x-axis of TCRS.    
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
    
elseif iTest==12
%% test input #12
% a system-wide tilt around y-axis of TCRS.     
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
    M1Motion = [dM1AZ*sin(angleRad)*1000 0 -dM1AZ*(1-cos(angleRad))*1000 0 angleAS 0];
    M2Motion = [dM2AZ*sin(angleRad)*1000 0 -dM2AZ*(1-cos(angleRad))*1000 0 angleAS 0];
    M3Motion = [dM3AZ*sin(angleRad)*1000 0 -dM3AZ*(1-cos(angleRad))*1000 0 angleAS 0];
    CamMotion = [dCamAZ*sin(angleRad)*1000 0 -dCamAZ*(1-cos(angleRad))*1000 0 angleAS 0];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
        
    %by convention, clockwise rotation around z is positive azimuth angle
    angleAS = -angleAS;
    
elseif iTest==22
    % azimuth rotation around z-axis
    alpha=0; %for tests
    beta=0;
    
    angleAS = 1; %angle in arcsec
    M1Motion = [0 0 0 0 0 angleAS];
    M2Motion = [0 0 0 0 0 angleAS];
    M3Motion = [0 0 0 0 0 angleAS];
    CamMotion = [0 0 0 0 0 angleAS];
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

elseif iTest==14
    
    % a azimuth rotation (around +z)
    alpha=0; %for tests
    beta=45;

    dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    dM1M3 = -233.8; %in mm
    dM1EL = -1895;
    dM2EL = dM1EL + dM1M2;
    dM3EL = dM1EL + dM1M3;
    dCamEL = dM1EL + dM1Cam;
    
    angleAS = 1; %*3600*3.5; %angle in arcsec
    angleDEG = angleAS/3600;
    angleRad= angleDEG/180*pi; %angle in radian
    M1Motion = [dM1EL*sin(angleRad)*1000*sqrt(2)/2 -dM1EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    M2Motion = [dM2EL*sin(angleRad)*1000*sqrt(2)/2 -dM2EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    M3Motion = [dM3EL*sin(angleRad)*1000*sqrt(2)/2 -dM3EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    CamMotion = [dCamEL*sin(angleRad)*1000*sqrt(2)/2 -dCamEL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    vx = [M1Motion M2Motion M3Motion CamMotion]';
    
    %by convention, clockwise rotation around z is positive azimuth angle
    angleAS = -angleAS*sqrt(2)/2;
        
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

    vx = shift_CG2Vtx(vxCG, alpha, beta);

    %% independently check FEA output, not really needed
    %     dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    %     dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    %     dM1M3 = -233.8; %in mm
    %     dM1EL = -1895;
    %     dM2EL = dM1EL + dM1M2;
    %     dM3EL = dM1EL + dM1M3;
    %     dCamEL = dM1EL + dM1Cam;
    %
    %     angleDEG = angleAS/3600;
    %     angleRad= angleDEG/180*pi; %angle in radian
    %     e = [1 0 0; 0 cos(beta/180*pi) -sin(beta/180*pi); 0 sin(beta/180*pi) cos(beta/180*pi)];
    %     M1MotionC = [(e*[0 -dM1EL*sin(angleRad)*1000 -dM1EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    %     M2MotionC = [(e*[0 -dM2EL*sin(angleRad)*1000 -dM2EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    %     M3MotionC = [(e*[0 -dM3EL*sin(angleRad)*1000 -dM3EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    %     CamMotionC = [(e*[0 -dCamEL*sin(angleRad)*1000 -dCamEL*(1-cos(angleRad))*1000]')' angleAS 0 0];
    %     plot(vx-[M1MotionC M2MotionC M3MotionC CamMotionC]');
    
    % another check
    %     a=cameraPos-m1m3Pos;
    %     a1=cameraPos+cameraTrans-(m1m3Pos+m1m3Trans);
    %     fprintf('%6.4f vs %6.4f\n',elevationRotAS, acos(sum(a.*a1)/sqrt(sum(a.^2)*sum(a1.^2)))/pi*180*3600);
    %     a=m2Pos-m1m3Pos;
    %     a1=m2Pos+m2Trans-(m1m3Pos+m1m3Trans);
    %     fprintf('%6.4f vs %6.4f\n',elevationRotAS, acos(sum(a.*a1)/sqrt(sum(a.^2)*sum(a1.^2)))/pi*180*3600);

elseif iTest ==15 % FEA output from Christoph

    load('FEAdata/deformedTelescopeAzimuthOnly');
    
    %Note: in the current FEA, the telescope pointing is a negative x-rotation from the zenith pointing.
    % so this is really -45 elevation angle, rather than +45 elevation angle.  
    alpha=asin(mirrorAxis(1))/pi*180;
    beta = -acos(mirrorAxis(3))/pi*180;
    
    angleAS = azimuthRotAS;

    %original: dx dy dz in meter, Rx,Ry,Rz in Rad
    % new: dx dy dz in um, Rx, Ry, Rz in arcsec

    M1Motion = [m1m3Trans*1e6 m1m3Rot/pi*180*3600]; 
    M2Motion = [m2Trans*1e6 m2Rot/pi*180*3600];
    CamMotion = [cameraTrans*1e6 cameraRot/pi*180*3600];
    vxCG = [M1Motion M2Motion CamMotion]';

    vx = shift_CG2Vtx(vxCG, alpha, beta);
     
    %% independently check FEA output, not really needed
    %     dM1M2 = 6628-472; % =6156, in mm (zemax has 6156.2006)
    %     dM1Cam = 5336-1938; % = 3398 in mm (zemax has 3398.6)
    %     dM1M3 = -233.8; %in mm
    %     dM1EL = -1895;
    %     dM2EL = dM1EL + dM1M2;
    %     dM3EL = dM1EL + dM1M3;
    %     dCamEL = dM1EL + dM1Cam;
    %
    %     angleDEG = -angleAS/3600;
    %     angleRad= angleDEG/180*pi; %angle in radian
    %     M1MotionC = [dM1EL*sin(angleRad)*1000*sqrt(2)/2 -dM1EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    %     M2MotionC = [dM2EL*sin(angleRad)*1000*sqrt(2)/2 -dM2EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    %     M3MotionC = [dM3EL*sin(angleRad)*1000*sqrt(2)/2 -dM3EL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];
    %     CamMotionC = [dCamEL*sin(angleRad)*1000*sqrt(2)/2 -dCamEL*(1-cos(angleRad))*1000*sqrt(2)/2 0 0 0 angleAS];

    %by convention, clockwise rotation around z is positive azimuth angle
    % FEA output already follows that convention
    angleAS = angleAS * sqrt(2)/2;

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
    
elseif iTest ==16 % FEA output from Christoph %same as iTest=5, but we use matrix form below
    load('FEAdata/deformedTelescopeAzimuthOnly');
    
    %Note: in the current FEA, the telescope pointing is a negative x-rotation from the zenith pointing.
    % so this is really -45 elevation angle, rather than +45 elevation angle.  
    alpha=asin(mirrorAxis(1))/pi*180;
    beta = -acos(mirrorAxis(3))/pi*180;
    
    angleAS = azimuthRotAS;

    %original: dx dy dz in meter, Rx,Ry,Rz in Rad
    % new: dx dy dz in um, Rx, Ry, Rz in arcsec

    M1Motion = [m1m3Trans*1e6 m1m3Rot/pi*180*3600]; 
    M2Motion = [m2Trans*1e6 m2Rot/pi*180*3600];
    CamMotion = [cameraTrans*1e6 cameraRot/pi*180*3600];
    vxCG = [M1Motion M2Motion CamMotion]';
    
    %by convention, clockwise rotation around z is positive azimuth angle
    % FEA output already follows that convention
    angleAS = angleAS * sqrt(2)/2;

end

if mod(iTest,10)<=5

    [QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);
    [LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);
else
    LOSM = LOS_matrix(alpha,beta);
    LOS = LOSM * vxCG;
    LOSx = LOS(1);
    LOSy = LOS(2);
end

fprintf('angle in arcsec: %6.4f\nLoS = (%6.4f, %6.4f) arcsec\n',angleAS, LOSx, LOSy);

end




