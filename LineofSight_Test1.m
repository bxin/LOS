function [] = LineofSight_Test1(elev,azTorque,elTorque, matrixForm, instant, outMatrix)

% author: Bo Xin (bxin@lsst.org)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

% 2/8/16  Tests with new data from Christoph, for 
% elev = 15, 45, 86.5, 
% azTorque = 0 or 1
% elTorque = 0 or 1
% matrix form = 1 or 0
% instant = 1 (work with 5001 instances) or 0 (work only with the last instance)
% outMatrix is the name of the desired filename for the LOS matrix output
%           use '' if no output is desired.

% example usage:
%  LineofSight_Test1(45,1,1,1,0,'simulink/LOSM_matrix45.mat')
%  LineofSight_Test1(15,1,1,1,0,'simulink/LOSM_matrix15.mat')
%  LineofSight_Test1(86.5,1,1,1,0,'simulink/LOSM_matrix86.5.mat')

% LineofSight_Test1(x,x,x,0,0) doesn't work, not yet implemented.

if nargin<6
    outMatrix = '';
end

azTorque = azTorque*1000;
elTorque = elTorque*1000;
if elev==86.5
    load(sprintf('testData/deg%4.1f/stepResponseTestData%4.1f El Torque %dNm Az Torque %dNm.mat', ...
        elev, elev, elTorque, azTorque));
else
    load(sprintf('testData/deg%d/stepResponseTestData%d El Torque %dNm Az Torque %dNm', ...
        elev, elev, elTorque, azTorque));
end

%by convention, clockwise rotation around z is positive azimuth angle
% FEA output already follows that convention
azimuthReading = azimuthReading*cos((elev)/180*pi);

%Note: in the current FEA, the telescope pointing is a negative x-rotation from the zenith pointing.
% this is obvious if we look at the undeformed FEA mat data files.
% so the 15, 45 and 86.5 deg angles are really negative elevation angle.
alpha=0;
beta = -(90-elev); %beta is zenith angle, = 90 - elevation angle
    
%original: dx dy dz in meter, Rx,Ry,Rz in Rad
% new: dx dy dz in um, Rx, Ry, Rz in arcsec
transIDX = [1 2 3 7 8 9 13 14 15];
angleIDX = [4 5 6 10 11 12 16 17 18];
if instant
    vxCG = testData(end,:);
    angleAS = vxCG(4)/pi*180*3600; % verified = elevationReading(end)
else
    vxCG = testData;
end
vxCG(:,transIDX) = vxCG(:,transIDX)*1e6;
vxCG(:,angleIDX) = vxCG(:,angleIDX)/pi*180*3600;

vxCG = transpose(vxCG);
if ~matrixForm
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
%     
%     e = [1 0 0; 0 cos(beta/180*pi) -sin(beta/180*pi); 0 sin(beta/180*pi) cos(beta/180*pi)];
%     M1MotionC = [(e*[0 -dM1EL*sin(angleRad)*1000 -dM1EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
%     M2MotionC = [(e*[0 -dM2EL*sin(angleRad)*1000 -dM2EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
%     M3MotionC = [(e*[0 -dM3EL*sin(angleRad)*1000 -dM3EL*(1-cos(angleRad))*1000]')' angleAS 0 0];
%     CamMotionC = [(e*[0 -dCamEL*sin(angleRad)*1000 -dCamEL*(1-cos(angleRad))*1000]')' angleAS 0 0];
%     plot(1:24,vx,'-r', 1:24,[M1MotionC M2MotionC M3MotionC CamMotionC]','-k');
%     
%     % another check
%     load('testData/undeformedTelescope15.mat');
%     a=cameraPos-m1m3Pos;
%     cameraTrans = testData(end,13:15);
%     m1m3Trans = testData(end,1:3);
%     a1=cameraPos+cameraTrans-(m1m3Pos+m1m3Trans);
%     fprintf('%6.4f vs %6.4f\n',angleAS, acos(sum(a.*a1)/sqrt(sum(a.^2)*sum(a1.^2)))/pi*180*3600);
%     m2Trans = testData(end,7:9);
%     a=m2Pos-m1m3Pos;
%     a1=m2Pos+m2Trans-(m1m3Pos+m1m3Trans);
%     fprintf('%6.4f vs %6.4f\n',angleAS, acos(sum(a.*a1)/sqrt(sum(a.^2)*sum(a1.^2)))/pi*180*3600);
%

[QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);
[LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);

else
    LOSM = LOS_matrix(alpha,beta);
    if ~isempty(outMatrix)
        save(outMatrix,'LOSM');
    end
    LOS = LOSM * vxCG;
    LOSx = LOS(1,:);
    LOSy = LOS(2,:);
end

if instant
    fprintf('angle in arcsec: %6.4f\nLoS = (%6.4f, %6.4f) arcsec\n',angleAS, LOSx, LOSy);

else
    clf;
    subplot(1,2,1);
    plot(time,azimuthReading,'-r',time,LOSx,'-k');
    legend({'encoder reading','LOS'},'location','best');
    text(0.2,0.5,'Azimuth motion','Units','Normalized');
    text(0.2,0.6,sprintf('%4.1f deg position',elev),'Units','Normalized');
    xlabel('time (s)'); ylabel('arcsec');
    subplot(1,2,2);
    plot(time,elevationReading,'-r',time,LOSy,'-k');
    legend({'encoder reading','LOS'},'location','best');
    text(0.2,0.5,'Elevation motion','Units','Normalized');
    text(0.2,0.6,sprintf('%4.1f deg position',elev),'Units','Normalized');
    xlabel('time (s)'); ylabel('arcsec');
end

end




