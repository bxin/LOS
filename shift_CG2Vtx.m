function vxnew = shift_CG2Vtx(vx, mysym)

%input is the coordinates of the CG (dummy mass) in the TCRS
%output is the cooridinates of the optical vertices (pivoting points as defined in ZCRS) in the TCRS

if nargin<2
    mysym = 0;
end
load('christoph/undeformedTelescope45');

% how much the optical vertex is above the CG mass, for M1,M2,M3,Cam
m1m3CGz = (m1m3Pos(3)-altitudeCenter(3))*sqrt(2)*1000; %in mm
m2CGz = (m2Pos(3)-altitudeCenter(3))*sqrt(2)*1000;
CamCGz = (cameraPos(3)-altitudeCenter(3))*sqrt(2)*1000;

stickLen = [-1895-m1m3CGz 6156-1895-m2CGz -1895-233.8-m1m3CGz 3398-1895-CamCGz]*1000;%in um
if mysym == 1
    vxnew = sym('vx',[length(vx) 1]);
else
    vxnew = zeros(length(vx),1);
end
vx = [vx(1:12); vx(1:6); vx(13:18)]; %insert m3motion, which is same as m1motion

for i=1:4 % for M1, M2, M3, Cam

    tRx = vx((i-1)*6+4); %t is for theta
    tRy = vx((i-1)*6+5);
    tRz = vx((i-1)*6+6);
    vxnew((i-1)*6+4)=tRx; %Rx is the same
    vxnew((i-1)*6+5)=tRy; %Ry is the same
    vxnew((i-1)*6+6)=tRz; %Rz is the same

    tRx = tRx/3600/180*pi; %convert arcsec into rad
    tRy = tRy/3600/180*pi;
    tRz = tRz/3600/180*pi;
    
    %     Rx = [1 0 0; 0 cos(tRx) -sin(tRx); 0 sin(tRx) cos(tRx)]; %3D rotation matrices
    %     Ry = [cos(tRy) 0 sin(tRy); 0 1 0; -sin(tRy) 0 cos(tRy)];
    %     Rz = [cos(tRz) -sin(tRz) 0; sin(tRz) cos(tRz) 0; 0 0 1];
    
    Rx = [1 0 0; 0 1 -tRx; 0 tRx 1]; %3D rotation matrices for small angles
    Ry = [1 0 tRy; 0 1 0; -tRy 0 1];
    Rz = [1 -tRz 0; tRz 1 0; 0 0 1];
        
    stick = [0 0 stickLen(i)]';
    vxnew((i-1)*6+1:(i-1)*6+3) = vx((i-1)*6+1:(i-1)*6+3)+(Rx*Ry*Rz-eye(3))*stick;
    
end

end
