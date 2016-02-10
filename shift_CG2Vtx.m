function vxnew = shift_CG2Vtx(vx, alpha, beta, mysym)

% author: Bo Xin (bxin@lsst.org)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

%input is the coordinates of the CG (dummy mass) in the TCRS
%output is the cooridinates of the optical vertices (pivoting points as defined in ZCRS) in the TCRS

if nargin<4
    mysym = 0;
end
load('FEAdata/undeformedTelescope45');

stickTopz = [-1895 6156-1895 -1895-233.8 3398-1895]*1000;%in um

ea = [cos(alpha/180*pi) -sin(alpha/180*pi) 0; sin(alpha/180*pi) cos(alpha/180*pi) 0; 0 0 1];
ez = [1 0 0; 0 cos(beta/180*pi) -sin(beta/180*pi); 0 sin(beta/180*pi) cos(beta/180*pi)];
    
%plug alpha = 0, beta = -45 into the above, so make everything independent of undeformedTelescope45.mat
% ea45 is a unit matrix, so ignore. To get ez45inv from ez45, reverse the sign of the angle
ez45inv = [1 0 0; 0 cos(pi/4) -sin(pi/4); 0 sin(pi/4) cos(pi/4)];
 
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
        
    if i==1 %M1
        stick = ea*ez*([0 0 stickTopz(i)]'-ez45inv*(m1m3Pos - altitudeCenter)'*1e6);
    elseif i==2 %M2
        stick = ea*ez*([0 0 stickTopz(i)]'-ez45inv*(m2Pos - altitudeCenter)'*1e6);
    elseif i==3 %M3
        stick = ea*ez*([0 0 stickTopz(i)]'-ez45inv*(m1m3Pos - altitudeCenter)'*1e6);
    elseif i==4 %Cam
        stick = ea*ez*([0 0 stickTopz(i)]'-ez45inv*(cameraPos - altitudeCenter)'*1e6);
    end
    
    vxnew((i-1)*6+1:(i-1)*6+3) = vx((i-1)*6+1:(i-1)*6+3)+(Rx*Ry*Rz-eye(3))*stick;
    
end

end
