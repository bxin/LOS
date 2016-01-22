clear all
close all
clc

%% Selection of analysis parameters

%Altitud position
position = 2; % 1 = 86.5º; 2 = 45º; 3 = 15º


%% Read Input File

IniGeometricData;

Jalt = 2.6e6; %La inercia de Noelia de 2.1e6 no es fiable porque es respecto a coordenadas globales Ixx %2.6e6; % Mikel tenía 2700000  % 2.6e6 en presentación PDR de Phase
Jaz = 9.06e6; %Noelia Izz 45deg %Mikel tenía 7000014  % 10.1e6 en presentación PDR de Phase

switch position
    case 1
        data = load('modal_disp_86_100modes_v4.mat');
        AltAng = 86.5 *pi/180;
    case 2
        data = load('modal_disp_45_100modes_v4.mat');
        AltAng = 45 *pi/180;
    case 3
        data = load('modal_disp_15_100modes_v4.mat');
        AltAng = 15 *pi/180;
end

datos_fem = data.mat_modelo;

[nfil,ncol]=size(datos_fem);

psi=1/100;  %amortiguamiento general para todos los modos.

for i=1:ncol
    stiff_maq(i,1)=(2*pi()*datos_fem(1,i))^2;
    damp_maq(i,1)=2*pi()*datos_fem(1,i)*2*psi;
    for j=2:nfil
        disp(j-1,i)=datos_fem(j,i);
    end
end

Phi=disp; %vel=disp;
PhiT=disp';

nModos=ncol;
nNodos= nfil-1;

%% Generate initial matrices for obtaining TFs

Mass = eye(nModos);
Stiff = diag(stiff_maq);
Damp = diag(damp_maq);

%Espacio de estado genérico para M, C y K en coordenadas modales q/Fq
%Para estas matrices: Nx=2*nModos; Nu=nModos; Ny=2*nModos
A=[zeros(nModos), eye(nModos); -inv(Mass)*Stiff, -inv(Mass)*Damp];
B = [zeros(nModos); inv(Mass)];
C = eye(2*nModos);
%D = zeros(2*nModos, nModos);
D = 0;
sis = ss (A, B, C, D);   %Solo para comprobar

%Espacio de estado en x/F
%Para estas matrices: Nx=2*nModos; Nu=nNodos; Ny=2*nNodos

Ax = A;
Bx = B*PhiT;
Cx = [Phi zeros(nNodos, nModos); zeros(nNodos, nModos) Phi]*C;
%Dx = [Phi zeros(nNodos, nModos); zeros(nNodos, nModos) Phi]*D*PhiT;
Dx = 0;
sisx = ss (Ax, Bx, Cx, Dx);

%% Coordinate transformation matrices for positions

CoordTransfMatrix33 = zeros(3,3);
CoordTransfMatrix33(1,1) = 1;
CoordTransfMatrix33(2,2) = sin(AltAng);
CoordTransfMatrix33(3,3) = sin(AltAng);
CoordTransfMatrix33(2,3) = -cos(AltAng);
CoordTransfMatrix33(3,2) = cos(AltAng);
InvCoordTransfMatrix33 = inv(CoordTransfMatrix33);

%% Obtain vectors for altitud torque input and theta output 

%Index for Altitud drives
IndexFirstInput = 31; %6*nPOI + 1
IndexLastInput = 66;  %6*nPOI + 2*nPAltDrives

k= 0;
InAltT = zeros (nNodos, 1);
for(j=IndexFirstInput:IndexLastInput)
    if mod(k,2) %Odd
        InAltT(j) = (1/R_alt_d)/(2*n_drives_alt); %Torque to force. 
    else 
        InAltT(j) = -(1/R_alt_d)/(2*n_drives_alt);
    end 
    k= k+1;
end

%Index for Altitud encoders
IndexFirstOutput = 163; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 1
IndexLastOutput = 170; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders
nEncAlt=4;

OutAltP = zeros(1, 2*nNodos);
k = 0;
% old incorrect
for(j=IndexFirstOutput:IndexLastOutput)
    if mod(k,2) %Odd
        OutAltP(j) = (1/nEncAlt)/R_alt_e; % dY %linear mov, to rot mov. 
    else 
        OutAltP(j) = (-1/nEncAlt)/R_alt_e; % dZ
    end 
    k= k+1;
end

% % new fixed
% multipliers = [3.1624 1.2383 3.1624 1.2383];
% ii = 1;
% IndexLastOutput = IndexFirstOutput + 2;
% nEncAlt = 1;
% for(j=IndexFirstOutput:IndexLastOutput)
%     if mod(k,2) %Odd
%         OutAltP(j) = 0; % ignore dY 
%     else 
%         OutAltP(j) = (-1/nEncAlt)/R_alt_e * multipliers(ii); % dZ
%         ii = ii + 1;
%     end 
%     k= k+1;
% end


%% Obtain vectors for azimuth torque input and theta output 

%Index for Azimuth drives
IndexFirstInput = 67; %6*nPOI + 2*nPAltDrives + 1
IndexLastInput = 162; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives

InAzT= zeros (nNodos, 1);
k = 0;
for(j=IndexFirstInput:IndexLastInput)
    if mod(k,2) %Odd
        InAzT(j) = -(1/R_az_d)/(2*n_drives_az); %Torque to force. 
    else %Even
        InAzT(j) = (1/R_az_d)/(2*n_drives_az);
    end 
    k= k+1;
end

%Index for Azimuth encoders
IndexFirstOutput = 171; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 1
IndexLastOutput = 178;  %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 2*nPAzEncoders
nEncAz=4;

OutAzP = zeros(1, 2*nNodos);
k= 0;
for(j=IndexFirstOutput:IndexLastOutput)
    if mod(k,2) %Odd
        OutAzP(j) = (-1/nEncAz)/R_az_e; %linear mov, to rot mov.
    else %Even
        OutAzP(j) = (1/nEncAz)/R_az_e;
    end 
    k= k+1;
end

%% fix problem with alitude encoders
% luckily it's a linear model so I can just do a linear compensation to fix
% it
OutAltP = OutAltP - 0.288288 * OutAzP;

%% Obtain vector input camera force input and displacement output
InCamY= zeros (nNodos, 1); InCamY(20)=1; %1*nPOI + 2 
InCamZ= zeros (nNodos, 1); InCamZ(21)=1;
% OutCam = zeros(1, 2*nNodos); OutCam(19)=1; %3*nPOI + 1

%% Obtain vectors for line of sight outputs 

%Step 1 Put the 18 values in a column vector
Step1LOS = zeros(3*6, 2*nNodos);

%Index for M1M3
IndexFirstOutput = 1; %
IndexLastOutput = 6;  %1*nPOI 
k = 1;
for(j=IndexFirstOutput:IndexLastOutput)
    Step1LOS(k,j) = 1;
    k = k+1;
end

%Index for M2 and Camera
IndexFirstOutput = 13; %2*nPOI+1 
IndexLastOutput = 24;  %2*nPOI
for(j=IndexFirstOutput:IndexLastOutput)
    Step1LOS(k, j) = 1;
    k = k+1;
end

%Step 2. Rearrange the 18 values column vector in a 24 values vector and 
%scale them
Step2LOS = zeros(4*6, 3*6);
%M1
Step2LOS(1,1) = 1e6; %um
Step2LOS(2,2) = 1e6;
Step2LOS(3,3) = 1e6;
Step2LOS(4,4) = (180/pi)*3600; %arcsec
Step2LOS(5,5) = (180/pi)*3600;
Step2LOS(6,6) = (180/pi)*3600;
%M2
Step2LOS(7,7) = 1e6;
Step2LOS(8,8) = 1e6;
Step2LOS(9,9) = 1e6;
Step2LOS(10,10) = (180/pi)*3600;
Step2LOS(11,11) = (180/pi)*3600;
Step2LOS(12,12) = (180/pi)*3600;
%M3=M1
Step2LOS(13,1) = 1e6;
Step2LOS(14,2) = 1e6;
Step2LOS(15,3) = 1e6;
Step2LOS(16,4) = (180/pi)*3600;
Step2LOS(17,5) = (180/pi)*3600;
Step2LOS(18,6) = (180/pi)*3600;
%Camera
Step2LOS(19,13) = 1e6;
Step2LOS(20,14) = 1e6;
Step2LOS(21,15) = 1e6;
Step2LOS(22,16) = (180/pi)*3600;
Step2LOS(23,17) = (180/pi)*3600;
Step2LOS(24,18) = (180/pi)*3600;

%Step 3 Apply the cooordinate transformation to the 24 values 
Step3LOS = zeros(4*5, 4*6);

CoordTransfMatrix = zeros(5,6);
CoordTransfMatrix(1,1) = 1;
CoordTransfMatrix(4,4) = 1;
CoordTransfMatrix(2,2) = sin(AltAng);
CoordTransfMatrix(3,3) = sin(AltAng);
CoordTransfMatrix(5,5) = sin(AltAng);
%CoordTransfMatrix(6,6) = sin(AltAng);
CoordTransfMatrix(2,3) = -cos(AltAng);
CoordTransfMatrix(3,2) = cos(AltAng);
CoordTransfMatrix(5,6) = -cos(AltAng);
%CoordTransfMatrix(6,5) = cos(AltAng);

Step3LOS = [CoordTransfMatrix, zeros(5,6), zeros(5,6), zeros(5,6);
    zeros(5,6), CoordTransfMatrix, zeros(5,6), zeros(5,6);
     zeros(5,6), zeros(5,6), CoordTransfMatrix, zeros(5,6);
     zeros(5,6), zeros(5,6), zeros(5,6), CoordTransfMatrix];

%Step 4 put the displacements in z, x, y sequence 
Step4LOS = zeros(4*5, 4*5);
CoordSwapMatrix=[0 0 1 0 0;1 0 0 0 0;0 1 0 0 0;0 0 0 1 0;0 0 0 0 1];
Step4LOS=[CoordSwapMatrix, zeros(5), zeros(5), zeros(5);
    zeros(5), CoordSwapMatrix, zeros(5), zeros(5);
    zeros(5), zeros(5), CoordSwapMatrix, zeros(5);
    zeros(5), zeros(5), zeros(5), CoordSwapMatrix];

%Matrices LOS
load LofS.mat;

OutLOS = Line_of_Sight * Step4LOS * Step3LOS * Step2LOS * Step1LOS;

%% new LOS matrix

% Step2LOS * Step1LOS * SSoutput = [m1x m1y m1z m1rx m1ry m1rz  m2 ...  m3 ...  cam ...]' 
% in microns and arc sec

% new LOS matrix components: 
%     M1Motion = [m1m3Trans*1e6 m1m3Rot/pi*180*3600]; 
%     M2Motion = [m2Trans*1e6 m2Rot/pi*180*3600];
%     CamMotion = [cameraTrans*1e6 cameraRot/pi*180*3600];
%     vxCG = [M1Motion M2Motion CamMotion]';
%     % 2 x 18 translations in micro meters rotations in arc sec
%     LOS = [m1x m1y m1z m1rx m1ry m1rz m2 ... cam ...
    
load('LOSM_matrix45.mat');
LOSMStretch = [LOSM(:,1:12) zeros(2,6) LOSM(:,13:end)];
OutLOS = LOSMStretch * Step2LOS * Step1LOS;

%% Obtain Inputs and outputs for dampers
%Index for Altitud drives
IndexFirstInput = 259; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 2*nPAzEncoders + 2*(nP_WF_TE+nP_WF_TR+nP_WF_CS) + 1;
IndexLastInput = 288; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 2*nPAzEncoders + 2*(nP_WF_TE+nP_WF_TR+nP_WF_CS)+ 3*nP_WF_DN;

InDamperNodes= zeros (nNodos, 3*10); %3*nP_WF_DN
OutDamperNodes = zeros(3*10, 2*nNodos); %3*nP_WF_DN

k = 1;
for(j=IndexFirstInput:IndexLastInput)
     InDamperNodes(j,k) = 1;
     OutDamperNodes(k,j) = 1;
     k= k+1;
end

%% Obtain outputs for raw azimuth encoder measurements
%Index for Azimuth encoders
IndexFirstOutput = 171; %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 1
IndexLastOutput = 178;  %6*nPOI + 2*nPAltDrives + 2*nPAzDrives + 2*nPAltEncoders + 2*nPAzEncoders

OutAzPRaw = zeros(2*nEncAz, 2*nNodos);
k = 1;
for(j=IndexFirstOutput:IndexLastOutput)
        OutAzPRaw(k, j) = 1;
        k=k+1;
end

%%  Obtain Inputs for Fx drives
%Index for Azimuth drives
IndexFirstInput = 66; %6*nPOI + 2*nPAltDrives (no hay que sumarle el +1)

IndexDrivesFx = (IndexFirstInput-1) + [1 24 25 48]*2;
Index_beta_d = [1 24 25 48];

nDrivesFx = length(IndexDrivesFx);

InDrivesFx= zeros (nNodos, 1); %Sólo meto fuerza en 4 puntos

for(j=1:nDrivesFx)
     InDrivesFx(IndexDrivesFx(j)) = (1/nDrivesFx)/(-sin(beta_d(Index_beta_d(j))));
     InDrivesFx(IndexDrivesFx(j)+1) = (-1/nDrivesFx)/(-sin(beta_d(Index_beta_d(j))));
end


%%  Obtain Inputs for Fy drives
%Index for Azimuth drives
IndexFirstInput = 66; %6*nPOI + 2*nPAltDrives (no hay que sumarle el +1)

IndexDrivesFy = (IndexFirstInput-1) + [12 13 36 37]*2;
Index_beta_d = [12 13 36 37];

nDrivesFy = length(IndexDrivesFy);

InDrivesFy= zeros (nNodos, 1); %Sólo meto fuerza en 4 puntos

for(j=1:nDrivesFy)
     InDrivesFy(IndexDrivesFy(j)) = (1/nDrivesFy)/cos(beta_d(Index_beta_d(j)));
     InDrivesFy(IndexDrivesFy(j)+1) = (-1/nDrivesFy)/cos(beta_d(Index_beta_d(j)));
end


%% Generate SS with all inputs and outputs

%Cambio matrices para sacar entradas y salidas deseadas. Para estas matrices:
%Nx=2*nModos; Nu=4; Ny=2+2

In = [InAltT InAzT, InCamY, InCamZ, InDamperNodes, InDrivesFx, InDrivesFy];
Out = [OutAltP; OutAzP; OutLOS; OutDamperNodes; OutAzPRaw];

ATF = Ax;
BTF = Bx*In;
CTF = Out*Cx;
%DTF = Out*Dx*In;
DTF = 0;
sisTF = ss (ATF, BTF, CTF, DTF);

%% Define position and speed control controllers

fAltSpeedLoopBandwidth =  0.8*11.5;
fAzSpeedLoopBandwidth =  0.8*11;

%Speed series Controller
KpAlt = fAltSpeedLoopBandwidth * 2*pi * Jalt; 
KpAz =  fAzSpeedLoopBandwidth * 2*pi * Jaz;  
SpeedPController = ss([], [], [],[KpAlt 0; 0 KpAz]);

%Integral Action
TiAlt = 3/(fAltSpeedLoopBandwidth* 2*pi);
TiAz = 3/(fAzSpeedLoopBandwidth* 2*pi);

%Position series Controller
KvAlt = 0.25 * fAltSpeedLoopBandwidth * 2*pi; 
KvAz =  0.25 * fAzSpeedLoopBandwidth * 2*pi;  
PosPController = ss([], [], [],[KvAlt 0; 0 KvAz]);

%% Dampers

%Dampers in central section
MDampCS = 750;
KDampCS = MDampCS* (6*2*pi)^2; %1066e3;
CDampCS = 2*0.05*sqrt(KDampCS/MDampCS)*MDampCS; %2827;

%Dampers in hexapod
MDampHex = 13.3;
KDampHex = MDampHex* (11.13*2*pi)^2; %75.8e3;
CDampHex = 2*0.05*sqrt(KDampHex/MDampHex)*MDampHex; %100.5; %

%Active Damper
wnMode4TransX = (7*2*pi);
wnMode3TransY = (6*2*pi);

MDampAct = 312e3;
KDampActX = MDampAct * wnMode4TransX^2;
CDampActX = 20*2*0.05*wnMode4TransX*MDampAct; %El 20 es un valor calculado del rootlocus
KDampActY = MDampAct * wnMode3TransY^2;
CDampActY = 5*2*0.05*wnMode3TransY*MDampAct; %El 5 es un valor calculado del rootlocus

%BandPass Filter with lead term (former lead lag filter: tf([1/(2*5*pi) 1], [1/(2*40*pi)])
NumBPX = [0 2*1*wnMode4TransX 0];
DenBPX = [1 2*1*wnMode4TransX wnMode4TransX^2];
TauLead = 1/(2*5*pi);
NumBPLeadX = conv(NumBPX, [TauLead 1]);
DenBPLeadX = DenBPX;

% P = bodeoptions;
% P.FreqUnits = 'Hz';
% figure; bodeplot(tf(NumBPLeadX, DenBPLeadX),P); grid on;
