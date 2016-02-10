function [SxFinal, SyFinal] = make_LOS_senM()

% author: Bo Xin (bxin@lsst.org)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

%% construct sensitivity matrix that can be multiplied with Qfinal
%load senM, unit is um/um or um/arcsec
load data/RB_sensitivities

% Calculate conversion from annular Zernike tilt (RMS um) to OPD tilt (arcsec)

R_CA_M1 = 4.18;             % meter
e = 0.61;
c = -(2e-6/sqrt(1+e^2))/R_CA_M1*180/pi*3600; %minus sign? check s2 s3, make physical sense of s2(5) & s3(4)

s23 = squeeze(senM(1,2:3,1:20)); %tip and tilt

Sxy_optics = zeros(4,4,4,2);
for k=1:2
    for i=1:4
        Sxy_optics(:,:,i,k)=[0 0  s23(k,(i-1)*5+5)/pi*180*3600*c s23(k,(i-1)*5+2)*c;
                             0 0 -s23(k,(i-1)*5+4)/pi*180*3600*c s23(k,(i-1)*5+3)*c;
                             0 0          0        s23(k,(i-1)*5+1)*c;
                             0 0          0            0];
    end
    if k==1
        SxFinal = blkdiag(squeeze(Sxy_optics(:,:,1,k)),squeeze(Sxy_optics(:,:,2,k)), ...
            squeeze(Sxy_optics(:,:,3,k)),squeeze(Sxy_optics(:,:,4,k)));
    else
        SyFinal = blkdiag(squeeze(Sxy_optics(:,:,1,k)),squeeze(Sxy_optics(:,:,2,k)), ...
            squeeze(Sxy_optics(:,:,3,k)),squeeze(Sxy_optics(:,:,4,k)));        
    end
    
end

end
