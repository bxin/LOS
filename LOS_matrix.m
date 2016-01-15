function LOSMnum = LOS_matrix(alphanum, betanum)

syms m1m3dx m1m3dy m1m3dz m1m3rx m1m3ry m1m3rz
syms m2dx m2dy m2dz m2rx m2ry m2rz
syms camdx camdy camdz camrx camry camrz
vxCG=transpose([m1m3dx m1m3dy m1m3dz m1m3rx m1m3ry m1m3rz ...
    m2dx m2dy m2dz m2rx m2ry m2rz ...
    camdx camdy camdz camrx camry camrz]);

vx = shift_CG2Vtx(vxCG,1);
[QM1, QM2, QM3, QCam] = make_4x4_TCRS(vx);

alpha=sym('alpha');
beta = sym('beta');

[LOSx, LOSy] = LineofSight(alpha, beta, QM1, QM2, QM3, QCam);

LOSM=sym('LOSM',[2, 18]);
for i=1:2
    if i==1
        bb = LOSx;
    elseif i==2
        bb = LOSy;
    end
    aa=(coeffs(vpa(bb),m1m3dx));LOSM(i,1)=aa(2);
    aa=(coeffs(vpa(bb),m1m3dy));LOSM(i,2)=aa(2);
    aa=(coeffs(vpa(bb),m1m3dz));LOSM(i,3)=aa(2);
    aa=(coeffs(vpa(bb),m1m3rx));LOSM(i,4)=aa(2);
    aa=(coeffs(vpa(bb),m1m3ry));LOSM(i,5)=aa(2);
    aa=(coeffs(vpa(bb),m1m3rz));LOSM(i,6)=aa(2);
    
    aa=(coeffs(vpa(bb),m2dx));LOSM(i,7)=aa(2);
    aa=(coeffs(vpa(bb),m2dy));LOSM(i,8)=aa(2);
    aa=(coeffs(vpa(bb),m2dz));if length(aa)<2, LOSM(i,9)=0; else LOSM(i,9)=aa(2);end
    aa=(coeffs(vpa(bb),m2rx));LOSM(i,10)=aa(2);
    aa=(coeffs(vpa(bb),m2ry));LOSM(i,11)=aa(2);
    aa=(coeffs(vpa(bb),m2rz));LOSM(i,12)=aa(2);
    
    aa=(coeffs(vpa(bb),camdx));LOSM(i,13)=aa(2);
    aa=(coeffs(vpa(bb),camdy));LOSM(i,14)=aa(2);
    aa=(coeffs(vpa(bb),camdz));if length(aa)<2, LOSM(i,15)=0; else LOSM(i,15)=aa(2);end
    aa=(coeffs(vpa(bb),camrx));LOSM(i,16)=aa(2);
    aa=(coeffs(vpa(bb),camry));LOSM(i,17)=aa(2);
    aa=(coeffs(vpa(bb),camrz));LOSM(i,18)=aa(2);
end

m1m3dx = 0;
m1m3dy = 0;
m1m3dz = 0;
m1m3rx = 0;
m1m3ry = 0;
m1m3rz = 0;

m2dx =0;
m2dy =0;
m2dz =0;
m2rx =0;
m2ry =0;
m2rz = 0;

camdx =0;
camdy =0;
camdz =0;
camrx =0;
camry =0;
camrz =0;

alpha = alphanum;
beta = betanum;

LOSMnum = subs(LOSM);

end
