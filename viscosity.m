%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [X_h,Y_h,N_h]=viscosity(u,v,r)

%r 角速度
%u y方向的速度
%v x方向的速度

%L 船长
%r0 静水阻力
%den 水的密度
%Xvv,Xvr,Xrr,Xvvvv 纵向水动力导数
%Nv,Nr,Nvvv,Nvvr,Nvrr,Nrrr 首向力矩导数
%Yv,Yr,Yvvv,Yvvr,Yvrr,Yrrr 横向水动力导数


L = 7;
den = 1000;
d = 0.46;
r0=0.0205;
Xvv=-0.04;Xvr=0.002;Xrr=0.011;Xvvvv=0.712;
Nv=-0.127;Nr=-0.049;Nvvv=-0.03;Nvvr=-0.294;Nvrr=0.049;Nrrr=-0.012;
Yv=-0.295;Yr=0.083;Yvvv=-1.682;Yvvr=0.379;Yvrr=-0.383;Yrrr=0.008;
%calculate the objective of vector velocity 

V=sqrt(u^2+v^2);
v=v./V;
r=r*L/V;
X_h = -r0+Xvv*v^2+Xvr*v*r+Xrr*r^2+Xvvvv*v^4;
Y_h = Yv*v+Yr*r+Yvvv*v^3+Yvvr*v^2*r+Yvrr*v*r^2+Yrrr*r^3;
N_h = Nv*v+Nr*r+Nvvv*v^3+Nvvr*v^2*r+Nvrr*v*r^2+Nrrr*r^3;
X_h=X_h.*(1/2*den*V^2*L*d);
Y_h=Y_h.*(1/2*den*V^2*L*d);
N_h=N_h.*(1/2*den*V^2*L^2*d);