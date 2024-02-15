%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
%np 螺旋桨转速
%sigma 舵机转角

%r 角速度
%beta 漂角
%r 船的角速度
%u y方向的速度
%v x方向的速度

%L 船长
%A_r 舵面积
%den 水的密度
%Dp 螺旋桨直径
%n 桨径和舵高的比值
%f_a 舵的升力梯度系数
%wp0 直航时的有效伴流分数
%kt0, kt1, kt2 泰勒展开系数

%m 船体质量大小
%I_z 船舶对ox轴的转动惯量
%mx 船在x方向上所增加的液体质量大小
%my 船在y方向上所增加的液体质量大小
%J_zz 船在转动方向上所增加的液体质量大小

%r0 静水阻力
%Xvv,Xvr,Xrr,Xvvvv 纵向水动力导数
%Nv,Nr,Nvvv,Nvvr,Nvrr,Nrrr 首向力矩导数
%Yv,Yr,Yvvv,Yvvr,Yvrr,Yrrr 横向水动力导数
clc,clear

L=7;
x_g = 0.25;
A_r=0.0539;
den=1000;
Dp=0.216;
n=0.216/0.345; 
f_a=2.747; 
wp0=0.40/0.35; 
kt0=0.2931; kt1=-0.2753; kt2=-0.1385; 
m = 3270; 
I_z=10008.25;
m_x=0.022;
m_y=0.223;
J_zz=0.011;
r0=0.0205;
Xvv=-0.04;Xvr=0.002;Xrr=0.011;Xvvvv=0.712;
Nv=-0.127;Nr=-0.049;Nvvv=-0.03;Nvvr=-0.294;Nvrr=0.049;Nrrr=-0.012;
Yv=-0.295;Yr=0.083;Yvvv=-1.682;Yvvr=0.379;Yvrr=-0.383;Yrrr=0.008;

syms  u_ v_ r_ np sigma u r v beta;
beta = (atan(-v/u));

% [X_h,Y_h,N_h]=viscosity(u,v,r,r0,den,L,Xvv,Xvr,Xrr,Xvvvv,Yv,Yr,Yvvv,Yvvr,Yvrr,Yrrr,Nv,Nr,Nvvv,Nvvr,Nvrr,Nrrr);
% [X_r, Y_r, N_r] = rudder(sigma, den, f_a, A_r, kt0, kt1, kt2, np, Dp, u, wp0, n, r, L, v, beta);
% [Xp] = propeller(tp, wp0, kt0, kt1, kt2, beta, r, L, V, u, np, Dp, den);
% [X_i, Y_i, N_i] = inertia(m,I_z,m_x,m_y,J_zz,u,v,r,u_,v_,r_,x_g);

[X_h,Y_h,N_h]=viscosity(u,v,r);
[X_r, Y_r, N_r] = rudder(sigma, np, u, r, v, beta);
[Xp] = propeller(beta, r, v , u, np);
[X_i, Y_i, N_i] = inertia(u,v,r,u_,v_,r_);

equation1 = X_i == X_h+X_r+Xp;
equation2 = Y_i == Y_h+Y_r;
equation3 = N_i == N_h+N_r;
global result
result = solve([equation1,equation2,equation3],[u_,v_,r_]);




