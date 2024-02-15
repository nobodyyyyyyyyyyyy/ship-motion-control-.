%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [X_r, Y_r, N_r] = rudder(sigma, np, u, r, v, beta)

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

L=7;
A_r=0.0539;
den=1000;
Dp=0.216;
n=0.216/0.345; 
f_a=2.747; 
wp0=0.40/0.35; 
kt0=0.2931; kt1=-0.2753; kt2=-0.1385; 

t_r = 0.387;
a_h = 0.312;
V = sqrt(v^2+u^2);
beta_p = beta - (-0.5 * r * L / V);
% if (beta_p >= 0)
%     C2 = 1.6;
% else
    C2 = 1.1;
% end
C1 = 2.0;
syms wp real positive
eq = (1 - wp) / (1 - wp0) == 1 + (1 - exp(-C1 * abs(beta_p))) * (C2 - 1);
sol = solve(eq, wp);
wp_value = sol;
Va = u * (1 - wp_value);
Jp = Va / (np * Dp);
Kt = kt0 + kt1 * Jp + kt2 * Jp^2;
U = sqrt(u^2 + v^2);
beta_r = beta + 0.71 * r * L / V;

% if (beta_r >= 0)
%     r_R = 0.64;
% else
    r_R = 0.395;
% end
u_r = 1.09 * u * (1 - wp_value) * sqrt(n * (1 + 0.5 * (sqrt(1 + 8 * Kt / (pi * Jp^2)) - 1))^2 + (1 - n));
v_r = U * r_R * beta_r;
U_r = sqrt(u_r^2 + v_r^2);
a_r = sigma - v_r / u_r;
F_n = 1/2 * den * f_a * A_r * U_r^2 * sind(a_r);
X_r = -(1 - t_r) * F_n * sind(sigma);
Y_r = -(1 + a_h) * F_n * cosd(sigma);
N_r = -(-0.5*L + a_h * (-0.464)*L)* F_n * cosd(sigma); 
end




% function [X_r, Y_r, N_r, U_r] = rudder(sigma, den, f_a, A_r, kt0, kt1, kt2, np, Dp, u, wp0, n, r, L, v, beta)
% 
% t_r = 0.387;
% a_h = 0.312;
% 
% V = dot(u, v);
% beta_p = beta - (-0.5 * r * L / V);
% if beta_p >= 0
%     C2 = 1.6;
% else
%     C2 = 1.1;
% end
% C1 = 2.0;
% syms wp real positive
% eq = (1 - wp) / (1 - wp0) == 1 + (1 - exp(-C1 * abs(beta_p))) * (C2 - 1);
% sol = solve(eq, wp);
% wp_value = double(sol); 
% Va = u * (1 - wp_value);
% Jp = Va / (np * Dp);
% Kt = kt0 + kt1 * Jp + kt2 * Jp^2;
% 
% F_n = 1/2 * den * f_a * A_r * U_r^2 * sin(a_r);
% 
% beta_r = beta + 0.71 * r * L / V;
% 
% if beta_r >= 0
%     r_R = 0.64;
% else
%     r_R = 0.395;
% end
% 
% X_r = -(1 - t_r) * F_n * sin(sigma);
% Y_r = -(1 + a_h) * F_n * cos(sigma);
% N_r = -(X_r + a_h * Y_r);
% 
% u_r = 1.09 * u * (1 - wp0) * sqrt(n * (1 + 0.5 * (sqrt(1 + 8 * Kt / pi * Jp^2) - 1))^2 + 1 - n);
% v_r = U * r_R * beta_r;
% 
% U = dot(u, v);
% U_r = dot(u_r, v_r);
% a_r = sigma - v_r / u_r;
% 
% end



