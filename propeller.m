%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [Xp] = propeller(beta, r, v , u, np)

    %np 螺旋桨转速
    %sigma 舵机转角

    %r 角速度
    %beta 漂角
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

    syms 
    L=7;
    tp = 0.220;
    A_r=0.0539;
    den=1000;
    Dp=0.216;
    n=0.216/0.345; 
    f_a=2.747; 
    wp0=0.40/0.35; 
    kt0=0.2931; kt1=-0.2753; kt2=-0.1385; 
    
    V=sqrt(u^2+v^2);
    beta_p = beta - (-0.5 * r * L / V);
%     if (beta_p >= 0)
%         C2 = 1.6;
%     else
        C2 = 1.1;
%     end
    C1 = 2.0;
    syms wp real positive
    eq = (1 - wp) / (1 - wp0) == 1 + (1 - exp(-C1 * abs(beta_p))) * (C2 - 1);
    sol = solve(eq, wp);
    wp_value = sol; 
    Va = u * (1 - wp_value);
    Jp = Va / (np * Dp);
    Kt = kt0 + kt1 * Jp + kt2 * Jp^2;
    T = den * np^2 * Dp^4 * Kt;
    Xp = (1 - tp) * T;
end



% 
% Xp=(1-tp)*T;
% T=den*np^2*Dp^4*Kt;
% Kt = kt0+kt1.*Jp+kt2.*Jp^2;
% Jp = Va./(np*Dp);
% (1-wp)/(1-wp0)=1+(1-exp(-C1*abs(beta_p)))*(C2-1);
% beta_p = beta-(-0.5*r*L/V);
% if beta_p>=0
%     C2=1.6;
% else 
%     C2=1.1;
% end
% C1=2.0;
% Va=u*(1-wp)