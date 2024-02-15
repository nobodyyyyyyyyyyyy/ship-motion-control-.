%% Hydrodynamic equations with Euler approximation.
function [x, y, ut, vt, rt] = hyddyn2(u,v,r,th)
%%Initialize.
%%

dt = 0.01;
t = 0:dt:100;
len = length(t);
ut = zeros(1,len);
vt = zeros(1,len);
rt = zeros(1,len);
th_t = zeros(1,len);
vxt = zeros(1,len);
vyt = zeros(1,len);

x = zeros(1,len);
y = zeros(1,len);
x(1) = 1;
y(1) = 1;

u = 1;
v = 1;
r = 0.03;
th = 0.05;

ax = 0.01;
ay = 0.01;

ut(1) = u;
vt(1) = v;
rt(1) = r;
th_t(1) = th;

vx = u*cos(th) + v*sin(th);
vy = v*cos(th) + u*sin(th);
vxt(1) = vx;
vyt(1) = vy;
for i = 1:length(t)-1

    if(mod(th_t(i),pi/4==0))
    th = th_t(i)+0.001;
    else
    th = th_t(i);
    end
A = [cos(th) sin(th);sin(th) cos(th)];
B = [ax*dt + ut(i)*cos(th) + vt(i)*sin(th) + rt(i)*dt*(vt(i)*cos(th)-ut(i)*sin(th)); ay*dt + rt(i)*dt*(ut(i)*cos(th)-vt(i)*sin(th)) + ut(i)*sin(th) + vt(i)*cos(th)];

X = linsolve(A,B);
ut(i+1) = X(1);
vt(i+1) = X(2);
vxt(i+1) = ut(i+1)*cos(th) + vt(i+1)*sin(th);
vyt(i+1) = vt(i+1)*cos(th) + ut(i+1)*sin(th);

ax = (vxt(i+1)-vxt(i))/dt;
ay = (vyt(i+1)-vyt(i))/dt;
% u_dot = 
% ax = 
th_t(i+1) = atan(vyt(i+1)/vxt(i+1));

rt(i+1) = (th_t(i+1)-th_t(i))/dt;
%th_t(i+1) = th_t(i) + r*dt;
x(i+1) = x(i) + vxt(i+1)*dt;
y(i+1) = y(i) + vyt(i+1)*dt;
end

figure()
subplot(2,2,1)
plot(t,ut)
xlabel('time, t')
ylabel('u(t)')
subplot(2,2,2)
plot(t,vt)
xlabel('time, t')
ylabel('v(t)')
subplot(2,2,3)
plot(t,th_t)
xlabel('time, t')
ylabel('\theta(t)')
subplot(2,2,4)
plot(x,y)
xlabel('x')
ylabel('y')
title('Trajectory')
pu = polyfit(t,ut,2);
pv = polyfit(t,vt,2);
pr = polyfit(t,rt,2);
end