clear;clc;close all;
syms q1 q2 q3 q4 qd1 qd2 qd3 qd4 qdd1 qdd2 qdd3 qdd4;
%% system information
m1=1.5; Pr=0.2; I1=0.5*(m1)*Pr^2; m2=0.5; m3=1; L=0.4; I3=(1/12)*m3*L^2; k1=75; d1=1; k2=50; d2=0.4;
cmpk1=-0.2; cmpk2=(m2+m3)*9.81/k2;
xG2=q1+q2*sin(q3); yG2=Pr-q2*cos(q3);
xdG2=qd1+qd2*sin(q3)+q2*qd3*cos(q3);  ydG2=-qd2*cos(q3)+q2*qd3*sin(q3);
xG3=q1+q2*sin(q3)+L/2*sin(q4); yG3=Pr-q2*cos(q3)-L/2*cos(q4);
xdG3=qd1+qd2*sin(q3)+q2*qd3*cos(q3)+L/2*qd4*cos(q4);
ydG3=-qd2*cos(q3)+q2*qd3*sin(q3)+L/2*qd4*sin(q4);
%% energy
T1=0.5*m1*(qd1)^2+0.5*I1*(qd1/Pr)^2; 
T2=0.5*m2*(xdG2^2+ydG2^2);
T3=0.5*m3*(xdG3^2+ydG3^2)+0.5*I3*qd4^2;
T=T1+T2+T3;
Vs=0.5*k1*(q1-d1)^2+0.5*k2*(q2-d2)^2;
Vg=m1*9.81*Pr+m2*9.81*yG2+m3*9.81*yG3;
V=Vs+Vg;
%% Lagrangian
L=T-V;

% %% ?????? A? ????? dL/dq ??????

% A(1,:)=jacobian(L,x);
% A(2,:)=jacobian(L,r);
% A(3,:)=jacobian(L,theta2);
% A(4,:)=jacobian(L,theta3);

%??????B? dL/dqdot ??????

% B(1,:)=jacobian(L,xd);
% B(2,:)=jacobian(L,rd);
% B(3,:)=jacobian(L,theta2d);
% B(4,:)=jacobian(L,theta3d);

%?????? C? ?? ????? ???? ?????? ?? d/dt(dL/dqdo)=dB/dqdot*qdoubledot ?????? 

% C(1,:)=jacobian(B(1),[x,xd,r,rd,theta2,theta2d,theta3,theta3d].')*[xd,xdd,rd,rdd,theta2d,theta2dd,theta3d,theta3dd].';
% C(2,:)=jacobian(B(2),[x,xd,r,rd,theta2,theta2d,theta3,theta3d].')*[xd,xdd,rd,rdd,theta2d,theta2dd,theta3d,theta3dd].';
% C(3,:)=jacobian(B(3),[x,xd,r,rd,theta2,theta2d,theta3,theta3d].')*[xd,xdd,rd,rdd,theta2d,theta2dd,theta3d,theta3dd].';
% C(4,:)=jacobian(B(4),[x,xd,r,rd,theta2,theta2d,theta3,theta3d].')*[xd,xdd,rd,rdd,theta2d,theta2dd,theta3d,theta3dd].';
% 
% D=C-A;
%% decoupling equations
dL_dq(1,:)=simplify(transpose(jacobian(L,q1)));
dL_dq(2,:)=simplify(transpose(jacobian(L,q2)));
dL_dq(3,:)=simplify(transpose(jacobian(L,q3)));
dL_dq(4,:)=simplify(transpose(jacobian(L,q4)));

dL_dqd(1,:)=simplify(transpose(jacobian(L,qd1)));
dL_dqd(2,:)=simplify(transpose(jacobian(L,qd2)));
dL_dqd(3,:)=simplify(transpose(jacobian(L,qd3)));
dL_dqd(4,:)=simplify(transpose(jacobian(L,qd4)));

A(1,:)=simplify(jacobian(dL_dqd,qd1));
A(2,:)=simplify(jacobian(dL_dqd,qd2));
A(3,:)=simplify(jacobian(dL_dqd,qd3));
A(4,:)=simplify(jacobian(dL_dqd,qd4));

B(1,:)=simplify(jacobian(dL_dqd,q1));
B(2,:)=simplify(jacobian(dL_dqd,q2));
B(3,:)=simplify(jacobian(dL_dqd,q3));
B(4,:)=simplify(jacobian(dL_dqd,q4));


%% solving equations and plotting result
x0=d1+cmpk1; xd0=0; r0=d2+cmpk2; rd0=0; theta2_0=0; theta2d_0=0; theta3_0=0; theta3d_0=0;
ICs=[x0,r0,theta2_0,theta3_0,xd0,rd0,theta2d_0,theta3d_0];
tspan=[0,10];
[t,x]=ode45(@decoupled_equation, tspan, ICs);

plot(t,x(:,1),'g','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('x (m)')
figure
plot(t,x(:,2),'b','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('r (m)')
figure
plot(t,x(:,3),'k','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('\theta_2 (rad)')
figure
plot(t,x(:,4),'m','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('\theta_3 (rad)')
figure

plot(t,x(:,5),'g','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('xdot (m/s)')
figure
plot(t,x(:,6),'b','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('rdot (m/s)')
figure
plot(t,x(:,7),'k','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('\theta_2 dot (rad/s)')
figure
plot(t,x(:,8),'m','linewidth',2);
grid on,
xlabel('Time (s)'),ylabel('\theta_3 dot (rad/s)')
