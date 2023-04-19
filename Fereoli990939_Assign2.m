%% Modeling and Simulation of Aerospace Systems (2022/2023)
% Assignment:     2
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1
clear; clc; close all;

%Initial guesses
X0=zeros(4,1);
k0=1;
b0=1;
t0=0;
tf=10;

%Initial integration with k,b guessed
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode113(@(t,X) rhs1(t,X, k0, b0),[t0, tf], X0,options1);

%Plot initial integration with k,b guessed
figure(1);
plot(tt,XX(:,1),'r','LineWidth',2);
hold on;
grid on;
grid minor;
plot(tt,XX(:,2),'b','LineWidth',2);
ylabel('$\theta$ [rad]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
legend('$\theta_1$','$\theta_1$','Interpreter','latex','Location',...
    'northwest');

figure(2);
plot(tt,XX(:,3),'r','LineWidth',2);
hold on;
grid on;
grid minor;
plot(tt,XX(:,4),'b','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\omega$ [rad/s]','Interpreter','latex');
legend('$\omega_1$','$\omega_1$','Interpreter','latex');

%Plot real Y from measurements
Samp=readmatrix('samples.txt');
tt_samp=Samp(:,1);
omegadot1_samp=Samp(:,2);
omegadot2_samp=Samp(:,3);
YYreal=[omegadot1_samp, omegadot2_samp];
figure(3);
plot(tt_samp,YYreal(:,1),'r','LineWidth',2);
hold on;
plot(tt_samp,YYreal(:,2),'b','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$\dot{\omega}_{Real}$ [$rad/s^2$]','Interpreter','latex');
legend('$\dot{\omega}_1$', '$\dot{\omega}_2$','Interpreter','latex');

%Finding k,b: parameter estimation with OE method
options2=optimset('Display','Iter','TolX',1e-12,'TolFun',1e-12);
sol=fminunc(@(par) costOE1(par), [k0; b0], options2);
k_true=sol(1);
b_true=sol(2);

%Plot model Y
[ttm, XXm]=ode113(@(t,X) rhs1(t,X,k_true,b_true),[t0, tf], X0, options1);
YYm=zeros(length(ttm),2);
for i=1:length(ttm)
    YYm_full=rhs1(ttm(i), (XXm(i,:))', k_true, b_true);
    YYm(i,:)=YYm_full(3:4);
end
figure(4);
plot(ttm,YYm(:,1),'r','LineWidth',2);
hold on;
plot(ttm,YYm(:,2),'b','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$\dot{\omega}_{Fitted}$ [$rad/s^2$]','Interpreter','latex');
legend('$\dot{\omega}_1$', '$\dot{\omega}_2$','Interpreter','latex');

%% Ex.2
clear; clc; close all;

%Initial guesses
X0=[1; 0];
t0=0;
tf=10;

%Integration without V(t) source in the model
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode113(@(t,X) rhs2A(t,X),[t0, tf], X0, options1);

%Plot Vc(t) without V(t)
figure(1);
plot(tt,XX(:,1),'r','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('Vc [V]','Interpreter','latex');
grid on;
grid minor;

%Integration with V(t) source in the model
[tt2, XX2]=ode113(@(t,X) rhs2B(t,X),[t0, tf], X0, options1);

%Plot Vc(t) with V(t)
figure(2);
plot(tt2,XX2(:,1),'B','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('Vc [V]','Interpreter','latex');
grid on;
grid minor;

%% Ex.3
clear; clc; close;

%Initial guesses with single-lumping
X0=20*ones(2,1);
t0=0;
tf=60;

%Integration with single-lumping
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode113(@(t,X) rhs3A(t,X),[t0, tf], X0, options1);

%Plot T2 and T4 with single-lumping
figure(1);
plot(tt,XX(:,1),'r','LineWidth',2);
hold on;
plot(tt,XX(:,2),'b','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$T_i [^{\circ} C]$','Interpreter','latex');
grid on;
grid minor;
legend('$T_{2}$','$T_{4}$','Interpreter','latex','Location','southeast');

%Initial guesses with multiple-lumping
X0=20*ones(4,1);
t0=0;
tf=60;

%Integration with multiple-lumping
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX2]=ode113(@(t,X) rhs3B(t,X),[t0, tf], X0, options1);

%Plot T2A, T2B, T4A, T4B with multiple-lumping
figure(2);
plot(tt,XX2(:,1),'r','LineWidth',2);
hold on;
plot(tt,XX2(:,2),'b','LineWidth',2);
plot(tt,XX2(:,3),'g','lineWidth',2);
plot(tt,XX2(:,4),'k','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$T_{i} [^{\circ} C]$','Interpreter','latex');
grid on;
grid minor;
legend('$T_{2A}$','$T_{2B}$','$T_{4A}$','$T_{4B}$','Interpreter',...
    'latex','Location','southeast');

%% Ex.4
clear; clc; close;

%Data
Km=20;
R=200;
v0=2;
omega=5;
beta=0.2;
L=2*1e-3;
J1=0.5;
J2=0.3;
b=0.1;
k=0.5;

%Matrix A
A=[-R/L, 0, 0, -Km/L, 0;...
    0, 0, 0, 1, 0;...
    0, 0, 0, 0, 1;...
    Km/J1, -k/J1, k/J1, -b/J1, b/J1;...
    0, k/J2, -k/J2, b/J2, -b/J2];

%Eigenvalues of A
eigA=eig(A);

%Integration
X0=zeros(5,1);
t0=0;
tf=30;
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode15s(@(t,X) rhs4(t,X),[t0, tf], X0, options1);

%Plots
figure(1);
plot(tt,XX(:,1),'r','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('i [A]','Interpreter','latex');
grid on;
grid minor;

figure(2);
subplot(2,2,1);
plot(tt,XX(:,2),'r','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\theta_1$ [rad]','Interpreter','latex');
grid on;
grid minor;
subplot(2,2,2);
plot(tt,XX(:,3),'b','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\theta_2$ [rad]','Interpreter','latex');
grid on;
grid minor;
subplot(2,2,3);
plot(tt,XX(:,4),'g','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\omega_1$ [$rad/s$]','Interpreter','latex');
grid on;
grid minor;
subplot(2,2,4);
plot(tt,XX(:,5),'k','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\omega_2$ [$rad/s$]','Interpreter','latex');
grid on;
grid minor;

%Finding Km, R: parameter identification with OE method
options2=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
sol=fminunc(@(par) costOE2(par), [Km; R], options2);
Km_true=sol(1);
R_true=sol(2);

%Unpacking Y measurements
Samp=readmatrix('profile.txt');
tt_samp=Samp(:,1);
omegadot2_samp=Samp(:,2);

%Integration with estimated Km, R
options=odeset('RelTol',2.24*1e-12, 'AbsTol',2.24*1e-12);
[~, XX_true]=ode15s(@(t,X) rhs4par(t,X,Km_true, R_true),tt_samp, X0,...
    options);

%Comparison fitted Km,R model and measurements
figure(3);
plot(tt_samp,omegadot2_samp(:,1),'r-','LineWidth',2);
hold on;
plot(tt_samp,XX_true(:,5),'b--','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$err_{\omega_2}$ [$rad/s$]','Interpreter','latex');
legend('Measurement','Model Output','Interpreter','latex');

%% Ex.5
clear; clc; close;

%Initial guesses
X0=[1; 0; 0; 320; 320];
t0=0;
tf=25;

%Integration
options1=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode15s(@(t,X) rhs5(t,X),[t0, tf], X0, options1);

%Plots of state dynamics
figure(1);
plot(tt,XX(:,1),'r','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$V_t$ [$m^3$]','Interpreter','latex');

figure(2);
subplot(1,2,1);
plot(tt,XX(:,2),'r','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$x_k$ [m]','Interpreter','latex');
subplot(1,2,2);
plot(tt,XX(:,3),'b','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$v_k$ [$m/s$]','Interpreter','latex');

figure(3);
plot(tt,XX(:,4),'r','LineWidth',2);
hold on;
plot(tt,XX(:,5),'b','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$T_i$ [K]','Interpreter','latex');
legend('$T_h$','$T_f$','Interpreter',...
    'latex','Location','northwest');

%Other quantities: pressures and volume flow rates
disc=length(tt);
Pvec=zeros(10,disc);
Pkvec=zeros(1,disc);
Qvec=zeros(2,disc);
for i=1:disc
    [P,Pk,Q]=response5(tt(i),XX(i,:));
    Pvec(:,i)=P;
    Pkvec(i)=Pk;
    Qvec(:,i)=Q;
end

%Other quantities plot: P(t) and Q(t)
figure(4);
subplot(1,2,1);
plot(tt,Pvec./1e6,'LineWidth',2);
grid on;
grid minor;
box on;
xlabel('Time [s]','Interpreter','latex');
ylabel('$P_i$ [MPa]','Interpreter','latex');
legend('$P_T$','$P_1$','$P_2$','$P_3$','$P_4$','$P_5$','$P_6$',...
    '$P_7$','$P_8$','$P_9$','Interpreter','latex');
subplot(1,2,2);
plot(tt,Pkvec/1e6,'b','LineWidth',2);
grid on;
grid minor;
box on;
xlabel('Time [s]','Interpreter','latex');
ylabel('$P_k$ [MPa]','Interpreter','latex');

figure(5);
plot(tt,Qvec(1,:),'b','LineWidth',2);
hold on;
plot(tt,Qvec(2,:),'r','LineWidth',2);
grid on;
grid minor;
xlabel('Time [s]','Interpreter','latex');
ylabel('$Q_i$ [$m^3/s$]','Interpreter','latex');
legend('$Q_2$','$Q_7$','Interpreter','latex');


%% Functions

%XDOT Ex.1
function dF=rhs1(~,X,k,b)
%Data 
J1=0.2;  %Kg m 
J2=0.1; 
T=0.1;

%Unpacking state
theta1=X(1);
theta2=X(2);
omega1=X(3);
omega2=X(4);

%Composition ODEs (x=[theta1,theta2,omega1,omega2])
dF=zeros(4,1); 
dF(1)=omega1;
dF(2)=omega2;
dF(3)=k*(theta2-theta1)/J1;
dF(4)=-sign(omega2)*b*omega2^2/J2-k*(theta2-theta1)/J2+T/J2;
end

%OE-method for Ex.1
function J=costOE1(par)
%Unpacking
k=par(1);
b=par(2);

%Read samples.txt
Samp=readmatrix('samples.txt');
tt_samp=Samp(:,1);
omegadot1_samp=Samp(:,2);
omegadot2_samp=Samp(:,3);
YYreal=[omegadot1_samp, omegadot2_samp];

%Integration with k,b from fminunc
X0=zeros(4,1);
options=odeset('RelTol',2.24*1e-14, 'AbsTol',2.24*1e-14);
[tt, XX]=ode113(@(t,X) rhs1(t,X,k, b),tt_samp, X0, options);

%Model output
YYm=zeros(length(tt_samp),2);
for i=1:length(tt_samp)
    YYm_full=rhs1(tt(i), (XX(i,:))', k, b);
    YYm(i,:)=YYm_full(3:4);
end

%Composition J
R=cov(YYreal);
err=YYreal-YYm;
J=0;
for i=1:length(err)
J=J+0.5*err(i,:)*R^-1*err(i,:)';
end

end

%XDOT-VersA Ex.2
function dF=rhs2A(~,X)
%Unpacking
Vc=X(1);
i2=X(2);

%Data
R1=100;
R2k=10;
C=1e-3;
L=10;

%Composition ODEs (x=[Vc, i2])
dF=[i2/C;...
    (1+2*R2k*i2/R1)^-1*(-i2/(C*R1)-(Vc+R2k*i2*i2)/L)];
end

%XDOT-VersB Ex.2
function dF=rhs2B(t,X)
%Unpacking
Vc=X(1);
i2=X(2);

%Data
R1=100;
R2k=10;
C=1e-3;
L=10;
f=5;
V=sin(2*pi*f*t)*atan(t);
Vdot=sin(2*pi*f*t)/(t^2+1)+2*pi*f*atan(t)*cos(2*pi*f*t);

%Composition ODEs (x=[Vc, i2])
dF=[i2/C;...
    (1+2*R2k*i2/R1)^-1*(-i2/(C*R1)-(Vc+R2k*i2*i2)/L)+(1+2*R2k*i2/R1)^-1*(V/L+Vdot/R1)];
end

%XDOT-VersA Ex.3
function dF=rhs3A(t,X)
%Unpacking state
T2=X(1);
T4=X(2);

%Inputs
To0=20;
Ti0=20;
Tif=1000;
To=To0;
if t<1
    Ti=(Tif-Ti0)*t+Ti0;
else 
    Ti=Tif;
end

%Data
A=1;        %m^2
rho2=7530;  %kg/m^3
rho4=8170;  
cp2=377.1;  %J/kg
cp4=431; 
k1=163.3;   %Tungsten, W/mK
k2=62.8;    %Al-Bronze
k4=12;      %Inconel750
k5=11.7;    %NIMONIC75
l1=5*1e-3;  %m
l2=1e-2;
l4=1e-2;
l5=5*1e-3;
C2=A*l2*rho2*cp2;            
C4=A*l4*rho4*cp4;
Rint=0.001;  %Metal-metal contact thermal resistance

%Useful definitions
R1=(l1/(A*k1)+l2/(A*k2*2));
R2=(Rint+l2/(A*k2*2)+l4/(A*k4*2));
R3=(l5/(A*k5)+l4/(A*k4*2));

%Composition of ODEs (x=[T2, T4])
dF=[(Ti-T2)/(R1*C2)-(T2-T4)/(R2*C2);...
    (T2-T4)/(R2*C4)-(T4-To)/(R3*C4)];

end

%XDOT-VersA Ex.3
function dF=rhs3B(t,X)
%Unpacking
T2A=X(1);
T2B=X(2);
T4A=X(3);
T4B=X(4);

%Inputs
To0=20;
Ti0=20;
Tif=1000;
To=To0;
if t<1
    Ti=(Tif-Ti0)*t+Ti0;
else 
    Ti=Tif;
end

%Data
A=1;       %m^2
rho2=7530; %kg/m^3
rho4=8170;  
cp2=377.1;  %J/kg
cp4=431; 
k1=163.3;   %Tungsten, W/mK
k2=62.8;    %Al-Bronze
k4=12;      %Inconel750
k5=11.7;    %NIMONIC75
l1=5*1e-3;  %m
l2=1e-2;
l4=1e-2;
l5=5*1e-3;
C2=A*l2*rho2*cp2;           
C4=A*l4*rho4*cp4;
Rint=0.001;  %Metal-metal contact thermal resistance 

%Useful definitions
R1=(l1/(A*k1)+l2/(A*k2*4));
R2=l2/(A*k2*2);
R3=(Rint+l2/(A*k2*4)+l4/(A*k4*4));
R4=l4/(A*k4*2);
R5=(l5/(A*k5)+l4/(A*k4*4));

%Composition of ODEs (x=[T2A, T2B, T4A, T4B])
dF=[(Ti-T2A)/(R1*C2*0.5)-(T2A-T2B)/(R2*C2*0.5);...
    (T2A-T2B)/(R2*C2*0.5)-(T2B-T4A)/(R3*C2*0.5);...
    (T2B-T4A)/(R3*C4*0.5)-(T4A-T4B)/(R4*C4*0.5);...
    (T4A-T4B)/(R4*C4*0.5)-(T4B-To)/(R5*C4*0.5)];
end

%XDOT Ex.4
function dF=rhs4(t,X)
%Data
Km=20;
R=200;
v0=2;
omega=5;
beta=0.2;
L=2*1e-3;
J1=0.5;
J2=0.3;
b=0.1;
k=0.5;

%Input
v=v0*cos(omega*t)*exp(-beta*t);

%Matrix A
A=[-R/L, 0, 0, -Km/L, 0;...
    0, 0, 0, 1, 0;...
    0, 0, 0, 0, 1;...
    Km/J1, -k/J1, k/J1, -b/J1, b/J1;...
    0, k/J2, -k/J2, b/J2, -b/J2];

%Composition of ODEs (x=[i, theta1, theta2, omega1, omega2])
dF=A*X+[v/L; 0; 0; 0; 0];
end

%XDOT for OE-method Ex.4
function dF=rhs4par(t,X,Km,R)
%Data
v0=2;
omega=5;
beta=0.2;
L=2*1e-3;
J1=0.5;
J2=0.3;
b=0.1;
k=0.5;

%Input
v=v0*cos(omega*t)*exp(-beta*t);

%Matrix A
A=[-R/L, 0, 0, -Km/L, 0;...
    0, 0, 0, 1, 0;...
    0, 0, 0, 0, 1;...
    Km/J1, -k/J1, k/J1, -b/J1, b/J1;...
    0, k/J2, -k/J2, b/J2, -b/J2];

%Composition of ODEs (x=[i, theta1, theta2, omega1, omega2])
dF=A*X+[v/L; 0; 0; 0; 0];
end

%OE-method for Ex.4
function J=costOE2(par)
%Unpacking
Km=par(1);
R=par(2);

%Read profile.txt
Samp=readmatrix('profile.txt');
tt_samp=Samp(:,1);
omegadot2_samp=Samp(:,2);
YYreal=omegadot2_samp;

%Initial integration with Km,R of fminunc
X0=zeros(5,1);
options=odeset('RelTol',2.24*1e-12, 'AbsTol',2.24*1e-12);
[~, XX]=ode15s(@(t,X) rhs4par(t,X,Km, R),tt_samp, X0, options);

%Model output
YYm=XX(:,end);

%Composition J
R=cov(YYreal);
err=YYreal-YYm;
J=0;
for i=1:length(err)
J=J+0.5*err(i,:)*R^-1*err(i,:)';
end

end

%XDOT for Ex.5
function dF=rhs5(t,X)
%Unpacking
xk=X(2);
vk=X(3);
Th=X(4);
Tf=X(5);

%Data fluid
rho=1000;   %km/m^3
Cf=4186;    %J/kg

%Data pump regulator
Pnom=5*101325;
mk=2;  %kg
F0=5;  %N
rk=1; 
Dk=1e-2; 
Ak=pi*Dk^2/4;
thetamax=deg2rad(20);
Lc=1e-1;
c=Lc*tan(thetamax);
h=(Pnom*Ak-F0)/c;
dk=1e-3;
Apr=pi*dk^2/4;
vp=vk*Ak/Apr;
Kp=2.5;
dP2k=0.5*Kp*rho*vp*abs(vp);   %dP between 2-K

%Data pump motor and Q1
n=66.66; %RPS
Np=9;       %Number of Pistons
Dpm=0.7*1e-2;
Apm=pi*Dpm^2/4;
dss=1.5*1e-2;
s=@(x) dss*(c-x)/Lc;
Q1=n*Np*Apm*s(xk);

%Filter 
etaf=2.5/100;   %Leakeage
Q7=(1-etaf)*Q1;

%Data pipes
Dp=20*1e-3;
Ap=pi*Dp^2/4;

%Distributor state
do=1e-2;
if t<2
    alpha=pi*t/2+pi;
    Aalpha=(do^2/8)*(alpha-sin(alpha));
else
    Aalpha=(do^2/8)*(pi*2-sin(pi*2));
end

%Data pressure losses
Kcv=2;
f34=0.032;
L34=1.5;
Kd=15;
f56=0.040;
L56=2.7;
Kf=35;                                                                      
f78=0.028;
L78=2.5;
f9t=0.032;
L9t=1;

%Pressure losses: concentrated and distributed
dP23=0.5*Kcv*rho*(Q1/Ap)*abs(Q1/Ap);
dP34=0.5*f34*(L34/Dp)*rho*(Q1/Ap)*abs(Q1/Ap); 
dP45=0.5*Kd*rho*(Q1/Aalpha)*abs(Q1/Aalpha);
dP56=0.5*f56*(L56/Dp)*rho*(Q1/Ap)*abs(Q1/Ap);
dP67=0.5*Kf*rho*(Q1/Ap)*abs(Q1/Ap);
dP78=0.5*f78*(L78/Dp)*rho*(Q7/Ap)*abs(Q7/Ap);
dP9t=0.5*f9t*(L9t/Dp)*rho*(Q7/Ap)*abs(Q7/Ap);

%Resistances Heat Exchanger
A=0.1;
t1=1e-2;        %R1
k1=395;
R1=t1/(k1*A);
t2=2.5*1e-2;    %R2
k2=310;
R2=t2/(k2*A*2);
t3=1e-2;        %R3
k3=125;
R3=t3/(k3*A);                                                           
h_conv=20;      %R4
R4=1/(h_conv*A);
Req1=R1+R2;
Req2=R2+R3+R4;

%Capacitance heat exchanger and fluid
rho_he=8620;
Cp_he=100;
M4C4=rho_he*t2*Cp_he*A;                
Le=0.5;
De=20*1e-3;
MfCf=rho*Cf*pi*Le*De^2/4;                                                 

%Data and other quantities heat exchanger
T0=350;                                                                     
kt=20;
omega=5;
Text=T0+kt*cos(omega*t);
Qdot2=(Th-Tf)/Req2;                                                        
k_ray=1000;

%Pressure relationships and computations
PT=0.1*1e6;
P9=PT+dP9t;
P8=P9/exp(-Qdot2/k_ray);
P7=P8+dP78;
P6=P7+dP67;
P5=P6+dP56;
P4=P5+dP45;
P3=P4+dP34;
P2=P3+dP23;
Pk=P2-dP2k;

%Composition ODEs (x=[Vt, xk, vk, Th, Tf])
dF=[-Q1+Q7;...
    vk;...
    (Pk*Ak-F0-h*xk-rk*vk)/mk;...
    (Text-Th)/(M4C4*Req1)-(Th-Tf)/(M4C4*Req2);...
    (Th-Tf)/(MfCf*Req2)];

end

%FULL RESPONSE for Ex.5
function [P,Pk,Q]=response5(t,X)
%Unpacking
xk=X(2);
vk=X(3);
Th=X(4);
Tf=X(5);

%Data Fluid
rho=1000; %kg/m^3

%Data pump regulator
Dk=1e-2;
Ak=pi*Dk^2/4;
thetamax=deg2rad(20);
Lc=1e-1;
c=Lc*tan(thetamax);
dk=1e-3;
Apr=pi*dk^2/4;
vp=vk*Ak/Apr;
Kp=2.5;
dP2k=0.5*Kp*rho*vp*abs(vp);

%Data pump motor and Q1
n=66.66;  %RPS
Np=9;
Dpm=0.7*1e-2;
Apm=pi*Dpm^2/4;
dss=1.5*1e-2;
s=@(x) dss*(c-x)/Lc;
Q1=n*Np*Apm*s(xk);

%Filter
etaf=2.5/100;   %Leakeage
Q7=(1-etaf)*Q1;

%Data pipes
Dp=20*1e-3;
Ap=pi*Dp^2/4;

%Distributor state
do=1e-2;
if t<2
    alpha=pi*t/2+pi;
    Aalpha=(do^2/8)*(alpha-sin(alpha));
else
    Aalpha=(do^2/8)*(pi*2-sin(pi*2));
end

%Data pressure losses
ft1=0.032;
Lt1=0.5;
Kcv=2;
f34=0.032;
L34=1.5;
Kd=15;
f56=0.040;
L56=2.7;
Kf=35;                                                                      
f78=0.028;
L78=2.5;
f9t=0.032;
L9t=1;

%Pressure losses: concentrated and distributed
dPt1=0.5*ft1*(Lt1/Dp)*rho*(Q1/Ap)*abs(Q1/Ap);
dP23=0.5*Kcv*rho*(Q1/Ap)*abs(Q1/Ap);
dP34=0.5*f34*(L34/Dp)*rho*(Q1/Ap)*abs(Q1/Ap); 
dP45=0.5*Kd*rho*(Q1/Aalpha)*abs(Q1/Aalpha);
dP56=0.5*f56*(L56/Dp)*rho*(Q1/Ap)*abs(Q1/Ap);
dP67=0.5*Kf*rho*(Q1/Ap)*abs(Q1/Ap);
dP78=0.5*f78*(L78/Dp)*rho*(Q7/Ap)*abs(Q7/Ap);
dP9t=0.5*f9t*(L9t/Dp)*rho*(Q7/Ap)*abs(Q7/Ap);

%Resistances Heat Exchanger
A=0.1;
t2=2.5*1e-2;     %R2
k2=310;
R2=t2/(k2*A*2);
t3=1e-2;         %R3
k3=125;
R3=t3/(k3*A);                                                           
h_conv=20;       %R4
R4=1/(h_conv*A);
Req2=R2+R3+R4;

%Data and other quantities heat exchanger
Qdot2=(Th-Tf)/Req2;                                                           
k_ray=1000; 

%Pressure relationships and computations
PT=0.1*1e6;
P9=PT+dP9t;
P8=P9/exp(-Qdot2/k_ray);
P7=P8+dP78;
P6=P7+dP67;
P5=P6+dP56;
P4=P5+dP45;
P3=P4+dP34;
P2=P3+dP23;
Pk=P2-dP2k;
P1=PT-dPt1;

%Composition of P and Q vectors
P=[PT; P1; P2; P3; P4; P5; P6; P7; P8; P9];
Q=[Q1; Q7];

end

