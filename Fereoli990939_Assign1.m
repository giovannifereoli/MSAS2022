%% Modeling and Simulation of Aerospace Systems (2022/2023)
% Assignment:     1
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1
clear; clc; close;

%Initialization
fun=@(x) sin(x)+x-1;
a=0; b=1;
n=8;
tol=10^-(n+1);

%Roots finding solutions and CPU time
[solbis, fevalbis]=bisection(fun, a, b, tol);
timebis=timeit(@() bisection(fun, a, b, tol));
[solsec, fevalsec]=secant(fun, a, b, tol); 
timesec=timeit(@() secant(fun, a, b, tol));
[solregf, fevalregf]=regulafalsi(fun,a,b,tol);
timeregf=timeit(@() regulafalsi(fun, a, b, tol));

%Plot function
xvec=-3:0.1:3;
fvec=fun(xvec);
figure(1);
plot(xvec,fvec,'b','LineWidth',2);
hold on;
plot(solbis,0,'r.','MarkerSize',20);
text(solbis,0,'Solution ','HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
plot(a, fun(a),'k.','MarkerSize',20);
plot(b, fun(b),'k.','MarkerSize',20);
text(a, fun(a),'f(a)<0   ','HorizontalAlignment','right');
text(b, fun(b),'f(b)>0   ','HorizontalAlignment', 'right');
xL=xlim;
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('x [-]');
ylabel('f [-]');
ylim([-3,3]);
grid on;

%Plots results
figure(2);
subplot(1,3,1);
scatter([1;2;3],[solbis,solsec,solregf],'filled');
xticks([1,2,3]);
xticklabels({'Bisection', 'Secant','Regula-Falsi'});
title('Solutions');
grid on;
subplot(1,3,2);
semilogy([1;2;3],[timebis; timesec; timeregf],'b-*','LineWidth',1);
xticks([1,2,3]);
xticklabels({'Bisection', 'Secant','Regula-Falsi'});
title('CPU Time');
grid on;
ylabel('CPU time [s]');
subplot(1,3,3);
plot([1;2;3],[fevalbis; fevalsec; fevalregf],'r-*','LineWidth',1);
xticks([1,2,3]);
xticklabels({'Bisection', 'Secant','Regula-Falsi'});
title('Number of Iterations');
grid on;
ylabel('FunEval [-]');
sgtitle('Numerical solutions f(x)=0');

%% Ex.2
clear; clc; close;

%Function, Jacobians definition (Exact and Finite Differences)
F=@(x) [x(1).^2+x(2)-5; x(2).^2-x(1)];                                      
JFexact=@(x) [2*x(1), 1; -1, 2*x(2)];
JFapprox=@(x) [(F(x+max(sqrt(eps),sqrt(eps)*norm(x))*[1;0])-...
        F(x))/max(sqrt(eps),sqrt(eps)*norm(x)), (F(x+max(sqrt(eps),...
        sqrt(eps)*norm(x))*[0;1])-F(x))/max(sqrt(eps),sqrt(eps)*norm(x))];

%Root-1 finding solution
tol=1e-12;
x0=[2;2];
[solex, itex]=newton(@(x) F(x), @(x) JFexact(x), x0, tol);
[solapp, itapp]=newton(@(x) F(x), @(x) JFapprox(x), x0, tol);

%Root-2 finding solution
tol=1e-12;
x02=[2;-2];
[solex2, itex2]=newton(@(x) F(x), @(x) JFexact(x), x02, tol);
[solapp2, itapp2]=newton(@(x) F(x), @(x) JFapprox(x), x02, tol);

%Plot, guesses
disc=100;
[X1,X2]=meshgrid(linspace(-5,+5,disc),linspace(-5,+5,disc));
Z=zeros(disc,disc);
for i=1:disc
    for j=1:disc
        Z(i,j)=norm(F([X1(i,j);X2(i,j)]));
    end
end
contourf(X1, X2, Z, 0:1:30);
contourcbar;
hold on;
plot(solapp(1), solapp(2),'r.','MarkerSize',15);
plot(solapp2(1), solapp2(2),'g.','MarkerSize',15);
legend('Surface Contour','First Solution','Second Solution');
xticks(-5:5);
xlabel('x [-]');
ylabel('y [-]');
zlabel('z [-]');

%Accuracy analysis
%Root-1
err_ex=norm(F(solex));
err_approx=norm(F(solapp));
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n']);
fprintf(['Roots finding error for FIRST solution Newton\''s'...
    ' method with exact Jacobian is:\n abs(f(x*))=%d in %d'...
    ' iterations.\n'],err_ex,itex);
fprintf(['Roots finding error for FIRST solutions Newton\''s'...
    ' method with approximated Jacobian is:\n abs(f(x*))=%d in %d'...
    ' iterations.\n'],err_approx,itapp);
%Root-2
err_ex2=norm(F(solex2));
err_approx2=norm(F(solapp2));
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n']);
fprintf(['Roots finding error for SECOND solution Newton\''s'...
    ' method with exact Jacobian is:\n abs(f(x*))=%d in %d'...
    ' iterations.\n'],err_ex2,itex2);
fprintf(['Roots finding error for SECOND solutions Newton\''s'...
    ' method with approximated Jacobian is:\n abs(f(x*))=%d in %d'...
    ' iterations.\n'],err_approx2,itapp2);
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n']);

%% Ex.3
clear; clc; close;

%Function definitions
xdot=@(x,t) -x*(t^2-1)/(t^2+1);                                             
xex=@(t) exp(2*atan(t)-t);
x0=1;

%Initialization
t0=0;                                                                      
tf=2;
h=[0.5; 0.2; 0.05; 0.01];

%Heuns integration and CPU time
[tspan1,xheun1]=RK2(@(x,t) xdot(x,t), t0, tf, h(1), x0);
theun1=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(1), x0));              
[tspan2,xheun2]=RK2(@(x,t) xdot(x,t), t0, tf, h(2), x0);
theun2=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(2), x0));  
[tspan3,xheun3]=RK2(@(x,t) xdot(x,t), t0, tf, h(3), x0);
theun3=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(3), x0));  
[tspan4,xheun4]=RK2(@(x,t) xdot(x,t), t0, tf, h(4), x0);
theun4=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(4), x0));  

%Heuns errors
err_heun1=abs(xheun1-xex(tspan1));                                 
err_heun2=abs(xheun2-xex(tspan2));                                  
err_heun3=abs(xheun3-xex(tspan3));                                  
err_heun4=abs(xheun4-xex(tspan4));    

%RK4 integration and CPU time
[~,xrk41]=RK4(@(x,t) xdot(x,t), t0, tf, h(1), x0);
trk41=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(1), x0));               
[~,xrk42]=RK4(@(x,t) xdot(x,t), t0, tf, h(2), x0);
trk42=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(2), x0));  
[~,xrk43]=RK4(@(x,t) xdot(x,t), t0, tf, h(3), x0);
trk43=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(3), x0));  
[~,xrk44]=RK4(@(x,t) xdot(x,t), t0, tf, h(4), x0);
trk44=timeit(@() RK2(@(x,t) xdot(x,t), t0, tf, h(4), x0));  

%Heuns errors
err_rk41=abs(xrk41-xex(tspan1));                                 
err_rk42=abs(xrk42-xex(tspan2));                                  
err_rk43=abs(xrk43-xex(tspan3));                                  
err_rk44=abs(xrk44-xex(tspan4)); 

%Integrators plot
figure(1);
subplot(2,2,1);
fplot(xex,'b','LineWidth',1);
hold on;
plot(tspan1,xheun1,'r-*','LineWidth',1);
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
xlim([t0, tf]);
title('h_1=0.5');
legend('Exact Solution','Location','southeast')
subplot(2,2,2);
plot(tspan2,xheun2,'g-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
xlim([t0, tf]);
title('h_2=0.2');
subplot(2,2,3);
plot(tspan3,xheun3,'c-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
xlim([t0, tf]);
title('h_3=0.05');
subplot(2,2,4);
plot(tspan4,xheun4,'m-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
xlim([t0, tf]);
title('h_4=0.01');

figure(2);
subplot(2,2,1);
fplot(xex,'b','LineWidth',1);
hold on;
plot(tspan1,xrk41,'r-*','LineWidth',1);
legend('Exact Solution','Location','southeast')
xlim([t0, tf]);
title('h_1=0.5');
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
subplot(2,2,2);
plot(tspan2,xrk42,'g-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlim([t0, tf]);
title('h_2=0.2');
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
subplot(2,2,3);
plot(tspan3,xrk43,'c-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlim([t0, tf]);
title('h_3=0.05');
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;
subplot(2,2,4);
plot(tspan4,xrk44,'m-*','LineWidth',1);
hold on;
fplot(xex,'b','LineWidth',1);
xlim([t0, tf]);
title('h_4=0.01');
xlabel('t [s]');
ylabel('x [-]');
ylim([1,1.8]);
grid on;

%Errors plot
figure(3);
subplot(1,2,1);
semilogy(tspan1,err_heun1,'r-*','LineWidth',1);
hold on;
semilogy(tspan2,err_heun2,'g-*','LineWidth',1);
semilogy(tspan3,err_heun3,'y-*','LineWidth',1);
semilogy(tspan4,err_heun4,'m-*','LineWidth',1);
xlim([t0, tf]);
legend('h_1=0.5', 'h_2=0.2','h_3=0.05', 'h_4=0.01',...
    'Location','southeast');
title('RK2')
xlabel('t [s]');
ylabel('Err [-]');
grid on;
subplot(1,2,2)
semilogy(tspan1,err_rk41,'r-*','LineWidth',1);
hold on;
semilogy(tspan2,err_rk42,'g-*','LineWidth',1);
semilogy(tspan3,err_rk43,'y-*','LineWidth',1);
semilogy(tspan4,err_rk44,'m-*','LineWidth',1);
xlim([t0, tf]);
legend('h_1=0.5', 'h_2=0.2','h_3=0.05', 'h_4=0.01',...
    'Location','southeast');
title('RK4')
xlabel('t [s]');
ylabel('x [-]');
grid on;

%CPU Time plots
figure(4);
subplot(1,2,1);
semilogy([1;2;3;4],[theun1; theun2; theun3; theun4],'b-*','LineWidth',1);
xticks([1,2,3,4]);
xticklabels({'h_1=0.5', 'h_1=0.2','h_1=0.05', 'h_1=0.01'});
title('CPU time with RK2');
grid on;
ylabel('CPU time [s]');
subplot(1,2,2);
semilogy([1;2;3;4],[trk41; trk42; trk43; trk44],'r-*','LineWidth',1);
xticks([1,2,3,4]);
xticklabels({'h_1=0.5', 'h_1=0.2','h_1=0.05', 'h_1=0.01'});
title('CPU time with RK4');
grid on;
ylabel('CPU time [s]');

%Trade-Off CPU time VS Error plot
figure(5);
loglog(theun1,max(err_heun1),'x','LineWidth',1,'MarkerSize',15);
hold on;
loglog(theun2,max(err_heun2),'x','LineWidth',1,'MarkerSize',15);
loglog(theun3,max(err_heun3),'x','LineWidth',1,'MarkerSize',15);
loglog(theun4,max(err_heun4),'x','LineWidth',1,'MarkerSize',15);
loglog(trk41,max(err_rk41),'.','LineWidth',1,'MarkerSize',15);
loglog(trk42,max(err_rk42),'.','LineWidth',1,'MarkerSize',15);
loglog(trk43,max(err_rk43),'.','LineWidth',1,'MarkerSize',15);
loglog(trk44,max(err_rk44),'.','LineWidth',1,'MarkerSize',15);
grid on;
ylabel('Err [-]');
xlabel('t_{CPU} [s]');
legend('RK2 with h=0.5', 'RK2 with h=0.1', 'RK2 with h=0.05',...
    'RK2 with h=0.01', 'RK4 with h=0.5', 'RK4 with h=0.1',...
    'RK4 with h=0.05', 'RK4 with h=0.05','Location','southwest');

%% Ex.4
clear; clc; close; 

%Eigenvalues of Ex.3 in (h,lambda)-plane
eigtime_ex3=-((0:0.1:2).^2-1)./((0:0.1:2).^2+1);
h_ex3=[0.5; 0.2; 0.05; 0.01];
hlam1=eigtime_ex3.*h_ex3(1); hlam2=eigtime_ex3.*h_ex3(2); 
hlam3=eigtime_ex3.*h_ex3(3); hlam4=eigtime_ex3.*h_ex3(4);

%Function A(alpha) and RK2/RK4 linear operators definition
A=@(a) [0, 1; -1, 2*cos(a)];
FRK2=@(h,a) eye(2)+h*A(a)+0.5*h^2*A(a)^2;
FRK4=@(h,a) eye(2)+h*A(a)+h^2*A(a)^2/2+h^3*A(a)^3/6+h^4*A(a)^4/24;

%Stability regions computations
[hregrk2, alphark2]=stabreg(FRK2,3);
[hregrk4, alphark4]=stabreg(FRK4,3);

%Plots
figure(1);
plot(hregrk2.*cos(alphark2), hregrk2.*sin(alphark2),'b', 'LineWidth', 2);
hold on;
plot(hregrk4.*cos(alphark4), hregrk4.*sin(alphark4),'r','LineWidth', 2);
plot([hlam1(1), hlam2(1), hlam3(1), hlam4(1)],[0, 0, 0, 0],'mx',...
    'LineWidth',2); 
plot([hlam1(end), hlam2(end), hlam3(end), hlam4(end)],[0, 0, 0, 0],'kx',...
    'LineWidth',2); 
legend('Marginal Stability RK2', 'Marginal Stability RK4',...
    'h_i\lambda(t_0) of Ex.3','h_i\lambda(t_f) of Ex.3',...
    'Location','northwest','NumColumns',1,'FontSize',8);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
xlim([-3,1]);
xticks(-3:1:1);
ylim([0,3]);
grid on;


%Solutions for alpha=pi
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%\n']);
fprintf(['RK2 problem solution h value for alpha = 3.14 is:' ...
    '\n %f \n '],hregrk2(end));
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%\n']);
fprintf(['RK4 problem solution h value for alpha = 3.14 is:' ...
    '\n %f \n '],hregrk4(end));
fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    '%%%%%%%%%%%%%%%%%%%\n']);

%% Ex.5
clear; clc; close; 

%Initialization
A=@(a) [0, 1; -1, 2*cos(a)];
x0=[1;1];
t0=0;
tf=1;
tol=[1e-3, 1e-4, 1e-5, 1e-6];

%Analytical solution
xex=@(t,a) expm(A(a)*t)*x0;

%RK1, RK2, RK4 linear operators definition
FRK1=@(h,a) eye(2)+h*A(a);
FRK2=@(h,a) eye(2)+h*A(a)+h^2*A(a)^2/2;
FRK4=@(h,a) eye(2)+h*A(a)+h^2*A(a)^2/2+h^3*A(a)^3/6+h^4*A(a)^4/24;

%Problem solution
disc=100;
alpha=linspace(0,2*pi,disc);
alpha=alpha';
hreg1=zeros(disc,length(tol));
hreg2=zeros(disc,length(tol));
hreg4=zeros(disc,length(tol));
options = optimset('Display','off','TolX',1e-22);

%Solving algorithm with RK1
for i=1:length(tol)
    for j=1:disc
        xrk1fh= @(h) FRK1(h,alpha(j))^((t0+fix((tf-t0)/h)*h)/h)*x0;
        xexfh=@(h) xex(t0+fix((tf-t0)/h)*h,alpha(j));
        Fh=@(h) norm(xexfh(h)-xrk1fh(h), Inf)-tol(i);
        [hsol,~,~,outputs]=fzero(@(h) Fh(h^2), sqrt(tol(i)), options);
        hreg1(j,i)=(hsol)^2;
    end
end

%Plot with RK1
figure(1);
subplot(2,2,1);
plot(hreg1(:,1).*cos(alpha), hreg1(:,1).*sin(alpha), 'b', 'LineWidth', 2);
title('Tol=10^{-3}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,2);
plot(hreg1(:,2).*cos(alpha),hreg1(:,2).*sin(alpha),'g', 'LineWidth', 2);
title('Tol=10^{-4}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,3);
plot(hreg1(:,3).*cos(alpha),hreg1(:,3).*sin(alpha),'r', 'LineWidth', 2);
title('Tol=10^{-5}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,4);
plot(hreg1(:,4).*cos(alpha),hreg1(:,4).*sin(alpha),'m', 'LineWidth', 2);
title('Tol=10^{-6}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;

% Solving algorithm with RK2
for i=1:length(tol)
    for j=1:disc
        xrk2fh= @(h) FRK2(h,alpha(j))^((t0+fix((tf-t0)/h)*h)/h)*x0;
        xexfh=@(h) xex(t0+fix((tf-t0)/h)*h,alpha(j));
        Fh=@(h) norm(xexfh(h)-xrk2fh(h), Inf)-tol(i);
        [hsol,~,~,outputs]=fzero(@(h) Fh(h^2), sqrt(tol(i)), options);
        hreg2(j,i)=(hsol)^2;
    end
end

%Plot with RK2
figure(2);
subplot(2,2,1);
plot(hreg2(:,1).*cos(alpha),hreg2(:,1).*sin(alpha),'b', 'LineWidth', 2);
title('Tol=10^{-3}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,2);
plot(hreg2(:,2).*cos(alpha),hreg2(:,2).*sin(alpha),'g', 'LineWidth', 2);
title('Tol=10^{-4}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,3);
plot(hreg2(:,3).*cos(alpha),hreg2(:,3).*sin(alpha),'r', 'LineWidth', 2);
title('Tol=10^{-5}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,4);
plot(hreg2(:,4).*cos(alpha),hreg2(:,4).*sin(alpha),'m', 'LineWidth', 2);
title('Tol=10^{-6}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;

%Solving algorithm with RK4
for i=1:length(tol)
    for j=1:disc
        xrk4fh= @(h) FRK4(h,alpha(j))^((t0+fix((tf-t0)/h)*h)/h)*x0;
        xexfh=@(h) xex(t0+fix((tf-t0)/h)*h,alpha(j));
        Fh=@(h) norm(xexfh(h)-xrk4fh(h), Inf)-tol(i);
        [hsol,~,~,outputs]=fzero(@(h) Fh(h^2), sqrt(tol(i)), options);
        hreg4(j,i)=(hsol)^2;
    end
end

%Plot with RK4
figure(3);
subplot(2,2,1);
plot(hreg4(:,1).*cos(alpha),hreg4(:,1).*sin(alpha),'b', 'LineWidth', 2);
title('Tol=10^{-3}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,2);
plot(hreg4(:,2).*cos(alpha),hreg4(:,2).*sin(alpha),'g', 'LineWidth', 2);
title('Tol=10^{-4}')
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,3);
plot(hreg4(:,3).*cos(alpha),hreg4(:,3).*sin(alpha),'r', 'LineWidth', 2);
title('Tol=10^{-5}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;
subplot(2,2,4);
plot(hreg4(:,4).*cos(alpha),hreg4(:,4).*sin(alpha),'m', 'LineWidth', 2);
title('Tol=10^{-6}');
xL=xlim;
yL=ylim;
line([0 0], yL,'Color','black','LineStyle','--','LineWidth',1);
line(xL, [0 0],'Color','black','LineStyle','--','LineWidth',1);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;

%Plot FunEval VS Tol in alpha=pi
figure(4);
loglog(tol,hreg1(50,1:4).^-1,'b*--','LineWidth',1.5);                       %RK1: 1 FunEval
hold on;
loglog(tol,2*hreg2(50,1:4).^-1,'r*--','LineWidth',1.5);                     %RK2: 2 FunEval
loglog(tol,4*hreg4(50,1:4).^-1,'k*--','LineWidth',1.5);                     %RK4: 4 Funval
grid on;
legend('RK1','RK2','RK4');
xlabel('Tol [-]');
ylabel('FunEval [-]');

%% Ex.6
clear; clc; close;

%Initialization 
A=@(a) [0, 1; -1, 2*cos(a)];
F=@(h,a,theta) (eye(2)-(1-theta)*h*A(a)+0.5*(1-theta)^2*h^2*A(a)^2)\...
    (eye(2)+theta*h*A(a)+0.5*theta^2*h^2*A(a)^2);
theta_vec=[0.4,0.1,0.3,0.7,0.9];
hreg=zeros(1000,length(theta_vec));

%Stability region computation
[hreg(:,1), alpha]=stabreg(@(h,a) F(h,a,theta_vec(1)),2.95);
[hreg(:,2), ~]=stabreg(@(h,a) F(h,a,theta_vec(2)),3);
[hreg(:,3), ~]=stabreg(@(h,a) F(h,a,theta_vec(3)),3.5);
[hreg(:,4), ~]=stabreg(@(h,a) F(h,a,theta_vec(4)),3.5);
[hreg(:,5), ~]=stabreg(@(h,a) F(h,a,theta_vec(5)),3.5);

%Plot
figure(1);
plot(hreg(:,1).*cos(alpha),hreg(:,1).*sin(alpha), 'LineWidth', 2);
hold on;
plot(hreg(:,2).*cos(alpha),hreg(:,2).*sin(alpha), 'LineWidth', 2);
plot(hreg(:,3).*cos(alpha),hreg(:,3).*sin(alpha), 'LineWidth', 2);
plot(hreg(:,4).*cos(alpha),hreg(:,4).*sin(alpha), 'LineWidth', 2);
plot(hreg(:,5).*cos(alpha),hreg(:,5).*sin(alpha), 'LineWidth', 2);
legend('\theta=0.4','\theta=0.1','\theta=0.3','\theta=0.7','\theta=0.9',...
    'Location','northwest');
ylim([0,6]);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;

%% Ex.7
clear; clc; close;

%Initialization
B=[-180.5, 219.5; 179.5, -220.5];
xdot=@(x,t) B*x;
x0=[1;1];
t0=0; tf=5; h=0.1;
eigB=eig(B); eigBh=h*eigB;

%Analytical solution
timevec=t0:h:tf;
len=length(timevec);
xes=zeros(2,len);
for i=1:len
    xes(:,i)=expm(B*timevec(i))*x0;
end

%RK4 solution and errors
[timerk4,xrk4]=RK4(@(x,t) xdot(x,t), t0, tf, h, x0);  
err_RK4=abs(xes-xrk4);

%IEX4 solution and errors
[~,xiex4]=IEX4(@(x,t) xdot(x,t), t0, tf, h, x0);  
err_IEX4=abs(xes-xiex4);

%Plots solutions
figure(1);
subplot(3,1,2);
plot(timerk4,xrk4(1,:),'b--*','LineWidth',1);
hold on;
plot(timerk4,xrk4(2,:),'r--*','LineWidth',1);
title('RK4');
xlabel('t [s]');
ylabel('x [-]');
legend('x_1','x_2','Location','northwest');
grid on;
subplot(3,1,3);
plot(timerk4,xiex4(1,:),'b--*','LineWidth',1);
hold on;
plot(timerk4,xiex4(2,:),'r--*','LineWidth',1);
title('IEX4');
xlabel('t [s]');
ylabel('x [-]');
legend('x_1','x_2','Location','northeast');
grid on;
subplot(3,1,1);
plot(timerk4,xes(1,:),'b','LineWidth',1);
hold on;
plot(timerk4,xes(2,:),'r','LineWidth',1);
title('EXACT');
xlabel('t [s]');
ylabel('x [-]');
legend('x_1','x_2','Location','northeast');
grid on;

%Plot errors
figure(2);
subplot(2,1,1);
semilogy(timerk4,err_RK4(1,:),'b--*','LineWidth',1);
hold on;
semilogy(timerk4,err_RK4(2,:),'r--*','LineWidth',1);
title('Exact - RK4');
xlabel('t [s]');
ylabel('Err [-]');
legend('err_{x_1}','err_{x_2}','Location','northwest');
grid on;
subplot(2,1,2);
semilogy(timerk4,err_IEX4(1,:),'b--*','LineWidth',1);
hold on;
semilogy(timerk4,err_IEX4(2,:),'r--*','LineWidth',1);
title('Exact - IEX4');
xlabel('t [s]');
ylabel('Err [-]');
legend('err_{x_1}','err_{x_2}','Location','northeast');
grid on;

%Function and operators definition
A=@(a) [0, 1; -1, 2*cos(a)];
FRK4=@(h,a) eye(2)+h*A(a)+h^2*A(a)^2/2+h^3*A(a)^3/6+h^4*A(a)^4/24;
FIEX4=@(h,a) -((eye(2)-h*A(a))^-1)/6+...
            4*(eye(2)-0.5*h*A(a))^-2-...
            27*((eye(2)-(h/3)*A(a))^-3)/2+...
            32*((eye(2)-0.25*h*A(a))^-4)/3;

%Stability regions computations
[hregrk4, alphark4]=stabreg(FRK4, 3);
[hregrk2, alphaiex4]=stabreg(FIEX4, 5.5);

%Plots stability regions
figure(3);
plot(hregrk4.*cos(alphark4),hregrk4.*sin(alphark4), 'LineWidth', 2);
hold on;
plot(hregrk2(2:end).*cos(alphaiex4(2:end)),...
    hregrk2(2:end).*sin(alphaiex4(2:end)), 'LineWidth', 2);
plot(eigBh,[0;0],'x','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6],'LineWidth',2);
legend('Marginal Stability RK4', 'Marginal Stability IEX4',...
    '\{\lambdah\} points','Location','northwest');
ylim([0,8]);
xlabel('Re\{\lambdah\}');
ylabel('Im\{\lambdah\}');
grid on;

%% Ex.8
clear; clc; close;

%Initialization
Fdot=@(x,y,t) [-(5/2)*(1+8*sin(t))*x; (1-x)*y+x];
X0=[1;1];
t0=0; tf=3;
h=0.1;

%Integrations
[tspan, xab3]=AB3(@(x,t) Fdot(x(1),x(2),t), t0, tf, h, X0);
[~, xam3]=AM3(@(x,t) Fdot(x(1),x(2),t), t0, tf, h, X0);
[~, xabm3]=ABM3(@(x,t) Fdot(x(1),x(2),t), t0, tf, h, X0);
[~, xbdf3]=BDF3(@(x,t) Fdot(x(1),x(2),t), t0, tf, h, X0);

%Plots
figure(1);
subplot(2,2,1);
plot(tspan,xab3(1,:),'b--*', 'LineWidth', 1);
hold on;
plot(tspan,xab3(2,:),'r--*', 'LineWidth', 1);
title('AB3');
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','y(t)','Location','northwest');
grid on;
hold on;
subplot(2,2,2)
plot(tspan,xam3(1,:),'b--*', 'LineWidth', 1);
hold on;
plot(tspan,xam3(2,:),'r--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','y(t)','Location','northwest');
grid on;
title('AM3');
subplot(2,2,3)
plot(tspan,xabm3(1,:),'b--*', 'LineWidth', 1);
hold on;
plot(tspan,xabm3(2,:),'r--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','y(t)','Location','northwest');
grid on;
title('ABM3');
subplot(2,2,4)
plot(tspan,xbdf3(1,:),'b--*', 'LineWidth', 1);
hold on;
plot(tspan,xbdf3(2,:),'r--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','y(t)','Location','northwest');
title('BDF3');
grid on;

%Discussion
eigA=@(t,x,y) [-2.5-2.5*8*sin(t), 0; -y+1, 1-x];
eigA_vec=zeros(2,length(tspan));
for i=1:length(tspan)
    eigA_vec(:,i)=eig(eigA(tspan(i),xbdf3(1,i),xbdf3(2,i)));  
end
figure(2)
plot(tspan, h*eigA_vec(1,:),'r*--','LineWidth',1);
hold on;
plot(tspan, h*eigA_vec(2,:),'b*--','LineWidth',1);
xlabel('t [s]');
ylabel('h\lambda [-]');
legend('h\lambda_1','h\lambda_2','Location','southwest');
grid on;

%ABM3 observation
figure(3);
subplot(2,2,1);
plot(tspan,xab3(1,:),'b--*', 'LineWidth', 1);
title('AB3');
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','Location','northwest');
grid on;
hold on;
subplot(2,2,2)
plot(tspan,xam3(1,:),'b--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','Location','northwest');
grid on;
title('AM3');
subplot(2,2,3)
plot(tspan,xabm3(1,:),'b--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','Location','northwest');
grid on;
title('ABM3');
subplot(2,2,4)
plot(tspan,xbdf3(1,:),'b--*', 'LineWidth', 1);
xlabel('t [s]');
ylabel('X [-]');
legend('x(t)','Location','northwest');
title('BDF3');
grid on;

%% Functions

%Bisection for Ex.1
function [sol, feval]=bisection(fun, a, b, tol)
    fa=fun(a);                                                             %Initialization
    fb=fun(b);
    feval=0;
    if fa*fb>0                                                             %Initial control
        fprintf('No roots in [a,b] interval!');
    elseif fa==0
        sol=a; return
    elseif fb==0
        sol=b; return
    end
    err=tol+1;
    xold=a;
    while err>=tol                                                          %Iterations
        feval=feval+2;
        xk=a+(b-a)*0.5;
        fk=fun(xk);
        if fa*fk<0
            b=xk;
            fb=fun(b);
        elseif fk*fb<0
            a=xk;
            fa=fun(a); 
        end                                                                      
        err=abs(xold-xk);
        xold=xk;
    end
    sol=xk;
end

%Secant for Ex.1
function [sol, feval]=secant(fun, x0, x1, tol)
    feval=0;                                                                %Initalization
    f0=fun(x0);
    f1=fun(x1);
    if f0==0                                                                %Initial control
        sol=x0; return
    elseif f1==1
        sol=x1; return
    end
    err=tol+1;
    xold=x0;
    while err>=tol                                                          %Iterations
        feval=feval+2;
        xk=x1-fun(x1)*((fun(x1)-fun(x0))/(x1-x0))^-1;
        x0=x1;
        x1=xk;   
        err=abs(xold-xk);
        xold=xk;
    end
    sol=xk;
end

%Regula-Falsi for Ex.1
function [sol, feval]=regulafalsi(fun,a,b,tol)
    fa=fun(a);                                                              %Initialization
    fb=fun(b);
    feval=0;
    if fa*fb>0                                                              %Initial control
        fprintf('No roots in [x0,x1] interval!');
    elseif fa==0
        sol=a; return
    elseif fb==0
        sol=b; return
    end
    err=tol+1;
    xold=a;
    while err>=tol                                                          %Iterations
        feval=feval+3;
        xk=b-fun(b)*((fun(b)-fun(a))/(b-a))^-1;
        fk=fun(xk);
        if fa*fk<0
            b=xk;
            fb=fun(b);
        elseif fk*fb<0
            a=xk;
            fa=fun(a);
        end
        err=abs(xold-xk);
        xold=xk;
    end
    sol=xk;
end

%Newton for Ex.2
function [sol, it]=newton(Floc, JFloc, x0, tol)
    err=tol+1;                                                              %Initialization
    it=0;
    xk=x0;
    while err>=tol                                                          %Iterations
        df=JFloc(xk);
        f=Floc(xk);
        delta=-df\f;
        xk=xk+delta;
        err=norm(delta);
        it=it+1;
    end
    sol=xk;
end

%RK2 for Ex.3
function [tspan,x]=RK2(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    x(:,1)=x0;
    for i=1:length(tspan)-1                                                 %Iterations
        xp=x(:,i)+h*rhs(x(:,i),tspan(i));
        xc=x(:,i)+0.5*h*(rhs(x(:,i),tspan(i))+rhs(xp,tspan(i+1)));
        x(:,i+1)=xc;
    end
end

%RK4 for Ex.3
function [tspan, x]=RK4(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    x(:,1)=x0;
    for i=1:length(tspan)-1                                                 %Iterations
        k1=rhs(x(:,i),tspan(i));
        k2=rhs(x(:,i)+0.5*h*k1,tspan(i)+0.5*h);
        k3=rhs(x(:,i)+0.5*h*k2,tspan(i)+0.5*h);
        k4=rhs(x(:,i)+h*k3,tspan(i)+h);
        xc=x(:,i)+h*(k1+2*k2+2*k3+k4)/6;
        x(:,i+1)=xc;
    end
end

%Stability region for Ex.4
function [hreg,alpha]=stabreg(Foperator,h0) 
    disc=1000;                                                              %Initialiation
    alpha=linspace(0,pi,disc);
    alpha=alpha';
    hreg=zeros(disc,1);
    options = optimset('Display','off', 'TolFun', 1e-22);
    for i=1:disc                                                            %Solve equation for each alpha
        Fh=@(h) max(abs(eig(Foperator(h,alpha(i)))))-1;
        hreg(i)=fsolve(@(h) Fh(h), h0, options);
    end
end

%BE for Ex.7
function [tspan,x]=BE(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initalization
    x=zeros(length(x0),length(tspan));
    x(:,1)=x0;
    options = optimset('Display','off', 'TolFun', 1e-22);
    for i=1:length(tspan)-1                                                 %Iterations
        fun=@(xnew) xnew-x(:,i)-h*rhs(xnew,tspan(i+1));
        x(:,i+1)=fsolve(@(xnew) fun(xnew), x(:,i),options);                 %Roots finding (i.e. explicit method)
    end
end

%IEX4 for Ex.7
function [tspan, x]=IEX4(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    x(:,1)=x0;
    for i=1:length(tspan)-1                                                 %Iterations
        [~,xvec1]=BE(rhs, tspan(i), tspan(i+1), h  , x(:,i));...
            xp1=xvec1(:,end);
        [~,xvec2]=BE(rhs, tspan(i), tspan(i+1), h/2, x(:,i));...
            xp2=xvec2(:,end);
        [~,xvec3]=BE(rhs, tspan(i), tspan(i+1), h/3, x(:,i));...
            xp3=xvec3(:,end);
        [~,xvec4]=BE(rhs, tspan(i), tspan(i+1), h/4, x(:,i));...
            xp4=xvec4(:,end);
        xc=-xp1/6+4*xp2-27*xp3/2+32*xp4/3;
        x(:,i+1)=xc;
    end
end

%AB3 for Ex.8
function [tspan, x]=AB3(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    [~, x(:,1:3)]=RK4(rhs, t0, t0+2*h, h, x0);                              %Initial single-step iterations
    for i=3:length(tspan)-1                                                 %Iterations
        x(:,i+1)=x(:,i)+h*(23*rhs(x(:,i),tspan(i))-...
            16*rhs(x(:,i-1),tspan(i-1))+5*rhs(x(:,i-2),tspan(i-2)))/12;
    end
end

%AM3 for Ex.8
function [tspan, x]=AM3(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    [~, x(:,1:2)]=RK4(rhs, t0, t0+h, h, x0);                                %Initial single-step iterations
    options = optimset('Display','off', 'TolFun', 1e-22);
    for i=2:length(tspan)-1                                                 %Iterations
        fun=@(xnew) xnew-x(:,i)-h*5*rhs(xnew,tspan(i+1))/12-...
            h*8*rhs(x(:,i),tspan(i))/12+h*rhs(x(:,i-1),tspan(i-1))/12;
        x(:,i+1)=fsolve(@(xnew) fun(xnew),x0,options);                      %Roots finding (i.e. explicit method)
    end
end

%ABM3 for Ex.8
function [tspan, x]=ABM3(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    [~, x(:,1:3)]=RK4(rhs, t0, t0+2*h, h, x0);                              %Initial single-step iterations
    for i=3:length(tspan)-1                                                 %Iterations
        xp=x(:,i)+h*(23*rhs(x(:,i),tspan(i))-16*rhs(x(:,i-1),tspan(i-1))+...
            5*rhs(x(:,i-2),tspan(i-2)))/12;
        xc=x(:,i)+h*(5*rhs(xp,tspan(i+1))+8*rhs(x(:,i),tspan(i))-...
            rhs(x(:,i-1),tspan(i-1)))/12;
        x(:,i+1)=xc;
    end
end

%BDF3 for Ex.8
function [tspan, x]=BDF3(rhs, t0, tf, h, x0)
    tspan=t0:h:tf;                                                          %Initialization
    x=zeros(length(x0),length(tspan));
    [~, x(:,1:3)]=RK4(rhs, t0, t0+2*h, h, x0);                              %Initial single-step iterations
    options = optimset('Display','off', 'TolFun', 1e-22);
    for i=3:length(tspan)-1                                                 %Iterations
        fun=@(xnew) xnew-18*x(:,i)/11+9*x(:,i-1)/11-2*x(:,i-2)/11-...
            h*6*rhs(xnew,tspan(i+1))/11;
        x(:,i+1)=fsolve(@(xnew) fun(xnew),x0,options);                      %Roots finding (i.e. explicit method)
    end
end