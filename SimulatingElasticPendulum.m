%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Project: The Elastic Pendulum
%
% Diana Zhen Zhang (805777341)
% Dr. Marcus Rüter, CEE 103
% Date: 08/27/2023
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clc;
clearvars;
close all;

%% Initialization
whichProblem = 1; 

% 1 - Part I: One pendulum
% 2 - Part II: Three pendulums

% Methods: 
% 1 - Heun's method
% 2 - Forward Euler Method

set(0,'DefaultFigureRenderer','painters'); % Create vector graphics

% Define conditions per method
if whichProblem == 1 % Part I
    whichMethod = 1;
    nPend = 1;
    t0 = 0;
    tf = 6;
    nSI = 2000; % number of sub intervals

elseif whichProblem == 2
    % Part II
    whichMethod = 2;
    nPend = 3;
    t0 = 0;
    tf = 6;
    nSI = 20000; % Define conditions per method
else
    error('Unknown problem!');
end

g = 9.796; % Gravity of Earth at UCLA [m/s^2]
L = 1.7; % Undeformed length of pendulum [m]
m = 0.6; % Mass of the sphere [kg]
r = 0.3; % Radius of the sphere [m]
data(1) = g; % Store g in data array
data(2) = L; % Store L in data array
data(3) = m; % Store m in data array
data(4) = r; % Store r in data array
Ve = 0; % Elastic potential energy
u0 = 0; % Initial value for u
theta0 = 175 * pi/180; % Initial value for theta (Part I)
v0 = 0; % Initial value for v

omega0 = 0; % Initial value for omega
nN = nSI + 1; % Number of nodes
h = (tf-t0)/nSI; % Subinterval length h
t = 0:h:tf; % Array of t-nodes
y = zeros(2, nN); % Approximate solution vector at nodes t_n

y(1,1,1) = u0; % Initial condition for u (pendulum 1)
y(2,1,1) = theta0; % Initial condition for theta (pendulum 1)
y(3,1,1) = v0; % Initial condition for v (pendulum 1)
y(4,1,1) = omega0; % Initial condition for omega (pendulum 1)
y(1,1,2) = u0; % Initial condition for u (pendulum 2)
y(2,1,2) = 176 * pi/180; % Initial condition for theta (pendulum 2)
y(3,1,2) = v0; % Initial condition for v (pendulum 2)
y(4,1,2) = omega0; % Initial condition for omega (pendulum 2)
y(1,1,3) = u0; % Initial condition for u (pendulum 3)
y(2,1,3) = 177 * pi/180; % Initial condition for theta (pendulum 3)
y(3,1,3) = v0; % Initial condition for v (pendulum 3)
y(4,1,3) = omega0; % Initial condition for omega (pendulum 3)

xSphere = zeros(1,nPend); % x-position of the sphere(s)
ySphere = zeros(1,nPend); % y-position of the sphere(s)
plotSpring = zeros(1,nPend); % Array for plotting each spring
plotSphere = zeros(1,nPend); % Array for plotting each sphere
dummyX = zeros(2,1); % Dummy x- and y-posit. of spring for plot.
dummyC = 0; % Dummy center of sphere for plotting

% Measured u-values of spring
u_i = [0.05 0.17 0.22 0.28 0.35 0.46 0.51 0.59 0.63 0.69 0.75 0.84 0.98 ...
1.12 1.25 1.33 1.47 1.64 1.82 2.11 2.32 2.64 2.92 3.27 3.65 4.31];

% Measured F-values of spring
F_i = [0.4 2.2 4.6 5.8 7.3 8.5 9.1 10.7 11.1 12.3 12.8 13.5 14.2 14.9 15.2 ...
15.8 16.4 16.3 16.7 17.1 17.5 18.4 18.9 20.7 21.9 24.1];
nDP = length(F_i); % Number of force displacement data points

polyOrd = 4; % Polynomial order for least squares
p = zeros(polyOrd,length(F_i)); % Vectors of basis functions
M = zeros(polyOrd,polyOrd); % M-matrix
b = zeros(polyOrd,1); % Right-hand side b-vector

uPlot = linspace(0,5,1000); % u-nodes for plotting
FPlot = zeros(1,length(uPlot)); % F-values (at u-nodes) for plotting

%% Part Ia) Least-squares Data Fitting
for currOrd = 1:polyOrd
    p(currOrd,:) = u_i .^ currOrd;
end

for fDisp = 1:nDP
    M = M + p(:,fDisp) * p(:,fDisp)';
    b = b + p(:,fDisp) * F_i(fDisp);
end

a = M \ b;
data(5:5 + polyOrd - 1) = a;

for currOrd = 1:polyOrd
    fPlot = a(currOrd) * uPlot.^(currOrd);
    FPlot = FPlot + fPlot;
end

%% Parts Ib) and IIb) Compute Numerical Results for Non-linear IVP
switch whichMethod % Compute result based on numerical method
    case 1 % .. method
        for i = 1:nPend % Loop over all pendulums
            y(:,:,i) = computeHeunSol(y(:,:,i),data,h,nSI);
        end
       
    case 2 % .. method
        for i = 1:nPend % Loop over all pendulums
            y(:,:,i) = computeFEulerSol(y(:,:,i),data,h,nSI);
        end
    
    otherwise
    error('Unknown numerical method!');
end

%% Part Ic) Compute Energy

if whichProblem == 1 % Only required for Part I
    KE = zeros(1, nN); % Initialize array for kinetic energy
    PE = zeros(1, nN); % Initialize array for potential energy
    E = zeros(1, nN);

for n = 1:nN
    % Calculate kinetic energy
    v_vector = [y(3, n, 1), (L + y(1, n, 1)) * y(4, n, 1)];
    KE(n) = 0.5 * m * dot(v_vector, v_vector);
    
    % Calculate gravitational potential energy
    h = L + y(1, n, 1);
    PE_gravity = m * g * h;
    
    % Calculate elastic potential energy
    PE_elastic = 0.5 * Ve * y(1, n, 1)^2;
    % Total energy at time node t_n
    E(n) = KE(n) + PE_gravity + PE_elastic;
end
end

%% Parts Id) and IIc) Create Pendulum Animation

figure(1); % Open figure 1
hold on; % Put hold to on

for i = 1:nPend % Loop over all pendulums
    plotSpring(i) = plot(dummyX,dummyX,'LineWidth',5); % Plot init. spring
    plotSphere(i) = plot(dummyC,dummyC,'.','MarkerSize',250); % Plot sph.
end

plot(0,0,'.','MarkerSize',50); % Plot pin
title('The Elastic Pendulum'); % Set title
xlabel('x','Interpreter','LaTeX'); % Set x-label
ylabel('y','Interpreter','LaTeX'); % Set y-label
xlim([-1.15*(L+max(y(1,:))) 1.15*(L+max(y(1,:)))]); % Set x-limits
ylim([-1.15*(L+max(y(1,:))) 1.15*(L+max(y(1,:)))]); % Set y-limits

grid on; % Turn grid on

axis square; % Use equal lengths in the x- and y-direc.

set(gcf,'Position',[50 50 900 900]); % Set position and size of figure
set(gca,'LineWidth',3,'FontSize',18); % Set axis widths and font size

for n = 1:nSI % Loop over all subintervals
    for i = 1:nPend % Loop over all pendulums
        % Calculate x- and y-positions of sphere (or spring) at t_n+1
        xSphere(i) = (L+y(1,n+1,i))*sin(y(2,n+1,i));
        ySphere(i) = -(L+y(1,n+1,i))*cos(y(2,n+1,i));
        
        % Copy positions of sphere into plotSpring and plotSphere
        set(plotSpring(i),'xdata',[0 xSphere(i)],'ydata',[0 ySphere(i)]);
        set(plotSphere(i),'xdata',xSphere(i),'ydata',ySphere(i));
    end

if n == 1
pause % Stop initial frame of animation
end

drawnow; % Make sure to update pendulum plot
end

%% Part Ie) Construct Interpolant
if whichProblem == 1 % Create interpolant only for Part I
    [tPlot,yPlot,EPlot] = constructPQuadInterpl(t,y,E,nSI,t0,tf);
end

%% Part If) Plot Results for Part I
if whichProblem == 1 % Create plots only for Part I

    figure(2); % Open figure 2 for results from Part Ia)
    hold on; % Put hold to on
    plot(u_i,F_i,'o','Markersize',10) % Plot (u_i,F_i) data points
    plot(uPlot,FPlot,'LineWidth',2); % Plot least-sq. F(u) interpolant
    
    title('Part Ia): Least-squares Approximant'); % Set title
    xlabel('$u$','Interpreter','LaTeX'); % Set x-label
    ylabel('$F(u)$','Interpreter','LaTeX'); % Set y-label
    xlim([0 5]); % Set x-limits
    ylim([0 25]); % Set y-limits
    grid on; % Turn grid on
    
    legend('Data points $(u_i,F_i)$','Least-squares approximant $F(u)$',...
    'Location','SE','Interpreter','LaTeX'); % Create legend
    set(gcf,'Position',[30 350 1200 750]); % Set position and size of fig.
    set(gca,'LineWidth',2,'FontSize',20); % Set axis widths and font size
    
    figure(3); % Open figure 3 for results from Part Ic)
    plot(t,E,'LineWidth',2); % Plot interpolant E(t)
    title('Part Ic): Energy'); % Set title
    xlabel('$t$','Interpreter','LaTeX'); % Set x-label
    ylabel('$E(t)$','Interpreter','LaTeX'); % Set y-label
    xlim([t0 tf]); % Set x-limits
    grid on; % Turn grid on
    
    legend('Approximate energy $E(t)$','Location','SE','Interpreter',...
    'LaTeX'); % Create legend
    set(gcf,'Position',[30 350 1200 750]); % Set position and size of fig.
    set(gca,'LineWidth',2,'FontSize',20); % Set axis widths and font size
    
    figure(4); % Open figure 4 for results from Part Ib)
    hold on; % Put hold to on
    plot(tPlot,yPlot(1,:),'LineWidth',2); % Plot interpolated u(t)
    plot(tPlot,yPlot(3,:),'LineWidth',2); % Plot interpolated v(t)
    
    title('Part Ib): Displacement and Speed'); % Set title
    xlabel('$t$','Interpreter','LaTeX'); % Set x-label
    ylabel('$u(t) / v(t)$','Interpreter','LaTeX'); % Set y-label
    xlim([t0 tf]); % Set x-limits
    grid on; % Turn grid on
    
    legend('Approximate solution $u(t)$','Approximate solution $v(t)$',...
    'Location','SW','Interpreter','LaTeX'); % Create legend
    set(gcf,'Position',[30 350 1200 750]); % Set position and size of fig.
    set(gca,'LineWidth',2,'FontSize',20); % Set axis widths and font size
    
    figure(5); % Open figure 5 for results from Part Ib)
    hold on; % Put hold to on
    plot(tPlot,yPlot(2,:),'LineWidth',2); % Plot interpolated theta(t)
    plot(tPlot,yPlot(4,:),'LineWidth',2); % Plot interpolated omega(t)
    
    title('Part Ib): Angular Position and Angular Velocity'); % Set title
    xlabel('$t$','Interpreter','LaTeX'); % Set x-l.
    ylabel('$\theta(t) / \omega(t)$','Interpreter','LaTeX'); % Set y-l.
    xlim([t0 tf]); % Set x-limits
    grid on; % Turn grid on
    
    legend('Approximate solution $\theta(t)$', ...
    'Approximate solution $\omega(t)$','Location','SW',...
    'Interpreter','LaTeX'); % Create legend
    set(gcf,'Position',[30 350 1200 750]); % Set position and size of fig.
    set(gca,'LineWidth',2,'FontSize',20); % Set axis widths and font size

end

%% Function computeHeunSol
function y = computeHeunSol(y, data, h, nSI)
    for n = 1:nSI
        k1 = f(y(:, n), data); % Calculate k1 using f and current state y
        k2 = f(y(:, n) + h * k1, data); % Calculate k2 using f and predicted state k1
        y(:, n + 1) = y(:, n) + 0.5 * h * (k1 + k2); % Update y for next time step
    end
end

%% Function computeFEulerSol
function y = computeFEulerSol(y,data,h,nSI)
    for n = 1:nSI % Loop over all subintervals
        y(:,n+1) = y(:,n) + h * f(y(:,n),data); % Compute y_n+1 using Forward Euler
    end
end

%% Function f
function fReturn = f(Y, data)
    % Define the ODEs based on the given system
    % Y is a vector containing [u, theta, v, omega]
    u = Y(1);
    theta = Y(2);
    v = Y(3);
    omega = Y(4);
    
    % Extract constants from the data vector
    g = data(1);
    L = data(2);
    m = data(3);
    r = data(4);
    a = data(5:end);
    
    % Calculate the sum term for the polynomial force F(u)
    sumTerm = 0;
    for ord = 1:length(a)
        sumTerm = sumTerm + a(ord) * u^ord;
    end
    
    % Compute the derivatives of the state variables
    du_dt = v;
    dtheta_dt = omega;
    dv_dt = (L+u) * omega^2 + g * cos(theta) - sumTerm / m;
    domega_dt = -((2*v*omega + g*sin(theta)) / ((L+u)^2 + (2/5)*r^2)) * (L+u);
    
    % Return the derivatives as a column vector
    fReturn = [du_dt; dtheta_dt; dv_dt; domega_dt];
end

%% Function constructPQuadInterpl
function [tPlot, yPlot, EPlot] = constructPQuadInterpl(t, y, E, nSI, t0, tf)

% Initialization
tPlot = linspace(t0, tf, length(t)); % t-nodes for plotting
yPlot = zeros(4, length(tPlot)); % y-values (at t-nodes) for plotting
EPlot = zeros(1, length(tPlot)); % E-values (at t-nodes) for plotting

% Define Lagrangean basis functions as anonymous functions
L_n = @(tPlotSI, n) (tPlotSI - t(n+1)) .* (tPlotSI - t(n+2)) / ((t(n) - t(n+1)) * (t(n) - t(n+2)));
L_np1 = @(tPlotSI, n) (tPlotSI - t(n)) .* (tPlotSI - t(n+2)) / ((t(n+1) - t(n)) * (t(n+1) - t(n+2)));
L_np2 = @(tPlotSI, n) (tPlotSI - t(n)) .* (tPlotSI - t(n+1)) / ((t(n+2) - t(n)) * (t(n+2) - t(n+1)));

for n = 1:2:nSI % Loop over every other subinterval
    indices = [n n+1 n+2]; % Indices of t-nodes in 2 subintervals
    tPlotSI = tPlot(indices); % t-nodes in 2 subintervals
    
    % Evaluate Lagrangean interpolant for y and E at tPlotSI
    yPlot(:, indices) = y(:, n) * L_n(tPlotSI, n) + y(:, n+1) * L_np1(tPlotSI, n) + y(:, n+2) * L_np2(tPlotSI, n);
    EPlot(indices) = E(n) * L_n(tPlotSI, n) + E(n+1) * L_np1(tPlotSI, n) + E(n+2) * L_np2(tPlotSI, n);
end
end
