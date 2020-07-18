%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations for Diamond Shape Airfoil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By: Gülce Topal (gulce.topal@gmail.com)

% Program to calculate and plot pressure coefficients, pressure values 
% and drag/lift coefficients over a symmetric double wedge airfoils
% at different angle of attacks and freestream Mach numbers.

% Assumptions: Inviscid, steady, compressible two dimensional flow for a
% calorically perfect gas. Airfoil is symmetric with a constant wedge angle
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% initialization

%parameters 
M1 = 2;       % freestream mach number
gamma = 1.4;  % isentropic expansion factor
alpha = 5;    % angle of attack [degrees]
wedge = 6;    % wedge half angle [degrees]
p1 = 1;       % freestream air pressure [atm]
c = 10;       % chord length [m]

%shock-expansion wave angle calculations
for INITIALS=1
syms p1;
x2 = abs(wedge + alpha);   %angle between flow and wedge surface for region 2
if x2 == 0
    x2 = x2 + 0.00000001;
end
x3 = abs(wedge - alpha);   %angle between flow and wedge surface for region 3
if x3 == 0
    x3 = x3 + 0.00000001;
end
x4 = 2*wedge;              %angle between flow and wedge surface for region 4
x5 = 2*wedge;              %angle between flow and wedge surface for region 5
initials = [["M1";"angle of attack";"wedge angle";"theta2";"theta3";"theta4";"theta5"],[M1; alpha; wedge; x2; x3; x4; x5]]

%if angle of attack is bigger than wedge angle, region three forms an
%expansion wave
expflag = 0; %for oblique shock
if alpha > wedge
    expflag = 1; %for expansion wave
end
end 

%% shock-expansion theory calculations
%oblique shock/expansion wave for region 3
if expflag == 0
    beta3 = thetabetaM(M1, x3, gamma);
    M1n = M1*sind(beta3);
    M3n = sqrt((1+((gamma-1)*M1n^2)/2)/(gamma*M1n^2 - (gamma-1)/2));
    M3 = M3n/sind(beta3-x3);
    p3p1 = 1 + ((2*gamma)/(gamma+1))*(M1n^2-1);
    p3 = p3p1 * p1;
end
if expflag == 1
    v1 = prandtlmeyer(M1, gamma, 'mu');                  
    v3 = x3 + v1;
    M3 = prandtlmeyer(v3, gamma, 'mach');                       
    p03p3 = totaltostatic(M3, gamma);
    p01p1 = totaltostatic(M1, gamma);
    p3p1 = p01p1/p03p3;
    p3 = p3p1 * p1;
end

%oblique shock for region 2
beta2 = thetabetaM(M1, x2, gamma);
M1n = M1*sind(beta2);
M2n = sqrt((1+((gamma-1)*M1n^2)/2)/(gamma*M1n^2 - (gamma-1)/2));
M2 = M2n/sind(beta2-x2);
p2p1 = 1 + ((2*gamma)/(gamma+1))*(M1n^2-1);
p2 = p2p1 * p1;

%expansion wave for region 5
v3 = prandtlmeyer(M3, gamma, 'mu');                  
v5 = x5 + v3;
M5 = prandtlmeyer(v5, gamma, 'mach');                       
p05p5 = totaltostatic(M5, gamma);
p03p3 = totaltostatic(M3, gamma);
p5p3 = p03p3/p05p5;
p5 = p5p3 * p3;

%expansion wave for region 4
v2 = prandtlmeyer(M2, gamma, 'mu');              
v4 = x4 + v2;
M4 = prandtlmeyer(v4, gamma, 'mach');
p04p4 = totaltostatic(M4, gamma);
p02p2 = totaltostatic(M2, gamma);
p4p2 = p02p2/p04p4;
p4 = p4p2 * p2;

%calculation of pressure coefficients
cp2 = double(pressurecoeff(p1, p2, M1, gamma));
cp3 = double(pressurecoeff(p1, p3, M1, gamma));
cp4 = double(pressurecoeff(p1, p4, M1, gamma));
cp5 = double(pressurecoeff(p1, p5, M1, gamma));

%calculation of drag and lift coefficients
cd = ((p2/p1-p5/p1)*sind(wedge+alpha)+(p3/p1-p4/p1)*sind(wedge-alpha))/(cosd(wedge)*1.4*4);
cl = ((p2/p1-p5/p1)*cosd(wedge+alpha)+(p4/p1-p3/p1)*cosd(wedge-alpha))/(cosd(wedge)*1.4*4);

%results
mach_numbers = [["M3";"M2"; "M4"; "M5"],[M3; M2; M4; M5]]
static_pressures = vpa([["p2p1"; "p3p1"; "p4p1"; "p5p1"],[vpa(p2p1,4); vpa(p3p1,4); vpa(p4p2*p2/p1,4); vpa(p5p3*p3/p1,4)]],4)
pressure_coeffs = vpa([["cp2"; "cp3"; "cp4"; "cp5"],[vpa(cp2,4); vpa(cp3,4); vpa(cp4,4); vpa(cp5,4)]],4)
coeffs_shockexpansion = vpa([["cd"; "cl"],[cd; cl]],4)

%% plotting
yv2 = [];yv3 = [];yv4 = [];yv5 = [];
l = (c/2)/cosd(6); %surface length [m]
h = (c/2)*tand(6); %height of the one half of the wedge [m]

%first and second halves of the chord length
xv1 = 0:0.01:c/2;
xv2 = c/2:0.01:c;

%wedge drawing
figure(1); subplot(2,1,1); hold on;
for x = 0:0.01:c
    if x<= c/2
    y3 = tand(6)*(x-c/2) + h;
    yv3 = [yv3 y3];
    y2 = -tand(6)*(x+c/2) + h;
    yv2 = [yv2 y2];
    end
    if x>c/2
    y5 = -tand(6)*(x-c);
    yv5 = [yv5 y5];
    y4 = tand(6)*(x-c);
    yv4 = [yv4 y4];
    end
end
yv5 = [yv5 0]; yv4=[yv4 0];
plot(xv1, yv3, '-k.'); 
plot(xv1, yv2 ,'-k.');
plot(xv2, yv5, '-k.');
plot(xv2, yv4, '-k.'); 

%pressure distribution
title('Pressure Distribution')
p3 = subs(p3,p1,1); p2 = subs(p2,p1,1); p4 = subs(p4,p1,1); p5 = subs(p5,p1,1);
pd3 = quiver(xv1,yv3, -ones(1,length(xv1)), p3*ones(1,length(xv1)),'color',[1 0 0],'ShowArrowHead',0,'AlignVertexCenters',1);
pd5 = quiver(xv2,yv5, ones(1,length(xv2)), p5*ones(1,length(xv2)),'color',[0 1 0],'ShowArrowHead',0,'AlignVertexCenters',1);
pd2 = quiver(xv1,yv2, -ones(1,length(xv1)), -p2*ones(1,length(xv1)),'color',[1 1 0],'ShowArrowHead',0,'AlignVertexCenters',1);
pd4 = quiver(xv2,yv4, ones(1,length(xv2)), -p4*ones(1,length(xv2)),'color',[0 0 1],'ShowArrowHead',0,'AlignVertexCenters',1);
hold off;

%pressure coefficient distribution
subplot(2,1,2); hold on;
title('Pressure Coefficient Distribution')
c3 = plot(xv1, cp3, '-r.'); 
c2 = plot(xv1, cp2 ,'-y.');
c4 = plot(xv2, cp5, '-g.');
c5 = plot(xv2, cp4, '-b.');  
legend([pd2 pd3 pd4 pd5],{'Region 2','Region 3','Region 4','Region 5'}, 'FontSize', 12); 
hold off;

%% functions

% function calculating the ratio total pressure to static pressure in an
% isentropic flow
function [p0p] = totaltostatic(M, gamma)
p0p = (1 + ((gamma-1)*M^2)/2)^(gamma/(gamma-1));
end

% function calculating the wave angle from given Mach and shock angle
function [B] = thetabetaM(M, theta, gamma)
options = optimset('Display','off');
theta = degtorad(theta);
f = @(b) 2*cot(b)*((M^2*(sin(b))^2-1)/(M^2*(gamma+cos(2*b))+2))-tan(theta); 
b = fsolve(f,theta,options); 
B = b/pi*180; 
end

% function calculating the Mach number or prandtl-meyer function value
% specified by machmu parameter
% 'mach' for mach number calculation from given function value (x)
% 'mu' for function value calculation from a given mach number (x)
function [out] = prandtlmeyer(x, gamma, machmu)
options = optimset('Display','off');
if strcmp(machmu,'mu')
        y = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(x^2-1)))-atand(sqrt(x^2-1));
end
if strcmp(machmu,'mach')
        f = @(M) sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atand(sqrt(M^2-1)) - x;
        y = fsolve(f, x, options);
end     
out = abs(y);
end

% function calculating the pressure coefficient
function [cp] = pressurecoeff(p1, p2, M, gamma)
cp = (2/(gamma*M^2))*((p2/p1)-1);
end