% Ashton Gray
% 11/9/2019
% ECE 202 Fall 2019
% MATLAB Project 2
% Phase 4, Exploring Results

clear
clf

% ----- User Input for C -----

C = input("Enter Value for C (between 0.2 and 0.3): "); % Dimensionless constant

% ----- Define Given Information -----

x0 = 0; y0 = 0; % start
v0mph = 112; % exit velocity, in mph
phi0deg = 32; % launch angle, in deg
g = 10; % % gravitational constant, in m/s^2
m = 0.145;  % mass of baseball, in kg

pair = 1.225; % density of air, in kg/m^3
A = pi * (0.074/2)^2; % Cross sectional areal of a baseball, (pi*radius^2), in m^2

% ----- Useful Variables -----

deg2rad = pi()/180; % Convert from degrees to radians
mph2mps = 1609/3600; % Convert from mph to feet per second
mps2mph = mph2mps^-1; % Convert from meters per second to miles per hour
ft2m = 0.3048; % Convert feet to meters, used for calculating
m2ft = 3.28; % Convert meters back to feet, used for plotting

phi0 = phi0deg * deg2rad; %launch angle in rad (no units in variable name)
v0 = v0mph * mph2mps; % initial velocity in mps

v0x = v0 * cos(phi0); % x-component of v0
v0y = v0 * sin(phi0); % y-component of v0

% ----- characteristics of trajectory -----

tH = v0y/g; % tome to reach max. height H, in s
tland = 2*tH; % time to land (assuming flat ground), in s

H = v0y^2/2/g * m2ft; % max. height, in ft
R = v0x * tland * m2ft; % range of baseball (where it lands) in ft

% ----- analytical solution -----

tmin = 0; tmax = tland; % stop when ball lands
N = 2000; % intervals
t = linspace(tmin,tmax,1+N);

xt = x0 + v0x*t; % x(t), ax = 0
yt = y0 + v0y*t - .5*g*t.^2; % y(t), ay = -g (no drag)

% ----- numerical solution -----

dt  = (tmax-tmin)/N;
y = zeros(1, 1+N); % initialize y(t)
x = zeros(1, 1+N); % initialize x(t)
y(1) = y0; % set first value
x(1) = x0; % set first value
vy = v0y; % initialize vy
vx = v0x; % initialize vx

fcf = (0.5)*C*pair*A; % drag force common factor

for n = 1:N
    
    v = sqrt(vx^2+vy^2); % changing v
    Fnetx = (-1)*fcf*v*vx; % get force in x and y direction
    Fnety = (-1)*fcf*v*vy -m*g;
    
    ax = Fnetx/m; % get acceleration from force
    ay = Fnety/m;
    
    x(n+1) = x(n) + vx*dt + .5*ax*dt^2; % get x and y position
    y(n+1) = y(n) + vy*dt + .5*ay*dt^2;
    
    vy = vy + ay*dt; % update velocity
    vx = vx + ax*dt; 
    
    if y(n)/y(n+1) < 0
        nFR = n; % used to calculate Flight and Range
        vF = v; % final speed, before ball hits ground
    end
end

% ----- plot trajectory -----

% Convert all units back into feet from meters
x = x*m2ft;
y = y*m2ft;
xt = xt*m2ft; 
yt = yt*m2ft;

plot(xt,yt,x,y,'LineWidth',2)
legend({"Without Drag (C=0)",...
    "With Drag (C = "+C+")"}, 'Location', 'northeastoutside', 'FontSize',18)
xlabel('x (ft)','FontSize',20)
ylabel('y (ft)','FontSize',20)
grid on
title({'ECE 202, Project 2, Phase 4, Exploring Results:',...
    'Trajectory of a Baseball, With and Without Drag'},'FontSize',25)
ylim([0 120])
set(gca, 'FontSize', 16)

% ----- Check Analytical vs Numerical Solution -----

checkSumX = sum(abs(x-xt)); % should return 0
checkSumY = sum(abs(y-yt)); % should return 0

% ----- Trajectory Traits -----

% From Graph
maxHeight_ft = max(y) % max height, in ft
tFlight_s = t(nFR) % flight time, in s
range_ft = x(nFR) % range of flight, in ft
speedf_mps = vF; % final speed of ball, in m/s
eLoss_J = 0.5*m*(speedf_mps^2-v0^2) % energy lost, in J
speedf_mph = speedf_mps * mps2mph % final speed, in mph

% Given Traits
givenHeight = 114; % given max height, in ft
givenRange = 446; % given range, in ft
giventFlight = 5.7; % given flight time, in s

% percent differences
maxHeight_pDiff = abs(maxHeight_ft-H)/(maxHeight_ft+H)/2*100
range_pDiff = abs(range_ft-R)/(range_ft+R)/2*100
tFlight_pDiff = abs(tFlight_s-tland)/(tFlight_s+tland)/2*100

% Question c)
% At a % difference of 7.1345 for max height, and a % difference of 3.7057
% for flight time, this is a significant difference. This large contrast is
% due to the fact it is hard to truly calculate the effect of air
% resistance and drag of anything in a real world setting using basic
% calculations.