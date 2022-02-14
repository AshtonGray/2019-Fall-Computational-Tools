% Ashton Gray
% 11/9/2019
% ECE 202 Fall 2019
% MATLAB Project 1
% Phase 6, Understanding the Taylor Series

clear
clf
format shortG

% ----- Define Variables and Inputs -----
w = 40; % omega, rad/s
A = 12; % amplitude
n0 = input('Enter the number of non zero terms to plot (positve integer): '); % number of nonzero terms

ymin = -1.2*A; ymax = 1.2*A; % min and max y value from amplitude
tmin = 0; 
tmax = input("Enter the maximum time value (in s): ");
x = 400; % number of points to plot

% Get ready for calculations
n = 0:2:2*(n0-1); % even numbers only (odd will equal to zero)
t = linspace(tmin,tmax,x+1); % t, in seconds
tms = t*1000; % t, in ms

% Define coefficents from handwork
an = A*(-1).^(n/2).*w.^(n)./factorial(n);
CoefficientTbl = transpose([n; an])

% ----- Calculations -----

% OLD METHOD: Get each function, add to previous because of Summation
f = an(1)*t.^n(1); 
f1 = f + an(2)*t.^n(2);
f2 = f1 + an(3)*t.^n(3);
f3 = f2 + an(4)*t.^n(4);
f4 = f3 + an(5)*t.^n(5);
f5 = f4 + an(6)*t.^n(6);

% Use for loop to create and plot functions
subtotal = zeros(1,x+1);
hold on
for i = 1:n0
    subtotal = subtotal + an(i)*t.^n(i); % create next
    if i<n0
        plot(tms,subtotal,'LineWidth',2) % plot
        
    end
end
plot(tms,subtotal,'LineWidth',5) % plot last line thicker
hold off

% Check old method vs new
checkf = sum(abs(f5-subtotal)) % should be 0

% ----- Graph Labeling -----

% Title and Legend
title({"ECE 202, Project 1, Phase 6, Understanding the Taylor Series, Expanding "+A+"cos("+w+"t):"...
    "Power Expansion Using the First "+n0+" nonzero Terms"...
    "Of Truncated Series About t = 0"},"FontSize",24)
legend("Up to n = " + n,'Location','northeastoutside')

% Axis stuff
ax = gca;
ax.FontSize = 15;
xlabel("time t (ms)",'FontSize', 20)
ylabel("f(t)", 'FontSize', 20)
ylim([ymin ymax])

% ----- Questions -----

% a) at n = 10 nonzero terms, the last function looks like the given
% function (12cos(40t)) for the first 200ms.

% b) During the next 200 ms, the last function does not show up and instead
% goes to negative infinity. I believe this happens because the series
% needs more terms in order to create a function that will continue as t
% approaches infinity the same the real function does.

% c) At n = 21 nonzero terms, the last function "looks like" the given
% function (12cos(40t)) with an ending time of 400ms.