% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 80; %80; default
j2 = 0.05;%0.05; default
j3 = 0.5;
j4 = 0.16;
j5 = 0:0.05:5;
j6 = 1;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
ta5 = 1;
% 
% timelength = 1000;
% strength = 125;
% stimFun = @(t) 0.003 + (t>=0).*strength - strength.*(t>30); % basal value
%         %of 0.003. Strength of square wave between time 0 and 30
% strength2 = 10;
% stimFun2 = @(t) (t>=300).*strength2.*subplus(sin(t./60)) - (t>=4000).*strength2.*subplus(sin(t./60));%.*subplus(sin(t./60));
% 
% initialGuess = ones(1,3).*1e-4;
% 
% %ODE solution for simple model. 
% [T,Y] = ode45(@(t,y) neuronFireODE(t,y,stimFun,j1,j2,j3,j4,j5,ta1,ta2,...
%     ta3,ta4,ta5), [0,timelength], initialGuess);
% disp('Trying to solve new ODE system');
% [T2,Y2] = ode45(@(t,y) neuronFireODENewTerm(t,y,2,stimFun,stimFun2,j1,j2,j3,j4,j5,j6,...
%     ta1,ta2,ta3,ta4,ta5), [0, timelength], [initialGuess 0]);
% 
% figure()
% plot(T, Y(:,1),'-.', 'LineWidth', 1.2);hold on;
% plot(T2, Y2(:,1));
% legend('Paper Model', 'OurModel');
% xlabel('Time (m)'); ylabel('[PKM]');



% make sure to set this to the parameter being varied. For plotting
% purposes
varying = j5;

% compute values and determine stability
[pkm, actin, rna, hs, stability, indexes] = computeSSNew(2,j1,j2,j3,j4,j5,j6);

stable = [];
stableIndex = [];
unstable = [];
unstableIndex = [];
for j = 1:length(pkm)
    for i = 1:length(pkm{j})
        if stability{j}(i) > 0
            stable = [stable, pkm{j}(i)];
            stableIndex = [stableIndex, varying(indexes{j}(i))];
        else
            unstable = [unstable, pkm{j}(i)];
            unstableIndex = [unstableIndex, varying(indexes{j}(i))];
        end
    end
end

figure()
scatter(stableIndex, stable, '.'); hold on;
scatter(unstableIndex, unstable, '.');
legend('Stable', 'Unstable');
ylabel('[PKM]'); xlabel('j');

