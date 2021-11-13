%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs simulations that show the effects of psychological 
% factors on the spread of a new technology:
% - experiment 1: effects of the strength of normative factors in the
% utility function;
% -experiment 2: effects of the strength of cognitive dissonance.
% This code uses mainfun.m for creating data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Parameters
bmin=0.5;                   % minimum benefit of using new technology
b0=1;                       % individual benefit of using old technology
bmax=1.5;                   % maximum benefit of using new technology
avw=0.25;                   % mean value of the foresight parameter
ave=0.5;                    % mean strength of normative factors in the utility function
s=0.1;                      % the speed of change in attitudes
f0=1;                       % parameter scaling the effects of cognitive dissonance
f1=1;                       % parameter scaling the effects of visibility of peers' actions
f2=1;                       % parameter scaling the effects of visibility of peers' attitudes
nu=0.1;                     % update probability
sw=0.025;                   % std of the foresight parameter
se=0.05;                    % std of the strength of normative factors
mb=0.1;                     % mean initial attitudes
sb=0.2;                     % std of initial attitudes
mu=15;                      % precision parameter
a0=0.05;                    % learning parameter
sa1=0.005;                  % std of learning parameter
N=1000;                     % number of individuals
T=200;                      % number of time steps
p0=0;                       % initial frequency of adopters
Runs=100;                   % number of runs
XX=21;                      % number of different values of a varying parameter in each experiment

% Experiments
exp=1;
if exp==1
    for parf3=1:3
        f3=0.25*2^(parf3-1);
        for parw=1:3
            avw=0.2+0.05*(parw-1);
            Out=NaN(XX*15,Runs);
            
            for step=1:XX
                ave=0.05*(step-1);
                [Q] = mainfun(b0,bmin,bmax,avw,ave,s,f0,f1,f2,f3,nu,sw,se,mb,sb,mu,a0,sa1);
                Out((step-1)*15+1:step*15,:)=Q(:,:);
            end
            dlmwrite(['Ps' num2str(exp) num2str(parf3) num2str(parw) '.txt'],Out)
        end
    end
end

if exp==2
    f3=1;
    for pare=1:3
        ave=0.25*2^(pare-1);
        for parw=1:3
            avw=0.2+0.05*(parw-1);
            Out=NaN(XX*15,Runs);
            
            for step=1:XX
                f0=0.1*(step-1);
                [Q] = mainfun(b0,bmin,bmax,avw,ave,s,f0,f1,f2,f3,nu,sw,se,mb,sb,mu,a0,sa1);
                Out((step-1)*15+1:step*15,:)=Q(:,:);
            end
            dlmwrite(['Ps' num2str(exp) num2str(pare) num2str(parw) '.txt'],Out)
        end
    end
end