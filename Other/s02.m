%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs simulations that show effects of attitudes on the
% adoption curves and produces Figure 6 in the main text.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

bmin=0.5;              % minimum benefit of using new technology
b0=1;                  % individual benefit of using old technology 
ave=0.5;               % mean strength of normative factors in the utility function
f0=1;                  % parameter scaling the effects of cognitive dissonance
f1=1;                  % parameter scaling the effects of visibility of peers' actions
f2=1;                  % parameter scaling the effects of visibility of peers' attitudes
nu=0.1;                % update probability
sw=0.025;              % std of the foresight parameter
se=0.05;               % std ot the strength of normative factors
mb=0.5;                % mean initial attitudes
sb=0.2;                % std of initial attitudes
mu=15;                 % precision parameter
a0=0.05;               % learning parameter
sa1=0.005;             % std of learning parameter
N=1000;                % number of individuals
T=200;                 % number of time steps
yG=1;                  % attitudes of the authority
p0=0;                  % initial number of adopters
bmax=1.5;              % maximum benefit of using new technology 
f3=1;                  % parameter scaling the effects of authority's effort
avw=0.25;              % mean value of the foresight parameter  

seedStat=0;
%Seed
if seedStat==0
    rng('shuffle','twister')
    seed=rng;
    %  save seed;
else
    %  load seed;
    if exist('seed')
        rng(seed);
    else
        rng('shuffle','twister');
        seed=rng;
    end
end

% Generate individual parameters
% a0,a1,c - in the individual learning rate
a1=0.05+sa1*randn(N,1);            %learning slope
a1(a1<0)=0;
a1(a1>1)=1;
c=0.1;                              %knowledge depreciation rate

% foresight parameter w
w=avw+sw*randn(N,1);                %foresight parameter
w(w<0)=0;
w(w>1)=1;

% strength of normative factors e
e=ave+se*randn(N,1);
e(e<0)=0;

%Initial conditions
b(:,1)=bmin*ones(N,1);
varb=sb^2;
ab=mb*(mb*(1-mb)/varb-1);
bb=(1-mb)*(mb*(1-mb)/varb-1);
y(:,1)=betarnd(ab,bb,N,1);
x(:,1)=zeros(N,1);
p(1)=sum(x(:,1))/N;
x(:,1)=zeros(N,1);

% k and beta - parameters in the utility function and attitudes
k=NaN(N,4);
beta=NaN(N,4);
for i=1:N
    k0=sort(rand(3,1));
    beta0=sort(rand(3,1));
    kvs=sort([k0(1); diff(k0); 1-k0(end)]);
    kvs=kvs(randperm(length(kvs)));
    k(i,:)=kvs';
    betavs=sort([beta0(1); diff(beta0); 1-beta0(end)]);
    betavs=betavs(randperm(length(betavs)));
    beta(i,:)=betavs';
end


figure
set(gcf, 'Position',  [100, 100, 800, 400])

% loop fot 3 cases:
% - case 1: dynamically changing attitudes;
% - case 2: fixed attitudes;
% - case 3: no attitudes.
for run=1:3
    if run==1
        s=1;
    else
        if run==2
            s=0;
        else
            s=0;
            for i=1:N
                k(i,:)=[0 k(i,2)/(k(i,2)+k(i,4)) 0 k(i,4)/(k(i,2)+k(i,4))];
                beta(i,:)=[0 beta(i,2)/(beta(i,2)+beta(i,4)) 0 beta(i,4)/(beta(i,2)+beta(i,4))];
            end
        end
    end
    
    
    
    %Simulations
    for t=1:T-1
        % The dynamics of the individual benefit
        for i=1:N
            b(i,t+1)=b(i,t)-c*(1-x(i,t))*(b(i,t)-bmin)+x(i,t)*(a0+a1(i)*p(t))*(bmax-b(i,t));
        end
        % Utility function and decision-making
        avy=sum(y(:,t))/N;
        for i=1:N
            z=rand;
            if z<nu
                dU=(1-w(i))*b(i,t+1)+w(i)*bmax-b0+e(i)*(f0*k(i,1)*(2*y(i,t)-1)+f1*k(i,2)*(2*p(t)-1)+f2*k(i,3)*(2*avy-1)+f3*k(i,4)*(2*yG-1));
                lU=1/(1+exp(mu*dU));
                if rand<lU
                    x(i,t+1)=0;
                else
                    x(i,t+1)=1;
                end
            else
                x(i,t+1)=x(i,t);
            end
        end
        p(t+1)=sum(x(:,t+1))/N;
        % The dynamics of attitudes
        for i=1:N
            y(i,t+1)=y(i,t)+s*(f0*beta(i,1)*(x(i,t+1)-y(i,t))+f1*beta(i,2)*(p(t+1)-y(i,t))+f2*beta(i,3)*(avy-y(i,t))+f3*beta(i,4)*(yG-y(i,t)));
        end
    end
    
    if run==1
        cl='black';
        wd=6;
    else
        if run==2
            cl='blue';
            wd=3;
        else
            cl='c';
            wd=3;
        end
    end
    
    %Figure section
    plot(linspace(1,T,T),p,'LineWidth',wd,'Color',cl)
    hold on
    set(gca,'FontSize',20)
    xlim([0 T])
    ylim([0,1])
    ylabel('p','Fontsize',25)
    xlabel('t','FontSize',25);
end
hold off
legend('changing attitudes','fixed attitudes','no attitudes','Fontsize',18,'Location','se')
print('At','-dpng')