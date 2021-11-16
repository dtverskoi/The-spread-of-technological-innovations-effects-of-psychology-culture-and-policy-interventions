%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs simulations that show values of K3 for initial 
% adopters, and produces Figure 5b in the main text.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Parameters
bmin=0.5;             % minimum benefit of using new technology
b0=1;                 % individual benefit of using old technology
ave=0.5;              % mean strength of normative factors in the utility function
s=0.1;                % the speed of change in attitudes
f0=1;                 % parameter scaling the effects of cognitive dissonance
f1=1;                 % parameter scaling the effects of visibility of peers' actions
f2=1;                 % parameter scaling the effects of visibility of peers' attitudes
nu=0.1;               % update probability
sw=0.025;             % std of the foresight parameter
se=0.05;              % std ot the strength of normative factors
mb=0.1;               % mean initial attitudes
sb=0.2;               % std of initial attitudes
mu=15;                % precision parameter
a0=0.05;              % learning parameter
sa1=0.005;            % std of learning parameter
N=1000;               % number of individuals
T=500;                % number of time steps
yG=1;                 % attitudes of the authority
p0=0;                 % initial frequency of adopters
bmax=1.5;             % maximum benefit of using new technology 
f3=1;                 % parameter scaling the effects of authority's effort
avw=0.25;             % mean value of the foresight parameter  

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

[a_sorted, a_order] = sort(k(:,4)+beta(:,4));
k = k(a_order,:);
beta=beta(a_order,:);

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
[Max,I]=maxk(y(:,1),p0);
x(I,1)=1;
p(1)=sum(x(:,1))/N;

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
    % The dynamics of attitudes
    p(t+1)=sum(x(:,t+1))/N;
    for i=1:N
        y(i,t+1)=y(i,t)+s*(f0*beta(i,1)*(x(i,t+1)-y(i,t))+f1*beta(i,2)*(p(t+1)-y(i,t))+f2*beta(i,3)*(avy-y(i,t))+f3*beta(i,4)*(yG-y(i,t)));
    end
end


%Figure section
figure
set(gcf, 'Position',  [100, 100, 850, 850])

subplot(3,1,1)
plot(linspace(1,T,T),p,'LineWidth',4)
set(gca,'FontSize',20)
xlim([0 T])
ylim([0,1])
ylabel('p','Fontsize',25)
xlabel('t','FontSize',25);

subplot(3,1,2)
hold on
ZZ=(k(:,4)+beta(:,4))/2;
for i=1:N
    plot(linspace(1,T,T),y(i,:),'LineWidth',2,'color',[ZZ(i) 0.2 1-ZZ(i)])
end
hold off
mycolors=NaN(1,3);
for i=1:N
mycolors = [mycolors; [ZZ(i) 0.2 1-ZZ(i)]];
end
colormap(mycolors);
c = colorbar;
c.Label.String = 'K_3';
set(gca,'FontSize',20)
xlim([0 T])
ylim([0,1])
ylabel('y','Fontsize',25)
xlabel('t','FontSize',25);

subplot(3,1,3)
hold on
for i=1:N
    plot(linspace(1,T,T),b(i,:),'LineWidth',2,'color',[ZZ(i) 0.2 1-ZZ(i)])
end
hold off
colormap(mycolors);
c = colorbar;
c.Label.String = 'K_3';
set(gca,'FontSize',20)
xlim([0 T])
ylim([bmin,bmax])
ylabel('b','Fontsize',25)
xlabel('t','FontSize',25);

print('Sw3','-dpng')