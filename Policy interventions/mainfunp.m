%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is a Matlab implementation of the model with a material
% support for early adopters. It authomatically simulates the model 100 
% runs and shows the results on an individual choice x, attitude y, and 
% benefit of using a new technology b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q] = mainfunp(b0,bmin,bmax,avw,ave,s,f0,f1,f2,f3,nu,sw,se,mb,sb,mu,a0,sa1,sub)

%Parameters
N=1000;                            % number of individuals
T=200;                             % time steps
yG=1;                              % attitudes of the authority
p0=0;                              % initial number of adopters
Runs=100;                          % number of runs

%Output
Q=zeros(15,Runs);

for run=1:Runs
    
    %Variables
    b=NaN(N,T);                    % material benefits of the new technology
    x=NaN(N,T);                    % actions
    y=NaN(N,T);                    % attitudes
    p=NaN(1,T);                    % frequency of adopters
    
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
    a1=0.05+sa1*randn(N,1);             %learning slope
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
    ft=zeros(1,N);
    
    %Simulations
    for t=1:T-1
        % The dynamics of the individual benefit
        for i=1:N
            b(i,t+1)=b(i,t)-c*(1-x(i,t))*(b(i,t)-bmin)+x(i,t)*(a0+a1(i)*p(t))*(bmax-b(i,t));
        end
        avy=sum(y(:,t))/N;
        % Utility function and decision-making
        for i=1:N
            z=rand;
            if z<nu
                if (ft(i)==0)&&(t<=5)
                    dU=(1-w(i))*b(i,t+1)+w(i)*bmax-b0+e(i)*(f0*k(i,1)*(2*y(i,t)-1)+f1*k(i,2)*(2*p(t)-1)+f2*k(i,3)*(2*avy-1)+f3*k(i,4)*(2*yG-1))+sub;
                else
                    dU=(1-w(i))*b(i,t+1)+w(i)*bmax-b0+e(i)*(f0*k(i,1)*(2*y(i,t)-1)+f1*k(i,2)*(2*p(t)-1)+f2*k(i,3)*(2*avy-1)+f3*k(i,4)*(2*yG-1));
                end
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
        ft(x(:,t+1)==1)=1;
        p(t+1)=sum(x(:,t+1))/N;
        % The dynamics of attitudes
        for i=1:N
            y(i,t+1)=y(i,t)+s*(f0*beta(i,1)*(x(i,t+1)-y(i,t))+f1*beta(i,2)*(p(t+1)-y(i,t))+f2*beta(i,3)*(avy-y(i,t))+f3*beta(i,4)*(yG-y(i,t)));
        end
    end
    
    % Output
    I1=find(x(:,50)==0);
    I2=find(x(:,50)==1);
    I3=find(x(:,100)==0);
    I4=find(x(:,100)==1);
    I5=find(x(:,200)==0);
    I6=find(x(:,200)==1);
    Q(1,run)=mean(b(I1,50));
    Q(2,run)=mean(b(I2,50));
    Q(3,run)=mean(y(I1,50));
    Q(4,run)=mean(y(I2,50));
    Q(5,run)=size(I2,1);
    Q(6,run)=mean(b(I3,100));
    Q(7,run)=mean(b(I4,100));
    Q(8,run)=mean(y(I3,100));
    Q(9,run)=mean(y(I4,100));
    Q(10,run)=size(I4,1);
    Q(11,run)=mean(b(I5,200));
    Q(12,run)=mean(b(I6,200));
    Q(13,run)=mean(y(I5,200));
    Q(14,run)=mean(y(I6,200));
    Q(15,run)=size(I6,1);
end

end
