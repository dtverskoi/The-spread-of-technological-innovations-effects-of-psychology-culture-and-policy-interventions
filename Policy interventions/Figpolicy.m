%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots all figures related to intervention strategies using data
% files created by the program file Mainpolicy.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
XX=21;
Runs=100;
N=1000;


for exp=1:5
    
    if exp==1
        for pare=1:3
            ave=0.25*2^(pare-1);
            for parw=1:3
                avw=0.2+0.05*(parw-1);
                figure
                set(gcf, 'Position',  [200, 0, 1400, 1000])
                for tt=1:3
                    Q=readmatrix(['P' num2str(exp) num2str(pare) num2str(parw) '.txt']);
                    p5=NaN(XX,Runs);
                    y05=NaN(XX,Runs);
                    y15=NaN(XX,Runs);
                    for st=1:XX
                        p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                        y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                        y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                    end
                    
                    subplot(2,3,tt)
                    hold on
                    bmax=zeros(XX,1);
                    for st=1:XX
                        bmax(st)=0.1*(st-1);
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                    end
                    
                    plot(bmax,mean(p5,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    %xlabel('f_3','Fontsize',35)
                    %if tt==1
                        ylabel({'\fontsize{30}p'})
                    %end
                    title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                    
                    
                    subplot(2,3,tt+3)
                    hold on
                    for st=1:XX
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                    end
                    
                    plot(bmax,nanmean(y05,2),'LineWidth',6,'Color','g')
                    plot(bmax,nanmean(y15,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    xlabel('f_3','Fontsize',35)
                    %if tt==1
                        ylabel({'\fontsize{30}y'})
                    %end
                end
                print(['Ae' num2str(pare) 'w' num2str(parw)],'-dpng')
            end
        end
    end
    
    
    if exp==2
        parf3=2;
        for pare=1:3
            ave=0.25*2^(pare-1);
            for parw=1:3
                avw=0.2+0.05*(parw-1);
                figure
                set(gcf, 'Position',  [200, 0, 1400, 1000])
                for tt=1:3
                    Q=readmatrix(['P' num2str(exp) num2str(parf3) num2str(pare) num2str(parw) '.txt']);
                    p5=NaN(XX,Runs);
                    y05=NaN(XX,Runs);
                    y15=NaN(XX,Runs);
                    for st=1:XX
                        p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                        y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                        y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                    end
                    
                    subplot(2,3,tt)
                    hold on
                    bmax=zeros(XX,1);
                    for st=1:XX
                        bmax(st)=0.1*(st-1);
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                    end
                    
                    plot(bmax,mean(p5,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    %xlabel('f_3','Fontsize',35)
                    %if tt==1
                        ylabel({'\fontsize{30}p'})
                    %end
                    title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)

                    subplot(2,3,tt+3)
                    hold on
                    for st=1:XX
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                        scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                    end
                    
                    plot(bmax,nanmean(y05,2),'LineWidth',6,'Color','g')
                    plot(bmax,nanmean(y15,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    xlabel('f_2','Fontsize',35)
                    %if tt==1
                        ylabel({'\fontsize{30}y'})
                    %end
                end
                print(['Dis' num2str(parf3) num2str(pare) num2str(parw)],'-dpng')
            end
        end
    end
    
    
    if exp==3
        parf3=2;
        for pare=1:3
            ave=0.25*2^(pare-1);
            for parw=1:3
                avw=0.2+0.05*(parw-1);
                figure
                set(gcf, 'Position',  [200, 0, 1400, 1000])
                for tt=1:3
                    Q=readmatrix(['P' num2str(exp) num2str(parf3) num2str(pare) num2str(parw) '.txt']);
                    p5=NaN(XX,Runs);
                    y05=NaN(XX,Runs);
                    y15=NaN(XX,Runs);
                    for st=1:XX
                        p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                        y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                        y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                    end
                    
                    subplot(2,3,tt)
                    hold on
                    bmax=zeros(XX,1);
                    for st=1:XX
                        bmax(st)=0.1*(st-1);
                        scatter(bmax(st)*ones(1,Runs),p5(st,:),'MarkerEdgeColor','black','LineWidth',3)
                    end
                    
                    plot(bmax,mean(p5,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    %xlabel('f_3','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}p'})
                    %end
                    title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                    
                    
                    subplot(2,3,tt+3)
                    hold on
                    for st=1:XX
                        scatter(bmax(st)*ones(1,Runs),y15(st,:),'MarkerEdgeColor','red','LineWidth',1)
                        scatter(bmax(st)*ones(1,Runs),y05(st,:),'MarkerEdgeColor','blue','LineWidth',1)
                    end
                    
                    plot(bmax,nanmean(y05,2),'LineWidth',6,'Color','g')
                    plot(bmax,nanmean(y15,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 2]);
                    ylim([0 1])
                    xlabel('f_1','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}y'})
                    %end
                end
                print(['C' num2str(parf3) num2str(pare) num2str(parw)],'-dpng')
            end
        end
    end
    
    
    if exp==4
        parf3=2;
        for pare=1:3
            ave=0.25*2^(pare-1);
            for parw=1:3
                avw=0.25+0.05*(parw-1);
                figure
                set(gcf, 'Position',  [200, 0, 1400, 1000])
                for tt=1:3
                    Q=readmatrix(['P' num2str(exp) num2str(parf3) num2str(pare) num2str(parw) '.txt']);
                    p5=NaN(XX,Runs);
                    y05=NaN(XX,Runs);
                    y15=NaN(XX,Runs);
                    for st=1:XX
                        p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                        y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                        y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                    end
                    
                    subplot(2,3,tt)
                    hold on
                    bmax=zeros(XX,1);
                    for st=1:XX
                        bmax(st)=0.025*(st-1);
                        scatter(bmax(st)*(ones(1,Runs)+0.005*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                    end
                    
                    plot(bmax,mean(p5,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 0.5]);
                    ylim([0 1])
                    %xlabel('f_3','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}p'})
                    %end
                    title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                    
                    
                    subplot(2,3,tt+3)
                    hold on
                    for st=1:XX
                        scatter(bmax(st)*(ones(1,Runs)+0.005*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                        scatter(bmax(st)*(ones(1,Runs)+0.005*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                    end
                    
                    plot(bmax,nanmean(y05,2),'LineWidth',6,'Color','g')
                    plot(bmax,nanmean(y15,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 0.5]);
                    ylim([0 1])
                    xlabel('b_s','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}y'})
                    %end
                end
                print(['S2' num2str(pare) 'w' num2str(parw)],'-dpng')
            end
        end
    end
    
    
    if exp==5
        parf3=2;
        for pare=1:3
            ave=0.25*2^(pare-1);
            for parw=1:3
                avw=0.25+0.05*(parw-1);
                figure
                set(gcf, 'Position',  [200, 0, 1400, 1000])
                for tt=1:3
                    Q=readmatrix(['P' num2str(exp) num2str(parf3) num2str(pare) num2str(parw) '.txt']);
                    p5=NaN(XX,Runs);
                    y05=NaN(XX,Runs);
                    y15=NaN(XX,Runs);
                    for st=1:XX
                        p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                        y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                        y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                    end
                    
                    subplot(2,3,tt)
                    hold on
                    bmax=zeros(XX,1);
                    for st=1:XX
                        bmax(st)=0.005*(st-1);
                        scatter(bmax(st)*(ones(1,Runs)+0.0001*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                    end
                    
                    plot(bmax,mean(p5,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 0.1]);
                    ylim([0 1])
                    %xlabel('f_3','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}p'})
                    %end
                    title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                    
                    
                    subplot(2,3,tt+3)
                    hold on
                    for st=1:XX
                        scatter(bmax(st)*(ones(1,Runs)+0.0001*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                        scatter(bmax(st)*(ones(1,Runs)+0.0001*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                    end
                    
                    plot(bmax,nanmean(y05,2),'LineWidth',6,'Color','g')
                    plot(bmax,nanmean(y15,2),'LineWidth',6,'Color','g')
                    hold off
                    set(gca,'FontSize',25)
                    xlim([0 0.1]);
                    ylim([0 1])
                    xlabel('p_0','Fontsize',35)
                    %if tt==1
                    ylabel({'\fontsize{30}y'})
                    %end
                end
                print(['ip2' num2str(pare) 'w' num2str(parw)],'-dpng')
            end
        end
    end
    
end