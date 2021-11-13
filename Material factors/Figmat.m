%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code visualizes the effects of material factors on the spread of a
% new technology using data files produced by Mainmat.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
XX=21;
Runs=100;
N=1000;

for exp=1:3
    
    if exp==1
        for parf3=2:2
            figure
            set(gcf, 'Position',  [200, 0, 500, 1000])
            for tt=2:2
                Q=readmatrix(['NEx' num2str(exp) num2str(parf3*10) '.txt']);
                p5=NaN(XX,Runs);
                y05=NaN(XX,Runs);
                y15=NaN(XX,Runs);
                for st=1:XX
                    p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                    y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                    y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                end
                
                subplot(2,1,1)
                hold on
                bmax=zeros(XX,1);
                for st=1:XX
                    bmax(st)=0.0+0.025*(st-1);
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                end
                
                plot(bmax,mean(p5,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([0.0 0.5])
                ylim([0 1])
                xlabel('$\bar{\omega}$','Interpreter','LaTex','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}p'})
                end
                %title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                
                
                subplot(2,1,2)
                hold on
                bmax=zeros(XX,1);
                for st=1:XX
                    bmax(st)=0.0+0.025*(st-1);
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                end
                
                plot(bmax,nanmean(y05,2),'LineWidth',4,'Color','g')
                plot(bmax,nanmean(y15,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([0.0 0.5])
                ylim([0 1])
                xlabel('$\bar{\omega}$','Interpreter','LaTex','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}y'})
                end
            end
            print('N2a','-dpng')
        end
    end
    
    
    if exp==2
        for parf3=2:2
            figure
            set(gcf, 'Position',  [200, 0, 500, 1000])
            for tt=2:2
                Runs=100;
                Q=readmatrix(['NEx' num2str(exp) num2str(parf3*10) '.txt']);
                p5=NaN(XX,Runs);
                y05=NaN(XX,Runs);
                y15=NaN(XX,Runs);
                for st=1:XX
                    p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                    y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                    y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                end
                
                subplot(2,1,1)
                hold on
                bmax=zeros(XX,1);
                for st=1:XX
                    bmax(st)=1+0.05*(st-1);
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                end
                
                plot(bmax,mean(p5,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([1 2])
                ylim([0 1])
                xlabel('b_{max}','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}p'})
                end
                %title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                
                
                subplot(2,1,2)
                hold on
                for st=1:XX
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                end
                
                plot(bmax,nanmean(y05,2),'LineWidth',4,'Color','g')
                plot(bmax,nanmean(y15,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([1 2])
                ylim([0 1])
                xlabel('b_{max}','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}y'})
                end
            end
            print('N2c','-dpng')
        end
    end
    
    
    if exp==3
        for parf3=2:2
            figure
            set(gcf, 'Position',  [200, 0, 500, 1000])
            for tt=2:2
                Q=readmatrix(['NEx' num2str(exp) num2str(parf3*10) '.txt']);
                p5=NaN(XX,Runs);
                y05=NaN(XX,Runs);
                y15=NaN(XX,Runs);
                for st=1:XX
                    p5(st,:)=Q((st-1)*15+5*tt,:)/N;
                    y05(st,:)=Q((st-1)*15+5*(tt-1)+3,:);
                    y15(st,:)=Q((st-1)*15+5*(tt-1)+4,:);
                end
                
                subplot(2,1,1)
                hold on
                bmax=zeros(XX,1);
                for st=1:XX
                    bmax(st)=0.05*(st-1);
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),p5(st,:),'MarkerEdgeColor','black','LineWidth',0.1)
                end
                
                plot(bmax,mean(p5,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([0 1])
                ylim([0 1])
                xlabel('b_{min}','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}p'})
                end
                %title(['T=', num2str(tt*50+25*(tt-1)*(tt-2))],'Fontsize',35)
                
                
                subplot(2,1,2)
                hold on
                for st=1:XX
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y15(st,:),'MarkerEdgeColor','red','LineWidth',0.1)
                    scatter(bmax(st)*(ones(1,Runs)+0.02*rand(1,Runs)),y05(st,:),'MarkerEdgeColor','blue','LineWidth',0.1)
                end
                
                plot(bmax,nanmean(y05,2),'LineWidth',4,'Color','g')
                plot(bmax,nanmean(y15,2),'LineWidth',4,'Color','g')
                hold off
                set(gca,'FontSize',25)
                xlim([0 1])
                ylim([0 1])
                xlabel('b_{min}','Fontsize',35)
                if tt==2
                    ylabel({'\fontsize{30}y'})
                end
            end
            print('N2e','-dpng')
        end
    end
    
end
