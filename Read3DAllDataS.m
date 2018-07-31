clear all;
close all;


pConv = 1334;

suffix = 'm2';
name3D = '~/Documents/lab/KD-project/Artificial/RCA/m2/AllData' ; 

% Set overall parameters

AllData = load(name3D);
nRCR = 1;
nCOR_l = 6;
nCOR_r = 4;
nCycle = 4;


StartPlot = 2000;
EndPlot = 4000 ;

endStep_3D = EndPlot ; %size(AllData, 1);
singleCycle_3D = 1000 ; %endStep_3D/nCycle;

%Print Last cycle
%StartPlot = endStep_3D - singleCycle_3D ;
%EndPlot = endStep_3D ;

%% Do not need to modify below this line

nUnknowns = 10 + nRCR + 2*(nCOR_l + nCOR_r);
nAux = (nRCR + nCOR_l + nCOR_r) + 15;
nFaces = nRCR + nCOR_l + nCOR_r;
t_ind = nUnknowns + 1;
rcr_st = nUnknowns + 5;
lcor_st = nUnknowns + 5 + nRCR;
rcor_st = lcor_st + nCOR_l;


aortic_valve = find(AllData(StartPlot:EndPlot,nUnknowns + nFaces + 13) ==1) ;
aortic_valve_opening = [AllData( StartPlot + aortic_valve(1),t_ind), AllData( StartPlot+ aortic_valve(end), t_ind) ] - AllData( StartPlot,t_ind) ;

figures = [];

figures = [figures figure('position',[ 500 500 350 600]) ];
subplot(2,1,1)
plot(AllData(StartPlot:EndPlot,t_ind)-AllData(StartPlot, t_ind), -AllData(StartPlot:EndPlot, nUnknowns + nFaces + 5 ), 'r','linewidth',1.5);
hold on
plot(AllData(StartPlot:EndPlot,t_ind)-AllData(StartPlot, t_ind), AllData(StartPlot:EndPlot, 9), 'k','linewidth',1.5);
xlabel('t [s]')
ylabel('Flow [ml/s]')
legend('Inlet flow to 3D model','Ejected flow from LV');
hold on
y1=get(gca,'ylim');
ylim(y1)
line(aortic_valve_opening(1)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
line(aortic_valve_opening(end)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
ax=gca;
ax.FontSize = 12;
hold off
subplot(2,1,2)
plot(AllData(StartPlot:EndPlot, t_ind)-AllData(StartPlot, t_ind), AllData(StartPlot:EndPlot, 10), 'b','linewidth',1.5);
hold on
plot(AllData(StartPlot:EndPlot, t_ind)-AllData(StartPlot, t_ind), AllData(StartPlot:EndPlot, nUnknowns + nFaces + 9), 'k','linewidth',1.5);
legend('Aortic pressure' , 'LV pressure')
xlabel('t [s]')
ylabel('Pressure [mmHg]')
hold on
y1=get(gca,'ylim');
ylim(y1)
line(aortic_valve_opening(1)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
line(aortic_valve_opening(end)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
ax=gca;
ax.FontSize = 12;


figures = [figures figure];
plot(AllData(StartPlot:EndPlot,8),AllData(StartPlot:EndPlot,nUnknowns + nFaces + 9),'b','linewidth',1.5)
title('Left ventricle - PV Loop')
xlabel('Volume [ml]')
ylabel('Pressure [mmHg]')
ax=gca;
ax.FontSize = 12;

figures = [figures figure];
plot(AllData(StartPlot:EndPlot,6),AllData(StartPlot:EndPlot,nUnknowns + nFaces + 8))
title('Right ventricle - PV Loop')

figures = [figures figure];
subplot(2,1,1)
plot(AllData(:,t_ind), AllData(:,nUnknowns  + 5 :  lcor_st - 1 ));
title('Aortic Outlets Flow')
subplot(2,1,2)
plot(AllData(:,t_ind), AllData(:, lcor_st :nUnknowns + 4 + nFaces));
title('Coronary Outlets Flow')


figures = [figures figure('position',[ 500 500 350 700])];
subplot(3,1,1)
plot(AllData(StartPlot:EndPlot,t_ind) - AllData(StartPlot,t_ind), AllData(StartPlot:EndPlot,10),'b','linewidth',1.5);
title('Aortic pressure')
xlabel('t [s]')
ylabel('Pressure [mmHg]')
hold on
y1=get(gca,'ylim');
ylim(y1)
line(aortic_valve_opening(1)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
line(aortic_valve_opening(end)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
ax=gca;
ax.FontSize = 12;
subplot(3,1,2)
plot(AllData(StartPlot:EndPlot,t_ind) - AllData(StartPlot,t_ind) , AllData(StartPlot:EndPlot,lcor_st :rcor_st-1 ),'color',[0 0.65 0.25],'linewidth',1.5);
title('Left Coronaries')
xlabel('t [s]')
ylabel('Flow [ml/s]')
y1=get(gca,'ylim');
ylim(y1)
line(aortic_valve_opening(1)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
line(aortic_valve_opening(end)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
ax=gca;
ax.FontSize = 12;
subplot(3,1,3)
plot(AllData(StartPlot:EndPlot,t_ind) - AllData(StartPlot,t_ind), AllData(StartPlot:EndPlot,rcor_st :nUnknowns + 4 + nFaces),'color',[0 0.65 0.25],'linewidth',1.5);
title('Right Coronaries')
xlabel('t [s]')
ylabel('Flow [ml/s]')
y1=get(gca,'ylim');
ylim(y1)
line(aortic_valve_opening(1)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
line(aortic_valve_opening(end)*ones(1,2), y1,'LineStyle','--','LineWidth',0.75,'color',[0,0,0])
ax=gca;
ax.FontSize = 12;

% save the produced figures to a set location 
for f = 1:numel(figures)
    fig = figures(f);
    filename = sprintf('Artificial/RCA/%s/Figure%02d_%s.pdf', suffix, f, suffix);
    print( fig, '-dpng', filename );
end


%Calculate Resuts for the last cycle

Q_rcr3D = 0;
Q_lcor3D = 0;
Q_rcor3D = 0;

    % SUM RCR FLUX
    for i=0:1:nRCR-1
      temp = trapz(AllData(endStep_3D-singleCycle_3D:endStep_3D,t_ind), ... 
          AllData(endStep_3D-singleCycle_3D:endStep_3D,rcr_st+i));
      Q_rcr3D = Q_rcr3D + temp; 
    end

    % SUM LEFT CORONARY FLUX
    for i=0:1:nCOR_l-1
      temp=trapz(AllData(endStep_3D-singleCycle_3D:endStep_3D,t_ind), ... 
          AllData(endStep_3D-singleCycle_3D:endStep_3D,lcor_st+i));
      Q_lcor3D = Q_lcor3D + temp; 
    end
    
    % SUM RIGHT CORONARY FLUX
    for i=0:1:nCOR_r-1
      temp=trapz(AllData(endStep_3D-singleCycle_3D:endStep_3D,t_ind), ... 
          AllData(endStep_3D-singleCycle_3D:endStep_3D,rcor_st+i));
      Q_rcor3D=Q_rcor3D + temp; 
    end
    
    Total_Q_out = Q_rcor3D + Q_lcor3D + Q_rcr3D;
    
    % Only calculate certain 3D results
    Qinlet3D = trapz(AllData(endStep_3D-singleCycle_3D:endStep_3D,t_ind), ... 
        AllData(endStep_3D-singleCycle_3D:endStep_3D,nUnknowns + nFaces + 5));
    SV = trapz(AllData(endStep_3D-singleCycle_3D:endStep_3D,t_ind), ... 
        AllData(endStep_3D-singleCycle_3D:endStep_3D,9));
    Aor_Cor_split3D = ((Q_lcor3D+Q_rcor3D)/(Q_lcor3D+Q_rcor3D+Q_rcr3D))*100.0;
    L_R_corsplit3D = (Q_lcor3D/(Q_lcor3D+Q_rcor3D))*100.0;
    Pao_max3D = max(AllData(endStep_3D-singleCycle_3D:endStep_3D,10));
    Pao_min3D = min(AllData(endStep_3D-singleCycle_3D:endStep_3D,10));
    %Pao_mean3D = mean(AllData(endStep_3D-singleCycle_3D:endStep_3D,10));
    
    fprintf('\n --- 3D RESULTS for %s --- \n', suffix);
    fprintf('Qinlet = %8.3f ml/cycle\n',abs(Qinlet3D));
    fprintf('SV = %8.3f ml \n',SV);
    fprintf('Total Q out = %8.3f ml/cycle\n',Total_Q_out);
    fprintf('Ao-Cor-Split = %8.3f %%\n',Aor_Cor_split3D);
    fprintf('L_R_Cor_Split = %8.3f %%\n',L_R_corsplit3D);
    fprintf('Pao_max = %8.3f mmHg\n',Pao_max3D);
    fprintf('Pao_min = %8.3f mmHg\n',Pao_min3D);
    %fprintf('Pao_mean = %8.3f mmHg\n',Pao_mean3D);


