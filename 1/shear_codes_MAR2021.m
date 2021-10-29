clear all;
close all;
clc;
% --- Written by Austin Weir --- March 2021
% ---
% ---       P2
% ---   P3      P1   ^+y  
% ---       P4       <+x 
% ---  
% ---   You stand here
% ---
% ---
% --------------  (Variables to EDIT) -----------------------------
delay = 1;     % For the phase delay between ATI and Capacitance
N = 12;         % Number of Taxels *** make sure to change for different sensors
% ---
% ---

%  ---------------- Initialization Stuff ---------------
colors={[1 0 0],[0 1 0],[0 0 1],[1 0 1],[1 0.5 0],[0 1 1],[0.5 0 0],[0 0.5 0],[0 0 0.5],[0.5 0.5 0],[0.5 0 0.5],[0 0.5 0.5],[0.5 0.5 0.5],[0.25 0 0],[0 0.25 0],[0 0 0.25]};
Data = '1\posY_1p6N_0p2to0p8mm.txt';
DataATI = '1\posY_1p6N_0p2to0p8mm.csv';
area=196; %square indenter in mm2
leg = ["P1", "P2", "P3", "P4"];

name = split(Data, '\');
name = name(end);         %this is so we can use the name to find out how to plot the x or y forces over the cap data

% ------------  ATI Load Cell Set Up -------------------
FreqATI = 62.5;  %This is the sampling rate of the ATI software
ForceATI=readmatrix(DataATI,'NumHeaderLines',1);

xF = -ForceATI(:,2);
yF =  ForceATI(:,1);
zF =  ForceATI(:,3);
xT = -ForceATI(:,5);
yT =  ForceATI(:,4);
zT =  ForceATI(:,6);

listA1=[];
listC1=[];
listA2=[];
listC2=[];
listA3=[];
listC3=[];
listA4=[];
listC4=[];

nLinesATI = numel(ForceATI)/6;


% ----------------- Sensor Data Set up ------------------

Freq = 50/N;  %This is the sample rate of our CDC (we have chosen this on the code, this can be changed - refer to datasheet)
              %Note, this frequency will have to be changed if the code
              %used to read the CDC is not interrupt driven. (sometimes the
              %arduino code will miss cycles) --- The interrupt driven STM
              %code does not miss RDY pin cycles so the exact freq will be
              %50Hz
M = readmatrix(Data);
nLines = numel(M)/N;

% Note, a,b,c,d are the shear addresses
a = M(:,1);
b = M(:,2);
c = M(:,3);
d = M(:,4);


% ----------------- Data Processing ---------------
xcap=linspace(0,nLines/Freq,nLines); %new cap time axis in seconds
xf=linspace(0,nLinesATI/FreqATI, nLinesATI); %new force time axis in second

interpxF =interp1(xf, xF, xcap); % interpolating Force values to match xcap
interpyF =interp1(xf, yF, xcap); % interpolating Force values to match xcap
interpzF =interp1(xf, zF, xcap); % interpolating Force values to match xcap
interpxT =interp1(xf, xT, xcap);
interpyT =interp1(xf, yT, xcap);
interpzT =interp1(xf, zT, xcap);

%Find (C - C0)/C0
deltaA=((a-mean(a(1:10)))/(mean(a(1:10))).*100);
deltaB=((b-mean(b(1:10)))/(mean(b(1:10))).*100);
deltaC=((c-mean(c(1:10)))/(mean(c(1:10))).*100);
deltaD=((d-mean(d(1:10)))/(mean(d(1:10))).*100);

% --------- Find Max and Min Values to Align Axis better ---------------
minValue = min([deltaA(:);deltaB(:);deltaC(:);deltaD(:)]);
maxValue = max([deltaA(:);deltaB(:);deltaC(:);deltaD(:)]);
maxValue = 1.1*maxValue; %make it 10% bigger
maxT = max([interpxT(:);interpyT(:);interpzT(:)])+1;
maxF = max([interpxF(:);interpyF(:);interpzF(:)])+0.2;
minT = min([interpxT(:);interpyT(:);interpzT(:)])-1;
minF = min([interpxF(:);interpyF(:);interpzF(:)])-0.2;


% ------------------------ Now we Plot ----------------------------------

%Plot Raw Cap
figure(1);
ylabel('Capacitance','fontweight','bold');
xlabel('Time(s)','fontweight','bold')
plot(xcap,M,'LineWidth',3)
xlim([0 xcap(end)]);
set(gca,'FontSize',20);
title('Capacitance of Taxels');


%Plot Delta C / C0 using xcap
figure(2);
ylabel('\DeltaC / C_0','fontweight','bold');
xlabel('Time (s)','fontweight','bold')
plot(xcap,deltaA,'g','LineWidth',3)
hold on
plot(xcap,deltaB,'b','LineWidth',3)
hold on
plot(xcap,deltaC,'r','LineWidth',3)
hold on
plot(xcap,deltaD,'m','LineWidth',3)
xlim([0 xcap(end)]);
set(gca,'FontSize',20);
title('Change in Capacitance');
legend(leg);

% SHEAR PLOTS
fig=figure(3)
left_color = [0 0 1];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

if name{1}(4) == 'X'
    subplot (2,1,1)
    yyaxis right
    plot(xcap+delay,interpxF,'k','LineWidth',3);
    ylim([minF maxF]);
    xlim([0 xcap(end)])
    ylabel('Force (N) - X dir','fontweight','bold');
    xlabel('Time (s)')
    hold on;
    yyaxis left
    plot(xcap,deltaB,'g','LineWidth',3)
    ylabel('\DeltaC/C0','fontweight','bold')
    ylim([minValue maxValue])
    hold on
    plot(xcap,deltaD,'-b','LineWidth',3)
    title('Taxels 1 & 3, Along Shear')
    set(gca,'FontSize',20)
    
    subplot(2,1,2)
    yyaxis right
    plot(xcap,interpyF,'k','LineWidth',3);
    xlim([0 xcap(end)])
    ylim([minF maxF]);
    ylabel('Force (N) - Y dir','fontweight','bold');
    xlabel('Time (s)')
    hold on;
    yyaxis left
    plot(xcap,deltaA,'r','LineWidth',3)
    ylim([-5 maxValue]);
    ylabel('\DeltaC/C0','fontweight','bold')
    hold on
    plot(xcap,deltaC,'-m','LineWidth',3)
    title('Taxels 2 & 4, Perpendicular Shear')
    set(gca,'FontSize',20)
    set(gcf, 'Position',  [11, -180, 1200, 850])
else
    subplot (2,1,1)
    yyaxis right
    plot(xcap+delay,interpyF,'k','LineWidth',3);
    ylim([minF maxF]);
    xlim([0 xcap(end)])
    ylabel('Force (N) - Y dir','fontweight','bold');
    xlabel('Time (s)')
    hold on;
    yyaxis left
    plot(xcap,deltaA,'g','LineWidth',3)
    ylabel('\DeltaC/C0','fontweight','bold')
    ylim([minValue maxValue])
    hold on
    plot(xcap,deltaC,'-b','LineWidth',3)
    title('Taxels 2 & 4, Along Shear')
    set(gca,'FontSize',20)
    
    subplot(2,1,2)
    yyaxis right
    plot(xcap,interpxF,'k','LineWidth',3);
    xlim([0 xcap(end)])
    ylim([minF maxF]);
    ylabel('Force (N) - X dir','fontweight','bold');
    xlabel('Time (s)')
    hold on;
    yyaxis left
    plot(xcap,deltaB,'r','LineWidth',3)
    ylim([-5 maxValue]);
    ylabel('\DeltaC/C0','fontweight','bold')
    hold on
    plot(xcap,deltaD,'-m','LineWidth',3)
    title('Taxels 1 & 3, Perpendicular Shear')
    set(gca,'FontSize',20)
    set(gcf, 'Position',  [11, -180, 1200, 850])
end



% ----------- MAIN SHEAR PLOT ---------------------- 
figure(4);
subplot (4,1,1) %plot forces

plot(xcap+delay,interpxF,'LineWidth',3)
hold on 
plot(xcap+delay,interpyF,'LineWidth',3)
hold on 
plot(xcap+delay,interpzF,'LineWidth',3)
ylim([minF maxF])
xlim([0 xcap(end)])
ylabel('Force (N)','fontweight','bold')
title('Forces Values')
legend('Fx', 'Fy', 'Fz','fontweight','bold')
set(gca,'FontSize',18)

subplot (4,1,2) %plot torques
plot(xcap+delay,interpxT,'LineWidth',3)
hold on 
plot(xcap+delay,interpyT,'LineWidth',3)
hold on 
plot(xcap+delay,interpzT,'LineWidth',3)
ylim([minT maxT])
xlim([0 xcap(end)])
ylabel('Torque (N.mm)','fontweight','bold')
title('Torques Values')
legend('Tx', 'Ty', 'Tz','fontweight','bold')
set(gca,'FontSize',18)

subplot (4,1,3) %plot cap data
plot(xcap,deltaA, 'g','LineWidth',3)
hold on
plot(xcap,deltaB,'b','LineWidth',3)
hold on
plot(xcap,deltaC,'r','LineWidth',3)
hold on
plot(xcap,deltaD,'m','LineWidth',3)
ylabel('\DeltaC/C0','fontweight','bold')
xlim([0 xcap(end)])
ylim([minValue maxValue])
title('Change in Capacitance Values')
set(gca,'FontSize',18)

ax=subplot (4,1,4) %plot delta C and corresponding force data
plot(xcap,deltaA, 'g','LineWidth',3)
hold on
plot(xcap,deltaB,'b','LineWidth',3)
hold on
plot(xcap,deltaC,'r','LineWidth',3)
hold on
plot(xcap,deltaD,'m','LineWidth',3)
hold on
ylabel('\DeltaC/C0','fontweight','bold')
ylim([minValue maxValue])
xlim([0 xcap(end)])
xlabel('Time (s)')

yyaxis right
set(gca,'FontSize',18)
set(ax,{'YColor'}, {'k'});

if name{1}(4) == 'X'
    plot(xcap+delay,interpxF,'k','LineWidth',3);
    title('Change in Capacitance Values & Fx');
    ylabel('Fx (N)','fontweight','bold');
else
    plot(xcap+delay,interpyF,'k','LineWidth',3);
    title('Change in Capacitance Values & Fy');
    ylabel('Fy (N)','fontweight','bold');
end
ylim([minF maxF])
xlim([0 xcap(end)])
set(gcf, 'Position',  [11, -180, 1555, 1000]) %[1, 41, 1500, 963]



%ML
k=0;
l=0;
m=0;
n=0;

for i = 1:length(xcap)
    if ((xcap(i)>16.21) && (xcap(i)<18.06))
    k=k+1;
 
    listA1(k)=(deltaA(i));
    listC1(k)=(deltaC(i));
    end
   
    if ((xcap(i)>26.49) && (xcap(i)<29.0134))
    l=l+1;
    listA2(l)=(deltaA(i));
    listC2(l)=(deltaC(i));
    end
    
    if ((xcap(i)>37.2) && (xcap(i)<40.65))
    m=m+1;
    listA3(m)=(deltaA(i));
    listC3(m)=(deltaC(i));
    end
    
    if ((xcap(i)>49.5) && (xcap(i)<52.258))
    n=n+1;
    listA4(n)=(deltaA(i));
    listC4(n)=(deltaC(i)); 
    end
end

listA1=sum(listA1)/length(listA1);
listC1=sum(listC1)/length(listC1);
listA2=sum(listA2)/length(listA2);
listC2=sum(listC2)/length(listC2);
listA3=sum(listA3)/length(listA3);
listC3=sum(listC3)/length(listC3);
listA4=sum(listA4)/length(listA4);
listC4=sum(listC4)/length(listC4);

Abase=[];
k=0;
Cbase=[];
for i = 1:length(xcap)
    if ((xcap(i)>10.1145) && (xcap(i)<13.9677))
        k=k+1;
        Abase(k)=deltaA(i);
        Cbase(k)=deltaC(i);
        
    end
end
Abase=sum(Abase)/length(Abase);
Cbase=sum(Cbase)/length(Cbase);

for i= 1:length(listA1)
    listA1(i)=Abase-listA1(i);
    listA2(i)=Abase-listA2(i);
    listA3(i)=Abase-listA3(i);
    listA4(i)=Abase-listA4(i);
    listC1(i)=abs(Cbase-listC1(i));
    listC2(i)=abs(Cbase-listC2(i));
    listC3(i)=abs(Cbase-listC3(i));
    listC4(i)=abs(Cbase-listC4(i));
end
%W

AC1=0.5*listA1+0.5*listC1
AC2=0.5*listA2+0.5*listC2
AC3=0.5*listA3+0.5*listC3
AC4=0.5*listA4+0.5*listC4

ACP1=4.0334*exp(0.0808*AC1);
ACP2=4.0334*exp(0.0808*AC2);
ACP3=4.0334*exp(0.0808*AC3);
ACP4=4.0334*exp(0.0808*AC4);

Afixed1=[];
Cfixed1=[];
k=0;
check=0
for i = 1:length(xcap)
    if ((xcap(i)>14.93) && (xcap(i)<19.50))
    k=k+1;
    Cfixed1(k)=Cbase+ACP1;
    check=1;
    end
   
    if ((xcap(i)>25.5) && (xcap(i)<30.34))
    k=k+1;
    Cfixed1(k)=Cbase+ACP2;
    check=1;
    end
    
    if ((xcap(i)>36.36) && (xcap(i)<41.903))
    k=k+1;
    Cfixed1(k)=Cbase+ACP3;
    check=1;
    end
    
    if ((xcap(i)>48.16) && (xcap(i)<53.45))
    k=k+1;
    Cfixed1(k)=Cbase+ACP4;
    check=1;
    end
    if(xcap(i)<9.63 ||xcap(i)> 64.54)
        k=k+1;
        Cfixed1(k)=0;
        check=1;
    end
    
    if(check==0)
        k=k+1;
        Cfixed1(k)=Cbase;
    end 
    check=0;
    
end



k=0;
check=0
for i = 1:length(xcap)
    if ((xcap(i)>14.93) && (xcap(i)<19.50))
    k=k+1;
    Afixed1(k)=Abase-ACP1;
    check=1;
    end
   
    if ((xcap(i)>25.5) && (xcap(i)<30.34))
    k=k+1;
    Afixed1(k)=Abase-ACP2;
    check=1;
    end
    
    if ((xcap(i)>36.36) && (xcap(i)<41.903))
    k=k+1;
    Afixed1(k)=Abase-ACP3;
    check=1;
    end
    
    if ((xcap(i)>48.16) && (xcap(i)<53.45))
    k=k+1;
    Afixed1(k)=Abase-ACP4;
    check=1;
    end
    if(xcap(i)<9.63 ||xcap(i)> 64.54)
        k=k+1;
        Afixed1(k)=0;
        check=1;
    end
    
    if(check==0)
        k=k+1;
        Afixed1(k)=Abase;
    end 
    check=0;
    
end
figure(5);
ylabel('\DeltaC / C_0 (original data from Characterization Setup)','fontweight','bold');
xlabel('Time (s)','fontweight','bold')

hold on
plot(xcap,deltaC,'r','LineWidth',3)
hold on
plot(xcap,deltaA,'m','LineWidth',3)
xlim([0 xcap(end)]);
set(legend("deltaC","deltaA"));

title('Change in Capacitance (Original)');
ylabel('\DeltaC / C_0 ','fontweight','bold');
xlabel('Time (s)','fontweight','bold')
legend("deltaC","deltaA");


figure(6);

plot(xcap,Cfixed1,'g','LineWidth',3)
hold on
plot(xcap,Afixed1,'b','LineWidth',3)
hold on
xlim([0 xcap(end)]);
set(gca,'FontSize',7);
title('Change in Capacitance ( fixed using ML (assuming 2 individual original baselines))');
ylabel('\DeltaC / C_0 ','fontweight','bold');
xlabel('Time (s)','fontweight','bold')
legend("deltaC","deltaA");


Afixed2=[];
Cfixed2=[];
ACBase=0.5*Abase+0.5*Cbase;

k=0;
check=0
for i = 1:length(xcap)
    if ((xcap(i)>14.93) && (xcap(i)<19.50))
    k=k+1;
    Cfixed2(k)=ACBase+ACP1;
    Afixed2(k)=ACBase-ACP1
    check=1;
    end
   
    if ((xcap(i)>25.5) && (xcap(i)<30.34))
    k=k+1;
    Cfixed2(k)=ACBase+ACP2;
    Afixed2(k)=ACBase-ACP2;
    check=1;
    end
    
    if ((xcap(i)>36.36) && (xcap(i)<41.903))
    k=k+1;
    Cfixed2(k)=ACBase+ACP3;
    Afixed2(k)=ACBase-ACP3;
    check=1;
    end
    
    if ((xcap(i)>48.16) && (xcap(i)<53.45))
    k=k+1;
    Cfixed2(k)=ACBase+ACP4;
    Afixed2(k)=ACBase-ACP4;
    check=1;
    end
    if(xcap(i)<9.63 ||xcap(i)> 64.54)
        k=k+1;
        Cfixed2(k)=0;
        Afixed2(k)=0;
        check=1;
    end
    
    if(check==0)
        k=k+1;
        Cfixed2(k)=ACBase;
        Afixed2(k)=ACBase;
    end 
    check=0;
    
end

figure(7);
plot(xcap,Cfixed2,'g','LineWidth',3)
hold on
plot(xcap,Afixed2,'b','LineWidth',3)
hold on
xlim([0 xcap(end)]);
set(gca,'FontSize',7);
title('Change in Capacitance ( fixed using ML (assuming equal Baseline (Average of Individual Baselines))');
ylabel('\DeltaC / C_0 ','fontweight','bold');
xlabel('Time (s)','fontweight','bold')
legend("deltaC","deltaA");   
    
    




































%%%%% BELOW HERE IS PREVIOUS CODE FROM PREVIOUS PEOPLE - austin 

% % Plot Capacitance vs Force or Displacement 
% % Here we want to find one cap value at each peak corresponding to a certain
% % displacement or force
% 
% displacement= [0 0.2 0.4 0.6 0.8]; %Input displacement matrix from on experiment
% startforce=0.023889 %where is Fx, Fy, or Fz at the start(be careful with signs), right before shear but after first initial pressure. Be careful if negative. Take the middle of the peak.
% force_shear=[0 0.248576-startforce 0.47336-startforce 0.689295-startforce 0.886005-startforce]; % Look at Fx or Fy peaks corresponding to each displacement. 
% pressure=((force_shear)./area).*1000;
% 
% % Raw capacitance to pick point
% figx = figure (9);
% plot(xcap-delay,interpA,'g')
% hold on
% plot(xcap-delay,interpB,'b')
% hold on
% plot(xcap-delay,interpC,'r')
% hold on
% plot(xcap-delay,interpD,'m')
% xlabel('Time [seconds]','FontSize',20,'fontweight','bold')
% ylabel('Capacitance [pF]','FontSize',30,'fontweight','bold')
% set(gca,'FontSize',18)
% set(gcf, 'Position',  [1, -180, 1100, 963])
% set(gca,'color','none','FontSize',20)
% set(gcf, 'color', 'none'); 
% 
% num_input = input(sprintf('Number of taxels?  ')); %How many capacitance values are we looking at?
% clearvars c_info
% peaks= input(sprintf('Please select in order: green (baseline then peaks), blue (baseline then peaks), red (baseline then peaks), magenta(baseline then peaks). Number of points to detect (include baseline)?  ')); % How many peaks do we want to plot?
% 
% %click on BASELINE (after pressure) FIRST THEN each peak of one taxel, then move on to second taxel, etc.
% for j=1:num_input
%       for i = 1:peaks
%               shg
%               dcm_obj = datacursormode(figx);
%               set(dcm_obj,'DisplayStyle','window',...
%                   'SnapToDataVertex','off','Enable','on')
%               waitforbuttonpress
%               c_info{i,j} = getCursorInfo(dcm_obj);
%              
%       end
% end
% 
% disp('Click line to display a data tip, then press Enter.')
% % Wait while the user does this.
% pause 
% % 
% 
% figure (10) 
% plot(xcap-delay,interpA,'g')
% hold on
% plot(xcap-delay,interpB,'b')
% hold on
% plot(xcap-delay,interpC,'r')
% hold on
% plot(xcap-delay,interpD,'m')
% hold on
% 
% for j=1:num_input
%       for i = 1:size(c_info,1)
%             position{i,j} = c_info{i,j}.Position; %Put the clicked values into arrays: positions{point#,taxel#)
%             plot(position{i,j}(1),position{i,j}(2),'r*');%plot these values and make sure they look correct graphically
%             select_cap(i,j)=(position{i,j}(2)-position{1,j}(2))/(position{1,j}(2)).*100; % create an array of points selected for each taxels in the order of green, blue, red, magenta
%             % First row is baseline, then points selected peak values DC over C0. 
%             % 1st to 4th columns are taxel numbers.
% %             select_cap(i,j)=(position{i,j}(2)-position{1,j}(2)); % create an array of points selected for each taxels in the order of green, blue, red, magenta
% %             % First row is baseline, then points selected peak values DC. 
% %             % 1st to 4th columns are taxel numbers.
%       end                   
% end
% 
% %Plot each taxel DC over C0 vs Displacement
% 
% figure(11)
% subplot(1,2,1)
% plot(displacement,select_cap(:,1),'--go','MarkerSize',12,'LineWidth',3)
% hold on 
% plot(displacement,select_cap(:,2),'--bo','MarkerSize',12,'LineWidth',3)
% hold on
% plot(displacement,select_cap(:,3),'--ro','MarkerSize',12,'LineWidth',3)
% hold on 
% plot(displacement,select_cap(:,4),'--mo','MarkerSize',12,'LineWidth',3)
% 
% ylim([-25 25])
% xlim([0 1])
% title('%\DeltaC/C0 vs Shear Displacement in x')
% %title('%\DeltaC/C0 vs Shear Force in -y')
% xlabel('Displacement (mm)','fontweight','bold')
% ylabel('% \DeltaC after initial pressure ','fontweight','bold')
% %ylabel('% \DeltaC/C0 ','fontweight','bold')
% set(gca,'color','none','FontSize',20)
% set(gcf, 'color', 'none');
% 
% % 
% % 
% % Plot each taxel DC over C0 vs Shear Force
% %Figure dimension: 1300 width, 963 height in pixels (in figure properties)
% % 
% subplot(1,2,2)
% plot(force_shear,select_cap(:,1),'--go','MarkerSize',12,'LineWidth',3)
% hold on 
% plot(force_shear,select_cap(:,2),'--bo','MarkerSize',12,'LineWidth',3)
% hold on
% plot(force_shear,select_cap(:,3),'--ro','MarkerSize',12,'LineWidth',3)
% hold on 
% plot(force_shear,select_cap(:,4),'--mo','MarkerSize',12,'LineWidth',3)
% ylim([-25 25])
% xlim([0 1])
% %title('% \DeltaC vs Pressure')
% %title('% \DeltaC vs Shear Force in x')
% title('% Capacitance Change vs Shear Force in x')
% 
% %xlabel('Pressure (kPa)','fontweight','bold')
% xlabel('Force (N)','fontweight','bold')
% ylabel('% \DeltaC after initial pressure','fontweight','bold')
% set(gca,'color','none','FontSize',20)
% set(gcf, 'color', 'none');
% 
% set(gcf, 'Position', [0,0, 400,1300])
% green=[position{1,1}(2);position{2,1}(2); position{3,1}(2) ;position{4,1}(2);position{5,1}(2)]
% blue=[position{1,2}(2);position{2,2}(2); position{3,2}(2) ;position{4,2}(2);position{5,2}(2)]
% red=[position{1,3}(2);position{2,3}(2); position{3,3}(2) ;position{4,3}(2);position{5,3}(2)]
% magenta=[position{1,4}(2);position{2,4}(2); position{3,4}(2) ;position{4,4}(2);position{5,4}(2)]
% fp=[0 force_shear(2)+startforce force_shear(3)+startforce force_shear(4)+startforce force_shear(5)+startforce];
% fp'


