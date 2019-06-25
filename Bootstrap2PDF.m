% This code takes the bootstrapped data matrix output by Python codes and
% rearranges it into the 2-column period-power format of NTK. This can then
% be loaded into codes to generate PSD PDFs

clear all

%Pre-load the Bootstrapped Matrix
load ANMO_STS6_May2019_Bootstrap.mat
load LHZ_freqs.txt



PSD_M = ANMO_STS6;
PSD_M = fliplr(PSD_M);
Periods = 1./LHZ_freqs;

b = length(Periods); 
Periods = flipud(Periods);

% Get the percentile statistics 

MedianPSD = prctile(PSD_M,50);
MeanPSD = mean(PSD_M);

p2 = prctile(PSD_M,2.5);
p97 = prctile(PSD_M,97.5);

% Separately Define Periods with Earthquakes as found from the 25 to 33s
% period band

Pmin = 25;
Pmax = 33.333;

SI = find(Periods > Pmin,1, 'first');
SE = find(Periods > Pmax,1, 'first');


% Take the average across the frequency band

Fband_Means = mean(PSD_M(:,SI:SE),2);

EQI = find(Fband_Means > -173);

PSD_EQF = PSD_M;

PSD_EQF(EQI,:) = [];

Median_EQF = prctile(PSD_EQF,50);
Mean_EQF = mean(PSD_EQF);

% Repeat for Only Earthquake PSDs

Median_EQ = prctile(PSD_M(EQI,:),50);






runs = size(PSD_M,1);

%%
%Powers = reshape(PSD_M,runs*b,1);
%Period_V = repelem(Periods,runs);







%%

%Chart = [Period_V, Powers];

histcent = [-200:.1:-80];
[counts] = hist(PSD_M(:,:), histcent);

%Periods = Periods(1:10:end);


% Make the Peterson curves
fs=250;
dlP=.05;
PSDTOL=15;
[LNMA,HNMA,lpd1,lpd2]=peterson_acc(dlP,fs);

%Smoothed Peterson curves for plotting
NMplotind=(0.001:dlP:10);
LNMAp=spline(lpd1,LNMA,NMplotind);
HNMAp=spline(lpd2,HNMA,NMplotind);

pd1 = 10.^(lpd1);
pd2 = 10.^(lpd2);


figure(20); clf
%bookfonts

Period_lines =  [5,10,50,100,500];
ticks = log10(Period_lines);

h = pcolor(Periods,histcent,log10(counts))

cmap = viridis;
% Make values 0-5 black:
cmap(1,:) = 0.3*ones(1,3);
colormap(cmap);
c=colorbar
xlim([2.5 500])
ylim([-200 -80])
caxis([0 4])
set(gca,'FontSize',20)


hold on



H5 = plot(Periods,MedianPSD,'k');
H8 = plot(Periods,Median_EQF,'w');
H9 = plot(Periods,Median_EQ,'w');

H2 = plot(10.^(lpd1),LNMA,'w:');
H3 = plot(10.^(lpd2),HNMA,'w:');
H6 = plot(Periods,p2,'k:');
H7 = plot(Periods,p97,'k:');


set(H2,'LineWidth',5.0);
set(H3,'LineWidth',5.0);
set(H5,'LineWidth',3.0);
set(H6,'LineWidth',3.0);
set(H7,'LineWidth',3.0);
set(H8,'LineWidth',1.0);
set(H9,'LineWidth',2.0);




lgd = legend([H5 H6 H8 H2], 'Median PSD', '2.5%/97.5% PSD', 'EQ/EQ removed Median PSDs', 'NHNM/NLNM');
lgd.Color = [0.7 0.7 0.7];



set(h, 'EdgeColor', 'none');


grid off



axis off
axbot = gca;
set(axbot, 'XScale', 'log', 'YScale', 'linear');
axtop = axes('Position',get(axbot,'Position'),'Color','none',...
            'Xlim',get(axbot,'XLim'), 'Ylim',get(axbot,'YLim'),...
            'XScale', 'log', 'YScale', 'linear' , ...
            'YMinorTick','off' , 'YMinorGrid','off'....
            ) ;
set(gca,'FontSize',20) 

%ticks = get(axtop,'XTickLabel')
%Xticklabels = cellstr(ticks, '10^%d');

%set(axtop,'Xticklabel',Xticklabels)
%XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));

axtop.LineWidth = 3;

        

xlabel('Period (s)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')
ylabel(c,'Log_{10}(Counts)') 




%% Repeat Using Only the EQ Free Data



p2F = prctile(PSD_EQF,2.5);
p97F = prctile(PSD_EQF,97.5);

figure(14); clf

histcent = [-200:.1:-80];
[counts2] = hist(PSD_EQF(:,:), histcent);

h = pcolor(Periods,histcent,log10(counts2))

cmap = jet;
% Make values 0-5 black:
cmap(1,:) = 0.3*ones(1,3);
colormap(cmap);
c=colorbar
xlim([2.5 300])
ylim([-200 -80])
caxis([0 3.7])
set(gca,'FontSize',20)


hold on




H8 = plot(Periods,Median_EQF,'k');

H2 = plot(10.^(lpd1),LNMA,'w:');
H3 = plot(10.^(lpd2),HNMA,'w:');
H6 = plot(Periods,p2F,'k:');
H7 = plot(Periods,p97F,'k:');


set(H2,'LineWidth',5.0);
set(H3,'LineWidth',5.0);
set(H5,'LineWidth',3.0);
set(H6,'LineWidth',2.0);
set(H7,'LineWidth',2.0);
set(H8,'LineWidth',3.0);





lgd = legend([H5 H6 H2], 'Median PSD', '2.5%/97.5% PSD', 'NHNM/NLNM');
lgd.Color = [0.7 0.7 0.7];



set(h, 'EdgeColor', 'none');


grid off



axis off
axbot = gca;
set(axbot, 'XScale', 'log', 'YScale', 'linear');
axtop = axes('Position',get(axbot,'Position'),'Color','none',...
            'Xlim',get(axbot,'XLim'), 'Ylim',get(axbot,'YLim'),...
            'XScale', 'log', 'YScale', 'linear' , ...
            'YMinorTick','off' , 'YMinorGrid','off'....
            ) ;
set(gca,'FontSize',20) 

%ticks = get(axtop,'XTickLabel')
%Xticklabels = cellstr(ticks, '10^%d');

%set(axtop,'Xticklabel',Xticklabels)
%XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));

axtop.LineWidth = 3;

        

xlabel('Period (s)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')
ylabel(c,'Log_{10}(Counts)') 





%% Differences

Median_Diff = MedianPSD - Median_EQF;
Mean_Diff = MeanPSD - Mean_EQF;

figure(15); clf
semilogx(Periods,Median_Diff,'b','LineWidth',2)
hold on
semilogx(Periods,Mean_Diff,'b--','LineWidth',2)
xlabel('Period (s)')
ylabel('dB Difference')


set(gca,'FontSize',20) 
xlim([2.5 500])






