% This code calculates the difference between 1 and 3-hour bootstrapped PSD
% estimates. Periods and powers where the 3 hour estimate is more likely
% (less)are coded in red (blue). 


clear all

%Load in all the data

load Difference_1hr_3Hr/histcent.mat 
load Difference_1hr_3Hr/Periods.mat 

load Difference_1hr_3Hr/Counts_1Hr.mat 
counts1 = counts;


load Difference_1hr_3Hr/Counts_3Hr.mat 
counts3 = counts;

counts_diff = counts3-counts1;

% Make a matrix of Polarities 

Pol = ones(size(counts_diff));

for kk = 1:size(Pol,2)
    NI = find(counts_diff(:,kk) < 0);
    
    for jj = 1:length(NI)
        Pol(NI(jj),kk) = -1;
    end
    
end


%% Make the Figure 


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


figure(19); clf

% Find negative counts

ncounts = Pol.*log10(abs(counts_diff));

ncounts(isinf(ncounts)) = 0;


h = pcolor(Periods,histcent,ncounts);

cmap = redblue;
%Make values 0-5 black:
%cmap(1,:) = 0.3*ones(1,3);
colormap(cmap);
c=colorbar
xlim([2.5 500])
ylim([-200 -80])
caxis([-4 4])
set(gca,'FontSize',20)


hold on


H2 = plot(10.^(lpd1),LNMA,'k:');
H3 = plot(10.^(lpd2),HNMA,'k:');

set(H2,'LineWidth',5.0);
set(H3,'LineWidth',5.0);


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

axtop.LineWidth = 3;

        

xlabel('Period (s)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')
ylabel(c,'Log_{10}(Counts)') 







