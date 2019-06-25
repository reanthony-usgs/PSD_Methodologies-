clear all 

% load in the data

% load in un-smoothed NTK without deconvolution filter
load TA.O18A..BHZ.2008-08-14T12:00:00.000000.3600.period_NoDeconCorr.txt
NTK_Raw = TA_O18A__BHZ_2008_08_14T12_00_00_000000_3600_period_NoDeconCorr;

% load in un-smoothed NTK
load TA.O18A..BHZ.2008-08-14T12:00:00.000000.3600.period_CSD_Seg5IR.txt
NTK_Decon=TA_O18A__BHZ_2008_08_14T12_00_00_000000_3600_period_CSD_Seg5IR; 



% load in Adam's PSD calculation 
load Adam_O18A_PSD_Seg5.txt
Adam = Adam_O18A_PSD_Seg5;

% Mustang Metrics - in frequency (Hz)
load Mustang_O18A_BHE_2008_08_14_12.dat
Must_E = Mustang_O18A_BHE_2008_08_14_12;
load Mustang_O18A_BHN_2008_08_14_12.dat
Must_N = Mustang_O18A_BHN_2008_08_14_12;
load Mustang_O18A_BHZ_2008_08_14_12.dat
Must_Z = Mustang_O18A_BHZ_2008_08_14_12;

% Keiths 
load 2008.08.14.12.eig.out
Keith = X2008_08_14_12_eig;

% IRIS NTK (1/16th smooth) - returns period (s)
load TA.O18A.--.BHEpsd.dat
NTK_E = TA_O18A____BHEpsd(1:96,7:8);
load TA.O18A.--.BHNpsd.dat
NTK_N = TA_O18A____BHNpsd(1:96,7:8);
load TA.O18A.--.BHZpsd.dat
NTK_Z = TA_O18A____BHZpsd(1:96,7:8);

% fix bad values

DI = find(NTK_Z(:,2) <= -190);
NTK_Z(DI,:) = [];


% IRIS NTK (1/16th smooth) - returns period (s)
load TA.O18A.--.BHEpsd_1.dat
NTK_E1 = TA_O18A____BHEpsd_1(1:96,7:8);
load TA.O18A.--.BHNpsd_1.dat
NTK_N1 = TA_O18A____BHNpsd_1(1:96,7:8);
load TA.O18A.--.BHZpsd_1.dat
NTK_Z1 = TA_O18A____BHZpsd_1(1:96,7:8);



%% Make the figure
%Hz_lines = [15,10,5,3,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001];

%ticks = (1./Hz_lines);
%HZ_label = (1./ticks);

%Make the Peterson curves
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

figure(1);clf 
H7 = semilogx(Adam(:,1),Adam(:,2), 'k-');
hold on
%H1 = semilogx(1./Must_Z(:,1),Must_Z(:,2),'b');
%H5 = semilogx(NTK_Z1(:,1),NTK_Z1(:,2),'m');
%H6 = semilogx(NTK_Z(:,1),NTK_Z(:,2),'c');
%H2 = semilogx(1./Keith(:,1),Keith(:,2),'r');
H3 = semilogx(pd1,LNMA,'k:');
H4 = semilogx(pd2,HNMA,'k:');
%H5 = semilogx(log10(Periods),Calm_PSD,'b');
%H6 = semilogx(log10(Periods),Wind_PSD,'r');
%legend('High Discharge Median','Low Discharge Median')

set(gca,'FontSize',30)
xlim([0.1 200])
ylim([-200 -90])
set(gca,'ydir','normal')

%set(H1,'LineWidth',3.0);
%set(H5,'LineWidth',3.0);
set(H4,'LineWidth',3.0);
set(H3,'LineWidth',3.0);
set(H7,'LineWidth',3.0);
%set(H2,'LineWidth',3.0);
%set(H6,'LineWidth',3.0);

%set(gca,'xtick',ticks)
%set(gca,'Xticklabel',HZ_label)
%set(gca,'xdir','reverse')
%legend('Welchs Method', 'IRIS MUSTANG', 'Koper and Burlacu (2015)', 'NHNM/NLNM')

xlabel('Period (s)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')


%% Make a plot of the difference

%NTK_Raw = NTK_Raw(2:end,:);

%NTK_Decon = NTK_Decon(2:end,:);

PSD_Diff = NTK_Decon(:,2) - Adam(:,2);


figure(2);clf 
semilogx(NTK_Decon(:,1),PSD_Diff, 'k.','linewidth',3);

set(gca,'FontSize',30)
xlim([0.05 1000])
ylim([-10 5])
set(gca,'ydir','normal')

set(gca,'xtick',ticks)
set(gca,'Xticklabel',HZ_label)
set(gca,'xdir','reverse')

xlabel('Frequency (Hz)')
ylabel('PSD Difference (from Deconvolution Fiter only) (dB)')



PSD_Diff = NTK_Decon(:,2) - Adam(:,2);


figure(3);clf 
semilogx(NTK_Decon(:,1),PSD_Diff, 'k-','linewidth',3);

set(gca,'FontSize',30)
xlim([0.05 1000])
%ylim([-3 3])
set(gca,'ydir','normal')

set(gca,'xtick',ticks)
set(gca,'Xticklabel',HZ_label)
set(gca,'xdir','reverse')

xlabel('Frequency (Hz)')
ylabel('PSD Difference (dB)')



