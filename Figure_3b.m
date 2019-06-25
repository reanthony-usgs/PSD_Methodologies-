% Quick histogram code to show slices through time

% %%%%%%%%%%%% User Defined parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Frequency band of interest 

clear all

fmin = 0.03;
fmax = 0.04;

% Define bin spacing (in dB)

BS = 0.1;


% if necessary load in the data

load ANMO_STS6_May2019_3Hr_Bootstrap.mat
load LHZ_freqs.txt

load ANMO_3Hr_STS6_May2019.mat

%load Hours_ANMO_2018_279_2_9.txt



%%
PSD_M = ANMO_STS6_3Hr;
Hour = ANMO_Hours;

% Cut down to the day of interest 

TI = find(Hour > 173 & Hour < 174);

Hour = Hour(TI);
PSD = PSD_M(TI,:);

% Convert to Hours

Hour = (Hour-173)*24;


SI = find(LHZ_freqs > fmin,1, 'first');
SE = find(LHZ_freqs > fmax,1, 'first');


%% Take the average across the frequency band

Fband_Means = mean(PSD(:,SI:SE),2);

EQ_I = find(Fband_Means > -173);


%% Make a figure showing temporal evolution of power



figure(8)
clf
plot(Hour,Fband_Means,'kx')
hold on
plot([0,24],[-173,-173],'r-')
xlabel('Hour of Day')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')
xlim([0,24])
ylim([-190, -120])
set(gca,'FontSize',20)


