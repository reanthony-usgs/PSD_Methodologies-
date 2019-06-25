% This code uses the earthquake culling methodology of Aster et al., 2010
% to count the number of bootstrapped ANMO PSDs that contain earthquakes.
% This number is then compared against the number determined to have
% earthquakes based on 25-33 s power. 

clear all


% Set defaut parameters 
EQ_Thres = 0.95;

% load in the data (3 Hour) data 


load ANMO_STS6_May2019_3Hr_Bootstrap.mat
load LHZ_freqs.txt



PSD_M = ANMO_STS6_3Hr;
PSD_M = fliplr(PSD_M);
Periods = 1./LHZ_freqs;

b = length(Periods); 
Periods = flipud(Periods);


% Get 80th percentile Statistics 
p80 = prctile(PSD_M,80);


%% Find Periods greater than 30s and cut down the PSD and 80th percentile

SI = find(Periods > 30, 1, 'first');

PSD_M_LP = PSD_M(:,SI:end);
p80 = p80(SI:end);


%% Now we find the number of PSD bins that are 95% above the 80th percentile statistic for periods > 30 s

EQbinNum = length(p80);
EQ_Index = [];

% Now actually pick out PSDs containing earthquakes 
EQ_Index_Counter = 0;
for kk = 1:size(PSD_M_LP,1)
    kk
    Bin_Counter = 0;
    for mm = 1:EQbinNum
        if p80(mm) < PSD_M_LP(kk,mm)
            Bin_Counter = Bin_Counter+1;
        end
    end
       
        if Bin_Counter/EQbinNum > EQ_Thres
            %display('Earthquake detected')
            EQ_Index(EQ_Index_Counter+1) = kk;
            EQ_Index_Counter = EQ_Index_Counter+1;
            
        end
        
end


