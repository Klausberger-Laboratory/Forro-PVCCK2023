function VRShuffling

% script to check shuffling statististics, check significant differences 
% between no-odour and odour trial activities or between cell types or 
% any type of two matrices

% copy in manually the matrices to be compared

%script shuffles the identities of the trials

% p-value of 0.01 or 0.05 can be set at the end

% values above the threshold line or below are considered significant

close all
clear all
HighGroup = 0;%% add manually the matrix of the trial type with more trials, meaning no-odour trials, each row is one trial, columns are spatial bins
LowGroup = 0;

clearvars -except HighGroup LowGroup

firingrates = cat(1,HighGroup,LowGroup);

%load firingrates.txt
%%%%  number of non odour and odour trials
cellsHigh = size(HighGroup,1);

cellsLow = size(LowGroup,1);

alln = cellsHigh + cellsLow;

OrigFiringHigh = nanmean(HighGroup,1);
OrigFiringLow = nanmean(LowGroup,1);
OrigDiff = OrigFiringHigh - OrigFiringLow;

figure
hold on
%plot(OrigFiringHigh,'b')
%plot(OrigFiringLow,'r')
plot(OrigDiff,'g')

for shff = 1:10000

%%% new version

randomizedIND = randperm(alln);
    
randomizedfiring = firingrates(randomizedIND,:);    

randomizedHighfiring = randomizedfiring(1:cellsHigh,:);
randomizedLowfiring = randomizedfiring((cellsHigh+1):end,:);

randmeanFiringHigh = nanmean(randomizedHighfiring,1);
randmeanFiringLow = nanmean(randomizedLowfiring,1);
randDiff(shff,:) = randmeanFiringHigh - randmeanFiringLow;

%%%% make a select version where same amount of cells are selected

randHighcells = randperm(cellsHigh);
%select first same as number of low group cells
indselHighGroup = randHighcells(1:cellsLow);

selHighGroup = HighGroup(indselHighGroup,:);
selcellsHigh = cellsLow;
selfiringrates = cat(1,selHighGroup,LowGroup);
selalln = 2*cellsLow; 

selrandomizedIND = randperm(selalln);
    
selrandomizedfiring = selfiringrates(selrandomizedIND,:);    

selrandomizedHighfiring = selrandomizedfiring(1:selcellsHigh,:);
selrandomizedLowfiring = selrandomizedfiring((selcellsHigh+1):end,:);

selrandmeanFiringHigh = nanmean(selrandomizedHighfiring,1);
selrandmeanFiringLow = nanmean(selrandomizedLowfiring,1);
selrandDiff(shff,:) = selrandmeanFiringHigh - selrandmeanFiringLow;

end

%%%% correct for 4 sectors equals to 0.6333%
% sortedDiff = sort(randDiff);
% correc = 10000 * 0.005
% 10000 - correc;
% 
% upper = sortedDiff(10000-correc,:)
% lower = sortedDiff(correc,:)
% 
% figure
% hold on
% plot(upper)
% plot(lower)
% plot(OrigDiff,'.g')


%%%% plot selected same group size
selsortedDiff = sort(selrandDiff); %each column sorted
%selcorrec = 10000 * 0.025
selcorrec = 10000 * 0.005 % this is for p-value for 0.01
10000 - selcorrec;

selupper = selsortedDiff(10000-selcorrec,:)
sellower = selsortedDiff(selcorrec,:)

figure
hold on
plot(selupper,'g')
plot(sellower,'b')
plot(OrigDiff,'.r','MarkerSize',10)


end