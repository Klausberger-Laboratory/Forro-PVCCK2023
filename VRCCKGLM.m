function VRCCKGLM
%%%% This script has been used for the two CCK cells, it check the firing rate and theta LFP in defined episodes (given as txt files)
%%%% has to be run in the respective folder for the two CCK cells as well
%%%% Different GLM models can be fitted and tested in this script

% TO RUN THIS SCRIPT YOU NEED TO GO INTO THE MATLAB GLM Folder for each CCK
% cell

clear all

    %%% ADJUST HERE FOR WHICH CELL TO PROCESS
    load('TF107b1.mat')
    %load('TF131b1.mat')

    load('VR.txt')

    %%% add speed column to VR

    for t1 = 1:(size(VR,1)-1)

        x1 = VR(:,2);
        y1 = VR(:,3);

        xdist = x1(t1+1) - x1(t1);
        ydist = y1(t1+1) - y1(t1);
        totaldist = sqrt(xdist*xdist + ydist*ydist);
        %speed(t1) = totaldist;
        speed(t1) = (totaldist/(1/20));

    end
    
    speed1 = speed(1);
    speed = [speed1;speed'];
    
    VRSpeed = [VR speed];



    spiketimes = TFjx107b1_Ch7.times;
    %spiketimes = TFjx131b1_Ch7.times;

    %define start and end of spikes to use
    startT = 0;
    %endT = 420;
    endT = 390;


    TF = 1;

    binNs = endT/0.4;


    %%% bin into 400 ms bins 
    currT = 0;
    for spkBins = 1:binNs

        clear spkInds

        currT = currT + 0.4;

        spkInds = find(spiketimes >= (currT-0.4) & spiketimes < currT);

        if isempty(spkInds) == 1
            FRbinned(spkBins) = 0;
        else
            FRbinned(spkBins) = length(spkInds)/0.4;
        end
    end


    VRSpeedTimes = VRSpeed(:,1);

    %%% do the same for speed
    currT = 0;
    for speedBins = 1:binNs

        clear spkInds currVRSpeed

        currT = currT + 0.4;

        spkInds = find(VRSpeedTimes >= (currT-0.4) & VRSpeedTimes < currT);

        if isempty(spkInds) == 1
            Speedbinned(speedBins) = 0;
        else
            currVRSpeed = VRSpeed(spkInds,4)
            Speedbinned(speedBins) = nanmean(currVRSpeed,1);
        end
    end



    %%%% next assign to each FRbin a expl. variable
    %1. VR Speed (contin)
    %2. Reward (categorical)
    %3. Odour- no odour (categorical)
    %3. Stuck in VR (not at reward)
    %4. Sitting in VR (not during reward)

    BinTimes = 0.4:0.4:endT;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%  load Rewardtimes
    load('Reward.txt')
    allRewID = [];
    for rewT = 1:size(Reward,1)
        currRew = Reward(rewT,:);
        currRewID = find(BinTimes >= currRew(1) & BinTimes < currRew(2));
        allRewID = [allRewID; currRewID'];
    end

    RewVar(1:binNs) = 0;
    RewVar(allRewID) = 1; RewVar = RewVar';

    %%%%%%%%%%%%%%%%%%%%%%  

    %%%%%% Load odour no odour stimulus times
    load('NoOdourST.txt')
    load('OdourST.txt')

    load('AllST.txt')
    NoOdourST = AllST; %%% IF WE DONT WANT DIFFERENTIATION


    allNoOdID = [];
    for NoT = 1:size(NoOdourST,1)
        currNO = NoOdourST(NoT,:);
        currNOD = find(BinTimes >= currNO(1) & BinTimes < currNO(2));
        allNoOdID = [allNoOdID; currNOD'];
    end

    STVar(1:binNs) = 0;
    STVar(allNoOdID) = 1; STVar = STVar';

    allOdID = [];
    for oT = 1:size(OdourST,1)
        currO = OdourST(oT,:);
        currOD = find(BinTimes >= currO(1) & BinTimes < currO(2));
        allOdID = [allOdID; currOD'];
    end

    %STVar(allOdID) = 2; uncomment if want to look at odour specifically


    %%%%%%%%%%%%%%%% NO VR RUN 
    load('NoVRRun.txt')
    allNoVRID = [];
    for NoVRT = 1:size(NoVRRun,1)
        currNoVR = NoVRRun(NoVRT,:);
        currNoVRID = find(BinTimes >= currNoVR(1) & BinTimes < currNoVR(2));
        allNoVRID = [allNoVRID; currNoVRID'];
    end

    NoVRVar(1:binNs) = 0;
    NoVRVar(allNoVRID) = 1; NoVRVar = NoVRVar';



    %%%%%%%%%%%%%%%%% VR locked run

    load('LockedVRRun.txt')
    allLockedVRID = [];
    for LockedVRT = 1:size(LockedVRRun,1)
        currLockedVR = LockedVRRun(LockedVRT,:);
        currLockedVRID = find(BinTimes >= currLockedVR(1) & BinTimes < currLockedVR(2));
        allLockedVRID = [allLockedVRID; currLockedVRID'];
    end

    LockedVRVar(1:binNs) = 0;
    LockedVRVar(allLockedVRID) = 1; LockedVRVar = LockedVRVar';


    %%%%%% sitting VR not REW

    load('SitVR.txt')
    allSitVRID = [];
    for SitVRT = 1:size(SitVR,1)
        currSitVR = SitVR(SitVRT,:);
        currSitVRID = find(BinTimes >= currSitVR(1) & BinTimes < currSitVR(2));
        allSitVRID = [allSitVRID; currSitVRID'];
    end

    SitVRVar(1:binNs) = 0;
    SitVRVar(allSitVRID) = 1; SitVRVar = SitVRVar';
   
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 






    %%%% extract mean firings of variables

    %mean reward
    allRewFR = FRbinned(allRewID);
    avrgRewFR = nanmean(allRewFR);
    stdallRewFR = nanstd(allRewFR);
    SERewFR = stdallRewFR/sqrt(length(allRewFR));


    %mean locked VR running
    allLockedVRFR = FRbinned(allLockedVRID);
    avrgallLockedVRFR = nanmean(allLockedVRFR);
    stdLockedVRFR = nanstd(allLockedVRFR);
    SELockedVR = stdLockedVRFR/sqrt(length(allLockedVRFR));


    %mean Sitting in VR
    allSitVRFR = FRbinned(allSitVRID);
    avrgallSitVRFR = nanmean(allSitVRFR);
    stdallSitVRFR = nanstd(allSitVRFR);
    SESitVR = stdallSitVRFR/sqrt(length(allSitVRFR));


    % mean no VR running
    allNoVRFR = FRbinned(allNoVRID);
    avrgallNoVRFR = nanmean(allNoVRFR);
    stdallNoVRFRFR = nanstd(allNoVRFR);
    SENoVR = stdallNoVRFRFR/sqrt(length(allNoVRFR));    

    % mean stimulus firing
    %NO ODOUR
    allNoOdFR = FRbinned(allNoOdID);
    avrgallNoOdFR = nanmean(allNoOdFR);
    stdallNoOdFR = nanstd(allNoOdFR );
    SENoOd = stdallNoOdFR/sqrt(length(allNoOdFR));    

    %ODOUR
    allOdFR = FRbinned(allOdID);
    avrgallOdFR = nanmean(allOdFR);
    stdallOdFR = nanstd(allOdFR );
    SEOd = stdallOdFR/sqrt(length(allOdFR));  


    %%% running VR non event episodes
    alleventIDs = [allRewID; allLockedVRID; allSitVRID; allNoVRID; allNoOdID]%; allOdID];

    FRbinnedNonEvent = FRbinned;
    
    FRbinnedNonEvent(alleventIDs) = [];

    avrgallNonEvFR = nanmean(FRbinnedNonEvent);
    stdallNonEvFR = nanstd(FRbinnedNonEvent);
    SENonEv = stdallNonEvFR /sqrt(length(FRbinnedNonEvent)); 


    %Make a locomotion variable
    SpeedbinnedNonEvent = Speedbinned;
    SpeedbinnedNonEvent(alleventIDs) = [];
    avrgallSpeedbinnedNonEvent = nanmean(SpeedbinnedNonEvent);
    
%     Locomotionbinned = Speedbinned;
%     Locomotionbinned(allLockedVRID) = avrgallSpeedbinnedNonEvent;
%     Locomotionbinned(allNoVRID) = avrgallSpeedbinnedNonEvent;

    Locomotionbinned(1:binNs) = 1; %one is moving
    Locomotionbinned(allRewID) = 0;
    Locomotionbinned(allSitVRID) = 0;

    FRbinnedRunning = FRbinned;
    FRbinnedRunning([allRewID; allSitVRID]) = [];
    nanmean(FRbinnedRunning)
    

    FRbinnedSitting = [FRbinned(allRewID) FRbinned(allSitVRID)];
    nanmean(FRbinnedSitting)


    %%%% ALSO LOAD LFP
    %%%% THETA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TF = 1;

    LFPtheta = TFjx107b1_Ch12.values;

    LFPthetaTimesFrames = 1:length(LFPtheta);
    LFPthetaTimes = LFPthetaTimesFrames/1000;
    
    LFPthetaPower = abs(hilbert(LFPtheta));

    LFPthetaTimesPower = [LFPthetaTimes' LFPthetaPower];

    %%% BIN THE POWER FIRST INTO 400ms bins
    
    %%% do the same for speed
    currT = 0;
    for lfpBins = 1:binNs

        clear spkInds currthPOW

        currT = currT + 0.4;

        lfpIDs = find(LFPthetaTimes >= (currT-0.4) & LFPthetaTimes < currT);

        if isempty(lfpIDs) == 1
            LFPthetaPowerBinned(lfpBins) = 0;
        else
            currthPOW = LFPthetaPower(lfpIDs);
            LFPthetaPowerBinned(lfpBins) = nanmean(currthPOW,1);
        end
    end

    allRewIDth  = allRewID(allRewID < 856);

    allRewthPOW = LFPthetaPowerBinned(allRewIDth);
    avrgRewthPOW = nanmean(allRewthPOW);
    stdallRewthPOW = nanstd(allRewthPOW);
    SERewthPOW = stdallRewthPOW/sqrt(length(allRewthPOW));


    %mean locked VR running

    allLockedVRIDth  = allLockedVRID(allLockedVRID < 856);

    allLockedVRthPOW = LFPthetaPowerBinned(allLockedVRIDth);
    avrgallLockedVRthPOW= nanmean(allLockedVRthPOW);
    stdLockedVRthPOW = nanstd(allLockedVRthPOW);
    SELockedVRthPOW = stdLockedVRthPOW/sqrt(length(allLockedVRthPOW));


    %mean Sitting in VR
    allSitVRIDth  = allSitVRID(allSitVRID < 856);

    allSitthPOW = LFPthetaPowerBinned(allSitVRIDth);
    avrgallSitthPOW = nanmean(allSitthPOW);
    stdallSitthPOW = nanstd(allSitthPOW);
    SESitthPOW = stdallSitthPOW/sqrt(length(allSitthPOW));


    % mean no VR running
    allNoVRIDth  = allNoVRID(allNoVRID < 856);

    allNoVRthPOW = LFPthetaPowerBinned(allNoVRIDth);
    avrgallNoVRthPOW = nanmean(allNoVRthPOW);
    stdallNoVRthPOW = nanstd(allNoVRthPOW);
    SENoVRthPOW = stdallNoVRthPOW/sqrt(length(allNoVRthPOW));    

    % mean stimulus firing
    %NO ODOUR
    allNoOdIDth  = allNoOdID(allNoOdID < 856);

    allNoOdthPOW= LFPthetaPowerBinned(allNoOdIDth);
    avrgallNoOdthPOW = nanmean(allNoOdthPOW);
    stdallNoOdthPOW = nanstd(allNoOdthPOW);
    SENoOdthPOW = stdallNoOdthPOW/sqrt(length(allNoOdthPOW));    

    %ODOUR
    allOdIDth  = allOdID(allOdID < 856);

    allOdthPOW= LFPthetaPowerBinned(allOdIDth);
    avrgallOdthPOW = nanmean(allOdthPOW);
    stdallOdthPOW= nanstd(allOdthPOW);
    SEOdthPOW = stdallOdthPOW/sqrt(length(allOdthPOW));  


    alleventIDsth  = alleventIDs(alleventIDs < 856);
    LFPthetaPowerBinnedNonEvent = LFPthetaPowerBinned;
    LFPthetaPowerBinnedNonEvent(alleventIDsth) = [];

    avrgallNonEvthPOW = nanmean(LFPthetaPowerBinnedNonEvent);
    stdallNonEvthPOW = nanstd(LFPthetaPowerBinnedNonEvent);
    SENonEvthPOW = stdallNonEvthPOW/sqrt(length(LFPthetaPowerBinnedNonEvent));



    all_LFPthetaPower = [avrgallNoVRthPOW avrgRewthPOW avrgallSitthPOW avrgallLockedVRthPOW avrgallNonEvthPOW avrgallNoOdthPOW avrgallOdthPOW];

    all_LFPthetaPowerSE = [SENoVRthPOW SERewthPOW SESitthPOW SELockedVRthPOW SENonEvthPOW SENoOdthPOW SEOdthPOW];









    % MAKE GLM
    %Speedbinned(allNoVRID) = NaN;
    Speedbinned = Speedbinned';
    Locomotionbinned = Locomotionbinned';
    FRbinned = FRbinned';

    %TEST DIFFERENT MODELS: 

    GLMTAB = table(Speedbinned,RewVar,STVar,NoVRVar,LockedVRVar,FRbinned,'VariableNames',{'Speedbinned','RewVar','STVar','NoVRVar','LockedVRVar','FRbinned'});

    modelspec = 'FRbinned ~ Speedbinned*RewVar + STVar + NoVRVar + LockedVRVar';
    testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4 5]);

%%%%%%%%%%%%%%
%     GLMTAB = table(Speedbinned,RewVar,STVar,NoVRVar,LockedVRVar,Locomotionbinned,FRbinned,'VariableNames',{'Speedbinned','RewVar','STVar','NoVRVar','LockedVRVar', 'Locomotionbinned','FRbinned'});
% 
%     modelspec = ['FRbinned ~ Speedbinned*RewVar + Locomotionbinned  + STVar + NoVRVar  + LockedVRVar'];
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4 5 6]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
% 
% 
% 
%     GLMTAB = table(Speedbinned,RewVar,STVar,NoVRVar,LockedVRVar,Locomotionbinned,FRbinned,'VariableNames',{'Speedbinned','RewVar','STVar','NoVRVar','LockedVRVar', 'Locomotionbinned','FRbinned'});
% 
%     modelspec = ['FRbinned ~ Speedbinned*RewVar + Locomotionbinned  + STVar + NoVRVar  + LockedVRVar'];
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4 5 6]);





 
%     
%     FRbinned(allNoVRID) = [];
%     RewVar(allNoVRID) = [];
%     LockedVRVar(allNoVRID) = [];
%     STVar(allNoVRID) = [];
%     Speedbinned(allNoVRID) = [];
% 
% 
%     GLMTAB = table(Speedbinned,RewVar,STVar,LockedVRVar,FRbinned','VariableNames',{'Speedbinned','RewVar','STVar','LockedVRVar','FRbinned'});
% 
%     modelspec = 'FRbinned ~ Speedbinned*RewVar + STVar + LockedVRVar';
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4]);
%  


% 
%     GLMTAB = table(Locomotionbinned,RewVar,STVar,NoVRVar,LockedVRVar,FRbinned,'VariableNames',{'Locomotionbinned','RewVar','STVar','NoVRVar','LockedVRVar','FRbinned'});
% 
%     modelspec = 'FRbinned ~ Locomotionbinned*RewVar + STVar + NoVRVar + LockedVRVar';
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4 5]);
%  
% 
% 
% 
% 
% 
% 
% 
%     GLMTAB = table(Locomotionbinned,RewVar,STVar,NoVRVar,LockedVRVar,FRbinned,'VariableNames',{'Locomotionbinned','RewVar','STVar','NoVRVar','LockedVRVar','FRbinned'});
% 
%     modelspec = 'FRbinned ~ Locomotionbinned + RewVar + STVar + NoVRVar + LockedVRVar';
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3 4 5]);
% 
% 
% 
% 
%     GLMTAB = table(Locomotionbinned,RewVar,STVar,FRbinned,'VariableNames',{'Locomotionbinned','RewVar','STVar','FRbinned'});
% 
%     modelspec = 'FRbinned ~ Locomotionbinned*RewVar + STVar';
%     testglm = fitglm(GLMTAB,modelspec,'Distribution','poisson','CategoricalVars',[2 3]);
% 
% 
% 
% 
%  

    [rho p] = corr(Speedbinned,FRbinned)

    [rho p] = corr(Locomotionbinned,FRbinned)



    average_values = [avrgallNoVRFR avrgRewFR avrgallSitVRFR avrgallLockedVRFR avrgallNonEvFR avrgallNoOdFR avrgallOdFR];
    SE_alls = [SENoVR SERewFR SESitVR SELockedVR SENonEv SENoOd SEOd];




    [p h] = ranksum(allNoVRFR,allRewFR)
    [p h] = ranksum(allNoVRFR,allLockedVRFR)

    putallFR = [allNoVRFR allRewFR(1:119) allLockedVRFR FRbinnedNonEvent allNoOdFR(1:88)];
    catIDNOVR(1:length(allNoVRFR)) = 1;
    catIDRewFR(1:length(allRewFR)) = 2;
    catIDLockedVRFR(1:length(allLockedVRFR)) = 3;
    catIDNonEvFR(1:length(FRbinnedNonEvent)) = 4;
    catIDNoOdFR(1:length(allNoOdFR)) = 5;
    catallIDs = [catIDNOVR catIDRewFR catIDLockedVRFR catIDNonEvFR catIDNoOdFR];

    [p,tableanova,stats] = anovan(putallFR',catallIDs')
    [c,m] = multcompare(stats)



    putallTH = [allNoVRthPOW allRewthPOW allLockedVRthPOW LFPthetaPowerBinnedNonEvent(1:603) allNoOdthPOW];
    catIDNOVRTH(1:length(allNoVRthPOW)) = 1;
    catIDRewTH(1:length(allRewthPOW)) = 2;
    catIDLockedVRTH(1:length(allLockedVRthPOW)) = 3;
    catIDNonEvTH(1:length(LFPthetaPowerBinnedNonEvent)) = 4;
    catIDNoOdTH(1:length(allNoOdthPOW)) = 5;
    catallIDsTH = [catIDNOVRTH catIDRewTH catIDLockedVRTH catIDNonEvTH catIDNoOdTH];

    [p,tableanova,stats] = anovan(putallTH',catallIDsTH')
    [c,m] = multcompare(stats)

   [rh p] = corr(putallFR',putallTH')




    TF = 1




end