
function VRneuronTask
%%%% main analysis script for the odour/no-odour association task, can be used for single
%%%% cells or larger populations, run the script in the according folders
%%%% with the databases and textfiles (for example folder for all pyramidal
%%%% cells, interneurons or all PV basket cells

%%%% NOTE THIS SCRIPT CALCULATES MANY DIFFERENT OUTPUTS WHICH CAN BE THEN
%%%% USED FOR OTHER ANALYSIS, SUCH AS RAW FIRING RATES, DIFFERENTIATION SCORE,
%%%% THETA COUPLING ETC

%%% Required: circular statistics toolbox by Philipp Berens

clear all

%%% load text files with information which cells and according channels are
%%% to be used


Cells = importdata('CellNames.txt')
spikeCh = importdata('spikeChNames.txt')
LFPCh = importdata('LFPChNames.txt')
VRNames = importdata('VRNames.txt')
REWTHENDNAMES = importdata('REWTHENDNAMES.txt')

%load('sigIntsIDs.mat') %option to load only some signficant cells


allIDs = [1:151];%define how many cells

diffZ = []; FRZ = [];

%example cells = 5,43,109

for all_cellnum = 1:151%adjust according to cell numbers


    %cellnum = sigIntsID(all_cellnum);
    cellnum = allIDs(all_cellnum);

    
clearvars -except cellnum Cells spikeCh LFPCh VRNames REWTHENDNAMES normfiringtogetherplot ...
           absfiringdiffplot firingdiffplot manUconfvalues CellAngles NORHO NOPVALSpeedF ...
           ORHO OPVALSpeedF popmod2NOAngPhasPrecess popmod2OAngPhasPrecess rewNOpval rewOpval ...
           allCellssummedNOSPINFO allCellssummedOSPINFO pyrsignificances combinedfiring ...
           avrgNOREW2 avrgOREW2 NumberofNoOdourtrials NumberofOdourtrials RHOspeed PVALSpeedF ...
           avrgRHOspeed avrgPVALSpeedF totalspikesIn totalspikesOut avrgOspktrials avrgNOspktrials ...
           Oskptrialnumber NOskptrialnumber diffZ all_cellnum sigIntsID FRZ allIDs totalTimeNOodourSlow ...
           totalTimeOodourSlow meanTimeNOodourSlow meanTimeOodourSlow%FiringCorr

%clear all

% load the data from the text files

%load VR coordinate information
VR = load(VRNames{cellnum});


%infor for theta
REWTHEND = load(REWTHENDNAMES{cellnum});

thetaend = REWTHEND(2);

rewAnalysis = REWTHEND(1);



%%%%%%%%%% spike loading: 


%load MW03b1Cl33

filename = Cells{cellnum};

spikechannel = spikeCh{cellnum};
LFPchannel = LFPCh{cellnum};


datafile = load(filename)


spiketimesorig = datafile.(spikechannel).times;
thetaLFP =  datafile.(LFPchannel).values;


%%%%% restrict theta in recordings with labelling attempts

%thetaend = 2000;

thetaPhases = angle(hilbert(thetaLFP));

%thetaPhases = downsample(orithetaPhases,10);

thetaPower = abs(hilbert(thetaLFP));

%thetaPower = downsample(orithetaPower,10);

%plot(thetaPhases(1:10000))

%%% convert to angle as we use it

%thetaPhasesA = 180 + thetaPhases .*180/pi;

%%% make matrix with phases and spike 2 times

indexthetaLFP = find(thetaPhases);

thetatimes = indexthetaLFP .* 0.001;

%thetatimes(1:1000)

%thetatimes = downsample(orithetatimes,10);

%thetaPhasesAt = thetaPhasesA(1:10000);

thetamatrix = cat(2,thetaPhases,thetatimes);
thetamatrixpower = cat(2,thetaPower,thetatimes);


%end

%%%%%%%%%%%%%%%%%%%%%%%%% find trials

VRlength = size(VR,1);

trialcount = 0
for v1 = 1:(VRlength-1)
    
    yposDiff = abs(VR(v1+1,3) - VR(v1,3));
    
    if yposDiff > 800 && VR(v1+1,3) <= -2000 
        
        trialcount = trialcount + 1;        
        starttimes(trialcount,:) = VR(v1+1,:);
        
    elseif yposDiff > 800 && VR(v1+1,3) >= 9500 && VR(v1+1,3) <= 10000 
        
        trialcount = trialcount + 1;        
        starttimes(trialcount,:) = VR(v1+1,:);
        
        
    end
    
end

%%% subselect only those spikes after experiment started and before ended

indVRspks = find(spiketimesorig > starttimes(1,1) & spiketimesorig <= VR(VRlength,1));

spiketimes = spiketimesorig(indVRspks);


%%%%%%%%%%%%%%%%%%%%%%%% add speed values to VR


for stp = 1:(length(VR)-1)
        
        
        ys1 = VR(stp,3);
        ys2 = VR(stp+1,3);
        xs1 = VR(stp,2);
        xs2 = VR(stp+1,2); 
        
        %distance covered 
        stpdist = sqrt((ys2-ys1)*(ys2-ys1) + (xs2-xs1)*(xs2-xs1));
        
        VR(stp+1,4) = stpdist;
end


%%%%%%&%%%%%%%%%%% calculate trials

trialn = size(starttimes,1)

for t1 = 1:(trialn)

    trials(t1,1) = starttimes(t1,1);       
    if t1 == trialn        
        trials(t1,2) = VR(VRlength,1);        
    else                
        trials(t1,2) = starttimes(t1+1,1);        
    end
end

%%%%%%%%%%%%%%%%%% assign with odour or without odour matrix

for t3 = 1:trialn

   if starttimes(t3,3) > 0
       trialind(t3) = 1;
   else
       trialind(t3) = 0;
   end
end

odourInd = find(starttimes(:,3) > 0);
NOodourInd = find(starttimes(:,3) < 0);

%%%%%% trials end times 

NOendtimes = trials(NOodourInd,2);
Oendtimes = trials(odourInd,2);


%%%%%%%%%%%%%%%%%%%% assign coordinates to all spikes

spknumb = length(spiketimes)

for spk = 1:spknumb
    
        
        %%% check if any spike has the same time as coordinate
        
        equalind = find(VR(:,1) == spiketimes(spk));

        %%% if yes then the coordinates can be directly assigned

        if isempty(equalind) == 0
            
            x = VR(equalind,2)
            y = VR(equalind,3)
            z = VR(equalind,4)
            
            spkVR(spk,:) = [spiketimes(spk) x y z];
            
        %%%%% but usually not, so coordinates have to be interpolated

        else
            
            %select all VRdata that is greater than the current spike time
            vrindposA = find(VR(:,1) > spiketimes(spk));
            %select the first value which is then index of the next VR
            %datapoint
            vrA = vrindposA(1);
            %%% do the same for the data point before the spike time
            vrindposB = find(VR(:,1) < spiketimes(spk));
            vrB = vrindposB(length(vrindposB));

            x1 = VR(vrB,2);
            x2 = VR(vrA,2);
            y1 = VR(vrB,3);
            y2 = VR(vrA,3);
            
            z1 = VR(vrB,4);
            z2 = VR(vrA,4);
       

            spktimediff = spiketimes(spk) - VR(vrB);

            spktimepropor = spktimediff/0.5;
            distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        
              %%%%%%%% if the mouse doesnt move the coordinates will not change
               %%%%%% and the distance will be 0
               %%%% then the coordinates of the previous time have to be assigned
     
               if distance == 0 || VR(vrA,3) == 9550 || VR(vrA,3) == -2450

                    spkVR(spk,:) = [spiketimes(spk) VR(vrB,2) VR(vrB,3) 0];

                else

                    spkdistance = distance * spktimepropor;
                    currangle = asind((y2-y1)/distance);

                        if (x2-x1) < 0  
                            xmult = -1;
                        else
                            xmult = 1;
                        end

                    spkxdist = cosd(currangle)*spkdistance*(xmult);
                    spkydist = sind(currangle)*spkdistance*(1);

                    spkz = z1 + abs(z2-z1)/2;

                    spkx = x1 + spkxdist;
                    spky = y1 + spkydist;

                    spkVR(spk,:) = [spiketimes(spk) spkx spky spkz];

                end
             
        end
    
end

%%%%% assign phase values to each spike

% delete spikes which occur at low speeds

%fastindspkVR = find(spkVR(:,4) > 1);
%spkVRfast = spkVR(fastindspkVR,:);

%spknumbfast = length(spkVRfast(:,1))


%spktime2 = [1:spknumb]

spkphaseval = [1:spknumb]

for spkph = 1:spknumb
   
    
    spktime2 = spkVR(spkph,1);
    
    rndspktime2 = round(spktime2*1000);
    
    %thetatimes(rndspktime2)
    %roundedtheta = round(spktime2)
    %restricttheta = thetatimes((roundedtheta-1):(roundedtheta+1))
    
    %indtht1 = find(thetatimes >= spktime2); %& thetatimes < spktime2 + 0.01);
    %firstv1 = indtht1(1);
    
    %spktime2 = spkVR(1,1)
    
    %abstest = abs(spktime2-thetatimes);
   % abstest(3333)
    %abstest(3334)
    % abstest(3332)
    
   % thetatimes(248)
     
     %min(abstest)
   % [min indmind] = min(abs(spktime2-thetatimes)) 
    
   % spkphaseval(spkph) = thetamatrix(indmind,1);
         
    %spkphaseval(spkph) = thetamatrix(firstv1,1);
        
    spkphaseval(spkph) = thetamatrix(rndspktime2,1);
    
end

spkphaseval2 = spkphaseval';


%%%%% addd phase values to matrix


spkVR = cat(2,spkVR,spkphaseval2);

%%% thetaend restriction
%%% selected of VR

%%indthetaspkVR = find(spkVR(:,1) < thetaend)

%%spkphasevalRestrict = spkVR(1:length(indthetaspkVR),5);


%%%%%% add powervalue to VR matrix

%%%%%% exclude teleportation spikes


%%%%%%%%%%%%%%%%%%%% assign spikes with coordinates to trials

for t2 = 1:trialn
    
    ind = find(spkVR(:,1) >= trials(t2,1) & spkVR(:,1) < trials(t2,2));
    
    spkttrials{t2} = spkVR(ind,:);
    
    
end   

%%%%%%%%%%%%%%%%%%%% assign VR to trials

for t3 = 1:trialn
    
    ind = find(VR(:,1) >= trials(t3,1) & VR(:,1) < trials(t3,2));
    
    VRtrials{t3} = VR(ind,:);
    
    
end   


%%%%%%%%%%%% analyze trials without odour

%%%% select trials without odour

clear meanrad

NOodourtrials = spkttrials(NOodourInd);
NOVRtrials = VRtrials(NOodourInd);



NOsec = -2500:450:8750;

NOsn = length(NOsec);


NOvectorAngles = [];

for t4 = 1:length(NOodourInd)
    
    
    NOcurrVR = NOVRtrials{t4};
    
    NOcurrspk = NOodourtrials{t4};
    
    
    clear sectind
    clear NosecVR
    clear totalsectime
    clear corrsecttime
    clear spksectind
    clear NOsecspk
    clear NOslowspeedind
    clear NOfastspkind
    clear NOfastspknumbs
    clear NOfastspeedind
    clear speedfiltNOsecVR
    clear sectspeed
    
    clear NOrewardperiodind
    clear NOrewardperiod
    clear NOspkrewardindp1
    clear NOspkrewardindp2
    clear NOspkrewardindp3
    
   
    
    clear NOtrialmeanangle
    clear NOtrialmeanrad
        
    clear NOtrialmeanangle
    clear NOtrialmeanrad
    clear NOrvalueSector
    
    NOTrialAngles = [];
    
   
    
    for sc = 1:(NOsn-1)
        
        %clear NOsecAngles
    
        %%%% select the VR data of each segment
        sectind = find(NOcurrVR(:,3) >= NOsec(sc) & NOcurrVR(:,3) < NOsec(sc+1));        
        
        %%% if sector jump in the beginning
        if isempty(sectind) == 1
            
            
            NOtrialmeanrad(sc) = NaN;
            
            sectspeed(sc) = NaN;
            NOfastspknumbs(sc) = NaN;
            
        else
            
            NOsecVR = NOcurrVR(sectind,:);
            
        
            %%%%% select the spikes of each segment

            spksectind = find(NOcurrspk(:,3) >= NOsec(sc) & NOcurrspk(:,3) < NOsec(sc+1));

            NOsecspk = NOcurrspk(spksectind,:);
            

            %%% calculate total time spent in the segment WITH REFINED TIME
            %%% CALCULATION
            
            %%%% first start times
            
            if sc == 1 
                
                NoVRstart = NOsecVR(1,1);
                
            else
            
                NosecBdiffs = NOsecVR(1,3)- NOsec(sc);
                NosecPdiffs = NOsecVR(1,3) - NOcurrVR(sectind(1)-1,3);
                NoVRratioS = NosecBdiffs/NosecPdiffs;
            
                if NoVRratioS == 0
                          
                          
                     NoVRstart = NOsecVR(1,1);
                
                else
                
                    NoVRstart = NOsecVR(1,1) - 0.05*NoVRratioS;
            
                end
                
            end
            
            
            %%%%% calculate end times
            
            if sc > 24
                
                NoVRend = NOsecVR(length(sectind),1);
                
            else 
                
                NosecBdiffe = NOsec(sc+1) - NOsecVR(length(sectind),3);
                NosecPdiffe = NOcurrVR(sectind(length(sectind))+1,3) - NOsecVR(1,3);
                NoVRratioE = NosecBdiffe/NosecPdiffe;
                
           
           
                if NoVRratioE == 0
                
                     NoVRend = NOsecVR(length(sectind),1);
                else 
                
                     NoVRend = NOsecVR(length(sectind),1) + 0.05*NoVRratioE;
                end
                
                
            end
            
            
            totalsectime(sc) = NoVRend - NoVRstart;

            
            
            %%%%% filter and correct the time for when he is not moving

            NOslowspeedind = find(NOsecVR(:,4) < 1);
            NOslowpoints = length(NOslowspeedind);
            
            corrsecttime(sc) = totalsectime(sc) - 0.05*NOslowpoints;

            %%%%%% filter spikes for the same speed

            NOfastspkind = find(NOsecspk(:,4) >= 1);
            NOfastspknumbs(sc) = length(NOfastspkind);


            %%% ALSO SAVE NUMBER OF SPIKES TAKEN OUT
            NOslowspkind = find(NOsecspk(:,4) < 1);
            NOslowspknumbs(sc) = length(NOslowspkind);

            %%% also save time taken out per bin

            slowtime(sc) = 0.05*NOslowpoints;






            %%%%% filter speed values

            NOfastspeedind = find(NOsecVR(:,4) >= 1);
            speedfiltNOsecVR = NOsecVR(NOfastspeedind,4);
            sectspeed(sc) = mean(speedfiltNOsecVR);
            
            %%%%%% calculate average phase angle for sector
            
            
            
            if isempty(NOsecspk(:,5)) == 1
                
                
                NOtrialmeanrad(sc) = NaN;
                NOrvalueSector(sc) = NaN;
                NOsecAngles = NaN;
                
            elseif NOsecspk(1,1) < thetaend
               
                %%% variable with all spike angles in that sector
                
                
                NOsecAngles = NOsecspk(:,5);
                                
                %%% mean angle calculated for all sectors of 1 trial
                NOtrialmeanrad(sc) = circ_mean(NOsecAngles);
                
                %%% r-values of 1 trial
                NOrvalueSector(sc) = circ_r(NOsecAngles);
                
                
                %%% NOspikeAngles gathers all angles of spikes in all
                %%% nonodour trials
                NOvectorAngles = cat(1,NOvectorAngles,NOsecAngles);
                
            else
                
                NOtrialmeanrad(sc) = NaN;
                NOrvalueSector(sc) = NaN;
                NOsecAngles = NaN;
                
            end
    
            %Speed selection? for what?
            indNOfastspeedindANGLE = find(NOsecspk(:,4) >= 1 & NOsecspk(:,4) < 70);
            
            NOsecAnglesNew = NOsecspk(indNOfastspeedindANGLE,[3 5]);
            
            %%%%%%%%% thetapowerextraction
            clear NOindthetapower
            clear currNOthetapower
            
            NOindthetapower = find(thetamatrixpower(:,2)> NoVRstart & thetamatrixpower(:,2)< NoVRend);
            currNOthetapower = thetamatrixpower(NOindthetapower,1);
            meancurrNOthetapower = mean(currNOthetapower);
            
        end
        
        trialNOthetapower(sc) = meancurrNOthetapower;
        
        
        NOTrialAngles = cat(1,NOTrialAngles,NOsecAnglesNew);
        
        
                
    end
    
    
    
    overallNOthpow(t4,:) = trialNOthetapower;
    
    NOSpkAngles{t4} = NOTrialAngles;
    
    %meanangle(t4,:) = trialmeanangle
    NOmeanrad(t4,:) = NOtrialmeanrad;
    NOrvalues(t4,:) = NOrvalueSector;
    
    
    %%% TO CALCULATE THE NUMBER OF SPIKES DO NOT DIVIDE BY TIME:
    %NOfiringrate(t4,:) = NOfastspknumbs./corrsecttime;
    
    %uncomment for counting spikes:
    NOfiringrate(t4,:) = NOfastspknumbs%./corrsecttime;

    %SAVE ALSO SLOW SPIKE NUMBERS
    NOfiringrateSLOW(t4,:) = NOslowspknumbs;

    %Save also times in slow 
    NOtimesSLOW(t4,:) = slowtime;


    NOspeed(t4,:) = sectspeed;
    
    
    %%%%%%% filter out reward periods
            
                    
            %%%%% endtimes minus 5 seconds is start time
            
            NOendtime = NOendtimes(t4);
            
            %%%% reward period1 times
            
            NOrewardperiod1start = NOendtime - 5;
            NOrewardperiod1end = NOrewardperiod1start + 1; 

                        
            %%%%%%% select spikes FOR rewardperiod 1
                
            NOspkrewardindp1 = find(NOsecspk(:,1) > NOrewardperiod1start & NOsecspk(:,1) < NOrewardperiod1end);
            NOspkrewardperiod1 = length(NOspkrewardindp1);
            
            
            %%%% reward period2 times
            
            NOrewardperiod2start =  NOrewardperiod1end;
            NOrewardperiod2end =  NOrewardperiod2start + 1;
            
            %%%%%%% select spikes FOR rewardperiod 2
                
            NOspkrewardindp2 = find(NOsecspk(:,1) > NOrewardperiod2start & NOsecspk(:,1) < NOrewardperiod2end);
            NOspkrewardperiod2 = length(NOspkrewardindp2);
            
            
            %%%% reward period3 times
            
            NOrewardperiod3start =  NOrewardperiod2end;
            NOrewardperiod3end =  NOrewardperiod3start + 3;
            
            %%%%%%% select spikes FOR rewardperiod 2
                
            NOspkrewardindp3 = find(NOsecspk(:,1) > NOrewardperiod3start & NOsecspk(:,1) < NOrewardperiod3end);
            NOspkrewardperiod3 = length(NOspkrewardindp3);
            
            
            
            
            %%%% calculate firing rate during reward periods
%             NOFRatep1(t4) = NOspkrewardperiod1/1;
%             NOFRatep2(t4) = NOspkrewardperiod2/1;  
%             NOFRatep3(t4) = NOspkrewardperiod3/3;


            %%% SPIKES INSTEAD OF FR FROM REWARD PERIOD
            NOFRatep1(t4) = NOspkrewardperiod1;%/1;
            NOFRatep2(t4) = NOspkrewardperiod2;%/1;  
            NOFRatep3(t4) = NOspkrewardperiod3;%/3;
    
    
end

NOFRatep1 = NOFRatep1';

NOFRatep2 = NOFRatep2';

NOFRatep3 = NOFRatep3';

NOOverallREWFR = cat(2,NOFRatep1,NOFRatep2,NOFRatep3);



%%% also calculate average NOrewfiring and select highest number
meanNOREWFR = mean(NOOverallREWFR,1);
maxmeanNOREWFR = max(meanNOREWFR);

maxNOIND = find(maxmeanNOREWFR == meanNOREWFR);

if length(maxNOIND) > 1 
    
    maxNOIND = maxNOIND(1);
    
end

maxNOREWFIRING = NOOverallREWFR(:,maxNOIND);

%avergaged trial firing
meandurNOfiring = mean(NOfiringrate(:,4:end),2);

%%%% check if max firing rew greater than during task firing

if meanNOREWFR(maxNOIND)> mean(meandurNOfiring)
    
    currNOpval = ranksum(maxNOREWFIRING,meandurNOfiring);
    if currNOpval < 0.05
        rewNOpval(cellnum) = 1;
    else
        rewNOpval(cellnum) = 0.5;
    end
else 
    rewNOpval(cellnum) = 0.1;
end




%%% calculate NOOdourReward 
 
%NORewP1FR = mean(NOFRatep1) 

%NORewP2FR = mean(NOFRatep2) 

%NORewP3FR = mean(NOFRatep3) 



%figure
%hold on

%overallNOthpow

%plot(overallNOthpow(5,:))



avrgNOfiringrate = mean(NOfiringrate,1);
avrgNOspeed = mean(NOspeed,1);

%figure
%hold on

%plot(avrgNOfiringrate(2:25),'r')
%plot(avrgNOspeed(2:25))




%%%%%%%  plot phase to position to check phase precession
%figure
%hold on

ANGLNOsec = -2500:250:8750;
ANGLNOsn = length(ANGLNOsec);
NoRadValues = {};

Noy2 = 1 

allNox = []
allNoy = []

for aglNo = 1:length(NOodourInd)

    clear Notrial
    clear Nox
    clear Noy
    
    Notrial = NOSpkAngles{aglNo};
    if isempty(Notrial) == 1
        
        Nox = 1
        Noy = -1
        
    else
        Nox = Notrial(:,1);
        Noy = Notrial(:,2);
    
    end
    
    %Nox = Notrial(:,1);
    %Noy = Notrial(:,2);
    %Noy = Noy + 10 *aglNo  
    Noy2 = Noy2+1;
    for noyn = 1:length(Noy);
        
        if Noy(noyn) > 0
            Noy(noyn) = Noy(noyn) - 6.283;
        end
    end
    %if aglNo > 6 && aglNo < 10   
    %plot(Nox,Noy,'LineStyle','none','Marker','.')
    %caxis([-2500 8500])
    
    %plot(Nox,Noy2,'LineStyle','none','Marker','.','MarkerSize',1,'Color','blue')
    %end
    
    %%%plot pSTH
    %line([Nox Nox], [Noy2-1 Noy2],'Color','blue')
    
    clear anglsecInd
    clear secNoyVals
    
    for AsecNO = 1:(ANGLNOsn-1) 
        
        %%%%% calculate mean phases in defined spatial sectors
        anglsecInd = find(Nox > (-2500 + 250*(AsecNO-1)) & Nox < (-2500 + 250*(AsecNO))) 
        
        if isempty(anglsecInd) == 1
            secNoyVals{AsecNO} = {};
            %Anglemeansecvals(AsecNO) = NaN
        else
                        
            secNoyVals{AsecNO} = Noy(anglsecInd);
            %meansecvals = circ_mean(secNoyVals)
            %Anglemeansecvals(AsecNO) = circ_rad2ang(meansecvals)
        
        end
    end
    
    NoRadValues = [NoRadValues; secNoyVals];
    
    allNox = cat(1,allNox,Nox)
    allNoy = cat(1,allNoy,Noy) 
    
end



for NOSA = 1:45   
 
 clear allNOval   
 clear allNOval2
 clear NOradVal
    
 allNOval1 = cat(1,NoRadValues{:,NOSA});
 
 
 if isempty(allNOval1) == 1
     
   allNOval2 = [];
   NOAnglePhasePrecVals(NOSA) = NaN;
 else
     
    if iscell(allNOval1) == 1
        allNOval2 = allNOval1{:};
     else 
         allNOval2 = allNOval1;
    end
    
    NOradVal = circ_mean(allNOval2);
    NOAnglePhasePrecVals(NOSA) = circ_rad2ang(NOradVal);
 end
 
 
 
 %NOAnglePhasePrecVals(NOSA) = circ_rad2ang(NOradVal)

end

modNOAnglesPhasePrec = NOAnglePhasePrecVals + 180;

%for modAng = 1:length(modNOAnglesPhasePrec)
    
    %if modNOAnglesPhasePrec(modAng) < 180 
        
       % mod2NOAnglesPhasePrec(modAng)= modNOAnglesPhasePrec(modAng) + 360
    %else 
       % mod2NOAnglesPhasePrec(modAng)= modNOAnglesPhasePrec(modAng)

    %end
    
%end

mod2NOAnglesPhasePrec = modNOAnglesPhasePrec;

for modAng = 2:(length(modNOAnglesPhasePrec)-1)
    
    if abs((modNOAnglesPhasePrec(modAng)+360) - mod2NOAnglesPhasePrec(modAng-1)) < abs((modNOAnglesPhasePrec(modAng)) - mod2NOAnglesPhasePrec(modAng-1)) 
        
        mod2NOAnglesPhasePrec(modAng)= modNOAnglesPhasePrec(modAng) + 360;
        
    else 
        
        mod2NOAnglesPhasePrec(modAng)= modNOAnglesPhasePrec(modAng);

    end
    
    %if mod2NOAnglesPhasePrec(modAng) > 540 
        
       % mod2NOAnglesPhasePrec(modAng) =  mod2NOAnglesPhasePrec(modAng)-360
        
    %elseif mod2NOAnglesPhasePrec(modAng) < 180
        
        %mod2NOAnglesPhasePrec(modAng) =  mod2NOAnglesPhasePrec(modAng)+360
   %end
   
end

popmod2NOAngPhasPrecess(cellnum,:) = mod2NOAnglesPhasePrec;


%%%%%%%%%%%% analyze trials with odour
clear avrgOfiringrate
clear Ofiringrate
clear Omeanrad
%%%% select trials withodour

odourtrials = spkttrials(odourInd);
oVRtrials = VRtrials(odourInd);

%spkttrials{2}
%odourtrials{2}
%oVRtrials{9}

Osec = 9500:450:18550;

Osn = length(Osec);

OvectorAngles = [];

%OTrialAngles = []

%%% length(odourInd)

for t5= 1:length(odourInd)
    
    
    OcurrVR = oVRtrials{t5};
    
    Ocurrspk = odourtrials{t5};
    
    
    clear Osectind
    clear OsecVR
    clear Ototalsectime
    clear Ocorrsecttime
    clear Ospksectind
    clear Osecspk
    clear Oslowspeedind
    clear Ofastspkind
    clear Ofastspknumbs
    clear Osectspeed
    clear speedfiltOsecVR
    clear Ofastspeedind
    
    clear Orewardperiodind
    clear Orewardperiod
    clear Ospkrewardindp1
    clear Ospkrewardindp2
        
    
    clear Otrialmeanangle
    clear Otrialmeanrad
    clear OrvalueSector
    
    OTrialAngles = []
    
    
    
    for osc = 1:(Osn-1)
    
        %%%% select the VR data of each segment
        Osectind = find(OcurrVR(:,3) >= Osec(osc) & OcurrVR(:,3) < Osec(osc+1)); 
        
        %%% if sector jump in the beginning
        if isempty(Osectind) == 1
            
            
            Otrialmeanrad(osc) = NaN;
            
            Osectspeed(osc) = NaN;
            Ofastspknumbs(osc) = NaN;
            Ocorrsecttime(osc) = NaN;
            
            
        else
        
            OsecVR = OcurrVR(Osectind,:);

            %%%%% select the spikes of each segment

            Ospksectind = find(Ocurrspk(:,3) >= Osec(osc) & Ocurrspk(:,3) < Osec(osc+1));

            Osecspk = Ocurrspk(Ospksectind,:);


            %%% calculate total time spent in the segment

            
            %%% calculate total time spent in the segment WITH REFINED TIME
            %%% CALCULATION
            
            %%%% first start times
            
            if osc == 1 
                
                oVRstart = OsecVR(1,1);
                
            else
            
                osecBdiffs = OsecVR(1,3)- Osec(osc);
                osecPdiffs = OsecVR(1,3) - OcurrVR(Osectind(1)-1,3);
                oVRratioS = osecBdiffs/osecPdiffs;
            
                if oVRratioS == 0
                          
                          
                     oVRstart = OsecVR(1,1);
                
                else
                
                    oVRstart = OsecVR(1,1) - 0.05*oVRratioS;
            
                end
                
            end
            
            
            %%%%% calculate end times
            
            if osc > 19
                
                oVRend = OsecVR(length(Osectind),1);
                
            else 
                
                osecBdiffe = Osec(osc+1) - OsecVR(length(Osectind),3);
                osecPdiffe = OcurrVR(Osectind(length(Osectind))+1,3) - OsecVR(1,3);
                oVRratioE = osecBdiffe/osecPdiffe;
                
           
           
                if oVRratioE == 0
                
                     oVRend = OsecVR(length(Osectind),1);
                else 
                
                     oVRend = OsecVR(length(Osectind),1) + 0.05*oVRratioE;
                end
                
                
            end
            
            
            Ototalsectime(osc) = oVRend - oVRstart;
            
            %Ototalsectime(osc) = OsecVR(length(Osectind),1)- OsecVR(1,1);

            %%%%% filter and correct the time for when he is not moving

            Oslowspeedind = find(OsecVR(:,4) < 1);
            Oslowpoints = length(Oslowspeedind);

            Ocorrsecttime(osc) = Ototalsectime(osc) - 0.05*Oslowpoints;

            Oslowtime(osc) = 0.05*Oslowpoints;

            %%%%%% filter spikes for the same speed

            Ofastspkind = find(Osecspk(:,4) >= 1);
            Ofastspknumbs(osc) = length(Ofastspkind);

            %ALSO GET SLOW SPIKES
            Oslowspkind = find(Osecspk(:,4) < 1);
            Oslowspknumbs(osc) = length(Oslowspkind);



            %%%%% filter speed values

            Ofastspeedind = find(OsecVR(:,4) >= 1);
            speedfiltOsecVR = OsecVR(Ofastspeedind,4);
            Osectspeed(osc) = mean(speedfiltOsecVR);
            
            %%%%%% calculate average phase angle for sector
            
            if isempty(Osecspk(:,5)) == 1
                
                Otrialmeanrad(osc) = NaN;
                OrvalueSector(osc) = NaN;
                OsecAngles = NaN;
                
            elseif Osecspk(1,1) < thetaend
            
                OsecAngles = Osecspk(:,5);
            
                
                %%% nonodour trials
                OvectorAngles = cat(1,OvectorAngles,OsecAngles);
            
            else
                
                Otrialmeanrad(osc) = NaN;
                OrvalueSector(osc) = NaN;
                OsecAngles = NaN;
            
            
            
            
            end
            
            %%%% select speed
            indOfastspeedindANGLE = find(Osecspk(:,4) >= 1 & Osecspk(:,4) < 70);
            
            OsecAnglesNew = Osecspk(indOfastspeedindANGLE,[3 5]);

            Otrialmeanrad(osc) = circ_mean(OsecAngles);
            OrvalueSector(osc) = circ_r(OsecAngles);
        
        end
                
        
        OTrialAngles = cat(1,OTrialAngles,OsecAnglesNew);
        
        
    end
    
    OSpkAngles{t5} = OTrialAngles;
    
    
    Omeanrad(t5,:) = Otrialmeanrad;
    Orvalues(t5,:) = OrvalueSector;
    
    %Ofiringrate(t5,:) = Ofastspknumbs./Ocorrsecttime;
    %uncommment for counting spikes:
    Ofiringrate(t5,:) = Ofastspknumbs;%./Ocorrsecttime;
    OfiringrateSLOW(t5,:) = Oslowspknumbs;



        %Save also times in slow 
    OtimesSLOW(t5,:) = Oslowtime;

    Ospeed(t5,:) = Osectspeed;
    
        
    %%%%%%% filter out reward periods
            
            
            %%%%% endtimes minus 5 seconds is start time
            
            Oendtime = Oendtimes(t5);
            
            %%%% reward period1 times
            
            Orewardperiod1start = Oendtime - 5;
            Orewardperiod1end = Orewardperiod1start + 1; 

           
            %%%% reward period2 times
            
            Orewardperiod2start =  Orewardperiod1end;
            Orewardperiod2end =  Orewardperiod2start + 1;
            
            
            %%%% reward period3 times
            
            Orewardperiod3start =  Orewardperiod2end;
            Orewardperiod3end =  Orewardperiod3start + 3;
            
            
            %%%%%%% select spikes FOR rewardperiod 1
                
            Ospkrewardindp1 = find(Osecspk(:,1) > Orewardperiod1start & Osecspk(:,1) < Orewardperiod1end);
            Ospkrewardperiod1 = length(Ospkrewardindp1);
            
            %%%%%%% select spikes FOR rewardperiod 2
                
            Ospkrewardindp2 = find(Osecspk(:,1) > Orewardperiod2start & Osecspk(:,1) < Orewardperiod2end);
            Ospkrewardperiod2 = length(Ospkrewardindp2);
            
            %%%%%%% select spikes FOR rewardperiod 3
                
            Ospkrewardindp3 = find(Osecspk(:,1) > Orewardperiod3start & Osecspk(:,1) < Orewardperiod3end);
            Ospkrewardperiod3 = length(Ospkrewardindp3);
            
            
            
            %%%% calculate firing rate during reward periods
            
%             OFRatep1(t5) = Ospkrewardperiod1/1;
%             OFRatep2(t5) = Ospkrewardperiod2/1;     
%             OFRatep3(t5) = Ospkrewardperiod3/3;   

            OFRatep1(t5) = Ospkrewardperiod1;%/1;
            OFRatep2(t5) = Ospkrewardperiod2;%/1;     
            OFRatep3(t5) = Ospkrewardperiod3;%/3;   
        
  
end

%ORewP1FR = mean(OFRatep1) 

%ORewP2FR = mean(OFRatep2) 

%ORewP3FR = mean(OFRatep3) 

OFRatep1 = OFRatep1';
OFRatep2 = OFRatep2';
OFRatep3 = OFRatep3';

OverallREWFR = cat(2,OFRatep1,OFRatep2,OFRatep3);




%%%% OPTION TO CALCULATE MEAN AND STD FOR ZSCORE VALUES!!!

%First put all into one big vector
AllRewFRVec = [OverallREWFR(:); NOOverallREWFR(:)];
OVSELFIRING = (Ofiringrate(:,4:end)); OVECFIRING = OVSELFIRING(:);
NOVSELFIRING = (NOfiringrate(:,4:end)); NOVECFIRING = NOVSELFIRING(:);

ALLVECFIRING = [AllRewFRVec; OVECFIRING; NOVECFIRING];


%use these values to do a zscore normalisation
meanFRall = nanmean(ALLVECFIRING);
stdFRall = nanstd(ALLVECFIRING);

OverallREWFR_Z = (OverallREWFR-meanFRall)/stdFRall;  
NOOverallREWFR_Z = (NOOverallREWFR-meanFRall)/stdFRall; 



%%% DO A PEAK NORMALIZSATION

    avrgcutOfiring = nanmean(OVSELFIRING,1);
    avrgcutNOfiring = nanmean(NOVSELFIRING,1); 
    avrgREWfiring = nanmean(OverallREWFR,1);

    AllavrgFiring = [avrgcutOfiring avrgcutNOfiring avrgREWfiring];
    maxAllFR = max(AllavrgFiring); %USE THIS VALUE TO NORMALIZE TO PEAK

%%%%%%%    








meanOREWFR = mean(OverallREWFR,1);
maxmeanOREWFR = max(meanOREWFR);

maxOIND = find(maxmeanOREWFR == meanOREWFR);

if length(maxOIND) > 1 
    
    maxOIND = maxOIND(1);
    
end

maxOREWFIRING = OverallREWFR(:,maxOIND);

%avergaged trial firing
meandurOfiring = mean(Ofiringrate(:,4:end),2);

%%%% check if max firing rew greater than during task firing

if meanOREWFR(maxOIND)> mean(meandurOfiring)
    
    currOpval = ranksum(maxOREWFIRING,meandurOfiring);
    if currOpval < 0.05
        rewOpval(cellnum) = 1;
       else
        rewOpval(cellnum) = 0.5;
    end
else 
    rewOpval(cellnum) = 0.1;
end







%%%%%%%  plot phase to position to check phase precession

%figure
%hold on

ANGLOsec = 9500:250:18550;
ANGLOsn = length(ANGLOsec);
oRadValues = {};

oy2 = 1


for aglo = 1:length(odourInd)

    clear otrial
    clear ox
    clear oy
    
    
    
    otrial = OSpkAngles{aglo};
    ox = otrial(:,1);
    oy = otrial(:,2);
    %oy = oy + 10 *aglo;    
    oy2 = oy2 + 1;
    for oyn = 1:length(oy)
        
        if oy(oyn) > 0
            
            oy(oyn) = oy(oyn) - 6.283;
            
        end
        
    end

   %plot(ox,oy,'LineStyle','none','Marker','.','Color','red')
    %line([ox ox], [oy-1 oy],'Color','red')
    
   %line([ox ox], [oy2-1 oy2],'Color','red')
        
    clear OanglsecInd
    clear secoyVals
    
    
        for AsecO = 1:(ANGLOsn-1) 
        
        %%%%% calculate mean phases in defined spatial sectors
        OanglsecInd = find(ox > (9500 + 250*(AsecO-1)) & ox < (9500 + 250*(AsecO))); 
        
        if isempty(OanglsecInd) == 1
            secoyVals{AsecO} = {};
            
        else
                        
            secoyVals{AsecO} = oy(OanglsecInd);
            
        end
    end
    
    oRadValues = [oRadValues; secoyVals];


end

for OSA = 1:36   
 
 clear allOval   
 clear allOval2
 clear OradVal
    
 allOval1 = cat(1,oRadValues{:,OSA});
 
 
 if isempty(allOval1) == 1
     
   allOval2 = [];
   OAnglePhasePrecVals(OSA) = NaN;
 else
     
    if iscell(allOval1) == 1
        allOval2 = allOval1{:};
     else 
         allOval2 = allOval1;
    end
    
    OradVal = circ_mean(allOval2);
    OAnglePhasePrecVals(OSA) = circ_rad2ang(OradVal);
 end

end

modOAnglesPhasePrec = OAnglePhasePrecVals + 180;

modRadiansPhasePrec = circ_ang2rad(modOAnglesPhasePrec);

%figure
%rho = 1:36
%polar(modRadiansPhasePrec,rho)



%for modAng = 1:length(modOAnglesPhasePrec)
    
    %if modOAnglesPhasePrec(modAng) < 180 
        
        %mod2OAnglesPhasePrec(modAng)= modOAnglesPhasePrec(modAng) + 360
    %else 
       % mod2OAnglesPhasePrec(modAng)= modOAnglesPhasePrec(modAng)

    %end
    
%end

mod2OAnglesPhasePrec = modOAnglesPhasePrec;

for modAng = 2:(length(modOAnglesPhasePrec)-1)
    
    if abs((modOAnglesPhasePrec(modAng)+360) - mod2OAnglesPhasePrec(modAng-1)) < abs((modOAnglesPhasePrec(modAng)) - mod2OAnglesPhasePrec(modAng-1)) 
        
        mod2OAnglesPhasePrec(modAng)= modOAnglesPhasePrec(modAng) + 360;
        
    else 
        
        mod2OAnglesPhasePrec(modAng)= modOAnglesPhasePrec(modAng);

    end
    
    %if mod2OAnglesPhasePrec(modAng) > 540 
        
        %mod2OAnglesPhasePrec(modAng) =  mod2OAnglesPhasePrec(modAng)-360
        
    %elseif mod2OAnglesPhasePrec(modAng) < 180
        
        %mod2OAnglesPhasePrec(modAng) =  mod2OAnglesPhasePrec(modAng)+ 360
        
        
        
    %end
    
end


popmod2OAngPhasPrecess(cellnum,:) = mod2OAnglesPhasePrec;








avrgOfiringrate = mean(Ofiringrate,1);
avrgOspeed = mean(Ospeed,1);



%figure
%hold on

%plot(avrgOfiringrate(2:20),'r')
%plot(avrgOspeed(2:20))


%figure
%hold on

%%plot(avrgOfiringrate(3:20),'r')
%%%plot(avrgNOfiringrate(3:25),'b')


%plot(avrgNOspeed(2:25),'g')
%plot(avrgOspeed(2:20),'m')

%%%%%%%% phase coupling histogramm values

%%spkphaseval3 = 180 + spkphasevalRestrict.*180/pi;

%%meanphase = circ_mean(spkphaseval2)
%%meanphaseangle = 180 + meanphase*180/pi


%%%%%%%%%%%%%%%%%%  calculate mean angles in sectors excluding NAN values

for NOang = 1:25
   
    NOmeanradcolsect = NOmeanrad(:,NOang);
    NOradsindsect = find(isnan(NOmeanradcolsect)==0);
    NOradselectedsect = NOmeanradcolsect(NOradsindsect);
    
    NOselectaverg(NOang) = 180 + circ_mean(NOradselectedsect)*180/pi;
    
end

for Oang = 1:20
   
    Omeanradcolsect = Omeanrad(:,Oang);
    Oradsindsect = find(isnan(Omeanradcolsect)==0);
    Oradselectedsect = Omeanradcolsect(Oradsindsect);
    
    Oselectaverg(Oang) = 180 + circ_mean(Oradselectedsect)*180/pi;
    
end

%figure
%hold on
%plot(NOselectaverg)
%plot(Oselectaverg,'r')




%calculate mean angles for NONodour and Odour from single spike values


%%% calculate the mean angles for NO Odour trials and Odour trials 
%%% and all trial meanangle

NOOdourmeanAngle = 180 + circ_mean(NOvectorAngles)*180/pi;



OOdourmeanAngle = 180 + circ_mean(OvectorAngles)*180/pi;


alltrialradAngles = cat(1,NOvectorAngles, OvectorAngles);

TRIALMEANANGLE = 180 + circ_mean(alltrialradAngles)*180/pi;
trialrvalue = circ_r(alltrialradAngles);


alltrialangleAngles = 180 + alltrialradAngles.*180/pi;


ast = -18;
aend = 0;

for binang = 1:20
    
    clear Aspkind
    
    ast = ast + 18;
    aend = aend + 18;
    Aspkind = find(alltrialangleAngles >= ast & alltrialangleAngles < aend);
    Aspknumb(binang) = length(Aspkind);
    
end

%%%%% calculated with weighted bins

mphbins = [9 27 45 63 81 99 117 135 153 171 189 207 225 243 261 279 297 315 333 351];

mphbinsrad = circ_ang2rad(mphbins)';

meanphase2  = circ_mean(mphbinsrad,Aspknumb');

binnrvalue = circ_r(mphbinsrad,Aspknumb');

binnedpvalue = circ_rtest(mphbinsrad,Aspknumb');

BINNEDTRIALMEANANGLE = meanphase2*180/pi;






%%% calculate averages in rows without NaNs

for NoR = 1:25
    
    NOROW = NOmeanrad(:,NoR);
    NOindROW = find(isnan(NOROW)==0);
    NOROWSELECTED = NOROW(NOindROW);
    NOWROWavr(NoR) = 180 + circ_mean(NOROWSELECTED)*180/pi;

end



for oR = 1:20
    
    OROW = Omeanrad(:,oR);
    OindROW = find(isnan(OROW)==0);
    OROWSELECTED = OROW(OindROW);
    OWROWavr(oR) = 180 + circ_mean(OROWSELECTED)*180/pi;

end

clear cutavrgNOFiringrate
clear cutavrgOFiringrate

cutavrgNOFiringrate = avrgNOfiringrate(:,4:25);
cutavrgOFiringrate = avrgOfiringrate(:,4:20);

%%%% add reward period firing, first avergae reward periods!

clear cutrunrewNOFiring
clear cutrunrewOFiring

if rewAnalysis == 2%changed to dummy 2 for time purposes
    
    avrgNOREW2(cellnum,:) = mean(NOOverallREWFR,1);
    avrgOREW2(cellnum,:) = mean(OverallREWFR,1);    
    
    avrgNOREW = mean(NOOverallREWFR,1);
    avrgOREW = mean(OverallREWFR,1); 
    
    cutrunrewNOFiring = cat(2,cutavrgNOFiringrate,avrgNOREW);
    cutrunrewOFiring = cat(2,cutavrgOFiringrate,avrgOREW);
    

else
    
    avrgNOREW2(cellnum,:) = -0.1;
    avrgOREW2(cellnum,:) = -0.1; 
    
    avrgNOREW(1:3) = -0.1;
    avrgOREW(1:3) =  -0.1;
    
    
    cutrunrewNOFiring = cat(2,cutavrgNOFiringrate,avrgNOREW);
    cutrunrewOFiring = cat(2,cutavrgOFiringrate,avrgOREW);

end
    
    

%%%% normalized average Firing rate sector 4:end

clear cutallFiring

cutallFiring = cat(2,cutrunrewNOFiring,cutrunrewOFiring);

maxFiring = max(cutallFiring);

normcutavrgNOFiringrate = cutrunrewNOFiring./maxFiring;
normcutavrgOFiringrate = cutrunrewOFiring./maxFiring;

if rewAnalysis == 1

    normcutavrgNOFiringrate(24:26) = normcutavrgNOFiringrate(23:25);
    normcutavrgNOFiringrate(23) = -0.1;

    normcutavrgOFiringrate(24:26) = normcutavrgOFiringrate(18:20);
    normcutavrgOFiringrate(18:23) = -0.1;

    normfiringtogether = cat(1,normcutavrgNOFiringrate,normcutavrgOFiringrate);

    
else 
    
    normcutavrgNOFiringrate(23:26) = -0.1;
    normcutavrgOFiringrate(18:26) = -0.1;
    normfiringtogether = cat(1,normcutavrgNOFiringrate,normcutavrgOFiringrate);

    
end


%%%% difference in firing (odour - no odour)

%%% for each cell the difference between odour and no-odour trials are 
%%% calculated 

firingdiff = normcutavrgOFiringrate(1:17) - normcutavrgNOFiringrate(1:17);

firingdiffplot(cellnum,:) = firingdiff;


absfiringdiff = abs(firingdiff);


normfiringtogetherplot((cellnum*2-1):(cellnum*2),:) = normfiringtogether;
absfiringdiffplot(cellnum,:) = absfiringdiff;


togetherfiringdiff = cat(1,firingdiff,absfiringdiff);


%save(['savedvar' filename],'cutrunrewNOFiring','cutrunrewOFiring','normcutavrgNOFiringrate','normcutavrgOFiringrate','togetherfiringdiff')






cutcurrNOfiring = NOfiringrate(:,4:25)
cutNOSpeed = NOspeed(:,4:25)
cutcurrOfiring = Ofiringrate(:,4:20)
cutOSpeed = Ospeed(:,4:20)

%%%% correlation of speed and firing analysis

preRAWvecNOFiring = cutcurrNOfiring';
RAWvecNOFiring = preRAWvecNOFiring(:);

preRAWvecNOSpeed = cutNOSpeed';
RAWvecNOSpeed = preRAWvecNOSpeed(:);

preRAWvecOFiring = cutcurrOfiring';
RAWvecOFiring = preRAWvecOFiring(:);

preRAWvecOSpeed = cutOSpeed';
RAWvecOSpeed = preRAWvecOSpeed(:);

allvecFiring = cat(1,RAWvecNOFiring,RAWvecOFiring)
allvecSpeed = cat(1,RAWvecNOSpeed,RAWvecOSpeed)

[RHOspeed(cellnum),PVALSpeedF(cellnum)] = corr(allvecFiring,allvecSpeed)





%%%% apply GLM if p-value significant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CONVERT TO ZSCORED VALUES FIRST

allvecFiringZ = (allvecFiring-meanFRall)/stdFRall;
cutcurrOfiringZ = (cutcurrOfiring-meanFRall)/stdFRall;
cutcurrNOfiringZ = (cutcurrNOfiring-meanFRall)/stdFRall;




if PVALSpeedF(cellnum) < 0.05%>1%< 0.05 % IF SPEED CORRELARTED DO THE FOLLOWING:

[b,dev,stats] = glmfit(allvecSpeed,allvecFiringZ,'normal','constant','on');

intercept = b(1);
coeff1 = b(2);
        
coefficienGLM = [intercept coeff1];

%%% calculate residuals of No-odour and odour firing

NOFRspeedCorrect = coefficienGLM(1) + cutNOSpeed*coefficienGLM(2);
        
%calculate residuals, meaning difference from model
NOresidSpeedCorrectZ = cutcurrNOfiringZ-NOFRspeedCorrect; 



OFRspeedCorrect = coefficienGLM(1) + cutOSpeed*coefficienGLM(2);
        
%calculate residuals, meaning difference from model
OresidSpeedCorrectZ = cutcurrOfiringZ-OFRspeedCorrect; 


TF = 1;

%%% take the difference

avrgNOresidSpeedCorrectZ = nanmean(NOresidSpeedCorrectZ(:,1:17),1);
avrgOresidSpeedCorrectZ = nanmean(OresidSpeedCorrectZ,1);
nanDiffVec(1:5) = NaN; 

overall_avrgFRVR = nanmean([nanmean(NOresidSpeedCorrectZ,1); [avrgOresidSpeedCorrectZ nanDiffVec]],1);
overall_avrgFFREW = nanmean([NOOverallREWFR_Z; OverallREWFR_Z],1);


currdiff = avrgOresidSpeedCorrectZ - avrgNOresidSpeedCorrectZ;
diffZ = [diffZ; currdiff];

curroverallFR = [overall_avrgFRVR overall_avrgFFREW];
FRZ = [FRZ; curroverallFR];


NOresidSpeedCorrectZandREW = [NOresidSpeedCorrectZ NOOverallREWFR_Z];
OresidSpeedCorrectZandREW = [OresidSpeedCorrectZ OverallREWFR_Z];

avrgNOresidSpeedCorrectZandREW = nanmean(NOresidSpeedCorrectZandREW);
avrgOresidSpeedCorrectZandREW = [OresidSpeedCorrectZ OverallREWFR_Z];

%(cellnum,:)

else %IF NOT SPEED CORRELATED DO THIS: 

    avrgNOFRZ = nanmean(cutcurrNOfiringZ(:,1:17),1);
    avrgOFRZ = nanmean(cutcurrOfiringZ,1);
    nanDiffVec(1:5) = NaN; 

    overall_avrgFRVR = nanmean([nanmean(cutcurrNOfiringZ,1); [avrgOFRZ nanDiffVec]],1);
    overall_avrgFFREW = nanmean([NOOverallREWFR_Z; OverallREWFR_Z],1);

    curroverallFR = [overall_avrgFRVR overall_avrgFFREW];
    FRZ = [FRZ; curroverallFR];

    currdiff = avrgOFRZ - avrgNOFRZ;
    diffZ = [diffZ; currdiff];
    %diffZ(cellnum,:) = avrgOFRZ - avrgNOFRZ;



end %%% end of if significant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



TF = 1;


%%% correlation of average values


vecNOFiring = avrgNOfiringrate(4:end);
vecNOSpeed = avrgNOspeed(4:end);

vecOFiring = avrgOfiringrate(4:end);
vecOSpeed = avrgOspeed(4:end);


avrgallvecFiring = cat(2,vecNOFiring,vecOFiring)
avrgallvecSpeed = cat(2,vecNOSpeed,vecOSpeed)

[avrgRHOspeed(cellnum),avrgPVALSpeedF(cellnum)] = corr(avrgallvecFiring',avrgallvecSpeed')





%figure
%plot(RAWvecNOFiring,RAWvecNOSpeed,'.')

%[NORHO(cellnum),NOPVALSpeedF(cellnum)] = corr(RAWvecNOFiring,RAWvecNOSpeed)



%preRAWvecOFiring = cutcurrOfiring';
%RAWvecOFiring = preRAWvecOFiring(:);

%preRAWvecOSpeed = cutOSpeed';
%RAWvecOSpeed = preRAWvecOSpeed(:);

%[ORHO(cellnum),OPVALSpeedF(cellnum)] = corr(RAWvecOFiring,RAWvecOSpeed)


%%%%%%%% correlation of noodour odour firing rates
%for corrs = 1:17

    %%% randomly sample from nonodour with same length as odour trials
    
    %for Nocorr = 1:50
        %NOodourrand = randperm(length(NOodourInd));
        %selNOodourrand = NOodourrand(1:length(odourInd));
        %NOcorrsectfiring = cutcurrNOfiring(selNOodourrand,corrs);
        %Ocorrsectfiring = cutcurrOfiring(1:length(odourInd),corrs);
        %currFiringCorr(Nocorr) = corr(NOcorrsectfiring,Ocorrsectfiring,'Type','Spearman');
    %end
    
    %FiringCorr(cellnum,corrs) = mean(currFiringCorr);
    
%end
%%%%%%%%%%%%%%%%%%%%%%%%

    for manU = 1:17
    
        
        NOdoursecfiring = cutcurrNOfiring(:,manU);
        Odoursecfiring = cutcurrOfiring(:,manU);
        
        if mean(NOdoursecfiring)> 0 && mean(Odoursecfiring)> 0
            
            manUpvalues(manU) = ranksum(NOdoursecfiring,Odoursecfiring);
        
        else
            
            manUpvalues(manU) = 0;    
            
        end
    
    
    end
    
    manUconfvalues(cellnum,:) = 1 - manUpvalues;
    
    CellAngles(cellnum) = BINNEDTRIALMEANANGLE;
    
    
    
    %%%%%%% calculate autocorrelogramm
    
    %%%%%thetaend
    %selspkind = find(spiketimes < thetaend);
    %selectspiketimes = spiketimes(selspkind);
    
    %binsize = 2;
    %corrwidth = 60;
    %bins = corrwidth/binsize;
    
    
    %for autospk = 1:length(selectspiketimes)
        
        %currtime = selectspiketimes(autospk);
        
        %for autobins = 1:bins
            
            %autoinds = find(selectspiketimes < (currtime + (binsize/1000)*autobins) & selectspiketimes > currtime + (binsize/1000)*(autobins-1)); 
            %abincorr(autospk,autobins) = length(autoinds);
            
        %end
            
    %end
    
    %autcorrhistright = sum(abincorr,1)
    %autcorrhistleft = fliplr(autcorrhistright)
    %combautocorr = cat(2,autcorrhistleft,autcorrhistright)
    
    
    
    %%%%%% COMBINE FIRING FOR SIMPLER PLOT
    
%     avrgcutcurrOfiring = mean(cutcurrOfiring);
%     avrgcutcurrNOfiring = mean(cutcurrNOfiring); 
%     
%     puttogethfiring = cat(1,avrgcutcurrNOfiring(1:17),avrgcutcurrOfiring);
%     avrgcombinedfiring = mean(puttogethfiring);
%     
%     combinedfiring(cellnum,:) = cat(2,avrgcombinedfiring,avrgcutcurrNOfiring(18:22));

    %DO NOT MEAN BUT ADD TOGETHER SPIKE NUMBERS

    avrgcutcurrOfiring = nansum(cutcurrOfiring,1); 
    avrgcutcurrNOfiring = nansum(cutcurrNOfiring,1); 
    
    puttogethfiring = cat(1,avrgcutcurrNOfiring(1:17),avrgcutcurrOfiring);
    avrgcombinedfiring = nansum(puttogethfiring,1);
    
    combinedfiring(cellnum,:) = cat(2,avrgcombinedfiring,avrgcutcurrNOfiring(18:22));

    currencombinedfiring = cat(2,avrgcombinedfiring,avrgcutcurrNOfiring(18:22));
    totalspikesIn(cellnum) = nansum(NOOverallREWFR,'all') + nansum(OverallREWFR,'all') + nansum(currencombinedfiring,'all');

    %sum of spikes taken out: 
    totalspikesOut(cellnum) = nansum(NOfiringrateSLOW(:,2:24),'all') + nansum(OfiringrateSLOW(:,2:19),'all');

    totalTimeNOodourSlow(cellnum) = nansum(NOtimesSLOW(:,2:22),'all')/size(NOtimesSLOW,1);
    totalTimeOodourSlow(cellnum) = nansum(OtimesSLOW(:,2:17),'all')/size(OtimesSLOW,1);

    meanTimeNOodourSlow(cellnum,:) = nanmean(NOtimesSLOW(:,4:22),1);
    meanTimeOodourSlow(cellnum,:) = nanmean(OtimesSLOW(:,4:17),1);

    %average spike per trial types
    OspkNumbstrials = nansum(cutcurrOfiring,2);
    NOspkNumbstrials = nansum(cutcurrNOfiring,2);

    NOspkNumbstrialsREW = nansum(NOOverallREWFR,2);
    OspkNumbstrialsREW = nansum(OverallREWFR,2);

    avrgOspktrials(cellnum) = mean([OspkNumbstrials + OspkNumbstrialsREW]);
    avrgNOspktrials(cellnum) = mean([NOspkNumbstrials + NOspkNumbstrialsREW]);

    Oskptrialnumber(cellnum) = length(OspkNumbstrials);
    NOskptrialnumber(cellnum) = length(NOspkNumbstrials);


%     NOfiringrateSLOW;
%     OfiringrateSLOW;

    

    
    %%%%%%%%%%%% TEST WITH SHUFFLING THE SIGNIFICANCE
    
    NoFiring = cutcurrNOfiring(:,1:17)
    oFiring = cutcurrOfiring
    
    firingrates = cat(1,NoFiring(:,1:size(oFiring,2)),oFiring)


    %load firingrates.txt
    %%%%  number of non odour and odour trials
    NOn = size(NoFiring,1)

    On = size(oFiring,1)

    alln = NOn + On

    OrigFiringNO = mean(firingrates((1:NOn),:));
    OrigFiringO = mean(firingrates(((NOn+1):alln),:));
    OrigDiff = OrigFiringNO - OrigFiringO;

    NONodourFiring = firingrates((1:NOn),:);
    OdourFiring = firingrates(((NOn+1):alln),:);

    %figure
    %hold on

    shff = 1;

    for shff = 1:10000

        %%% new version

        randomizedIND = randperm(alln);
    
        randomizedfiring = firingrates(randomizedIND,:);    

        randomizedNOfiring = randomizedfiring(1:NOn,:);
        randomizedOfiring = randomizedfiring((NOn+1):end,:);

        randmeanFiringNO = mean(randomizedNOfiring,1);
        randmeanFiringO = mean(randomizedOfiring,1);
        randDiff(shff,:) = randmeanFiringNO - randmeanFiringO;

    end

    %%%% correct for 4 sectors equals to 0.6333%
    sortedDiff = sort(randDiff);
    correc = 10000 * 0.005
    10000 - correc;

    upper = sortedDiff(10000-correc,:)
    lower = sortedDiff(correc,:)

    
        
    
    
    %figure
    %hold on
    %plot(upper)
    %plot(lower)
    %plot(OrigDiff,'.g')
    
    %%%%% test if there is a significant point with difference
    
    %if the firing rate is very low do not consider significant, that is
    %why maxfiring checked - difference should be greater than 10time diff
    maxavrgFiring = max(avrgcombinedfiring)
    sigIndices = find(OrigDiff > upper | OrigDiff < lower & abs(OrigDiff) > maxavrgFiring/10)
    if isempty(sigIndices) == 1
        
        pyrsignificances{cellnum} = 0
    else
        pyrsignificances{cellnum} = sigIndices
        
    end
    
    
    
    
    %%%% do different shuffling to test spatial modulation
    
    alln2 = 17
    
    OrigMean = mean(firingrates,1);

    shff2 = 1;

    for shff2 = 1:10000

    %%% new version

    for tn = 1:(NOn+On)
    
    randomizedIND2 = randperm(alln2);
    randomizedfiring2(tn,:) = firingrates(tn,randomizedIND2);    

    end

    randmeanFiring2(shff2,:) = mean(randomizedfiring2,1);

    end

    %%%% correct for 4 sectors equals to 0.6333%

    %%% sort firing in each bin
    for s1 = 1:17
    sortedFiring2(:,s1) = sort(randmeanFiring2(:,s1));

    end

    correc2 = 10000 * 0.01
    10000 - correc2;

    upper = sortedFiring2(10000-correc2,:)
    %lower = sortedFiring2(correc2,:)
    
    %figure
    %hold on
    %plot(upper)
    %plot(lower)
    %plot(OrigMean,'.g')
    
    
    sigIndices2 = find(OrigMean > upper) %| OrigMean < lower)
    if isempty(sigIndices2) == 1
        
        pyrsignificances2{cellnum} = 0
    else
        pyrsignificances2{cellnum} = sigIndices2
        
    end
    
    
    
    
    
    
    
    %%%%% Calculate spatial information for each cell for odour and no
    %%%%% odour trials
    
    %%%% first no odour trials
    
    
    
    
    avrgcutcurrOfiring = mean(cutcurrOfiring);
    avrgcutcurrNOfiring = mean(cutcurrNOfiring); 
    
    summedNOSPINFO = 0;
    summedOSPINFO = 0;
    
    
    
    for NOB1 = 1:22
               
        NORATIO = avrgcutcurrNOfiring(NOB1)/mean(avrgcutcurrNOfiring);
        
        if NORATIO == 0 
            summedNOSPINFO = summedNOSPINFO;
            allCellssummedNOSPINFO(cellnum) = summedNOSPINFO;
        else             
            summedNOSPINFO = summedNOSPINFO + 1/25*NORATIO*log2(NORATIO); 
            allCellssummedNOSPINFO(cellnum) = summedNOSPINFO;
            
        end
    end
    
    
    for OB1 = 1:17
        
        ORATIO = avrgcutcurrOfiring(OB1)/mean(avrgcutcurrOfiring);
               
        if ORATIO == 0  
            summedOSPINFO = summedOSPINFO;
            allCellssummedOSPINFO(cellnum) = summedOSPINFO;
            
        else
            ORATIO = avrgcutcurrOfiring(OB1)/mean(avrgcutcurrOfiring);
            summedOSPINFO = summedOSPINFO + 1/17*ORATIO*log2(ORATIO); 
            allCellssummedOSPINFO (cellnum) = summedOSPINFO;
            
        end
    end
    
    NumberofNoOdourtrials(cellnum) = length(NOodourInd);
    NumberofOdourtrials(cellnum) = length(odourInd);
    
       
    
end %end of loop for all cells

%
%
%
%
%
%
%
%
%
%
%

%
%
%


%%%%%%%%%%% ADDITIONAL ANALYSIS THAT CAN BE DONE

TF = 1; 

%  rewNOpval = rewNOpval'
%  rewOpval= rewOpval'

save('diffZ.mat','diffZ')

 diffZ = 0;
 FRZ = 0;

 sigmaxVals = max(diffZ(:,1:13),[],2);
 sigminVals = min(diffZ(:,1:13),[],2); sigminValsABS = abs(sigminVals);

 checkdiff = sigminValsABS > sigmaxVals;
 minGreaterID = find(checkdiff == 1);

 %save IDs
% 
%  save('minGreaterID.mat','minGreaterID')
% 
%  load('minGreaterID.mat')


%matchPop1 = 0;
%matchPop2 = 0;
vec1 = 0;
vec2 = 0;

vec1L = length(vec1);

vec2L = length(vec2);



AllmatchIDs = []
if vec1L < vec2L

    for i = 1:vec1L
        clear indmatch
        currvec1 = vec1(i);

        indmatch = find(vec2 == currvec1);
        %matchIDs = vec2(indmatch);

        AllmatchIDs = [AllmatchIDs indmatch];
    end

else

    for i = 1:vec2L
        clear indmatch
        currvec2 = vec2(i);

        indmatch = find(vec1 == currvec2);
        %matchIDs = vec1(indmatch);

        AllmatchIDs = [AllmatchIDs indmatch];
    end

end

d = 1:151;

d(d ~= AllmatchIDs')

nonSIG = find(d ~= AllmatchIDs');

sigIDs = 0;

d(sigIDs) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffZ = 0;


minGreaterID = [45	5	47	55	20	38	19	1	18	25	28	42	23	30	64	2	24	29	31 ...
                32	39	44	52	56	27	35	26	33	34	37	40	4	36	41	43	50	58	63];


maxGreaterID = [48	46	62	15	16	53	57	60	6	7	8	10	13	49	9	11	12	14	61	17	22	51	59	54	21	3];




minGreaterID = [84	6	88	109	34	61	33	1	32	44	47	79	42	49	151	2	43	48	51	52	65	82	104	110	46	57	45	55	56	59	69	5	58	74	81	98	119	150];

maxGreaterID = [91	87	142	27	28	106	111	121	8	11	15	20	25	96	19	22	23	26	131	29	41	103	120	107	39	3];

nondiffID = [4	7	9	10	12	14	16	17	18	21	24	30	35	36	37	38	40	50	53	54	60	62	63	64	66	67	68	70	71	72	73	75	76	77	78	80	86	89	90 ...
            92	93	94	95	97	99	100	101	102	105	108	112	113	114	115	116	117	118	122	123	127	133	134	135	136	137	138	139	140	141	143	144	145	146	147	148	149];

rewIDs = [13 31	83	85	124	125	126	128	129	130	132];




 Noodourdifferences = diffZ(rewIDs,:)

 %Noodourdifferences = 0;
 
 [sigminVals2 minPos] = min(Noodourdifferences(:,1:13),[],2);
 [sigminVals2 minPos] = max(Noodourdifferences(:,1:13),[],2);
 
 [sortOUT sortedNOODSIG_IDs] = sort(minPos);

 sortedNoodourdifferences = Noodourdifferences(sortedNOODSIG_IDs,:);

 NoodourFRZs = FRZ(rewIDs,1:22);
 sortedNoodourFRZs = NoodourFRZs(sortedNOODSIG_IDs,:);

clear currFiring smoothedFiring y

x = (1:22)'
for pyrs = 1:11
    
    currFiring = NoodourFRZs(pyrs,:);
    y = fit(x,currFiring','smoothingspline');
    smoothedFiring(pyrs,:) = y(1:0.25:22);
    
    
end

smoothedFiring = 0;


    figure

    %load diffColorMap
diffColorMap = 0;

    figure
    map = diffColorMap
    
    imagesc(smoothedFiring)
    colormap(map)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([-1.6 1.6])



    figure
    map = jet
    
    imagesc(smoothedFiring)
    colormap(map)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([-1 1])









 imagesc(sortedNoodourdifferences)





 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 %
 

%clear all

CellAngles = CellAngles'

%%%%%%% check angles and convert - Angles to our system

for as1 = 1:length(CellAngles)
    
    if CellAngles(as1) <= 0        
        corrCellAngles(as1) = 360+CellAngles(as1);
        
    else        
        corrCellAngles(as1) = CellAngles(as1);        
    end
        
end


orginalpvalues = 1 - manUconfvalues;

%%%% the firingdiffplot has all firing difference values, from this the
%%%% diffscore is calculated

diffscore = firingdiffplot.*manUconfvalues;



%%%%%%%% 
end





