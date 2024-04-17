function plotsWaffle()

clear;
close all;

dataDir = sprintf('%s/data', pwd);
subjectID="waffles_pilot001";
sessionNum=1;
datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
date="240412_1707";
datafile = sprintf('%s/%s_s%d_predv2_s%s.mat', datafileDir, subjectID, sessionNum,date);
file = load(datafile);
%% check that staircase is staircasing
trialCount=length(file.d.stairIdx);
x=[1:1:trialCount];
figure();
plot(x,file.d.stairIdx);
ylabel("stair ID");
xlabel("trial ID");

%% plot contrast judgement task from waffles

waffleContrastsIdx=file.d.plaidContrast;
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];
waffleContrasts=[];
for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
end

precueValidities=file.d.precueValidity;
precueValidities(skipNotWaffle)=[];

waffleOrientation=file.d.plaidOrientation;
waffleOrientation(skipNotWaffle)=[];

responses= file.d.targetResponseKey;
responses(skipNotWaffle)=[];

correctResponse=[];
for i=1:length(responses)
    if responses(i)==1 && waffleOrientation(i)==1
        correctResponse=[correctResponse 1];
    elseif responses(i)==2 && waffleOrientation(i)==2
        correctResponse=[correctResponse 1];
    else
        correctResponse=[correctResponse 0];
    end
end
contrastJudgement= cat(4,waffleContrasts,precueValidities,waffleOrientation,correctResponse);
nTrialsWaffle=length(waffleContrastsIdx);
contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,4]);
%sort based on contrast levels
contrastJudgement=sortrows(contrastJudgement,1); 
% sort based on precue validity
contrastJudgement=sortrows(contrastJudgement,2);

correctsEx=[];
correctsUnex=[];
meanCorrectPerConEx=[];
meanCorrectPerConUnex=[];
for i=1:length(responses)
    if contrastJudgement(i,2)==1%% expected
        if i==1
            correctsEx=[correctsEx contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            correctsEx=[correctsEx contrastJudgement(i,4)];
            if i==length(responses)
                correctsEx=[correctsEx contrastJudgement(i,4)];
                meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            end

        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsEx=[];
            correctsEx=[correctsEx contrastJudgement(i,4)];
        end
    elseif contrastJudgement(i,2)==2 
        if contrastJudgement(i-1,2)==1
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsUnex=[correctsUnex contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            if i==length(responses)
                correctsUnex=[correctsUnex contrastJudgement(i,4)];
                meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            end
            correctsUnex=[correctsUnex contrastJudgement(i,4)];


        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            correctsUnex=[];
            correctsUnex=[contrastJudgement(i,4)];
        end
    end
end
contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);

allconsEx=unique(contrastJudgementEx);
allconsUnex=unique(contrastJudgementUnex);
figure();
plot(allconsEx,meanCorrectPerConEx,'blue-v');
hold on
plot(allconsUnex,meanCorrectPerConUnex,'red-v');
xlabel("contrasts")
ylabel("p(correct stronger)");
title("p(correct stronger) for expected vs unexpected ",subjectID);
legend("expected","unexpected");

%% p (correct tilt discrimination)

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];

precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];

response= file.d.targetResponseKey;
response(skipWaffle)=[];

diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];

static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);
correct_exp=correct_r==1 & precueValidities==1;
correct_unexp= correct_r==1 & precueValidities==2;

mean(correct_exp);
mean(correct_unexp);
orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
nTrialsGrating=length(stair);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
%sort based on stair levels
orientationJudgement=sortrows(orientationJudgement,1); 


corrects=[];
meanCorrectPerCon=[];
for i=1:length(response)
    if i==1
        corrects=[corrects orientationJudgement(i,3)];

    elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
        corrects=[corrects orientationJudgement(i,3)];
        if i==length(response)
            corrects=[corrects orientationJudgement(i,3)];
            meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
        end

    elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
        meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
        corrects=[];
        corrects=[corrects orientationJudgement(i,3)];
        if i==length(response)
           meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
        end
   
    end
end

allStairs=unique(orientationJudgement(:,1));

allstairsId=[];
for i=1:length(allStairs)
    allstairsId=[allstairsId find(file.p.stairs==allStairs(i))];
end


figure();
plot(allStairs,meanCorrectPerCon,'blue-v');
xlabel("stair value")
ylabel("p (correct tilt discrimination)");
title("p (correct tilt discrimination) ",subjectID);

figure();
plot(allstairsId,meanCorrectPerCon,'blue-v');
xlabel("stair id")
ylabel("p (correct tilt discrimination)");
title("p (correct tilt discrimination) for expected vs unexpected ",subjectID);

mean(correct_r)

%% p (correct tilt discrimination)

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];

precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];

response= file.d.targetResponseKey;
response(skipWaffle)=[];

diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];

static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);

orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
nTrialsGrating=length(stair);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
%sort based on stair levels
orientationJudgement=sortrows(orientationJudgement,1); 
% sort based on precue validity
orientationJudgement=sortrows(orientationJudgement,2);

correctsEx=[];
correctsUnex=[];
meanCorrectPerConEx=[];
meanCorrectPerConUnex=[];
for i=1:length(response)
    if orientationJudgement(i,2)==1%% expected
        if i==1
            correctsEx=[correctsEx orientationJudgement(i,3)];

        elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
            correctsEx=[correctsEx orientationJudgement(i,3)];
            if i==length(response)
                correctsEx=[correctsEx orientationJudgement(i,3)];
                meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            end

        elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsEx=[];
            correctsEx=[correctsEx orientationJudgement(i,3)];
            if i==length(response)
               meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            end
        end
    elseif orientationJudgement(i,2)==2 
        if orientationJudgement(i-1,2)==1
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsUnex=[correctsUnex orientationJudgement(i,3)];

        elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
            if i==length(response)
                correctsUnex=[correctsUnex orientationJudgement(i,3)];
                meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            end
            correctsUnex=[correctsUnex orientationJudgement(i,3)];


        elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
            meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            correctsUnex=[];
            correctsUnex=[orientationJudgement(i,3)];
            if i==length(response)
               meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            end
        end
    end
end

orientationJudgementEx=orientationJudgement(orientationJudgement(:,2)==1);
orientationJudgementUnex=orientationJudgement(orientationJudgement(:,2)==2);

allstairsEx=unique(orientationJudgementEx);
allstairsUnex=unique(orientationJudgementUnex);
allstairsIdEx=[];
allstairsIdUnex=[];
for i=1:length(allstairsEx)
    allstairsIdEx=[allstairsIdEx find(file.p.stairs==allstairsEx(i))];
end

for i=1:length(allstairsUnex)
    allstairsIdUnex=[allstairsIdUnex find(file.p.stairs==allstairsUnex(i))];
end


figure();
plot(allstairsEx,meanCorrectPerConEx,'blue-v');
hold on
plot(allstairsUnex,meanCorrectPerConUnex,'red-v');
xlabel("stair value")
ylabel("p (correct tilt discrimination)");
title("p (correct tilt discrimination) for expected vs unexpected ",subjectID);
legend("expected","unexpected");

figure();
plot(allstairsIdEx,meanCorrectPerConEx,'blue-v');
hold on
plot(allstairsIdUnex,meanCorrectPerConUnex,'red-v');
xlabel("stair id")
ylabel("p (correct tilt discrimination)");
title("p (correct tilt discrimination) for expected vs unexpected ",subjectID);
legend("expected","unexpected");


%% p (correct contrast discrimination with respect to difference)

waffleContrastsIdx=file.d.plaidContrast;
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];

waffleContrasts=[];
difference=[];
waffleOrientation=file.d.plaidOrientation; 
waffleOrientation(skipNotWaffle)=[];

for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
    if waffleOrientation(i)==2 %if orientation is +45 then difference is positive
        difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
    elseif waffleOrientation(i)==1 %if orientation is -45 then difference is positive, 
        difference(i)=file.p.plaidContrasts2-waffleContrasts(i); %so difference based on
    end
end

precueValidities=file.d.precueValidity;
precueValidities(skipNotWaffle)=[];

response= file.d.targetResponseKey;
response(skipNotWaffle)=[];

correct_r=(response==1 & waffleOrientation==1) | (response==2 & waffleOrientation==2)

contrastJudgement= cat(4,difference,precueValidities,correct_r,response);
nTrialsWaffle=length(waffleContrastsIdx);
contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,4]);
%sort based on difference values
contrastJudgement=sortrows(contrastJudgement,1); 
% sort based on precue validity
contrastJudgement=sortrows(contrastJudgement,2);

correctsEx=[];
correctsUnex=[];
meanCorrectPerConEx=[];
meanCorrectPerConUnex=[];
for i=1:length(response)
    if contrastJudgement(i,2)==1%% expected
        if i==1
            correctsEx=[correctsEx contrastJudgement(i,3)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            correctsEx=[correctsEx contrastJudgement(i,3)];
            if i==length(response)
                correctsEx=[correctsEx contrastJudgement(i,3)];
                meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            end

        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsEx=[];
            correctsEx=[correctsEx contrastJudgement(i,3)];
            if i==length(response)
               meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            end
        end
    elseif contrastJudgement(i,2)==2 
        if contrastJudgement(i-1,2)==1
            meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
            correctsUnex=[correctsUnex contrastJudgement(i,3)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            if i==length(response)
                correctsUnex=[correctsUnex contrastJudgement(i,3)];
                meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            end
            correctsUnex=[correctsUnex contrastJudgement(i,3)];


        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            correctsUnex=[];
            correctsUnex=[contrastJudgement(i,3)];
            if i==length(response)
               meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
            end
        end
    end
end

contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);

alldiffEx=unique(contrastJudgementEx);
alldiffUnex=unique(contrastJudgementUnex);


figure();
plot(alldiffEx,meanCorrectPerConEx,'blue-v');
hold on
plot(alldiffUnex,meanCorrectPerConUnex,'red-v');
xlabel("difference value (with respect to +45)")
ylabel("p (correct contrast discrimination with respect to difference)");
title("p (correct contrast discrimination with respect to difference) for expected vs unexpected ",subjectID);
legend("expected","unexpected");

%% p (responses with +45/-45 with respect to difference)

waffleContrastsIdx=file.d.plaidContrast;
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];

waffleContrasts=[];
difference=[];
waffleOrientation=file.d.plaidOrientation; 
waffleOrientation(skipNotWaffle)=[];

for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
    if waffleOrientation(i)==2 %if orientation is +45 then difference is positive
        difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
    elseif waffleOrientation(i)==1 %if orientation is -45 then difference is positive, 
        difference(i)=file.p.plaidContrasts2-waffleContrasts(i); %so difference based on
    end
end

precueValidities=file.d.precueValidity;
precueValidities(skipNotWaffle)=[];

response= file.d.targetResponseKey;
response(skipNotWaffle)=[];

response_45=(response==2);
response_135=(response==1);

correct_r=(response==1 & waffleOrientation==1) | (response==2 & waffleOrientation==2)

contrastJudgement= cat(7,difference,precueValidities,correct_r,response_45,response_135,response,waffleContrasts);
nTrialsWaffle=length(waffleContrastsIdx);
contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,7]);
%sort based on difference values
contrastJudgement=sortrows(contrastJudgement,1); 
% sort based on precue validity
contrastJudgement=sortrows(contrastJudgement,2);

responseEx=[];
responseUnex=[];
meanResponsePerConEx=[];
meanResponsePerConUnex=[];
for i=1:length(response)
    if contrastJudgement(i,2)==1%% expected
        if i==1
            responseEx=[responseEx contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            responseEx=[responseEx contrastJudgement(i,4)];
            if i==length(response)
                responseEx=[responseEx contrastJudgement(i,4)];
                meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
            end

        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
            responseEx=[];
            responseEx=[responseEx contrastJudgement(i,4)];
            if i==length(response)
               meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
            end
        end
    elseif contrastJudgement(i,2)==2 
        if contrastJudgement(i-1,2)==1
            meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
            responsesUnex=[responseUnex contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            if i==length(response)
                responseUnex=[responseUnex contrastJudgement(i,4)];
                meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
            end
            responseUnex=[responseUnex contrastJudgement(i,4)];


        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
            responseUnex=[];
            responseUnex=[contrastJudgement(i,4)];
            if i==length(response)
               meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
            end
        end
    end
end

responseEx135=[];
responsesUnex135=[];
meanResponsePerConEx_135=[];
meanResponsePerConUnex_135=[];

%get hits: responsePresent (they responsed with +45) and orientation
%present (orientation was in fact +45)
%get false alarms: they responsed with +45 and orientation was NOT +45 (was
%-45) 
% get d'
    % hitsE=responsePresent==1 & orientationPresent==1;
    % hit_rate_ex=mean(hitsE);
    % hit_rates_ex(i)=hit_rate_ex;
  % cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
  %   cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
    % orientationAbsent=(responseEx(prop_exp*i-prop_exp+1:prop_exp*i,3)==2);
    % fa_ex=responsePresent==1 & orientationAbsent==1;
    % fa_rate_ex=mean(fa_ex);
    % fa_rates_ex(i)=fa_rate_ex;


for i=1:length(response)
    if contrastJudgement(i,2)==1%% expected
        if i==1
            responseEx135=[responseEx135 contrastJudgement(i,5)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            responseEx135=[responseEx135 contrastJudgement(i,5)];
            if i==length(response)
                responseEx135=[responseEx135 contrastJudgement(i,5)];
                meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
            end

        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
            responseEx135=[];
            responseEx135=[responseEx135 contrastJudgement(i,5)];
            if i==length(response)
               meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
            end
        end
    elseif contrastJudgement(i,2)==2 
        if contrastJudgement(i-1,2)==1
            meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
            responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            if i==length(response)
                responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];
                meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
            end
            responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];


        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
            responsesUnex135=[];
            responsesUnex135=[contrastJudgement(i,5)];
            if i==length(response)
               meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
            end
        end
    end
end

contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);

alldiffEx=unique(contrastJudgementEx);
alldiffUnex=unique(contrastJudgementUnex);

% x-axis is difference between grating contrasts inside plaid with respect
% to +45
figure();
yyaxis left
plot(alldiffEx,meanResponsePerConEx,'-v');
hold on
plot(alldiffUnex,meanResponsePerConUnex,'--v');
ylabel("p (reported +45 with respect to difference)");

yyaxis right
plot(alldiffEx,meanResponsePerConEx_135,'-v');
hold on
plot(alldiffEx,meanResponsePerConUnex_135,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported -45 with respect to difference)");
title("p  (reported -45/+45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected 45","unexpected 45","expected -45","unexpected -45");

figure();
plot(alldiffEx,meanResponsePerConEx,'-v');
hold on
plot(alldiffUnex,meanResponsePerConUnex,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p  (reported +45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected 45","unexpected 45");

figure();
plot(alldiffEx,meanResponsePerConEx_135,'-v');
hold on
plot(alldiffEx,meanResponsePerConUnex_135,'--v');
xlabel("difference value (with respect to -45)")
ylabel("p (reported -45 with respect to difference)");
title("p (reported -45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected -45","unexpected -45");

% x-axis is strength of plaid 

allStrength=[flip(file.p.plaidContrasts1) file.p.plaidContrasts1 ];

figure();
yyaxis left
plot(allStrength,meanResponsePerConEx,'-v');
hold on
plot(allStrength,meanResponsePerConUnex,'--v');
ylabel("p (reported +45 with respect to difference)");

yyaxis right
plot(allStrength,meanResponsePerConEx_135,'-v');
hold on
plot(allStrength,meanResponsePerConUnex_135,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported -45 with respect to difference)");
title("p (reported -45/+45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected 45","unexpected 45","expected -45","unexpected -45");

figure();
plot(allStrength,meanResponsePerConEx,'-v');
hold on
plot(allStrength,meanResponsePerConUnex,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p (reported +45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected 45","unexpected 45");

figure();
plot(allStrength,meanResponsePerConEx_135,'-v');
hold on
plot(allStrength,meanResponsePerConUnex_135,'--v');
xlabel("difference value (with respect to -45)")
ylabel("p (reported -45 with respect to difference)");
title("p (reported -45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected -45","unexpected -45");

%% d' and criterion for expected vs unexpected

% %signal present= -45 (1)
% %signal absent = +45 (2)
signalPresent=1;
signalAbsent=2;

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

% get all stair values and remove waffle trials
stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];
% get all validities and remove waffle trials
precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];
% get all responses and remove waffle trials
response= file.d.targetResponseKey;
response(skipWaffle(1:length(response)))=[];
% get all difference multipliers and remove waffle trials
diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];
% get all static grating values and remove waffle trials
static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

% a correct response occurs when the response is 9 and static is standard
% and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% response is 0 and static is test and multiplier is 1, response is 0 and
% static is standard and multiplier is 2
correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);

response_45=(response==2); %response was 0 ie +45
response_135=(response==1); %response was 9 ie -45 or 135

orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
nTrialsGrating=length(gratingContrastsIdx);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
%sort based on stair values
orientationJudgement=sortrows(orientationJudgement,1); 
% sort based on precue validity
orientationJudgement=sortrows(orientationJudgement,2);

 
% Expected
%response present is +45
responsePresent=response_45;
%orientation present is +45
orientationPresent=(static==1 & diffMult==1) | (static==2 & diffMult==-1);
%precue validity is valid ie 1
precueValid=(precueValidities==1);
% get hits for expected
hitsE=responsePresent==1 & orientationPresent==1 & precueValid==1;

%orientation absent: is -45
orientationAbsent=(static==2 & diffMult==1) | (static==1 & diffMult==-1);
% get false alarms for expected
fa_ex=responsePresent==1 & orientationAbsent==1 & precueValid==1;
  


% response absent is -45
responseAbsent=response_135;
% misses for expected
miss_ex1= responseAbsent==1 & orientationPresent==1 & precueValid==1;
%correct rejection
cr_ex1=responseAbsent==1 & orientationAbsent==1;


% Unexpected
%precue validity is valid ie 1
precueInvalid=(precueValidities==2);
% get hits for unexpected
hitsU=responsePresent==1 & orientationPresent==1 & precueInvalid==1;
% get false alarms for unexpected
fa_unex=responsePresent==1 & orientationAbsent==1  & precueInvalid==1;
% get misses  for unexpected
miss_unex1= responseAbsent==1 & orientationPresent==1 & precueInvalid==1;    
% get correct rejection for unexpected
cr_unex1=responseAbsent==1 & orientationAbsent==1  & precueInvalid==1;    


orientationJudgementEx=orientationJudgement(orientationJudgement(:,2)==1);
orientationJudgementUnex=orientationJudgement(orientationJudgement(:,2)==2);

allstairsEx=unique(orientationJudgementEx);
allstairsUnex=unique(orientationJudgementUnex);
allstairsIdEx=[];
allstairsIdUnex=[];

for i=1:length(allstairsEx)
    allstairsIdEx=[allstairsIdEx find(file.p.stairs==allstairsEx(i))];
end

for i=1:length(allstairsUnex)
    allstairsIdUnex=[allstairsIdUnex find(file.p.stairs==allstairsUnex(i))];
end

hit_rates_ex=[];
fa_rates_ex=[];
cons_dprime_ex=[];
cons_c_ex=[];
miss_ex=[];
cr_ex=[];


for i=1:length(allstairsEx)

    % Expected
    hitsEstair=(hitsE==1) & (stair==allstairsEx(i));
    hit_rate_ex=mean(hitsEstair);
    hit_rates_ex(i)=hit_rate_ex;

    faEstair=(fa_ex==1) & (stair==allstairsEx(i));
    fa_rate_ex=mean(faEstair);
    fa_rates_ex(i)=fa_rate_ex;

    %   get dprime and criterion    
    if fa_rate_ex==0
        fa_rate_ex=1/(length(orientationJudgementEx));
    elseif fa_rate_ex==1 
        fa_rate_ex=(length(orientationJudgementEx)-1)/(length(orientationJudgementEx));
    end
    if hit_rate_ex==0 
        hit_rate_ex=1/(length(orientationJudgementEx));
    elseif hit_rate_ex==1
        hit_rate_ex=(length(orientationJudgementEx)-1)/(length(orientationJudgementEx));
    end

    cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
    cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
    
    missEstair=miss_ex1==1 & (stair==allstairsEx(i));
    miss_ex(i)=mean(missEstair);
    crEstair=cr_ex1==1 & (stair==allstairsEx(i));
    cr_ex(i)=mean(crEstair);

end

hit_rates_unex=[];
fa_rates_unex=[];
cons_dprime_unex=[];
cons_c_unex=[];
miss_unex=[];
cr_unex=[];
% Unexpected    
for i=1:length(allstairsUnex)

    hitsU_stair= hitsU==1 & (stair==allstairsUnex(i));
    hit_rate_unex=mean(hitsU_stair);
    hit_rates_unex(i)=hit_rate_unex;

    faU_stair=fa_unex==1 & (stair==allstairsUnex(i));
    fa_rate_unex=mean(faU_stair);
    fa_rates_unex(i)=fa_rate_unex;

    if fa_rate_unex==0
        fa_rate_unex=1/(length(orientationJudgementUnex));
    elseif fa_rate_unex==1 
        fa_rate_unex=(length(orientationJudgementUnex)-1)/(length(orientationJudgementUnex));
    end
    if hit_rate_unex==0 
        hit_rate_unex=1/(length(orientationJudgementUnex));
    elseif hit_rate_unex==1
        hit_rate_unex=(length(orientationJudgementUnex)-1)/(length(orientationJudgementUnex));
    end

    %   get dprime and criterion
    cons_dprime_unex(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
    cons_c_unex(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
    missUstair=miss_unex1==1 & (stair==allstairsUnex(i));
    miss_unex(i)=mean(missUstair);
    crUstair= cr_unex1==1 & (stair==allstairsUnex(i));
    cr_unex(i)=mean(crUstair);
end
   
% stair ids
figure();
subplot(2,1,1);
plot(allstairsIdEx,cons_dprime_ex,'blue-v');
hold on
plot(allstairsIdUnex,cons_dprime_unex,'red-v');
xlabel("stairs")
ylabel("d'");
title("d' expected vs unexpected");
legend("expected","unexpected");

subplot(2,1,2);
plot(allstairsIdEx,cons_c_ex,'blue-v');
hold on
plot(allstairsIdUnex,cons_c_unex,'red-v');
xlabel("stairs")
ylabel("criterion");
title("criterion expected vs unexpected, ",subjectID);
legend("expected","unexpected");

% stair values

figure();
subplot(2,1,1);
plot(allstairsEx,cons_dprime_ex,'blue-v');
hold on
plot(allstairsUnex,cons_dprime_unex,'red-v');
xlabel("stairs")
ylabel("d'");
title("d' expected vs unexpected");
legend("expected","unexpected");

subplot(2,1,2);
plot(allstairsEx,cons_c_ex,'blue-v');
hold on
plot(allstairsUnex,cons_c_unex,'red-v');
xlabel("stairs")
ylabel("criterion");
title("criterion expected vs unexpected, ",subjectID);
legend("expected","unexpected");


%% New adjustments

% adjust for ceiling or floor values
% nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% nh(nh==0) = 1;
% nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% nfa(nfa==0) = 1;

% %signal present= -45 (1)
% %signal absent = +45 (2)
signalPresent=1;
signalAbsent=2;

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

% get all stair values and remove waffle trials
stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];
% get all validities and remove waffle trials
precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];
% get all responses and remove waffle trials
response= file.d.targetResponseKey;
response(skipWaffle(1:length(response)))=[];
% get all difference multipliers and remove waffle trials
diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];
% get all static grating values and remove waffle trials
static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

% a correct response occurs when the response is 9 and static is standard
% and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% response is 0 and static is test and multiplier is 1, response is 0 and
% static is standard and multiplier is 2
correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);

response_45=(response==2); %response was 0 ie +45
response_135=(response==1); %response was 9 ie -45 or 135

orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
nTrialsGrating=length(gratingContrastsIdx);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
%sort based on stair values
orientationJudgement=sortrows(orientationJudgement,1); 
% sort based on precue validity
orientationJudgement=sortrows(orientationJudgement,2);

 
% Expected
%response present is +45
responsePresent=response_45;
%orientation present is +45
orientationPresent=(static==1 & diffMult==1) | (static==2 & diffMult==-1);
%precue validity is valid ie 1
precueValid=(precueValidities==1);
% get hits for expected
hitsE=responsePresent==1 & orientationPresent==1 & precueValid==1;

%orientation absent: is -45
orientationAbsent=(static==2 & diffMult==1) | (static==1 & diffMult==-1);
% get false alarms for expected
fa_ex=responsePresent==1 & orientationAbsent==1 & precueValid==1;
 
% response absent is -45
responseAbsent=response_135;
% misses for expected
miss_ex1= responseAbsent==1 & orientationPresent==1 & precueValid==1;
%correct rejection
cr_ex1=responseAbsent==1 & orientationAbsent==1;

% Unexpected
%precue validity is valid ie 1
precueInvalid=(precueValidities==2);
% get hits for unexpected
hitsU=responsePresent==1 & orientationPresent==1 & precueInvalid==1;
% get false alarms for unexpected
fa_unex=responsePresent==1 & orientationAbsent==1  & precueInvalid==1;
% get misses  for unexpected
miss_unex1= responseAbsent==1 & orientationPresent==1 & precueInvalid==1;    
% get correct rejection for unexpected
cr_unex1=responseAbsent==1 & orientationAbsent==1  & precueInvalid==1;    

orientationJudgementEx=orientationJudgement(orientationJudgement(:,2)==1);
orientationJudgementUnex=orientationJudgement(orientationJudgement(:,2)==2);

allstairsEx=unique(orientationJudgementEx);
allstairsUnex=unique(orientationJudgementUnex);
allstairsIdEx=[];
allstairsIdUnex=[];

for i=1:length(allstairsEx)
    allstairsIdEx=[allstairsIdEx find(file.p.stairs==allstairsEx(i))];
end

for i=1:length(allstairsUnex)
    allstairsIdUnex=[allstairsIdUnex find(file.p.stairs==allstairsUnex(i))];
end

hit_rates_ex=[];
fa_rates_ex=[];
cons_dprime_ex=[];
cons_c_ex=[];
miss_ex=[];
cr_ex=[];


for i=1:length(allstairsEx)

    % Expected
    hitsEstair=(hitsE==1) & (stair==allstairsEx(i));
    hit_rate_ex=mean(hitsEstair);
    hit_rates_ex(i)=hit_rate_ex;

    faEstair=(fa_ex==1) & (stair==allstairsEx(i));
    fa_rate_ex=mean(faEstair);
    fa_rates_ex(i)=fa_rate_ex;

    %   get dprime and criterion  
    nh=sum(hitsEstair);
    nfa=sum(faEstair);
    signal=orientationPresent==1 & (stair==allstairsEx(i));
    nsignal=sum(signal);
    noise=orientationAbsent==1 & (stair==allstairsEx(i));
    nnoise=sum(noise);
    if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
           nh = nh + 0.5;
           nfa = nfa + 0.5;
           nnoise = nnoise + 1;
           nsignal = nsignal + 1;
           % proportions
           hit_rate_ex = nh./nsignal;
           fa_rate_ex = nfa./nnoise;
    end

    cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
    cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
    
    missEstair=miss_ex1==1 & (stair==allstairsEx(i));
    miss_ex(i)=mean(missEstair);
    crEstair=cr_ex1==1 & (stair==allstairsEx(i));
    cr_ex(i)=mean(crEstair);

end

hit_rates_unex=[];
fa_rates_unex=[];
cons_dprime_unex=[];
cons_c_unex=[];
miss_unex=[];
cr_unex=[];
% Unexpected    
for i=1:length(allstairsUnex)

    hitsU_stair= hitsU==1 & (stair==allstairsUnex(i));
    hit_rate_unex=mean(hitsU_stair);
    hit_rates_unex(i)=hit_rate_unex;

    faU_stair=fa_unex==1 & (stair==allstairsUnex(i));
    fa_rate_unex=mean(faU_stair);
    fa_rates_unex(i)=fa_rate_unex;

     %   get dprime and criterion  
    nh=sum(hitsU_stair);
    nfa=sum(faU_stair);
    signal=orientationPresent==1 & (stair==allstairsUnex(i));
    nsignal=sum(signal);
    noise=orientationAbsent==1 & (stair==allstairsUnex(i));
    nnoise=sum(noise);
    if fa_rate_unex==0 || fa_rate_unex==1  || hit_rate_unex==0  || hit_rate_unex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
           nh = nh + 0.5;
           nfa = nfa + 0.5;
           nnoise = nnoise + 1;
           nsignal = nsignal + 1;
           % proportions
           hit_rate_unex = nh./nsignal;
           fa_rate_unex = nfa./nnoise;
    end

    %   get dprime and criterion
    cons_dprime_unex(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
    cons_c_unex(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
    missUstair=miss_unex1==1 & (stair==allstairsUnex(i));
    miss_unex(i)=mean(missUstair);
    crUstair= cr_unex1==1 & (stair==allstairsUnex(i));
    cr_unex(i)=mean(crUstair);
end
   
% stair ids
figure();
subplot(2,1,1);
plot(allstairsIdEx,cons_dprime_ex,'blue-v');
hold on
plot(allstairsIdUnex,cons_dprime_unex,'red-v');
xlabel("stairs")
ylabel("d'");
title("d' expected vs unexpected");
legend("expected","unexpected");

subplot(2,1,2);
plot(allstairsIdEx,cons_c_ex,'blue-v');
hold on
plot(allstairsIdUnex,cons_c_unex,'red-v');
xlabel("stairs")
ylabel("criterion");
title("criterion expected vs unexpected, ",subjectID);
legend("expected","unexpected");

% stair values

figure();
subplot(2,1,1);
plot(allstairsEx,cons_dprime_ex,'blue-v');
hold on
plot(allstairsUnex,cons_dprime_unex,'red-v');
xlabel("stairs")
ylabel("d'");
title("d' expected vs unexpected");
legend("expected","unexpected");

subplot(2,1,2);
plot(allstairsEx,cons_c_ex,'blue-v');
hold on
plot(allstairsUnex,cons_c_unex,'red-v');
xlabel("stairs")
ylabel("criterion");
title("criterion expected vs unexpected, ",subjectID);
legend("expected","unexpected");


%% New adjustments + Signal is expected

% adjust for ceiling or floor values
% nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% nh(nh==0) = 1;
% nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% nfa(nfa==0) = 1;

% %signal present= expected
% %signal absent = unexpected
signalPresent=1;
signalAbsent=2;

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

% get all stair values and remove waffle trials
stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];
% get all validities and remove waffle trials
precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];
% get all responses and remove waffle trials
response= file.d.targetResponseKey;
response(skipWaffle(1:length(response)))=[];
% get all difference multipliers and remove waffle trials
diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];
% get all static grating values and remove waffle trials
static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

% a correct response occurs when the response is 9 and static is standard
% and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% response is 0 and static is test and multiplier is 1, response is 0 and
% static is standard and multiplier is 2
correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);

response_45=(response==2); %response was 0 ie +45
response_135=(response==1); %response was 9 ie -45 or 135

orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
nTrialsGrating=length(gratingContrastsIdx);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
%sort based on stair values
orientationJudgement=sortrows(orientationJudgement,1); 
% sort based on precue validity
orientationJudgement=sortrows(orientationJudgement,2);

 
% Expected
%response present is as expected (i.e. responsed 45 when it was 45 and -45
%when it was -45)
%observer responded with expected orientation 
responsePresent=(response_45==1 & diffMult==1 & static==1) | (response_45==1 & static==2 & diffMult==-1) | (response_135==1 & static==1 & diffMult==-1) | (response_135==1 & static==2 & diffMult==1);
% signal: expected orientation
validPresent=precueValidities==1;

% get hits 
hits=responsePresent==1 & validPresent==1;
%signal absent: was unexpected
validAbsent=precueValidities==2;
% get false alarms for expected
fa_ex=responsePresent==1 & validAbsent==1;
 
% response absent is responsed -45 when +45 and +45 when -45
responseAbsent=(response_135==1 & diffMult==1 & static==1) | (response_135==1 & static==2 & diffMult==-1) | (response_45==1 & static==1 & diffMult==-1) | (response_45==1 & static==2 & diffMult==1);
% misses  
miss_ex1= responseAbsent==1 & validPresent==1;
%correct rejection
cr_ex1=responseAbsent==1 & validAbsent==1;

allstairs=unique(stair);
allstairsId=[];

for i=1:length(allstairs)
    allstairsId=[allstairsId find(file.p.stairs==allstairs(i))];
end


hit_rates_ex=[];
fa_rates_ex=[];
cons_dprime_ex=[];
cons_c_ex=[];
miss_ex=[];
cr_ex=[];


for i=1:length(allstairs)

    % Expected
    hitsEstair=(hits==1) & (stair==allstairs(i));
    hit_rate_ex=mean(hitsEstair);
    hit_rates_ex(i)=hit_rate_ex;

    faEstair=(fa_ex==1) & (stair==allstairs(i));
    fa_rate_ex=mean(faEstair);
    fa_rates_ex(i)=fa_rate_ex;

    %   get dprime and criterion  
    nh=sum(hitsEstair);
    nfa=sum(faEstair);
    signal=validPresent==1 & (stair==allstairs(i));
    nsignal=sum(signal);
    noise=validAbsent==1 & (stair==allstairs(i));
    nnoise=sum(noise);
    if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
           nh = nh + 0.5;
           nfa = nfa + 0.5;
           nnoise = nnoise + 1;
           nsignal = nsignal + 1;
           % proportions
           hit_rate_ex = nh./nsignal;
           fa_rate_ex = nfa./nnoise;
    end

    cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
    cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
    
    missEstair=miss_ex1==1 & (stair==allstairs(i));
    miss_ex(i)=mean(missEstair);
    crEstair=cr_ex1==1 & (stair==allstairs(i));
    cr_ex(i)=mean(crEstair);

end

hit_rates_unex=[];
fa_rates_unex=[];
cons_dprime_unex=[];
cons_c_unex=[];
miss_unex=[];
cr_unex=[];
% Unexpected    
for i=1:length(allstairsUnex)

    hitsU_stair= hitsU==1 & (stair==allstairsUnex(i));
    hit_rate_unex=mean(hitsU_stair);
    hit_rates_unex(i)=hit_rate_unex;

    faU_stair=fa_unex==1 & (stair==allstairsUnex(i));
    fa_rate_unex=mean(faU_stair);
    fa_rates_unex(i)=fa_rate_unex;

     %   get dprime and criterion  
    nh=sum(hitsU_stair);
    nfa=sum(faU_stair);
    signal=orientationPresent==1 & (stair==allstairsUnex(i));
    nsignal=sum(signal);
    noise=orientationAbsent==1 & (stair==allstairsUnex(i));
    nnoise=sum(noise);
    if fa_rate_unex==0 || fa_rate_unex==1  || hit_rate_unex==0  || hit_rate_unex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
           nh = nh + 0.5;
           nfa = nfa + 0.5;
           nnoise = nnoise + 1;
           nsignal = nsignal + 1;
           % proportions
           hit_rate_unex = nh./nsignal;
           fa_rate_unex = nfa./nnoise;
    end

    %   get dprime and criterion
    cons_dprime_unex(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
    cons_c_unex(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
    missUstair=miss_unex1==1 & (stair==allstairsUnex(i));
    miss_unex(i)=mean(missUstair);
    crUstair= cr_unex1==1 & (stair==allstairsUnex(i));
    cr_unex(i)=mean(crUstair);
end
   
% stair ids
figure();
subplot(2,1,1);
plot(allstairsId,cons_dprime_ex,'blue-v');
xlabel("stair IDs")
ylabel("d'");
title("d'");

subplot(2,1,2);
plot(allstairsId,cons_c_ex,'blue-v');
xlabel("stair IDs")
ylabel("criterion");
title("criterion, ",subjectID);

% stair values

figure();
subplot(2,1,1);
plot(allstairs,cons_dprime_ex,'blue-v');
xlabel("stairs")
ylabel("d'");
title("d'");

subplot(2,1,2);
plot(allstairs,cons_c_ex,'blue-v');
xlabel("stairs")
ylabel("criterion");
title("criterion, ",subjectID);

%% d' criterion signal: expected, split by orientation (needs to be completed)

% %signal present= expected
% %signal absent = unexpected
signalPresent=1;
signalAbsent=2;

gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

% get all stair values and remove waffle trials
stair=file.d.stair;
stair(skipWaffle(1:length(stair)))=[];
% get all validities and remove waffle trials
precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];
% get all responses and remove waffle trials
response= file.d.targetResponseKey;
response(skipWaffle(1:length(response)))=[];
% get all difference multipliers and remove waffle trials
diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];
% get all static grating values and remove waffle trials
static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

% a correct response occurs when the response is 9 and static is standard
% and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% response is 0 and static is test and multiplier is 1, response is 0 and
% static is standard and multiplier is 2
correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1) | (response==2 & static==2 & diffMult==-1);

response_45=(response==2); %response was 0 ie +45
response_135=(response==1); %response was 9 ie -45 or 135

orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
nTrialsGrating=length(gratingContrastsIdx);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
%sort based on stair values
orientationJudgement=sortrows(orientationJudgement,1); 
% sort based on precue validity
orientationJudgement=sortrows(orientationJudgement,2);

 
% Expected
%response present is as expected (i.e. responsed 45 when it was 45 and -45
%when it was -45)
%observer responded with expected orientation 
responsePresent=(response_45==1 & diffMult==1 & static==1) | (response_45==1 & static==2 & diffMult==-1) | (response_135==1 & static==1 & diffMult==-1) | (response_135==1 & static==2 & diffMult==1);
% signal: expected orientation
validPresent=precueValidities==1;

% get hits 
hits=responsePresent==1 & validPresent==1;
%signal absent: was unexpected
validAbsent=precueValidities==2;
% get false alarms for expected
fa_ex=responsePresent==1 & validAbsent==1;
 
% response absent is responsed -45 when +45 and +45 when -45
responseAbsent=(response_135==1 & diffMult==1 & static==1) | (response_135==1 & static==2 & diffMult==-1) | (response_45==1 & static==1 & diffMult==-1) | (response_45==1 & static==2 & diffMult==1);
% misses  
miss_ex1= responseAbsent==1 & validPresent==1;
%correct rejection
cr_ex1=responseAbsent==1 & validAbsent==1;

allstairs=unique(stair);
allstairsId=[];

for i=1:length(allstairs)
    allstairsId=[allstairsId find(file.p.stairs==allstairs(i))];
end


hit_rates_ex=[];
fa_rates_ex=[];
cons_dprime_ex=[];
cons_c_ex=[];
miss_ex=[];
cr_ex=[];


for i=1:length(allstairs)

    % Expected
    hitsEstair=(hits==1) & (stair==allstairs(i));
    hit_rate_ex=mean(hitsEstair);
    hit_rates_ex(i)=hit_rate_ex;

    faEstair=(fa_ex==1) & (stair==allstairs(i));
    fa_rate_ex=mean(faEstair);
    fa_rates_ex(i)=fa_rate_ex;

    %   get dprime and criterion  
    nh=sum(hitsEstair);
    nfa=sum(faEstair);
    signal=validPresent==1 & (stair==allstairs(i));
    nsignal=sum(signal);
    noise=validAbsent==1 & (stair==allstairs(i));
    nnoise=sum(noise);
    if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
           nh = nh + 0.5;
           nfa = nfa + 0.5;
           nnoise = nnoise + 1;
           nsignal = nsignal + 1;
           % proportions
           hit_rate_ex = nh./nsignal;
           fa_rate_ex = nfa./nnoise;
    end

    cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
    cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
    
    missEstair=miss_ex1==1 & (stair==allstairs(i));
    miss_ex(i)=mean(missEstair);
    crEstair=cr_ex1==1 & (stair==allstairs(i));
    cr_ex(i)=mean(crEstair);

end

   
% stair ids
figure();
subplot(2,1,1);
plot(allstairsId,cons_dprime_ex,'blue-v');
xlabel("stair IDs")
ylabel("d'");
title("d'");

subplot(2,1,2);
plot(allstairsId,cons_c_ex,'blue-v');
xlabel("stair IDs")
ylabel("criterion");
title("criterion, ",subjectID);

% stair values

figure();
subplot(2,1,1);
plot(allstairs,cons_dprime_ex,'blue-v');
xlabel("stairs")
ylabel("d'");
title("d'");

subplot(2,1,2);
plot(allstairs,cons_c_ex,'blue-v');
xlabel("stairs")
ylabel("criterion");
title("criterion, ",subjectID);

% responses= cat(4,testContrasts,responses,testOri,precueValidities);
% nTrials=size(testContrasts,2);
% responses=reshape(responses,[nTrials,4]);
% responses=sortrows(responses,1); 


% responseEx=[];
% responseUnex=[];
% for i = 1:length(orientationJudgement) 
%     if orientationJudgement(i,2)==1
%         responseEx(i,:)=orientationJudgement(i,:);
%     elseif orientationJudgement(i,2)==2
%         responseUnex(i,:)=orientationJudgement(i,:);
%     end
% end
% 
% skipRowsEx=find(responseEx(:,1)==0);
% responseEx(skipRowsEx,:)=[];
% 
% skipRowsUnex=find(responseUnex(:,1)==0);
% responseUnex(skipRowsUnex,:)=[];

%%Commented Out
% %figure();
% %plot(trial_count,)
% 
% 
% % %% p(correct)
% % 
% % gratingContrastsIdx=file.d.gratingContrast;
% % skipWaffle=isnan(gratingContrastsIdx);
% % gratingContrastsIdx(skipWaffle)=[];
% % 
% % response= file.d.targetResponseKey;
% % response(skipWaffle)=[];
% % 
% % stairIdx=file.d.stairIdx;
% % stairIdx(skipWaffle)=[];
% % 
% % static=file.d.staticGrating;
% % static(skipWaffle(1:length(static)))=[];
% % 
% % precueValidities=file.d.precueValidity;
% % precueValidities(skipWaffle)=[];
% % 
% % correctness= cat(3,response,precueValidities,stairIdx, static);
% % 
% % nTrialsGrating=length(gratingContrastsIdx);
% % correctness=reshape(correctness,[nTrialsGrating,4]);
% % 
% % correctness=sortrows(correctness,3);  %sort based on stairIdx
% % correctness=sortrows(correctness,2); %sort expected vs unexpected
% % 
% % prop_exp=(nTrialsGrating/length(file.p.stairs))*(mean(file.p.precueValidities(1,:)==1));
% % prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
% % 
% % 
% 
% %% subjective report plot
% 
% 
% waffleContrastsIdx=file.d.plaidContrast;
% skipNotWaffle=isnan(waffleContrastsIdx);
% waffleContrastsIdx(skipNotWaffle)=[];
% waffleContrasts=[];
% difference=[];
% for i=1:length(waffleContrastsIdx)
%     index=waffleContrastsIdx(i);
%     waffleContrasts(i)=file.p.plaidContrasts1(index);
%     difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
% end
% 
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipNotWaffle)=[];
% 
% waffleOrientation=file.d.plaidOrientation;
% waffleOrientation(skipNotWaffle)=[];
% 
% responses= file.d.targetResponseKey;
% responses(skipNotWaffle)=[];
% 
% dataMatrix=cat(4,responses,precueValidities,waffleOrientation,difference);
% nTrialsWaffle=length(waffleContrastsIdx);
% dataMatrix=reshape(dataMatrix,[nTrialsWaffle,4]);
% % sort based on precue validity
% dataMatrix=sortrows(dataMatrix,2); 
% % sort based on difference
% dataMatrix=sortrows(dataMatrix,4);
% % sort based on orientation
% dataMatrix=sortrows(dataMatrix,3); 
% 
% allDiff=unique(difference);
% nDiff=length(allDiff);
% 
% prop_con=nTrialsWaffle/length(allDiff);
% cons_45pex=[];
% cons_45nex=[];
% cons_45punex=[];
% cons_45nunex=[];
% cons_135pex=[];
% cons_135nex=[];
% cons_135punex=[];
% cons_135nunex=[];
% % cons_45=[];
% % cons_135=[];
% cons_45p=[];
% cons_45n=[];
% cons_135p=[];
% cons_135n=[];
% 
% for i=1:length(allDiff)
%     cons_45p(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2); 
%     cons_45n(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1);
%     cons_135p(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1);
%     cons_135n(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2);
% 
%     cons_45pex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==1);
%     cons_45nex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==1);
%     cons_45punex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==2);
%     cons_45nunex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==2);
% 
%     cons_135pex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==1);
%     cons_135nex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==1);
%     cons_135punex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==2);
%     cons_135nunex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==2);
% 
% end
% 
% cons_45=[flip(cons_45n) cons_45p];
% cons_135=[cons_135p flip(cons_135n)];
% allDiffFull=[allDiff*(-1) allDiff];
% allDiffFull=sort(allDiffFull);
% figure();
% yyaxis left
% plot(allDiffFull,cons_45,"-v");
% xlabel("contrast difference") %grating(45)-grating(135)
% ylabel("p(reported +45 stronger)");
% yyaxis right
% plot(allDiffFull,cons_135,"-v");
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger), ",subjectID);
% legend("+45","-45");
% 
% cons_45ex=[flip(cons_45nex) cons_45pex];
% cons_135ex=[cons_135pex flip(cons_135nex)];
% cons_45unex=[flip(cons_45nunex) cons_45punex];
% cons_135unex=[cons_135punex flip(cons_135nunex)];
% allDiffFull=[allDiff*(-1) allDiff];
% allDiffFull=sort(allDiffFull);
% figure();
% yyaxis left
% plot(allDiffFull,cons_45ex,"-v");
% hold on
% plot(allDiffFull,cons_45unex,"--v");
% xlabel("contrast difference") %grating(45)-grating(135)
% ylabel("p(reported +45 stronger)");
% yyaxis right
% plot(allDiffFull,cons_135ex,"-v");
% hold on
% plot(allDiffFull,cons_135unex,"--v");
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger), ",subjectID);
% legend("+45 expected","+45 unexpected","-45 expected","-45 unexpected");
% 
% %% subjective report plot attempt 2
% 
% waffleContrastsIdx=file.d.plaidContrast;
% skipNotWaffle=isnan(waffleContrastsIdx);
% waffleContrastsIdx(skipNotWaffle)=[];
% waffleContrasts=[];
% difference=[];
% waffleOrientation=file.d.plaidOrientation;
% waffleOrientation(skipNotWaffle)=[];
% 
% for i=1:length(waffleContrastsIdx)
%     index=waffleContrastsIdx(i);
%     waffleContrasts(i)=file.p.plaidContrasts1(index);
%     if waffleOrientation(i)==2
%         difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
%     elseif waffleOrientation(i)==1
%         difference(i)=file.p.plaidContrasts2-waffleContrasts(i);
%     end
% end
% 
% allDiff=unique(difference);
% nDiff=length(allDiff);
% 
% responses= file.d.targetResponseKey;
% responses(skipNotWaffle)=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipNotWaffle)=[];
% 
% cons_45=[];
% cons_135=[];
% cons_45ex=[];
% cons_45unex=[];
% cons_135ex=[];
% cons_135unex=[];
% for i=1:length(allDiff)
%     cons_45(i)=mean((difference==allDiff(i))&(responses==2)); %for X difference they responsed that +45 was stronger
%     cons_135(i)=mean((difference==allDiff(i))&(responses==1)); %for X difference they responsed that -45 was stronger
% 
%     cons_45ex(i)=mean((difference==allDiff(i))&(responses==2)&(precueValidities==1));
%     cons_45unex(i)=mean((difference==allDiff(i))&(responses==2)&(precueValidities==2));
% 
%     cons_135ex(i)=mean((difference==allDiff(i))&(responses==1)&(precueValidities==1));
%     cons_135unex(i)=mean((difference==allDiff(i))&(responses==1)&(precueValidities==2));
% end
% 
% figure();
% yyaxis left
% plot(allDiff,cons_45,"-v");
% xlabel("contrast difference") %grating(45)-grating(135)
% ylabel("p(reported +45 stronger)");
% yyaxis right
% plot(allDiff,cons_135,"-v");
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger), ",subjectID);
% legend("+45","-45");
% 
% figure();
% yyaxis left
% plot(allDiff,cons_45ex,"-v");
% hold on
% plot(allDiff,cons_45unex,"--v");
% xlabel("contrast difference") %grating(45)-grating(135)
% ylabel("p(reported +45 stronger)");
% yyaxis right
% plot(allDiff,cons_135ex,"-v");
% hold on
% plot(allDiff,cons_135unex,"--v");
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger), ",subjectID);
% legend("+45 expected","+45 unexpected","-45 expected","-45 unexpected");
% 
% 
% %% p(correct orientation discrimination)
% 
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% stairs=file.d.stair;
% stairs(skipWaffle(1:length(stairs)))=[];
% 
% allStairs=file.p.stairs;
% nStairs=length(allStairs);
% 
% 
% correct=file.d.correct;
% correct(skipWaffle(1:length(correct)))=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% 
% cons_correct=[];
% cons_correct_exp=[];
% cons_correct_unexp=[];
% 
% for i=1:length(allStairs)
%     cons_correct(i)=mean((stairs==allStairs(i))& (correct==1)); %for X difference they responsed that +45 was stronger
%     cons_correct_exp(i)=mean((stairs==allStairs(i))&(correct==1)&(precueValidities==1));
%     cons_correct_unexp(i)=mean((stairs==allStairs(i))&(correct==1)& (precueValidities==2));
% end
% 
% figure();
% plot(allStairs,cons_correct,"-v");
% xlabel("difference in orientation (stair)")
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) ",subjectID);
% 
% figure();
% plot(allStairs,cons_correct_exp,"-v");
% hold on
% plot(allStairs,cons_correct_unexp,"-v");
% xlabel("difference in orientation (stair)") 
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) for expected vs unexpected, ",subjectID);
% legend("expected","unexpected");
% 
% figure();
% plot([1:20],cons_correct,"-v");
% xlabel("stair id") 
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) ",subjectID);
% 
% figure();
% plot([1:20],cons_correct_exp,"-v");
% hold on
% plot([1:20],cons_correct_unexp,"-v");
% xlabel("stair id") 
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) for expected vs unexpected, ",subjectID);
% % legend("expected","unexpected");
% 
% %% p(correct orientation discrimination)
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% stairs=file.d.stair;
% stairs(skipWaffle(1:length(stairs)))=[];
% 
% allStairs=file.p.stairs;
% nStairs=length(allStairs);
% 
% correct=file.d.correct;
% correct(skipWaffle(1:length(correct)))=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% 
% response=file.d.targetResponseKey;
% response(skipWaffle(1:length(response)))=[];
% 
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% 
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% cons_correct=[];
% cons_correct_exp=[];
% cons_correct_unexp=[];
% correct_r=[];
% 
% correct_r=(response==1 & static==1 & diffMult==-1) | (response==1 & static==2 & diffMult==1) | (response==2 & static==1 & diffMult==1)| (response==2 & static==2 & diffMult==-1);
% 
% for i=1:length(allStairs)
%     cons_correct(i)=mean((stairs==allStairs(i))& (correct_r==1)); %for X difference they responsed that +45 was stronger
%     cons_correct_exp(i)=mean((stairs==allStairs(i))&(correct_r==1)&(precueValidities==1));
%     cons_correct_unexp(i)=mean((stairs==allStairs(i))&(correct_r==1)& (precueValidities==2));
% end
% 
% figure();
% plot(allStairs,cons_correct,"-v");
% xlabel("difference in orientation (stair)")
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) ",subjectID);
% 
% figure();
% plot(allStairs,cons_correct_exp,"-v");
% hold on
% plot(allStairs,cons_correct_unexp,"-v");
% xlabel("difference in orientation (stair)") 
% ylabel("p(correct orientation discrimination)");
% title("p(correct orientation discrimination) for expected vs unexpected, ",subjectID);
% legend("expected","unexpected");



end
