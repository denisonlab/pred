function plotsWaffle()

clear;
close all;

subjectID = input('Enter subject ID:  ', 's');
dataDir = sprintf('%s/data', pwd);

% Prompt user for single file or all sessions
%choice = input('Analyze single file (1) or all sessions (2): ');

choice=1;

if choice==1
    sessionNum = input('Enter session number (1,2...) ');
    datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
    datafile = sprintf('%s/%s_s%d_predv2_s240411_1813.mat', datafileDir, subjectID, sessionNum);
    file = load(datafile);

end

%% check that staircase is staircasing
trialCount=length(file.d.stairIdx);
x=[1:1:trialCount];
figure();
plot(x,file.d.stairIdx);

%% plot contrast judgement task from waffles

waffleContrastsIdx=file.d.plaidContrast;
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];
waffleContrasts=[];
for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
end

%skipNotWaffleIdx=find(isnan(waffleContrastsIdx));
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
%prop_exp=(nTrials/length(file.p.plaidContrasts1))*(mean(file.p.precueValidities(1,:)==1));
%prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));

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
                correctsUnex=[correctsEx contrastJudgement(i,4)];
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


figure();
%plot(trial_count,)
end