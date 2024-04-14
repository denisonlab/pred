function plotsWaffle()

clear;
close all;
% 
% subjectID = input('Enter subject ID:  ', 's');
% dataDir = sprintf('%s/data', pwd);
% 
% % Prompt user for single file or all sessions
% %choice = input('Analyze single file (1) or all sessions (2): ');
% 
% choice=1;
% 
% if choice==1
%     sessionNum = input('Enter session number (1,2...) ');
%     datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
%     datafile = sprintf('%s/%s_s%d_predv2_s240412_1707.mat', datafileDir, subjectID, sessionNum);
%     file = load(datafile);
% 
% end
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


%figure();
%plot(trial_count,)


% %% p(correct)
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle)=[];
% 
% response= file.d.targetResponseKey;
% response(skipWaffle)=[];
% 
% stairIdx=file.d.stairIdx;
% stairIdx(skipWaffle)=[];
% 
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% 
% correctness= cat(3,response,precueValidities,stairIdx, static);
% 
% nTrialsGrating=length(gratingContrastsIdx);
% correctness=reshape(correctness,[nTrialsGrating,4]);
% 
% correctness=sortrows(correctness,3);  %sort based on stairIdx
% correctness=sortrows(correctness,2); %sort expected vs unexpected
% 
% prop_exp=(nTrialsGrating/length(file.p.stairs))*(mean(file.p.precueValidities(1,:)==1));
% prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
% 
% 

%% subjective report plot

% +45 stronger
% hits = responsePresent & orientationPresent; they said +45 and was
% actually +45
% fa= responsePresent & orientationAbsent;

waffleContrastsIdx=file.d.plaidContrast;
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];
waffleContrasts=[];
difference=[];
for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
    difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
end


precueValidities=file.d.precueValidity;
precueValidities(skipNotWaffle)=[];

waffleOrientation=file.d.plaidOrientation;
waffleOrientation(skipNotWaffle)=[];

responses= file.d.targetResponseKey;
responses(skipNotWaffle)=[];

dataMatrix=cat(4,responses,precueValidities,waffleOrientation,difference);
nTrialsWaffle=length(waffleContrastsIdx);
dataMatrix=reshape(dataMatrix,[nTrialsWaffle,4]);
% sort based on precue validity
dataMatrix=sortrows(dataMatrix,2); 
% sort based on difference
dataMatrix=sortrows(dataMatrix,4);
% sort based on orientation
dataMatrix=sortrows(dataMatrix,3); 

allDiff=unique(difference);
nDiff=length(allDiff);

diff45s=[];
diff135s=[];
diff45sEx=[];
diff135sEx=[];
diff45sUnex=[];
diff135sUnex=[];

for i=1:length(allDiff)
    diff45s{i}=(difference==allDiff(i)& waffleOrientation==2); % positive differences
    diff135s{i}=((difference==allDiff(i))& waffleOrientation==1); %negative differences
    
    diff45sEx{i}=((difference==allDiff(i))& waffleOrientation==2 & precueValidities==1);
    diff135sEx{i}=((difference==allDiff(i))& waffleOrientation==1 & precueValidities==1);
    
    diff45sUnex{i}=((difference==allDiff(i))& waffleOrientation==2 & precueValidities==1);
    diff135sUnex{i}=((difference==allDiff(i))& waffleOrientation==1 & precueValidities==1);
end

% prop_exp=(nTrialsWaffle/length(allDiff))*(mean(file.p.precueValidities==1));
% prop_unexp=(nTrialsWaffle/length(allDiff))*(mean(file.p.precueValidities==2));
prop_con=nTrialsWaffle/length(allDiff);
cons_45pex=[];
cons_45nex=[];
cons_45punex=[];
cons_45nunex=[];
cons_135pex=[];
cons_135nex=[];
cons_135punex=[];
cons_135nunex=[];
% cons_45=[];
% cons_135=[];
cons_45p=[];
cons_45n=[];
cons_135p=[];
cons_135n=[];

for i=1:length(allDiff)
    %diff45s{i}=diff45s{i}.';
    %diff135s{i}=diff135s{i}.';
    cons_45p(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2);
    cons_45n(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1);
    cons_135p(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1);
    cons_135n(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2);
    
    cons_45pex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==1);
    cons_45nex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==1);
    cons_45punex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==2);
    cons_45nunex(i)=mean(dataMatrix(:,1)==2 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==2);

    cons_135pex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==1);
    cons_135nex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==1);
    cons_135punex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==1 & dataMatrix(:,2)==2);
    cons_135nunex(i)=mean(dataMatrix(:,1)==1 & dataMatrix(:,4)==allDiff(i) & dataMatrix(:,3)==2 & dataMatrix(:,2)==2);

end

cons_45=[flip(cons_45n) cons_45p];
cons_135=[cons_135p flip(cons_135n)];
allDiffFull=[allDiff*(-1) allDiff];
allDiffFull=sort(allDiffFull);
figure();
yyaxis left
plot(allDiffFull,cons_45,"-v");
xlabel("contrast difference") %grating(45)-grating(135)
ylabel("p(reported +45 stronger)");
yyaxis right
plot(allDiffFull,cons_135,"-v");
xlabel("contrasts")
ylabel("p(reported -45 stronger)");
title("p(reported +45/-45 stronger), ",subjectID);
legend("+45","-45");

cons_45ex=[flip(cons_45nex) cons_45pex];
cons_135ex=[cons_135pex flip(cons_135nex)];
cons_45unex=[flip(cons_45nunex) cons_45punex];
cons_135unex=[cons_135punex flip(cons_135nunex)];
allDiffFull=[allDiff*(-1) allDiff];
allDiffFull=sort(allDiffFull);
figure();
yyaxis left
plot(allDiffFull,cons_45ex,"-v");
hold on
plot(allDiffFull,cons_45unex,"--v");
xlabel("contrast difference") %grating(45)-grating(135)
ylabel("p(reported +45 stronger)");
yyaxis right
plot(allDiffFull,cons_135ex,"-v");
hold on
plot(allDiffFull,cons_135unex,"--v");
xlabel("contrasts")
ylabel("p(reported -45 stronger)");
title("p(reported +45/-45 stronger), ",subjectID);
legend("+45 expected","+45 unexpected","-45 expected","-45 unexpected");
% prop_exp=(nTrialsWaffle/length(allDiff))*(mean(file.p.precueValidities==1));
% prop_unexp=(nTrialsWaffle/length(allDiff))*(mean(file.p.precueValidities==2));
% prop_con=nTrialsWaffle/length(allDiff);
% 
% 
% cons_45ex=[];
% cons_45unex=[];
% cons_135ex=[];
% cons_135unex=[];
% cons_45=[];
% cons_135=[];
% 
% for i=1:length(allDiff)
% 
%     cons_45(i)=mean((dataMatrix(prop_con*i-prop_con+1:prop_con*i,1)==1);
%     cons_135(i)=mean(dataMatrix(prop_con*i-prop_con+1:prop_con*i,1)==2);
% 
%     cons_45ex(i)=mean(dataMatrix(prop_exp*i-prop_exp+1:prop_exp*i,1)==1);
%     cons_135ex(i)=mean(dataMatrix(prop_exp*i-prop_exp+1:prop_exp*i,1)==2);
%     start=prop_exp*length(allDiff)+1;
%     cons_45unex(i)=mean(dataMatrix((start+prop_unexp*i)-prop_unexp:start+prop_unexp*i-1,1)==1);
%     cons_135unex(i)=mean(dataMatrix((start+prop_unexp*i)-prop_unexp:start+prop_unexp*i-1,1)==2);
% end
% 
% figure();
% yyaxis left
% plot(allDiff,cons_45);
% xlabel("contrasts")
% ylabel("p(reported +45 stronger)");
% yyaxis right
% plot(allDiff,cons_135);
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger), ",subjectID);
% legend("+45","-45");
% 
% figure();
% yyaxis left
% plot(allDiff,cons_45ex);
% hold on
% plot(allDiff,cons_45unex);
% xlabel("contrasts")
% ylabel("p(reported +45 stronger)");
% 
% yyaxis right
% plot(allDiff,cons_135ex);
% hold on
% plot(allDiff,cons_135unex);
% xlabel("contrasts")
% ylabel("p(reported -45 stronger)");
% title("p(reported +45/-45 stronger) for expected vs unexpected, ",subjectID);
% legend("expected 45","unexpected 45","expected -45","unexpcted -45");
% 



end
