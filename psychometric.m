function psychometric()

% subjectID = input('Enter subject ID:  ','s');
% total=input('Enter number of sessions (1,2...) ');
% sessionNum = input('Enter session number (1,2...) ');
% dataDir= sprintf('%s/data',pwd);
% datafile=sprintf('%s/%s_s%d_pred.mat',dataDir,subjectID,sessionNum);
% file=load(datafile);

% melissa2 has 2 sessions, each of 448 trials for contrast range 2.5 to 50%
% and no masks
% testBEH has 2 session, 896 trials with mask contrast = 0.8 and range 20%
% to 80%

subjectID = input('Enter subject ID:  ', 's');
dataDir = sprintf('%s/data', pwd);

% Prompt user for single file or all sessions
choice = input('Analyze single file (1) or all sessions (2): ');

% Load data based on user choice
if choice == 1
    % Get session number for single file
    sessionNum = input('Enter session number (1,2...) ');
    datafile = sprintf('%s/%s_s%d_pred.mat', dataDir, subjectID, sessionNum);
    file = load(datafile);
else
    % Find all session files for the subject
    sessionFiles = dir(fullfile(dataDir, [subjectID, '_s*.mat']));
    allDataP = [];
    allDataD=[];

    % Loop through each session file
    for i = 1:length(sessionFiles)
        filename = sessionFiles(i).name;
        filepath = fullfile(dataDir, filename);
        data = load(filepath);

        % Concatenate fields from 'p' structures
        fieldNames = fieldnames(data.p);
        for j = 1:length(fieldNames)
            fieldName = fieldNames{j};
            if ~isfield(allDataP, fieldName)
                allDataP.(fieldName) = [];  % Add missing field if needed
            end
            allDataP.(fieldName) = [allDataP.(fieldName); data.p.(fieldName)];
        end

        % Concatenate fields from 'd' structures (similar logic)
        fieldNames = fieldnames(data.d);
        for j = 1:length(fieldNames)
            fieldName = fieldNames{j};
            if ~isfield(allDataD, fieldName)
                allDataD.(fieldName) = [];  % Add missing field if needed
            end

            allDataD.(fieldName) = [allDataD.(fieldName), data.d.(fieldName)];
        end
    end

    % Assign combined data to 'file' variable
    file.p = allDataP;
    file.d = allDataD;
    % % Load data from each session file
    % for i = 1:length(sessionFiles)
    %     filename = sessionFiles(i).name;
    %     filepath = fullfile(dataDir, filename);
    %     data = load(filepath);
    %     allDataP = [allDataP; data.p];  % Concatenate 'p' structures
    %     allDataD = [allDataD; data.d];  % Concatenate 'd' structures
    % end
    %
    % % Assign the combined data to 'file' variable
    % file = allData;
end

% Accessing specific variables within 'p' and 'd' structures
% Assuming 'allcons' needs the first 7 elements from 'file.p.gratingContrasts'
allcons = file.p.gratingContrasts(1,1:7);

%% P(test stronger) for unexpected vs expected

% first we get all responses
phigher= file.d.targetResponseKey ;

%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(phigher==0);
% delete skipped trials based on indices
phigher(skipRows)=[];

% NOT ANYMORE: if it has a zero, it means there was no response and that trial was
% skipped due to lack of fixation, so make the zeros NaN so they do not
% impact the calculations
%phigher(phigher == 0) = NaN;

% if the response was key 1 or key 0 i.e. 1 or 4, participants perceived
% the test as higher contrast/stronger so make those values 1 and all other
% values 0
phigher= phigher == (1 || 4);

%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];

precueValidities=file.d.precueValidity;
precueValidities(skipRows)=[];

phigher= cat(3,testContrasts,phigher,precueValidities);

nTrials=size(testContrasts,2);

phigher=reshape(phigher,[nTrials,3]);
%phigher = phigher(~any(isnan(phigher), 2), :);

%sort based on contrast levels
phigher=sortrows(phigher,1); 
% sort based on precue validity
phigher=sortrows(phigher,3);
%nTrials=length(phigher);
prop_exp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
cons_highex=[];
cons_highunex=[];
for i = 1:length(file.p.gratingContrasts(1,:))
    cons_highex(i)=mean(phigher(prop_exp*i-prop_exp+1:prop_exp*i,2)==1);
    start=prop_exp*length(file.p.gratingContrasts(1,:))+1;
    cons_highunex(i)=mean(phigher((start+prop_unexp*i)-prop_unexp:start+prop_unexp*i-1,2)==1);
end

figure();
plot(allcons,cons_highex,'blue-v');
hold on
plot(allcons,cons_highunex,'red-v');
xlabel("contrasts")
ylabel("percent that answered 0 or 1 aka perceived test as higher contrast");
title("p(test stronger) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;
%% P(correct) for unexpected vs expected

% no need to do this because skipRowsCorrect and skipRows should be the
% same but just a dummy check

corrects=file.d.correct;
skipRowsCorrect=find(isnan(corrects));

% delete skipped trials based on indices
corrects(skipRowsCorrect)=[];


%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRowsCorrect)=[];

precueValidities=file.d.precueValidity;
precueValidities(skipRowsCorrect)=[];


corrects=cat(3,testContrasts,corrects,precueValidities);
nTrials=size(testContrasts,2);
corrects=reshape(corrects,[nTrials,3]);
%corrects = corrects(~any(isnan(corrects), 2), :);
corrects=sortrows(corrects,1);
corrects=sortrows(corrects,3);
%nTrials=length(corrects);
prop_exp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
cons_ex=[];
cons_unex=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    cons_ex(i)=mean(corrects(prop_exp*i-prop_exp+1:prop_exp*i,2)==1);
    start=prop_exp*length(file.p.gratingContrasts(1,:))+1;
    cons_unex(i)=mean(corrects((start+prop_unexp*i)-prop_unexp:start+prop_unexp*i-1,2)==1);
end

figure();
plot(allcons(2:7),cons_ex(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex(2:7),'red-v');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;
%% P(correct)

corrects=file.d.correct;
skipRowsCorrect=find(isnan(corrects));

% delete skipped trials based on indices
corrects(skipRowsCorrect)=[];


%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRowsCorrect)=[];

corrects=cat(2,testContrasts,corrects);
nTrials=size(testContrasts,2);

corrects=reshape(corrects,[nTrials,2]);
%corrects = corrects(~any(isnan(corrects), 2), :);
corrects=sortrows(corrects,1);
%nTrials=length(corrects);
prop_con=nTrials/length(file.p.gratingContrasts(1,:));
%prop_exp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==1));
%prop_unexp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==2));
cons_cor=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    cons_cor(i)=mean(corrects(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons(2:7),cons_cor(2:7),'blue-v');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct)");
xscale log;
%% P(test stronger)

% first we get all responses
phigher= file.d.targetResponseKey ;

skipRows=find(phigher==0);
phigher(skipRows)=[];



phigher= phigher == (1 || 4);

%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];

phigher= cat(2,testContrasts,phigher);

nTrials=size(testContrasts,2);

%phigher = phigher(~any(isnan(phigher), 2), :);
%nTrials=size(file.d.testStatus,2)-sum(d.stopThisTrial==1);
phigher=reshape(phigher,[nTrials,2]);
phigher=sortrows(phigher,1);
%nTrials=length(phigher);
prop_con=nTrials/length(file.p.gratingContrasts(1,:));
cons_high=[];
for i = 1:length(file.p.gratingContrasts(1,:))
    cons_high(i)=mean(phigher(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons,cons_high,'blue-v');
xlabel("contrasts")
ylabel("percent that answered 0 or 1 aka perceived test as higher contrast");
title("p(test stronger)");
xscale log;


%% Does performance effect go away?

% first we get all responses
phigher= file.d.targetResponseKey ;

%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(phigher==0);
% delete skipped trials based on indices
phigher(skipRows)=[];
phigher= phigher == (1 || 4);

%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];
precueValidities=file.d.precueValidity;
precueValidities(skipRows)=[];

phigher= cat(3,testContrasts,phigher,precueValidities);
nTrials=size(testContrasts,2);
phigher=reshape(phigher,[nTrials,3]);
phigher1=phigher(1:nTrials/2,:,:);
phigher2=phigher((nTrials/2+1):nTrials,:,:);
phigher1=sortrows(phigher1,1); 
phigher1=sortrows(phigher1,3);
phigher2=sortrows(phigher2,1); 
phigher2=sortrows(phigher2,3);

prop_exp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
prop_exp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));

cons_ex1=[];
cons_unex1=[];
cons_ex2=[];
cons_unex2=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    cons_ex1(i)=mean(phigher1(prop_exp1*i-prop_exp1+1:prop_exp1*i,2)==1);
    start1=prop_exp1*length(file.p.gratingContrasts(1,:))+1;
    cons_unex1(i)=mean(phigher1((start1+prop_unexp1*i)-prop_unexp1:start1+prop_unexp1*i-1,2)==1);
    cons_ex2(i)=mean(phigher2(prop_exp2*i-prop_exp2+1:prop_exp2*i,2)==1);
    start2=prop_exp2*length(file.p.gratingContrasts(1,:))+1;
    cons_unex2(i)=mean(phigher2((start2+prop_unexp2*i)-prop_unexp2:start2+prop_unexp2*i-1,2)==1);
end

figure();
subplot(2,1,1);
plot(allcons(2:7),cons_ex1(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex1(2:7),'red-v');
xlabel("contrasts")
ylabel("perceived test as higher contrast");
title("p(test stronger) 1st half for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;

subplot(2,1,2);
plot(allcons(2:7),cons_ex2(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex2(2:7),'red-v');
xlabel("contrasts")
ylabel("perceived test as higher contrast");
title("p(test stronger) 2nd half for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;

%% Does performance effect go away?

corrects=file.d.correct;
skipRowsCorrect=find(isnan(corrects));
% delete skipped trials based on indices
corrects(skipRowsCorrect)=[];
%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRowsCorrect)=[];

precueValidities=file.d.precueValidity;
precueValidities(skipRows)=[];

corrects= cat(3,testContrasts,corrects,precueValidities);
nTrials=size(testContrasts,2);
corrects=reshape(corrects,[nTrials,3]);
corrects1=corrects(1:nTrials/2,:,:);
corrects2=corrects((nTrials/2+1):nTrials,:,:);
corrects1=sortrows(corrects1,1); 
corrects1=sortrows(corrects1,3);
corrects2=sortrows(corrects2,1); 
corrects2=sortrows(corrects2,3);

prop_exp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
prop_exp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));

cons_ex1=[];
cons_unex1=[];
cons_ex2=[];
cons_unex2=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    cons_ex1(i)=mean(corrects1(prop_exp1*i-prop_exp1+1:prop_exp1*i,2)==1);
    start1=prop_exp1*length(file.p.gratingContrasts(1,:))+1;
    cons_unex1(i)=mean(corrects1((start1+prop_unexp1*i)-prop_unexp1:start1+prop_unexp1*i-1,2)==1);
    cons_ex2(i)=mean(corrects2(prop_exp2*i-prop_exp2+1:prop_exp2*i,2)==1);
    start2=prop_exp2*length(file.p.gratingContrasts(1,:))+1;
    cons_unex2(i)=mean(corrects2((start2+prop_unexp2*i)-prop_unexp2:start2+prop_unexp2*i-1,2)==1);
end

figure();
subplot(2,1,1);
plot(allcons(2:7),cons_ex1(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex1(2:7),'red-v');
xlabel("contrasts")
ylabel("p(correct)");
title("p(correct) 1st half for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;
subplot(2,1,2);
plot(allcons(2:7),cons_ex2(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex2(2:7),'red-v');
xlabel("contrasts")
ylabel("p(correct)");
title("p(correct) 2nd half for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;

%% D' and criterion (general)

%signal present= +45 (2)
%signal absent = -45 (1)
signalPresent=1;
signalAbsent=2;


response=file.d.targetResponseKey;
%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(response==0);
response(skipRows)=[];

testOri=file.d.testOrientation;
testOri(skipRows)=[];

testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];

%precueValidities=file.d.precueValidity;
%precueValidities(skipRows)=[];

responses= cat(3,testContrasts,response,testOri);
nTrials=size(testContrasts,2);
responses=reshape(responses,[nTrials,3]);
responses=sortrows(responses,1); 

prop_con=nTrials/length(file.p.gratingContrasts(1,:));
cons_dprime=[];
cons_c=[];
hit_rates=[];
fa_rates=[];
miss=[];
CR=[];
for i = 1:length(file.p.gratingContrasts(1,:))
    hits=(responses(prop_con*i-prop_con+1:prop_con*i,2)==(1||2))==(responses(prop_con*i-prop_con+1:prop_con*i,3)==1);
    hit_rate=mean(hits);
    hit_rates(i)=hit_rate;
    fa=(responses(prop_con*i-prop_con+1:prop_con*i,2)==(1||2))==(responses(prop_con*i-prop_con+1:prop_con*i,3)==2);
    fa_rate=mean(fa);
    fa_rates(i)=fa_rate;
    cons_dprime(i)=norminv(hit_rate) - norminv(fa_rate);
    cons_c(i)=-(norminv(hit_rate) + norminv(fa_rate))/2;
    miss1=(responses(prop_con*i-prop_con+1:prop_con*i,2)==(3||4))==(responses(prop_con*i-prop_con+1:prop_con*i,3)==1);
    miss(i)=mean(miss1);
    CR1=(responses(prop_con*i-prop_con+1:prop_con*i,2)==(3||4))==(responses(prop_con*i-prop_con+1:prop_con*i,3)==2);
    CR(i)=mean(CR1);
end

figure();
plot(allcons,hit_rates,'blue-v');
xlabel("contrasts")
ylabel("hits");
title("hits");
xscale log;
figure();
plot(allcons,fa_rates,'blue-v');
xlabel("contrasts")
ylabel("fa");
title("fa");
xscale log;

figure();
plot(allcons,miss,'blue-v');
xlabel("contrasts")
ylabel("miss'");
title("miss");
xscale log;


figure();
plot(allcons,CR,'blue-v');
xlabel("contrasts")
ylabel("CR'");
title("CR");
xscale log;

figure();
subplot(2,1,1);
plot(allcons(1:7),cons_dprime(1:7),'blue-v');
xlabel("contrasts")
ylabel("d'");
title("d'");
xscale log;

subplot(2,1,2);
plot(allcons(1:7),cons_c(1:7),'blue-v');
xlabel("contrasts")
ylabel("criterion");
title("criterion");
xscale log;


%% Correctness based on perceived contrast actually stronger

% first we get all responses
responses= file.d.targetResponseKey ;

%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(responses==0);
% delete skipped trials based on indices
responses(skipRows)=[];

testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];

precueValidities=file.d.precueValidity;
precueValidities(skipRows)=[];

responses= cat(3,testContrasts,responses,precueValidities);
nTrials=size(testContrasts,2);
responses=reshape(responses,[nTrials,3]);


%sort based on contrast levels
phigher=sortrows(phigher,1); 
% sort based on precue validity
phigher=sortrows(phigher,3);
%nTrials=length(phigher);
prop_exp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
cons_highex=[];
cons_highunex=[];
for i = 1:length(file.p.gratingContrasts(1,:))
    cons_highex(i)=mean((responses(prop_exp*i-prop_exp+1:prop_exp*i,2)==(1||4))==(responses(prop_exp*i-prop_exp+1:prop_exp*i,3)==(5||6||7)));
    start=prop_exp*length(file.p.gratingContrasts(1,:))+1;
    cons_highunex(i)=mean((responses(prop_unexp*i-prop_unexp+1:prop_unexp*i,2)==(1||4))==(responses(prop_unexp*i-prop_unexp+1:prop_unexp*i,3)==(5||6||7)));
end



figure();
plot(allcons,cons_highex,'blue-v');
hold on
plot(allcons,cons_highunex,'red-v');
xlabel("contrasts")
ylabel("p(test stronger correct)");
title("p(test stronger correct) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");
xscale log;
end
