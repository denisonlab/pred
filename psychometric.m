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
% pilot002 has 2 session, 896 trials with mask contrast = 0.8 and range 20%
% to 80%
clear;
close all;

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

plotdir = sprintf('%s/plots',pwd);
if ~exist(plotdir, 'dir')
    mkdir(plotdir)
end

plotdir_subjectID = sprintf('%s/plots/%s',pwd,subjectID);
if ~exist(plotdir_subjectID, 'dir')
    mkdir(plotdir_subjectID)
end

if choice==1
    plotdir_sub = sprintf('%s/plots/%s/session_%d',pwd,subjectID,sessionNum);
    if ~exist(plotdir_sub, 'dir')
        mkdir(plotdir_sub)
    end
else
    total=length(sessionFiles);
    plotdir_sub = sprintf('%s/plots/%s/%d_combined',pwd,subjectID,total);
    if ~exist(plotdir_sub, 'dir')
        mkdir(plotdir_sub)
    end
end



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
phigher= ismember(phigher,[1 4]);
%phigher= phigher == (1 || 4);
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
ylabel("p(test stronger)");
title("p(test stronger) for expected vs unexpected ",subjectID);
legend("expected","unexpected");
saveas(figure(1),[plotdir_sub,'/pteststronger_cond_' subjectID '.jpg']);
%xscale log;
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
title("p(correct) for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
saveas(figure(2),[plotdir_sub,'/pcorrect_cond_' subjectID '.jpg']);
%xscale log;
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
saveas(figure(3),[plotdir_sub,'/pcorrect_' subjectID '.jpg']);

%xscale log;
%% P(test stronger)

% first we get all responses
phigher= file.d.targetResponseKey ;

skipRows=find(phigher==0);
phigher(skipRows)=[];



phigher= ismember(phigher,[1 4]);
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
ylabel("p(test stronger)");
title("p(test stronger), ",subjectID);
saveas(figure(4),[plotdir_sub,'/pteststronger_' subjectID '.jpg']);

%xscale log;


%% Does performance effect go away?

% first we get all responses
phigher= file.d.targetResponseKey ;

%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(phigher==0);
% delete skipped trials based on indices
phigher(skipRows)=[];
phigher= ismember(phigher,[1 4]);

%clean all used variables by removing skipped trials
testContrasts=file.d.testContrast;
testContrasts(skipRows)=[];
precueValidities=file.d.precueValidity;
precueValidities(skipRows)=[];

phigher= cat(3,testContrasts,phigher,precueValidities);
nTrials=size(testContrasts,2);
phigher=reshape(phigher,[nTrials,3]);

phigher1=[];
phigher2=[];
if choice==1
    phigher1=phigher(1:nTrials/2,:,:);
    phigher2=phigher((nTrials/2+1):nTrials,:,:);

elseif choice==2
    n=length(sessionFiles);
    start=0;
    for i= 1:n
        phigher1=[phigher1; phigher(1+start:start+nTrials/(n*2),:,:);];
        phigher2=[phigher2 ;phigher(start+nTrials/(n*2)+1:start+nTrials/n,:,:)];
        start=start+nTrials/n;
    end
end


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
plot(allcons,cons_ex1,'blue-v');
hold on
plot(allcons,cons_unex1,'red-v');
xlabel("contrasts")
ylabel("perceived test as higher contrast");
title("p(test stronger) 1st half for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
%xscale log;

subplot(2,1,2);
plot(allcons,cons_ex2,'blue-v');
hold on
plot(allcons,cons_unex2,'red-v');
xlabel("contrasts")
ylabel("p(test stronger)");
title("p(test stronger) 2nd half for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
%xscale log;
saveas(figure(5),[plotdir_sub,'/pteststronger_split' subjectID '.jpg']);

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
corrects1=[];
corrects2=[];
if choice==1
    corrects1=corrects(1:nTrials/2,:,:);
    corrects2=corrects((nTrials/2+1):nTrials,:,:);

elseif choice==2
    n=length(sessionFiles);
    start=0;
    for i= 1:n
        corrects1=[corrects1; corrects(1+start:start+nTrials/(n*2),:,:);];
        corrects2=[corrects2 ;corrects(start+nTrials/(n*2)+1:start+nTrials/n,:,:)];
        start=start+nTrials/n;
    end
end

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
title("p(correct) 1st half for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
%xscale log;
subplot(2,1,2);
plot(allcons(2:7),cons_ex2(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_unex2(2:7),'red-v');
xlabel("contrasts")
ylabel("p(correct)");
title("p(correct) 2nd half for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
%xscale log;
saveas(figure(6),[plotdir_sub,'/pcorrect_split_' subjectID '.jpg']);
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
responses=sortrows(responses,1); 
% sort based on precue validity
responses=sortrows(responses,3);
%nTrials=length(phigher);
prop_exp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));

cons_highex=[];
cons_highunex=[];
gratingContrasts=file.p.gratingContrasts(1,:);
stand=find(gratingContrasts==file.p.standardContrast(1,:));

for i = 1:length(file.p.gratingContrasts(1,:))
    if i<stand
        ansExW=ismember(responses(prop_exp*i-prop_exp+1:prop_exp*i,2),[2 3]);
        cons_highex(i)=mean(ansExW);
        start=prop_exp*length(file.p.gratingContrasts(1,:))+1;
        ansUnexW=ismember(responses(prop_unexp*i-prop_unexp+1:prop_unexp*i,2),[2 3]);
        cons_highunex(i)=mean(ansUnexW);
    elseif i>=stand
        ansExW=ismember(responses(prop_exp*i-prop_exp+1:prop_exp*i,2),[1 4]);
        cons_highex(i)=mean(ansExW);
        start=prop_exp*length(file.p.gratingContrasts(1,:))+1;
        ansUnexW=ismember(responses(prop_unexp*i-prop_unexp+1:prop_unexp*i,2),[1 4]);
        cons_highunex(i)=mean(ansUnexW);
        
    end
end
figure();
plot(allcons,cons_highex,'blue-v');
hold on
plot(allcons,cons_highunex,'red-v');
xlabel("contrasts")
ylabel("p(correct strength comparison judgement)");
title("p(correct strength comparison judgement) for expected vs unexpected, ",subjectID);
legend("expected","unexpected");
%xscale log;
saveas(figure(7),[plotdir_sub,'/pcorrectstrengthjudgement_cond_' subjectID '.jpg']);

%% d' and criterion general

% %signal present= -45 (1)
% %signal absent = +45 (1)
signalPresent=1;
signalAbsent=2;

% first we get all responses
responses= file.d.targetResponseKey ;
%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(responses==0);
% delete skipped trials based on indices
responses(skipRows)=[];

% get test orientations
testOri=file.d.testOrientation;
% delete skipped trials based on indices
testOri(skipRows)=[];

% get test contrasts
testContrasts=file.d.testContrast;
% delete skipped trials based on indices
testContrasts(skipRows)=[];


responses= cat(3,testContrasts,responses,testOri);
nTrials=size(testContrasts,2);
responses=reshape(responses,[nTrials,3]);
responses=sortrows(responses,1); 
prop_con=nTrials/length(file.p.gratingContrasts(1,:));


cons_dprime=[];
cons_c=[];
hit_rates=[];
fa_rates=[];
miss=[];
cr=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    
    responsePresent=ismember(responses(prop_con*i-prop_con+1:prop_con*i,2),[1 2]);
    orientationPresent=(responses(prop_con*i-prop_con+1:prop_con*i,3)==1);
%   get hits
    hits=responsePresent==1 & orientationPresent==1;
    hit_rate=mean(hits);
    hit_rates(i)=hit_rate;
%   get false alarms
    orientationAbsent=(responses(prop_con*i-prop_con+1:prop_con*i,3)==2);
    fa=responsePresent==1 & orientationAbsent==1;
    fa_rate=mean(fa);
    fa_rates(i)=fa_rate;
%   get dprime and criterion    
    if fa_rate==0
        fa_rate=1/(length(responses));
    elseif fa_rate==1 
        fa_rate=(length(responses)-1)/(length(responses));
    elseif hit_rate==0 
        hit_rate=1/(length(responses));
    elseif hit_rate==1
        hit_rate=(length(responses)-1)/(length(responses));
    end

    cons_dprime(i)=norminv(hit_rate) - norminv(fa_rate);
    cons_c(i)=-(norminv(hit_rate) + norminv(fa_rate))/2;
   
%   get misses
    responseAbsent=ismember(responses(prop_con*i-prop_con+1:prop_con*i,2),[3 4]);
    miss_1= responseAbsent==1 & orientationPresent==1;
    miss(i)=mean(miss_1);
%   get correct rejection
    cr_1=responseAbsent==1 & orientationAbsent==1;
    cr(i)=mean(cr_1);
end
figure();
subplot(2,1,1);
plot(allcons,cons_dprime,'blue-v');

xlabel("contrasts")
ylabel("d'");
title("d'");
%xscale log;

subplot(2,1,2);
plot(allcons,cons_c,'blue-v');

xlabel("contrasts")
ylabel("criterion");
title("criterion");
saveas(figure(8),[plotdir_sub,'/dprime_criterion_' subjectID '.jpg']);


%xscale log;
figure();
plot(allcons,hit_rates,'blue-v');
xlabel("contrasts")
ylabel("p");
title("hits, fa, miss, CR, ",subjectID);
%xscale log;
hold on
plot(allcons,fa_rates,'red-v');
hold on
plot(allcons,miss,'-v');
hold on
plot(allcons,cr,'-v');
legend("hit","fa","miss","CR");
saveas(figure(9),[plotdir_sub,'/hit_miss_fa_cr_' subjectID '.jpg']);

% figure();
% plot(fa_rates,hit_rates,'blue-v');
% xlabel("fa")
% ylabel("hit");
% title("ROC");


%% d' and criterion for expected vs unexpected

% %signal present= -45 (1)
% %signal absent = +45 (1)
signalPresent=1;
signalAbsent=2;

% first we get all responses
responses= file.d.targetResponseKey ;
%find all indices of elements that were skipped (have zero instead of a response key)
skipRows=find(responses==0);
% delete skipped trials based on indices
responses(skipRows)=[];

% get test orientations
testOri=file.d.testOrientation;
% delete skipped trials based on indices
testOri(skipRows)=[];

% get test contrasts
testContrasts=file.d.testContrast;
% delete skipped trials based on indices
testContrasts(skipRows)=[];

% get precue validities
precueValidities=file.d.precueValidity;
% delete skipped trials based on indices
precueValidities(skipRows)=[];

responses= cat(4,testContrasts,responses,testOri,precueValidities);
nTrials=size(testContrasts,2);
responses=reshape(responses,[nTrials,4]);
responses=sortrows(responses,1); 
responseEx=[];
responseUnex=[];
for i = 1:length(responses) 
    if responses(i,4)==1
        responseEx(i,:)=responses(i,:);
    elseif responses(i,4)==2
        responseUnex(i,:)=responses(i,:);
    end
end

skipRowsEx=find(responseEx(:,1)==0);
responseEx(skipRowsEx,:)=[];

skipRowsUnex=find(responseUnex(:,1)==0);
responseUnex(skipRowsUnex,:)=[];

prop_exp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));

cons_dprime_ex=[];
cons_dprime_unex=[];
cons_c_ex=[];
cons_c_unex=[];
hit_rates_ex=[];
hit_rates_unex=[];
fa_rates_ex=[];
fa_rates_unex=[];
miss_ex=[];
miss_unex=[];
cr_ex=[];
cr_unex=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    
    responsePresent=ismember(responseEx(prop_exp*i-prop_exp+1:prop_exp*i,2),[1 2]);
    orientationPresent=(responseEx(prop_exp*i-prop_exp+1:prop_exp*i,3)==1);
%   get hits
    hitsE=responsePresent==1 & orientationPresent==1;
    hit_rate_ex=mean(hitsE);
    hit_rates_ex(i)=hit_rate_ex;
%   get false alarms
    orientationAbsent=(responseEx(prop_exp*i-prop_exp+1:prop_exp*i,3)==2);
    fa_ex=responsePresent==1 & orientationAbsent==1;
    fa_rate_ex=mean(fa_ex);
    fa_rates_ex(i)=fa_rate_ex;
%   get dprime and criterion    
    if fa_rate_ex==0
        fa_rate_ex=1/(length(responseEx));
    elseif fa_rate_ex==1 
        fa_rate_ex=(length(responseEx)-1)/(length(responseEx));
    elseif hit_rate_ex==0 
        hit_rate_ex=1/(length(responseEx));
    elseif hit_rate_ex==1
        hit_rate_ex=(length(responseEx)-1)/(length(responseEx));
    end

    cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
    cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
   
%   get misses
    responseAbsent=ismember(responseEx(prop_exp*i-prop_exp+1:prop_exp*i,2),[3 4]);
    miss_ex1= responseAbsent==1 & orientationPresent==1;
    miss_ex(i)=mean(miss_ex1);
%   get correct rejection
    cr_ex1=responseAbsent==1 & orientationAbsent==1;
    cr_ex(i)=mean(cr_ex1);

    

    responsePresentUn=ismember(responseUnex(prop_unexp*i-prop_unexp+1:prop_unexp*i,2),[1 2]);
    orientationPresentUn=(responseUnex(prop_unexp*i-prop_unexp+1:prop_unexp*i,3)==1);
%   get hits
    hitsU=responsePresentUn==1 & orientationPresentUn==1;
    hit_rate_unex=mean(hitsU);
    hit_rates_unex(i)=hit_rate_unex;
%   get false alarms
    orientationAbsentUn=(responseUnex(prop_unexp*i-prop_unexp+1:prop_unexp*i,3)==2);
    fa_unex=responsePresentUn==1 & orientationAbsentUn==1;
    fa_rate_unex=mean(fa_unex);
    fa_rates_unex(i)=fa_rate_unex;

    if fa_rate_unex==0
        fa_rate_unex=1/(length(responseUnex));
    elseif fa_rate_unex==1 
        fa_rate_unex=(length(responseUnex)-1)/(length(responseUnex));
    elseif hit_rate_unex==0 
        hit_rate_unex=1/(length(responseUnex));
    elseif hit_rate_unex==1
        hit_rate_unex=(length(responseUnex)-1)/(length(responseUnex));
    end
%   get dprime and criterion
    cons_dprime_unex(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
    cons_c_unex(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
%   get misses
    responseAbsentUn=ismember(responseUnex(prop_unexp*i-prop_unexp+1:prop_unexp*i,2),[3 4]);
    miss_unex1= responseAbsentUn==1 & orientationPresentUn==1;
    miss_unex(i)=mean(miss_unex1);
%   get correct rejection
    cr_unex1=responseAbsentUn==1 & orientationAbsentUn==1;
    cr_unex(i)=mean(cr_unex1);

end
figure();
subplot(2,1,1);
plot(allcons(2:7),cons_dprime_ex(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_dprime_unex(2:7),'red-v');
xlabel("contrasts")
ylabel("d'");
title("d' expected vs unexpected");
legend("expected","unexpected");
%xscale log;

subplot(2,1,2);
plot(allcons(2:7),cons_c_ex(2:7),'blue-v');
hold on
plot(allcons(2:7),cons_c_unex(2:7),'red-v');
xlabel("contrasts")
ylabel("criterion");
title("criterion expected vs unexpected, ",subjectID);
legend("expected","unexpected");
saveas(figure(10),[plotdir_sub,'/dprime_criterion_split_' subjectID '.jpg']);
%xscale log;

% %% split by orientation
% 
% corrects=file.d.correct;
% skipRowsCorrect=find(isnan(corrects));
% % delete skipped trials based on indices
% corrects(skipRowsCorrect)=[];
% %clean all used variables by removing skipped trials
% testContrasts=file.d.testContrast;
% testContrasts(skipRowsCorrect)=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipRows)=[];
% 
% 
% corrects= cat(3,testContrasts,corrects,precueValidities);
% nTrials=size(testContrasts,2);
% corrects=reshape(corrects,[nTrials,3]);
% corrects1=[];
% corrects2=[];
% if choice==1
%     corrects1=corrects(1:nTrials/2,:,:);
%     corrects2=corrects((nTrials/2+1):nTrials,:,:);
% 
% elseif choice==2
%     n=length(sessionFiles);
%     start=0;
%     for i= 1:n
%         corrects1=[corrects1; corrects(1+start:start+nTrials/(n*2),:,:);];
%         corrects2=[corrects2 ;corrects(start+nTrials/(n*2)+1:start+nTrials/n,:,:)];
%         start=start+nTrials/n;
%     end
% end
% 
% corrects1=sortrows(corrects1,1); 
% corrects1=sortrows(corrects1,3);
% corrects2=sortrows(corrects2,1); 
% corrects2=sortrows(corrects2,3);
% 
% prop_exp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
% prop_unexp1=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
% prop_exp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==1));
% prop_unexp2=((nTrials/2)/length(file.p.gratingContrasts(1,:)))*(mean(file.p.precueValidities(1,:)==2));
% 
% cons_ex1=[];
% cons_unex1=[];
% cons_ex2=[];
% cons_unex2=[];
% 
% for i = 1:length(file.p.gratingContrasts(1,:))
%     cons_ex1(i)=mean(corrects1(prop_exp1*i-prop_exp1+1:prop_exp1*i,2)==1);
%     start1=prop_exp1*length(file.p.gratingContrasts(1,:))+1;
%     cons_unex1(i)=mean(corrects1((start1+prop_unexp1*i)-prop_unexp1:start1+prop_unexp1*i-1,2)==1);
%     cons_ex2(i)=mean(corrects2(prop_exp2*i-prop_exp2+1:prop_exp2*i,2)==1);
%     start2=prop_exp2*length(file.p.gratingContrasts(1,:))+1;
%     cons_unex2(i)=mean(corrects2((start2+prop_unexp2*i)-prop_unexp2:start2+prop_unexp2*i-1,2)==1);
% end
% 
% figure();
% subplot(2,1,1);
% plot(allcons(2:7),cons_ex1(2:7),'blue-v');
% hold on
% plot(allcons(2:7),cons_unex1(2:7),'red-v');
% xlabel("contrasts")
% ylabel("p(correct)");
% title("p(correct) 1st half for expected vs unexpected, ",subjectID);
% legend("expected","unexpected");
% %xscale log;
% subplot(2,1,2);
% plot(allcons(2:7),cons_ex2(2:7),'blue-v');
% hold on
% plot(allcons(2:7),cons_unex2(2:7),'red-v');
% xlabel("contrasts")
% ylabel("p(correct)");
% title("p(correct) 2nd half for expected vs unexpected, ",subjectID);
% legend("expected","unexpected");
% %xscale log;
% saveas(figure(6),[plotdir_sub,'/pcorrect_split_ori_' subjectID '.jpg']);
% 

end
