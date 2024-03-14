function psychometric()

% subjectID = input('Enter subject ID:  ','s');
% total=input('Enter number of sessions (1,2...) ');
% sessionNum = input('Enter session number (1,2...) ');
% dataDir= sprintf('%s/data',pwd);
% datafile=sprintf('%s/%s_s%d_pred.mat',dataDir,subjectID,sessionNum);
% file=load(datafile);

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
phigher= file.d.targetResponseKey 
phigher(phigher == 0) = NaN;
phigher= phigher == (1 || 4);
phigher= cat(3,file.d.testContrast,phigher,file.d.precueValidity);
nTrials=size(file.d.testStatus,2);
phigher=reshape(phigher,[nTrials,3]);
phigher = phigher(~any(isnan(phigher), 2), :);
phigher=sortrows(phigher,1);
phigher=sortrows(phigher,3);
nTrials=length(phigher);
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
plot(allcons,cons_highex,'blue');
hold on
plot(allcons,cons_highunex,'red');
xlabel("contrasts")
ylabel("percent that answered 0 or 1 aka perceived test as higher contrast");
title("p(test stronger) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");

%% P(correct) for unexpected vs expected

corrects=cat(3,file.d.testContrast,file.d.correct,file.d.precueValidity);
nTrials=size(file.d.testStatus,2);
corrects=reshape(corrects,[nTrials,3]);
corrects = corrects(~any(isnan(corrects), 2), :);
corrects=sortrows(corrects,1);
corrects=sortrows(corrects,3);
nTrials=length(corrects);
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
plot(allcons(2:7),cons_ex(2:7),'blue');
hold on
plot(allcons(2:7),cons_unex(2:7),'red');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");

%% P(correct)

corrects=cat(2,file.d.testContrast,file.d.correct);
nTrials=size(file.d.testStatus,2);
corrects=reshape(corrects,[nTrials,2]);
corrects = corrects(~any(isnan(corrects), 2), :);
corrects=sortrows(corrects,1);
nTrials=length(corrects);
prop_con=nTrials/length(file.p.gratingContrasts(1,:));
%prop_exp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==1));
%prop_unexp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==2));
cons_cor=[];

for i = 1:length(file.p.gratingContrasts(1,:))
    cons_cor(i)=mean(corrects(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons(2:7),cons_cor(2:7),'blue');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct)");

%% P(test stronger)
phigher= file.d.targetResponseKey 
phigher(phigher == 0) = NaN;
phigher= phigher == (1 || 4);
phigher= cat(2,file.d.testContrast,phigher);
phigher = phigher(~any(isnan(phigher), 2), :);
nTrials=size(file.d.testStatus,2);
phigher=reshape(phigher,[nTrials,2]);
phigher=sortrows(phigher,1);
nTrials=length(phigher);
prop_con=nTrials/length(file.p.gratingContrasts(1,:));
cons_high=[];
for i = 1:length(file.p.gratingContrasts(1,:))
    cons_high(i)=mean(phigher(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons,cons_high,'blue');
xlabel("contrasts")
ylabel("percent that answered 0 or 1 aka perceived test as higher contrast");
title("p(test stronger)");
end
