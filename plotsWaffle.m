function plotsWaffle()

clear;
close all;

dataDir = sprintf('%s/data', pwd);
subjectID="S0120_cleanupcode";
choice = 1; % a single session (1) or many (2)
sessionNum=1;
responseKeyIds=[1 2];

% Load data based on user choice
if choice == 1
    % Get session number for single file
    datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
    date="240503_1628";
    datafile = sprintf('%s/%s_s%d_predv2_s%s.mat', datafileDir, subjectID, sessionNum,date);
    file = load(datafile);
else
    % Find all session files for the subject
    dataDir=sprintf('%s/%s/',dataDir,subjectID);
    dircontent=dir(dataDir);
    subfolder=sum([dircontent.isdir])-2;
    dircontent(1)=[];
    dircontent(1)=[];
    sessionFiles=[];
    dataSubs=[];
    for i=1:subfolder
        dataSub=sprintf('%s%s/', dataDir,dircontent(i).name);
        dataSubs=[dataSubs string(dataSub)]
        sessionFiles = [sessionFiles dir(fullfile(dataSub, '*.mat'))];
    end

    allDataP = [];
    allDataD=[];

    % Loop through each session file
    for i = 1:length(sessionFiles)
        filename = sessionFiles(i).name;
        dataSubDir=dataSubs(i);
        filepath = fullfile(dataSubDir, filename);
        data = load(filepath);
        
        if i==1
            % Concatenate fields from 'p' structures
            fieldNames = fieldnames(data.p);
            for j = 1:length(fieldNames)
                fieldName = fieldNames{j};
                if ~isfield(allDataP, fieldName)
                    allDataP.(fieldName) = [];  % Add missing field if needed
                end
                allDataP.(fieldName) = [allDataP.(fieldName); data.p.(fieldName)];
            end
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

    % Assign combined data to 'file' variale
    file.p = allDataP;
    file.d = allDataD;
end

% Accessing specific variables within 'p' and 'd' structures
% Assuming 'allcons' needs the first 7 elements from 'file.p.gratingContrasts'

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

f=1;

%% check that staircase is staircasing
if choice==1 
    figure();
    stair=file.d.stairIdx;
    skipEye=(file.d.stairIdx==0);
    stair(skipEye(1:length(stair)))=[];
    skipNotStair=isnan(stair);
    stair(skipNotStair)=[];
    trialCount=length(stair);
    x=[1:1:trialCount];
    plot(x,stair);
    ylabel("stair ID");
    xlabel("trial ID");
    saveas(figure(f),[plotdir_sub,'/staircasing.jpg']);
    f=f+1;
    figure();
    stairExp=file.d.stairIdxExp;
    skipEyeEx=(file.d.stairIdxExp==0);
    stairExp(skipEyeEx(1:length(stairExp)))=[];
    skipNotStairExp=isnan(stairExp);
    stairExp(skipNotStairExp(1:length(stairExp)))=[];
    trialCount=length(stairExp);
    x=[1:1:trialCount];
    stairUn=file.d.stairIdxUn;
    skipEyeUn=(file.d.stairIdxUn==0);
    stairUn(skipEyeUn(1:length(stairUn)))=[];
    skipNotStairUn=isnan(stairUn);
    stairUn(skipNotStairUn)=[];
    trialCount2=length(stairUn);
    x2=[1:1:trialCount2];
    plot(x,stairExp,'o-');
    hold on
    plot(x2,stairUn,'o-');
    ylabel("stair ID");
    xlabel("trial ID");
    legend("expected","unexpected");
    saveas(figure(f),[plotdir_sub,'/staircasing_split.jpg']);
    f=f+1;
else
    for i=1:2
        figure();
        stair=file.d.stairIdx(i*768-768+1:i*768);
        skipNotStair=isnan(file.d.stairIdx(i*768-768+1:i*768));
        stair(skipNotStair)=[];
        trialCount=length(stair);
        x=1:1:trialCount;
        plot(x,stair,'-o');
        ylabel("stair ID");
        xlabel("trial ID");
        saveas(figure(f),[plotdir_sub,sprintf('/staircasing_session%d.jpg',i)]);
        f=f+1;
        
        figure();
        stairExp=file.d.stairIdxExp(i*768-768+1:i*768);
        skipNotStairExp=isnan(file.d.stairIdxExp(i*768-768+1:i*768));
        stairExp(skipNotStairExp)=[];
        trialCount=length(stairExp);
        x=1:1:trialCount;
        stairUn=file.d.stairIdxUn(i*768-768+1:i*768);
        skipNotStairUn=isnan(file.d.stairIdxUn(i*768-768+1:i*768));
        stairUn(skipNotStairUn)=[];
        trialCount2=length(stairUn);
        x2=1:1:trialCount2;
        %plot(x,stairExp,'o-blue','MarkerSize', 8,'LineWidth', 1);
        plot(x,stairExp,'-*','MarkerIndices',1:10:length(stairExp))
        hold on
        %plot(x2,stairUn,'o-red','MarkerSize', 8,'LineWidth', 1);
        plot(x2,stairUn,'-*','MarkerIndices',1:10:length(stairUn))
        ylabel("stair ID");
        xlabel("trial ID");
        legend("expected","unexpected");
        saveas(figure(f),[plotdir_sub,sprintf('/staircasing_split_session%d.jpg',i)]);
        f=f+1;
    end
end
%% Kok Relevant Plots: Mean Tilt + Percent Correct
gratingContrastsIdx=file.d.gratingContrast;
skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

precueValidities=file.d.precueValidity;
precueValidities(skipWaffle)=[];

response= file.d.targetResponseKey;
response(skipWaffle)=[];

diffMult=file.d.differenceMultiplier;
diffMult(skipWaffle(1:length(diffMult)))=[];

static=file.d.staticGrating;
static(skipWaffle(1:length(static)))=[];

stair=file.d.stair;
stair(stair==0)=NaN;
stair(skipWaffle(1:length(stair)))=[];
stair(isnan(stair))=[];

correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
correct_exp=correct_r==1 & precueValidities==1;
correct_unexp= correct_r==1 & precueValidities==2;

orientationJudgement= cat(3,precueValidities,correct_r,response);
nTrialsGrating=length(precueValidities);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,3]);
%sort based on stair levels
orientationJudgement=sortrows(orientationJudgement,1); 

sum(precueValidities==2)
sum(correct_unexp)
mean_unexp=sum(correct_unexp)/sum(precueValidities==2)
sum(precueValidities==1)
sum(correct_exp)
mean_exp=sum(correct_exp)/sum(precueValidities==1)

figure()
x = ["Expected" "Unexpected"];
bar(x,[mean_exp mean_unexp],'FaceColor',[0 .5 .5])
title('Percent Correct from Kok Task')
ylabel('Percent Correct')
saveas(figure(f),[plotdir_sub,'/percentcorrect_split.jpg']);

f=f+1;

meantilt_exp=sum(orientationJudgement(orientationJudgement(:,1)==1))/sum(orientationJudgement(:,1)==1);
meantilt_unexp=sum(orientationJudgement(orientationJudgement(:,1)==2))/sum(orientationJudgement(:,1)==2);

figure()
bar(x,[meantilt_exp meantilt_unexp],'FaceColor',[0 .5 .7])
title('Mean Tilt from Kok Task')
ylabel('Mean Tilt')
saveas(figure(f),[plotdir_sub,'/meantilt_split.jpg']);


%% Waffle Relevant Plots: P(test stronger) 

f=f+1;

waffleContrastsIdx=file.d.plaidContrast;
skipZero=file.d.stopThisTrial==1;
waffleContrastsIdx(skipZero)=[];
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];

waffleContrasts=[];
for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
end

precueValidities=file.d.precueValidity;
precueValidities(skipZero)=[];
precueValidities(skipNotWaffle)=[];

waffleOrientation=file.d.plaidOrientation;
waffleOrientation(skipZero)=[];
waffleOrientation(skipNotWaffle)=[];


responses= file.d.targetResponseKey;
responses(skipZero)=[];
responses(skipNotWaffle)=[];

testStronger=[];
for i=1:length(responses)
    % if they respond with -45 and orientation is -45, then they chose the
    % test as stronger
    if responses(i)==responseKeyIds(1) && waffleOrientation(i)==1 
        testStronger=[testStronger 1];

    % if they respond with +45 and orientation is +45, then they chose the
    % test as stronger
    elseif responses(i)==responseKeyIds(2) && waffleOrientation(i)==2
        testStronger=[testStronger 1];
    else
        testStronger=[testStronger 0];
    end
end
contrastJudgement= cat(4,waffleContrasts,precueValidities,waffleOrientation,testStronger);
nTrialsWaffle=length(waffleContrastsIdx);
contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,4]);
%sort based on contrast levels
contrastJudgement=sortrows(contrastJudgement,1); 
% sort based on precue validity
contrastJudgement=sortrows(contrastJudgement,2);

strongerEx=[];
strongerUnex=[];
meanStrongerPerConEx=[];
meanStrongerPerConUnex=[];
for i=1:length(responses)
    if contrastJudgement(i,2)==1%% expected
        if i==1
            strongerEx=[strongerEx contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            strongerEx=[strongerEx contrastJudgement(i,4)];
            if i==length(responses)
                strongerEx=[strongerEx contrastJudgement(i,4)];
                meanStrongerPerConEx=[meanStrongerPerConEx mean(strongerEx)];
            end

        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanStrongerPerConEx=[meanStrongerPerConEx mean(strongerEx)];
            strongerEx=[];
            strongerEx=[strongerEx contrastJudgement(i,4)];
        end
    elseif contrastJudgement(i,2)==2 
        if contrastJudgement(i-1,2)==1
            meanStrongerPerConEx=[meanStrongerPerConEx mean(strongerEx)];
            strongerUnex=[strongerUnex contrastJudgement(i,4)];

        elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
            if i==length(responses)
                strongerUnex=[strongerUnex contrastJudgement(i,4)];
                meanStrongerPerConUnex=[meanStrongerPerConUnex mean(strongerUnex)];
            end
            strongerUnex=[strongerUnex contrastJudgement(i,4)];


        elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
            meanStrongerPerConUnex=[meanStrongerPerConUnex mean(strongerUnex)];
            strongerUnex=[];
            strongerUnex=[contrastJudgement(i,4)];
        end
    end
end

contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);

allconsEx=unique(contrastJudgementEx);
allconsUnex=unique(contrastJudgementUnex);
figure();
plot(allconsEx,meanStrongerPerConEx,'blue-v');
hold on
plot(allconsUnex,meanStrongerPerConUnex,'red-v');
xlabel("contrasts")
ylabel("p(test stronger)");
title("p(test stronger) for expected vs unexpected ",subjectID);
legend("expected","unexpected");
saveas(figure(f),[plotdir_sub,'/pteststronger_split.jpg']);

%% P(reported +45)

f=f+1;
waffleContrastsIdx=file.d.plaidContrast;
skipZero=file.d.stopThisTrial==1;

waffleContrastsIdx(skipZero)=[];

skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];

waffleContrasts=[];
difference=[];
waffleOrientation=file.d.plaidOrientation; 
waffleOrientation(skipZero)=[];
waffleOrientation(skipNotWaffle)=[];


for i=1:length(waffleContrastsIdx)
    index=waffleContrastsIdx(i);
    waffleContrasts(i)=file.p.plaidContrasts1(index);
    if waffleContrasts(i)==file.p.plaidContrasts1(4)
            difference(i)=0;
    elseif waffleOrientation(i)==1   %if orientation is -45 then difference is positive, 
            difference(i)=file.p.plaidContrasts2-waffleContrasts(i); %so difference based on
           
    elseif waffleOrientation(i)==2 %if orientation is +45 then difference is positive
            difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
    end
end

precueValidities=file.d.precueValidity;
precueValidities(skipZero)=[];
precueValidities(skipNotWaffle)=[];

response= file.d.targetResponseKey;
response(skipZero)=[];
response(skipNotWaffle)=[];

response_45=(response==responseKeyIds(2));
response_135=(response==responseKeyIds(1));

correct_r=(response==responseKeyIds(1) & waffleOrientation==1) | (response==responseKeyIds(2) & waffleOrientation==2);

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

contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);

alldiffEx=unique(contrastJudgementEx(:,1));
alldiffUnex=unique(contrastJudgementUnex(:,1));

figure();
plot(alldiffEx,meanResponsePerConEx,'-v');
hold on
plot(alldiffUnex,meanResponsePerConUnex,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p  (reported +45 with respect to difference) for expected vs unexpected ",subjectID);
legend("expected 45","unexpected 45");
saveas(figure(f),[plotdir_sub,'/preported45_difference_split.jpg']); 

end
