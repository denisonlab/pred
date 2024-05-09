function plotsWaffle()

clear;
close all;

dataDir = sprintf('%s/data', pwd);
subjectID="S0004";
choice = 1; % a single session (1) or many (2)
sessionNum=1;
responseKeyIds=[1 2];

% Load data based on user choice
if choice == 1
    % Get session number for single file
    datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
    date="240506_0943";
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
    skipEye=file.d.stopThisTrial==1;
    stair(skipEye(1:length(stair)))=[];
    skipNotStair=isnan(stair);
    stair(skipNotStair)=[];
% split blocks 1-8 and 9-12 and 13-14
    stair1=stair(1:length(stair)/2);
    stair2=stair(length(stair)/2+1:12*length(stair)/16);
    stair3= stair(12*length(stair)/16 +1 :length(stair));
    
    count1=length(stair1);
    count2=length(stair2);
    count3=length(stair3);
    
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
skipEye=(file.d.stopThisTrial==1);
gratingContrastsIdx(skipEye)=[];

skipWaffle=isnan(gratingContrastsIdx);
gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];

precueValidities=file.d.precueValidity;
precueValidities(skipEye)=[];

precueValidities(skipWaffle)=[];

response= file.d.targetResponseKey;
response(skipEye)=[];

response(skipWaffle)=[];

diffMult=file.d.differenceMultiplier;
diffMult(skipEye)=[];

diffMult(skipWaffle(1:length(diffMult)))=[];

static=file.d.staticGrating;
static(skipEye)=[];

static(skipWaffle(1:length(static)))=[];

stair=file.d.stair;
stair(skipEye)=[];
stair(skipWaffle(1:length(stair)))=[];
stair(isnan(stair))=[];

correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
correct_exp=correct_r==1 & precueValidities==1;
correct_unexp= correct_r==1 & precueValidities==2;

orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
nTrialsGrating=length(precueValidities);
orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
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

meantilt_exp=sum(orientationJudgement(orientationJudgement(:,2)==1))/sum(orientationJudgement(:,2)==1);
meantilt_unexp=sum(orientationJudgement(orientationJudgement(:,2)==2))/sum(orientationJudgement(:,2)==2);

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

testStronger=file.d.testStronger;
testStronger(skipZero)=[];
testStronger(skipNotWaffle)=[];

% for i=1:length(responses)
%     % if they respond with -45 and orientation is -45, then they chose the
%     % test as stronger
%     if responses(i)==responseKeyIds(1) && waffleOrientation(i)==1 
%         testStronger=[testStronger 1];
% 
%     % if they respond with +45 and orientation is +45, then they chose the
%     % test as stronger
%     elseif responses(i)==responseKeyIds(2) && waffleOrientation(i)==2
%         testStronger=[testStronger 1];
%     else
%         testStronger=[testStronger 0];
%     end
% end
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
allconsEx=[1:5];
allconsUnex=[1:5];
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


%% plots split by blocks 1-8, 9-12, and 13-14
f=f+1;
if choice==1 
    figure();
    stair=file.d.stairIdx;
    stair(file.d.stopThisTrial==1)=[];
    % split blocks 1-8 and 9-12 and 13-14
    stair1=stair(1:length(stair)/2);
    stair2=stair(length(stair)/2+1:12*length(stair)/16);
    stair3= stair(12*length(stair)/16 +1 :length(stair));
    
    skipNotStair1=isnan(stair1);
    stair1(skipNotStair1(1:length(stair1)))=[];
    
    skipNotStair2=isnan(stair2);
    stair2(skipNotStair2(1:length(stair2)))=[];

    skipNotStair3=isnan(stair3);
    stair3(skipNotStair3(1:length(stair3)))=[];

    count1=[1:1:length(stair1)];
    count2=[1:1:length(stair2)];
    count3=[1:1:length(stair3)];

    subplot(3,1,1);
    plot(count1,stair1);
    ylabel("stair ID");
    xlabel("trial ID");

    subplot(3,1,2);
    plot(count2,stair2);
    ylabel("stair ID");
    xlabel("trial ID");
    subplot(3,1,3);
    plot(count3,stair3);
    ylabel("stair ID");
    xlabel("trial ID");

    saveas(figure(f),[plotdir_sub,'/staircasing_blocks.jpg']);
    f=f+1;

        stairExp=file.d.stairIdxExp;
    skipEyeEx=(file.d.stopThisTrial==1);

    stairExp(skipEyeEx(1:length(stairExp)))=[];
    stairExp1=stairExp(1:length(stairExp)/2);
    stairExp2=stairExp(length(stairExp)/2+1:12*length(stairExp)/16);
    stairExp3= stairExp(12*length(stairExp)/16 +1 :length(stairExp));

    skipNotStairExp1=isnan(stairExp1);
    stairExp1(skipNotStairExp1(1:length(stairExp1)))=[];
    
    skipNotStairExp2=isnan(stairExp2);
    stairExp2(skipNotStairExp2(1:length(stairExp2)))=[];

    skipNotStairExp3=isnan(stairExp3);
    stairExp3(skipNotStairExp3(1:length(stairExp3)))=[];
    
    count1=[1:1:length(stairExp1)];
    count2=[1:1:length(stairExp2)];
    count3=[1:1:length(stairExp3)];

    stairUn=file.d.stairIdxUn;
    skipEyeUn=(file.d.stopThisTrial==1);
    stairUn(skipEyeUn(1:length(stairUn)))=[];

    stairUn1=stairUn(1:length(stairUn)/2);
    stairUn2=stairUn(length(stairUn)/2+1:12*length(stairUn)/16);
    stairUn3= stairUn(12*length(stairUn)/16 +1 :length(stairUn));

    skipNotStairUn1=isnan(stairUn1);
    stairUn1(skipNotStairUn1(1:length(stairUn1)))=[];
    
    skipNotStairUn2=isnan(stairUn2);
    stairUn2(skipNotStairUn2(1:length(stairUn2)))=[];

    skipNotStairUn3=isnan(stairUn3);
    stairUn3(skipNotStairUn3(1:length(stairUn3)))=[];

    count1u=[1:1:length(stairUn1)];
    count2u=[1:1:length(stairUn2)];
    count3u=[1:1:length(stairUn3)];

    figure();
    plot(count1,stairExp1,'o-');
    hold on
    plot(count1u,stairUn1,'o-');
    
    ylabel("stair ID");
    xlabel("trial ID");
    title("blocks 1-8");
    legend("expected","unexpected");
    saveas(figure(f),[plotdir_sub,'/staircasing_split_blocks18.jpg']);

    f=f+1;

    figure();
    plot(count2,stairExp2,'o-');
    hold on
    plot(count2u,stairUn2,'o-');

    ylabel("stair ID");
    xlabel("trial ID");
    title("blocks 9-12");
    legend("expected","unexpected");
    saveas(figure(f),[plotdir_sub,'/staircasing_split_blocks912.jpg']);

    f=f+1;

    figure();
    plot(count3,stairExp3,'o-');
    hold on
    plot(count3u,stairUn3,'o-');

    ylabel("stair ID");
    xlabel("trial ID");
    title("blocks 13-16");

    legend("expected","unexpected");
    saveas(figure(f),[plotdir_sub,'/staircasing_split_blocks1316.jpg']);
    f=f+1;
end


%% Kok Relevant Plots: Mean Tilt + Percent Correct by blocks 1-8, 9-12, and 13-14
gratingContrastsIdx=file.d.gratingContrast;
gratingContrastsIdx(file.d.stopThisTrial==1)=[];
gratingContrastsIdx1=gratingContrastsIdx(1:length(gratingContrastsIdx)/2);
gratingContrastsIdx2=gratingContrastsIdx(length(gratingContrastsIdx)/2+1:12*length(gratingContrastsIdx)/16);
gratingContrastsIdx3= gratingContrastsIdx(12*length(gratingContrastsIdx)/16 +1 :length(gratingContrastsIdx));


skipWaffle1=isnan(gratingContrastsIdx1);
gratingContrastsIdx1(skipWaffle1(1:length(gratingContrastsIdx1)))=[];

skipWaffle2=isnan(gratingContrastsIdx2);
gratingContrastsIdx2(skipWaffle2(1:length(gratingContrastsIdx2)))=[];

skipWaffle3=isnan(gratingContrastsIdx3);
gratingContrastsIdx3(skipWaffle3(1:length(gratingContrastsIdx3)))=[];

precueValidities=file.d.precueValidity;
precueValidities(file.d.stopThisTrial==1)=[];

precueValidities1=precueValidities(1:length(precueValidities)/2);
precueValidities2=precueValidities(length(precueValidities)/2+1:12*length(precueValidities)/16);
precueValidities3= precueValidities(12*length(precueValidities)/16 +1 :length(precueValidities));

precueValidities1(skipWaffle1)=[];
precueValidities2(skipWaffle2)=[];
precueValidities3(skipWaffle3)=[];

response= file.d.targetResponseKey;
response(file.d.stopThisTrial==1)=[];

response1=response(1:length(response)/2);
response2=response(length(response)/2+1:12*length(response)/16);
response3= response(12*length(response)/16 +1 :length(response));

response1(skipWaffle1)=[];
response2(skipWaffle2)=[];
response3(skipWaffle3)=[];

diffMult=file.d.differenceMultiplier;
diffMult(file.d.stopThisTrial==1)=[];

diffMult1=diffMult(1:length(diffMult)/2);
diffMult2=diffMult(length(diffMult)/2+1:12*length(diffMult)/16);
diffMult3= diffMult(12*length(diffMult)/16 +1 :length(diffMult));

diffMult1(skipWaffle1(1:length(diffMult1)))=[];
diffMult2(skipWaffle2(1:length(diffMult2)))=[];
diffMult3(skipWaffle3(1:length(diffMult3)))=[];

static=file.d.staticGrating;
static(file.d.stopThisTrial==1)=[];

static1=static(1:length(static)/2);
static2=static(length(static)/2+1:12*length(static)/16);
static3= static(12*length(static)/16 +1 :length(static));
static1(skipWaffle1(1:length(static1)))=[];
static2(skipWaffle2(1:length(static2)))=[];
static3(skipWaffle3(1:length(static3)))=[];

stair=file.d.stair;
stair(file.d.stopThisTrial==1)=[];

stair1=stair(1:length(stair)/2);
stair2=stair(length(stair)/2+1:12*length(stair)/16);
stair3= stair(12*length(stair)/16 +1 :length(stair));

stair1(skipWaffle1(1:length(stair1)))=[];
stair1(isnan(stair1))=[];
stair2(skipWaffle2(1:length(stair2)))=[];
stair2(isnan(stair2))=[];
stair3(skipWaffle3(1:length(stair3)))=[];
stair3(isnan(stair3))=[];

correct_r1=(response1==responseKeyIds(1) & static1==1 & diffMult1==-1) | (response1==responseKeyIds(1) & static1==2 & diffMult1==1) | (response1==responseKeyIds(2) & static1==1 & diffMult1==1) | (response1==responseKeyIds(2) & static1==2 & diffMult1==-1);
correct_r2=(response2==responseKeyIds(1) & static2==1 & diffMult2==-1) | (response2==responseKeyIds(1) & static2==2 & diffMult2==1) | (response2==responseKeyIds(2) & static2==1 & diffMult2==1) | (response2==responseKeyIds(2) & static2==2 & diffMult2==-1);
correct_r3=(response3==responseKeyIds(1) & static3==1 & diffMult3==-1) | (response3==responseKeyIds(1) & static3==2 & diffMult3==1) | (response3==responseKeyIds(2) & static3==1 & diffMult3==1) | (response3==responseKeyIds(2) & static3==2 & diffMult3==-1);


correct_exp1=correct_r1==1 & precueValidities1==1;
correct_unexp1= correct_r1==1 & precueValidities1==2;

correct_exp2=correct_r2==1 & precueValidities2==1;
correct_unexp2= correct_r2==1 & precueValidities2==2;
correct_exp3=correct_r3==1 & precueValidities3==1;
correct_unexp3= correct_r3==1 & precueValidities3==2;

orientationJudgement1= cat(4,stair1,precueValidities1,correct_r1,response1);
nTrialsGrating1=length(precueValidities1);
orientationJudgement1=reshape(orientationJudgement1,[nTrialsGrating1,4]);
%sort based on stair levels
orientationJudgement1=sortrows(orientationJudgement1,1); 


orientationJudgement2= cat(4,stair2,precueValidities2,correct_r2,response2);
nTrialsGrating2=length(precueValidities2);
orientationJudgement2=reshape(orientationJudgement2,[nTrialsGrating2,4]);
%sort based on stair levels
orientationJudgement2=sortrows(orientationJudgement2,1); 

orientationJudgement3= cat(4,stair3,precueValidities3,correct_r3,response3);
nTrialsGrating3=length(precueValidities3);
orientationJudgement3=reshape(orientationJudgement3,[nTrialsGrating3,4]);
%sort based on stair levels
orientationJudgement3=sortrows(orientationJudgement3,1); 


sum(precueValidities1==2)
sum(correct_unexp1)
mean_unexp1=sum(correct_unexp1)/sum(precueValidities1==2)
sum(precueValidities1==1)
sum(correct_exp1)
mean_exp1=sum(correct_exp1)/sum(precueValidities1==1)

sum(precueValidities2==2)
sum(correct_unexp2)
mean_unexp2=sum(correct_unexp2)/sum(precueValidities2==2)
sum(precueValidities2==1)
sum(correct_exp2)
mean_exp2=sum(correct_exp2)/sum(precueValidities2==1)

sum(precueValidities3==2)
sum(correct_unexp3)
mean_unexp3=sum(correct_unexp3)/sum(precueValidities3==2)
sum(precueValidities3==1)
sum(correct_exp3)
mean_exp3=sum(correct_exp3)/sum(precueValidities3==1)


figure()
subplot(3,1,1)
x = ["Expected" "Unexpected"];
bar(x,[mean_exp1 mean_unexp1],'FaceColor',[0 .5 .5])
title('Percent Correct from Kok Task Blocks 1-8')
ylabel('Percent Correct')

subplot(3,1,2)
x = ["Expected" "Unexpected"];
bar(x,[mean_exp2 mean_unexp2],'FaceColor',[0 .5 .5])
title('Percent Correct from Kok Task Blocks 9-12')
ylabel('Percent Correct')

subplot(3,1,3)
x = ["Expected" "Unexpected"];
bar(x,[mean_exp3 mean_unexp3],'FaceColor',[0 .5 .5])
title('Percent Correct from Kok Task Blocks 13-16')
ylabel('Percent Correct')

saveas(figure(f),[plotdir_sub,'/percentcorrect_split_blocks.jpg']);

f=f+1;

meantilt_exp1=sum(orientationJudgement1(orientationJudgement1(:,2)==1))/sum(orientationJudgement1(:,2)==1);
meantilt_unexp1=sum(orientationJudgement1(orientationJudgement1(:,2)==2))/sum(orientationJudgement1(:,2)==2);

meantilt_exp2=sum(orientationJudgement2(orientationJudgement2(:,2)==1))/sum(orientationJudgement2(:,2)==1);
meantilt_unexp2=sum(orientationJudgement2(orientationJudgement2(:,2)==2))/sum(orientationJudgement2(:,2)==2);

meantilt_exp3=sum(orientationJudgement3(orientationJudgement3(:,2)==1))/sum(orientationJudgement3(:,2)==1);
meantilt_unexp3=sum(orientationJudgement3(orientationJudgement3(:,2)==2))/sum(orientationJudgement3(:,2)==2);


figure()
subplot(3,1,1)
bar(x,[meantilt_exp1 meantilt_unexp1],'FaceColor',[0 .5 .7])
title('Mean Tilt from Kok Task Blocks 1-8')
ylabel('Mean Tilt')
subplot(3,1,2)
bar(x,[meantilt_exp2 meantilt_unexp2],'FaceColor',[0 .5 .7])
title('Mean Tilt from Kok Task Blocks 9-12')
ylabel('Mean Tilt')
subplot(3,1,3)
bar(x,[meantilt_exp3 meantilt_unexp3],'FaceColor',[0 .5 .7])
title('Mean Tilt from Kok Task Blocks 13-16')
ylabel('Mean Tilt')

saveas(figure(f),[plotdir_sub,'/meantilt_split_blocks.jpg']);


%% Waffle Relevant Plots: P(test stronger) by blocks 1-8, 9-12, and 13-14
 
f=f+1;

waffleContrastsIdx=file.d.plaidContrast;
skipZero=file.d.stopThisTrial==1;
waffleContrastsIdx(skipZero)=[];

waffleContrastsIdx1=waffleContrastsIdx(1:length(waffleContrastsIdx)/2);
waffleContrastsIdx2=waffleContrastsIdx(length(waffleContrastsIdx)/2+1:12*length(waffleContrastsIdx)/16);
waffleContrastsIdx3= waffleContrastsIdx(12*length(waffleContrastsIdx)/16 +1 :length(waffleContrastsIdx));

skipNotWaffle1=isnan(waffleContrastsIdx1);
waffleContrastsIdx1(skipNotWaffle1)=[];

skipNotWaffle2=isnan(waffleContrastsIdx2);
waffleContrastsIdx2(skipNotWaffle2)=[];

skipNotWaffle3=isnan(waffleContrastsIdx3);
waffleContrastsIdx3(skipNotWaffle3)=[];

waffleContrasts1=[];
waffleContrasts2=[];
waffleContrasts3=[];
for i=1:length(waffleContrastsIdx1)
    index=waffleContrastsIdx1(i);
    waffleContrasts1(i)=file.p.plaidContrasts1(index);
end
for i=1:length(waffleContrastsIdx2)
    index=waffleContrastsIdx2(i);
    waffleContrasts2(i)=file.p.plaidContrasts1(index);
end
for i=1:length(waffleContrastsIdx3)
    index=waffleContrastsIdx3(i);
    waffleContrasts3(i)=file.p.plaidContrasts1(index);
end

precueValidities=file.d.precueValidity;
precueValidities(skipZero)=[];
precueValidities1=precueValidities(1:length(precueValidities)/2);
precueValidities2=precueValidities(length(precueValidities)/2+1:12*length(precueValidities)/16);
precueValidities3= precueValidities(12*length(precueValidities)/16 +1 :length(precueValidities));

precueValidities1(skipNotWaffle1)=[];
precueValidities2(skipNotWaffle2)=[];
precueValidities3(skipNotWaffle3)=[];



waffleOrientation=file.d.plaidOrientation;
waffleOrientation(skipZero)=[];

waffleOrientation1=waffleOrientation(1:length(waffleOrientation)/2);
waffleOrientation2=waffleOrientation(length(waffleOrientation)/2+1:12*length(waffleOrientation)/16);
waffleOrientation3= waffleOrientation(12*length(waffleOrientation)/16 +1 :length(waffleOrientation));

waffleOrientation1(skipNotWaffle1)=[];
waffleOrientation2(skipNotWaffle2)=[];
waffleOrientation3(skipNotWaffle3)=[];


responses= file.d.targetResponseKey;
responses(skipZero)=[];

responses1=responses(1:length(responses)/2);
responses2=responses(length(responses)/2+1:12*length(responses)/16);
responses3= responses(12*length(responses)/16 +1 :length(responses));

responses1(skipNotWaffle1)=[];
responses2(skipNotWaffle2)=[];
responses3(skipNotWaffle3)=[];

testStronger1=[];
testStronger2=[];
testStronger3=[];

for i=1:length(responses1)
    % if they respond with -45 and orientation is -45, then they chose the
    % test as stronger
    if responses1(i)==responseKeyIds(1) && waffleOrientation1(i)==1 
        testStronger1=[testStronger1 1];

    % if they respond with +45 and orientation is +45, then they chose the
    % test as stronger
    elseif responses1(i)==responseKeyIds(2) && waffleOrientation1(i)==2
        testStronger1=[testStronger1 1];
    else
        testStronger1=[testStronger1 0];
    end
end

for i=1:length(responses2)
    % if they respond with -45 and orientation is -45, then they chose the
    % test as stronger
    if responses2(i)==responseKeyIds(1) && waffleOrientation2(i)==1 
        testStronger2=[testStronger2 1];

    % if they respond with +45 and orientation is +45, then they chose the
    % test as stronger
    elseif responses2(i)==responseKeyIds(2) && waffleOrientation2(i)==2
        testStronger2=[testStronger2 1];
    else
        testStronger2=[testStronger2 0];
    end
end

for i=1:length(responses3)
    % if they respond with -45 and orientation is -45, then they chose the
    % test as stronger
    if responses3(i)==responseKeyIds(1) && waffleOrientation3(i)==1 
        testStronger3=[testStronger3 1];

    % if they respond with +45 and orientation is +45, then they chose the
    % test as stronger
    elseif responses3(i)==responseKeyIds(2) && waffleOrientation3(i)==2
        testStronger3=[testStronger3 1];
    else
        testStronger3=[testStronger3 0];
    end
end

contrastJudgement1= cat(4,waffleContrasts1,precueValidities1,waffleOrientation1,testStronger1);
nTrialsWaffle1=length(waffleContrastsIdx1);
contrastJudgement1=reshape(contrastJudgement1,[nTrialsWaffle1,4]);
%sort based on contrast levels
contrastJudgement1=sortrows(contrastJudgement1,1); 
% sort based on precue validity
contrastJudgement1=sortrows(contrastJudgement1,2);

contrastJudgement2= cat(4,waffleContrasts2,precueValidities2,waffleOrientation2,testStronger2);
nTrialsWaffle2=length(waffleContrastsIdx2);
contrastJudgement2=reshape(contrastJudgement2,[nTrialsWaffle2,4]);
%sort based on contrast levels
contrastJudgement2=sortrows(contrastJudgement2,1); 
% sort based on precue validity
contrastJudgement2=sortrows(contrastJudgement2,2);

contrastJudgement3= cat(4,waffleContrasts3,precueValidities3,waffleOrientation3,testStronger3);
nTrialsWaffle3=length(waffleContrastsIdx3);
contrastJudgement3=reshape(contrastJudgement3,[nTrialsWaffle3,4]);
%sort based on contrast levels
contrastJudgement3=sortrows(contrastJudgement3,1); 
% sort based on precue validity
contrastJudgement3=sortrows(contrastJudgement3,2);

strongerEx1=[];
strongerUnex1=[];
meanStrongerPerConEx1=[];
meanStrongerPerConUnex1=[];
for i=1:length(responses1)
    if contrastJudgement1(i,2)==1%% expected
        if i==1
            strongerEx1=[strongerEx1 contrastJudgement1(i,4)];

        elseif contrastJudgement1(i-1,1)==contrastJudgement1(i,1) 
            strongerEx1=[strongerEx1 contrastJudgement1(i,4)];
            if i==length(responses1)
                strongerEx1=[strongerEx1 contrastJudgement1(i,4)];
                meanStrongerPerConEx1=[meanStrongerPerConEx1 mean(strongerEx1)];
            end

        elseif  contrastJudgement1(i-1,1)~=contrastJudgement1(i,1)
            meanStrongerPerConEx1=[meanStrongerPerConEx1 mean(strongerEx1)];
            strongerEx1=[];
            strongerEx1=[strongerEx1 contrastJudgement1(i,4)];
        end
    elseif contrastJudgement1(i,2)==2 
        if contrastJudgement1(i-1,2)==1
            meanStrongerPerConEx1=[meanStrongerPerConEx1 mean(strongerEx1)];
            strongerUnex1=[strongerUnex1 contrastJudgement1(i,4)];

        elseif contrastJudgement1(i-1,1)==contrastJudgement1(i,1) 
            if i==length(responses1)
                strongerUnex1=[strongerUnex1 contrastJudgement1(i,4)];
                meanStrongerPerConUnex1=[meanStrongerPerConUnex1 mean(strongerUnex1)];
            end
            strongerUnex1=[strongerUnex1 contrastJudgement1(i,4)];


        elseif  contrastJudgement1(i-1,1)~=contrastJudgement1(i,1)
            meanStrongerPerConUnex1=[meanStrongerPerConUnex1 mean(strongerUnex1)];
            strongerUnex1=[];
            strongerUnex1=[contrastJudgement1(i,4)];
        end
    end
end

contrastJudgementEx1=contrastJudgement1(contrastJudgement1(:,2)==1);
contrastJudgementUnex1=contrastJudgement1(contrastJudgement1(:,2)==2);

strongerEx2=[];
strongerUnex2=[];
meanStrongerPerConEx2=[];
meanStrongerPerConUnex2=[];
for i=1:length(responses2)
    if contrastJudgement2(i,2)==1%% expected
        if i==1
            strongerEx2=[strongerEx2 contrastJudgement2(i,4)];

        elseif contrastJudgement2(i-1,1)==contrastJudgement2(i,1) 
            strongerEx2=[strongerEx2 contrastJudgement2(i,4)];
            if i==length(responses1)
                strongerEx1=[strongerEx2 contrastJudgement2(i,4)];
                meanStrongerPerConEx2=[meanStrongerPerConEx2 mean(strongerEx2)];
            end

        elseif  contrastJudgement2(i-1,1)~=contrastJudgement2(i,1)
            meanStrongerPerConEx2=[meanStrongerPerConEx2 mean(strongerEx2)];
            strongerEx2=[];
            strongerEx2=[strongerEx2 contrastJudgement2(i,4)];
        end
    elseif contrastJudgement2(i,2)==2 
        if contrastJudgement2(i-1,2)==1
            meanStrongerPerConEx2=[meanStrongerPerConEx2 mean(strongerEx2)];
            strongerUnex2=[strongerUnex2 contrastJudgement2(i,4)];

        elseif contrastJudgement2(i-1,1)==contrastJudgement2(i,1) 
            if i==length(responses2)
                strongerUnex2=[strongerUnex2 contrastJudgement2(i,4)];
                meanStrongerPerConUnex2=[meanStrongerPerConUnex2 mean(strongerUnex2)];
            end
            strongerUnex2=[strongerUnex2 contrastJudgement2(i,4)];


        elseif  contrastJudgement2(i-1,1)~=contrastJudgement2(i,1)
            meanStrongerPerConUnex2=[meanStrongerPerConUnex2 mean(strongerUnex2)];
            strongerUnex2=[];
            strongerUnex2=[contrastJudgement2(i,4)];
        end
    end
end

contrastJudgementEx2=contrastJudgement2(contrastJudgement2(:,2)==1);
contrastJudgementUnex2=contrastJudgement2(contrastJudgement2(:,2)==2);


strongerEx3=[];
strongerUnex3=[];
meanStrongerPerConEx3=[];
meanStrongerPerConUnex3=[];
for i=1:length(responses3)
    if contrastJudgement3(i,2)==1%% expected
        if i==1
            strongerEx3=[strongerEx3 contrastJudgement3(i,4)];

        elseif contrastJudgement3(i-1,1)==contrastJudgement3(i,1) 
            strongerEx3=[strongerEx3 contrastJudgement3(i,4)];
            if i==length(responses1)
                strongerEx3=[strongerEx3 contrastJudgement3(i,4)];
                meanStrongerPerConEx3=[meanStrongerPerConEx3 mean(strongerEx3)];
            end

        elseif  contrastJudgement3(i-1,1)~=contrastJudgement3(i,1)
            meanStrongerPerConEx3=[meanStrongerPerConEx3 mean(strongerEx3)];
            strongerEx3=[];
            strongerEx3=[strongerEx3 contrastJudgement3(i,4)];
        end
    elseif contrastJudgement3(i,2)==2 
        if contrastJudgement3(i-1,2)==1
            meanStrongerPerConEx3=[meanStrongerPerConEx3 mean(strongerEx3)];
            strongerUnex3=[strongerUnex3 contrastJudgement3(i,4)];

        elseif contrastJudgement3(i-1,1)==contrastJudgement3(i,1) 
            if i==length(responses2)
                strongerUnex3=[strongerUnex3 contrastJudgement3(i,4)];
                meanStrongerPerConUnex3=[meanStrongerPerConUnex3 mean(strongerUnex3)];
            end
            strongerUnex3=[strongerUnex3 contrastJudgement3(i,4)];


        elseif  contrastJudgement3(i-1,1)~=contrastJudgement3(i,1)
            meanStrongerPerConUnex3=[meanStrongerPerConUnex3 mean(strongerUnex3)];
            strongerUnex3=[];
            strongerUnex3=[contrastJudgement3(i,4)];
        end
    end
end

contrastJudgementEx3=contrastJudgement3(contrastJudgement3(:,2)==1);
contrastJudgementUnex3=contrastJudgement3(contrastJudgement3(:,2)==2);

allconsEx1=unique(contrastJudgementEx1);
allconsUnex1=unique(contrastJudgementUnex1);
allconsEx2=unique(contrastJudgementEx2);
allconsUnex2=unique(contrastJudgementUnex2);
allconsEx3=unique(contrastJudgementEx3);
allconsUnex3=unique(contrastJudgementUnex3);
figure();
plot(allconsEx1,meanStrongerPerConEx1,'blue-v');
hold on
plot(allconsUnex1,meanStrongerPerConUnex1,'red-v');
xlabel("contrasts")
ylabel("p(test stronger)");
title("p(test stronger) for expected vs unexpected, blocks 1-8",subjectID);
saveas(figure(f),[plotdir_sub,'/pteststronger_split_blocks18.jpg']);

f=f+1;
figure();
plot(allconsEx2,meanStrongerPerConEx2,'blue-v');
hold on
plot(allconsUnex2,meanStrongerPerConUnex2,'red-v');
xlabel("contrasts")
ylabel("p(test stronger)");
title("p(test stronger) for expected vs unexpected, blocks 9-12",subjectID);
saveas(figure(f),[plotdir_sub,'/pteststronger_split_blocks912.jpg']);

f=f+1;
figure();
plot(allconsEx3,meanStrongerPerConEx3,'blue-v');
hold on
plot(allconsUnex3,meanStrongerPerConUnex3,'red-v');
xlabel("contrasts")
ylabel("p(test stronger)");
title("p(test stronger) for expected vs unexpected, blocks 13-16 ",subjectID);
legend("expected","unexpected");
saveas(figure(f),[plotdir_sub,'/pteststronger_split_blocks1316.jpg']);


%% P(reported +45) for blocks 1-8,9-12,13-14

f=f+1;
waffleContrastsIdx=file.d.plaidContrast;
skipZero=file.d.stopThisTrial==1;

waffleContrastsIdx(skipZero)=[];
waffleContrastsIdx1=waffleContrastsIdx(1:length(waffleContrastsIdx)/2);
waffleContrastsIdx2=waffleContrastsIdx(length(waffleContrastsIdx)/2+1:12*length(waffleContrastsIdx)/16);
waffleContrastsIdx3= waffleContrastsIdx(12*length(waffleContrastsIdx)/16 +1 :length(waffleContrastsIdx));

skipNotWaffle1=isnan(waffleContrastsIdx1);
skipNotWaffle2=isnan(waffleContrastsIdx2);
skipNotWaffle3=isnan(waffleContrastsIdx3);

waffleContrastsIdx1(skipNotWaffle1)=[];
waffleContrastsIdx2(skipNotWaffle2)=[];
waffleContrastsIdx3(skipNotWaffle3)=[];

waffleContrasts1=[];
difference1=[];
waffleContrasts2=[];
difference2=[];
waffleContrasts3=[];
difference3=[];

waffleOrientation=file.d.plaidOrientation; 
waffleOrientation(skipZero)=[];
waffleOrientation1=waffleOrientation(1:length(waffleOrientation)/2);
waffleOrientation2=waffleOrientation(length(waffleOrientation)/2+1:12*length(waffleOrientation)/16);
waffleOrientation3= waffleOrientation(12*length(waffleOrientation)/16 +1 :length(waffleOrientation));

waffleOrientation1(skipNotWaffle1)=[];
waffleOrientation2(skipNotWaffle2)=[];
waffleOrientation3(skipNotWaffle3)=[];


for i=1:length(waffleContrastsIdx1)
    index=waffleContrastsIdx1(i);
    waffleContrasts1(i)=file.p.plaidContrasts1(index);
    if waffleContrasts1(i)==file.p.plaidContrasts1(4)
            difference1(i)=0;
    elseif waffleOrientation1(i)==1   %if orientation is -45 then difference is positive, 
            difference1(i)=file.p.plaidContrasts2-waffleContrasts1(i); %so difference based on
           
    elseif waffleOrientation1(i)==2 %if orientation is +45 then difference is positive
            difference1(i)=waffleContrasts1(i)-file.p.plaidContrasts2;
    end
end

for i=1:length(waffleContrastsIdx2)
    index=waffleContrastsIdx2(i);
    waffleContrasts2(i)=file.p.plaidContrasts1(index);
    if waffleContrasts2(i)==file.p.plaidContrasts1(4)
            difference2(i)=0;
    elseif waffleOrientation2(i)==1   %if orientation is -45 then difference is positive, 
            difference2(i)=file.p.plaidContrasts2-waffleContrasts2(i); %so difference based on
           
    elseif waffleOrientation2(i)==2 %if orientation is +45 then difference is positive
            difference2(i)=waffleContrasts2(i)-file.p.plaidContrasts2;
    end
end

for i=1:length(waffleContrastsIdx3)
    index=waffleContrastsIdx3(i);
    waffleContrasts3(i)=file.p.plaidContrasts1(index);
    if waffleContrasts3(i)==file.p.plaidContrasts1(4)
            difference3(i)=0;
    elseif waffleOrientation3(i)==1   %if orientation is -45 then difference is positive, 
            difference3(i)=file.p.plaidContrasts2-waffleContrasts3(i); %so difference based on
           
    elseif waffleOrientation3(i)==2 %if orientation is +45 then difference is positive
            difference3(i)=waffleContrasts3(i)-file.p.plaidContrasts2;
    end
end


precueValidities=file.d.precueValidity;
precueValidities(skipZero)=[];

precueValidities1=precueValidities(1:length(precueValidities)/2);
precueValidities2=precueValidities(length(precueValidities)/2+1:12*length(precueValidities)/16);
precueValidities3= precueValidities(12*length(precueValidities)/16 +1 :length(precueValidities));

precueValidities1(skipNotWaffle1)=[];
precueValidities2(skipNotWaffle2)=[];
precueValidities3(skipNotWaffle3)=[];

responses= file.d.targetResponseKey;
responses(skipZero)=[];
responses1=responses(1:length(responses)/2);
responses2=responses(length(responses)/2+1:12*length(responses)/16);
responses3= responses(12*length(responses)/16 +1 :length(responses));

responses1(skipNotWaffle1)=[];
responses2(skipNotWaffle2)=[];
responses3(skipNotWaffle3)=[];


response1_45=(responses1==responseKeyIds(2));
response1_135=(responses1==responseKeyIds(1));
response2_45=(responses2==responseKeyIds(2));
response2_135=(responses2==responseKeyIds(1));
response3_45=(responses3==responseKeyIds(2));
response3_135=(responses3==responseKeyIds(1));

correct_r1=(responses1==responseKeyIds(1) & waffleOrientation1==1) | (responses1==responseKeyIds(2) & waffleOrientation1==2);
correct_r2=(responses2==responseKeyIds(1) & waffleOrientation2==1) | (responses2==responseKeyIds(2) & waffleOrientation2==2);
correct_r3=(responses3==responseKeyIds(1) & waffleOrientation3==1) | (responses3==responseKeyIds(2) & waffleOrientation3==2);

contrastJudgement1= cat(7,difference1,precueValidities1,correct_r1,response1_45,response1_135,responses1,waffleContrasts1);
nTrialsWaffle1=length(waffleContrastsIdx1);
contrastJudgement1=reshape(contrastJudgement1,[nTrialsWaffle1,7]);
%sort based on difference values
contrastJudgement1=sortrows(contrastJudgement1,1); 
% sort based on precue validity
contrastJudgement1=sortrows(contrastJudgement1,2);

responseEx1=[];
responseUnex1=[];
meanResponsePerConEx1=[];
meanResponsePerConUnex1=[];
for i=1:length(responses1)
    if contrastJudgement1(i,2)==1%% expected
        if i==1
            responseEx1=[responseEx1 contrastJudgement1(i,4)];

        elseif contrastJudgement1(i-1,1)==contrastJudgement1(i,1) 
            responseEx1=[responseEx1 contrastJudgement1(i,4)];
            if i==length(responses1)
                responseEx1=[responseEx1 contrastJudgement1(i,4)];
                meanResponsePerConEx1=[meanResponsePerConEx1 mean(responseEx1)];
            end

        elseif  contrastJudgement1(i-1,1)~=contrastJudgement1(i,1)
            meanResponsePerConEx1=[meanResponsePerConEx1 mean(responseEx1)];
            responseEx1=[];
            responseEx1=[responseEx1 contrastJudgement1(i,4)];
            if i==length(responses1)
               meanResponsePerConEx1=[meanResponsePerConEx1 mean(responseEx1)];
            end
        end
    elseif contrastJudgement1(i,2)==2 
        if contrastJudgement1(i-1,2)==1
            meanResponsePerConEx1=[meanResponsePerConEx1 mean(responseEx1)];
            responseUnex1=[responseUnex1 contrastJudgement1(i,4)];

        elseif contrastJudgement1(i-1,1)==contrastJudgement1(i,1) 
            if i==length(responses1)
                responseUnex1=[responseUnex1 contrastJudgement1(i,4)];
                meanResponsePerConUnex1=[meanResponsePerConUnex1 mean(responseUnex1)];
            end
            responseUnex1=[responseUnex1 contrastJudgement1(i,4)];


        elseif  contrastJudgement1(i-1,1)~=contrastJudgement1(i,1)
            meanResponsePerConUnex1=[meanResponsePerConUnex1 mean(responseUnex1)];
            responseUnex1=[];
            responseUnex1=[contrastJudgement1(i,4)];
            if i==length(responses1)
               meanResponsePerConUnex1=[meanResponsePerConUnex1 mean(responseUnex1)];
            end
        end
    end
end


contrastJudgement2= cat(7,difference2,precueValidities2,correct_r2,response2_45,response2_135,responses2,waffleContrasts2);
nTrialsWaffle2=length(waffleContrastsIdx2);
contrastJudgement2=reshape(contrastJudgement2,[nTrialsWaffle2,7]);
%sort based on difference values
contrastJudgement2=sortrows(contrastJudgement2,1); 
% sort based on precue validity
contrastJudgement2=sortrows(contrastJudgement2,2);

responseEx2=[];
responseUnex2=[];
meanResponsePerConEx2=[];
meanResponsePerConUnex2=[];
for i=1:length(responses2)
    if contrastJudgement2(i,2)==1%% expected
        if i==1
            responseEx2=[responseEx2 contrastJudgement2(i,4)];

        elseif contrastJudgement2(i-1,1)==contrastJudgement2(i,1) 
            responseEx2=[responseEx2 contrastJudgement2(i,4)];
            if i==length(responses2)
                responseEx2=[responseEx2 contrastJudgement2(i,4)];
                meanResponsePerConEx2=[meanResponsePerConEx2 mean(responseEx2)];
            end

        elseif  contrastJudgement2(i-1,1)~=contrastJudgement2(i,1)
            meanResponsePerConEx2=[meanResponsePerConEx2 mean(responseEx2)];
            responseEx2=[];
            responseEx2=[responseEx2 contrastJudgement2(i,4)];
            if i==length(responses2)
               meanResponsePerConEx2=[meanResponsePerConEx2 mean(responseEx2)];
            end
        end
    elseif contrastJudgement2(i,2)==2 
        if contrastJudgement2(i-1,2)==1
            meanResponsePerConEx2=[meanResponsePerConEx2 mean(responseEx2)];
            responseUnex2=[responseUnex2 contrastJudgement2(i,4)];

        elseif contrastJudgement2(i-1,1)==contrastJudgement2(i,1) 
            if i==length(responses2)
                responseUnex2=[responseUnex2 contrastJudgement2(i,4)];
                meanResponsePerConUnex2=[meanResponsePerConUnex2 mean(responseUnex2)];
            end
            responseUnex2=[responseUnex2 contrastJudgement2(i,4)];


        elseif  contrastJudgement2(i-1,1)~=contrastJudgement2(i,1)
            meanResponsePerConUnex2=[meanResponsePerConUnex2 mean(responseUnex2)];
            responseUnex2=[];
            responseUnex2=[contrastJudgement2(i,4)];
            if i==length(responses2)
               meanResponsePerConUnex2=[meanResponsePerConUnex2 mean(responseUnex2)];
            end
        end
    end
end

contrastJudgement3= cat(7,difference3,precueValidities3,correct_r3,response3_45,response3_135,responses3,waffleContrasts3);
nTrialsWaffle3=length(waffleContrastsIdx3);
contrastJudgement3=reshape(contrastJudgement3,[nTrialsWaffle3,7]);
%sort based on difference values
contrastJudgement3=sortrows(contrastJudgement3,1); 
% sort based on precue validity
contrastJudgement3=sortrows(contrastJudgement3,2);

responseEx3=[];
responseUnex3=[];
meanResponsePerConEx3=[];
meanResponsePerConUnex3=[];
for i=1:length(responses3)
    if contrastJudgement3(i,2)==1%% expected
        if i==1
            responseEx3=[responseEx3 contrastJudgement3(i,4)];

        elseif contrastJudgement3(i-1,1)==contrastJudgement3(i,1) 
            responseEx3=[responseEx3 contrastJudgement3(i,4)];
            if i==length(responses3)
                responseEx3=[responseEx3 contrastJudgement3(i,4)];
                meanResponsePerConEx3=[meanResponsePerConEx3 mean(responseEx3)];
            end

        elseif  contrastJudgement3(i-1,1)~=contrastJudgement3(i,1)
            meanResponsePerConEx3=[meanResponsePerConEx3 mean(responseEx3)];
            responseEx3=[];
            responseEx3=[responseEx3 contrastJudgement3(i,4)];
            if i==length(responses3)
               meanResponsePerConEx3=[meanResponsePerConEx3 mean(responseEx3)];
            end
        end
    elseif contrastJudgement3(i,2)==2 
        if contrastJudgement3(i-1,2)==1
            meanResponsePerConEx3=[meanResponsePerConEx3 mean(responseEx3)];
            responseUnex3=[responseUnex3 contrastJudgement3(i,4)];

        elseif contrastJudgement3(i-1,1)==contrastJudgement3(i,1) 
            if i==length(responses3)
                responseUnex3=[responseUnex3 contrastJudgement3(i,4)];
                meanResponsePerConUnex3=[meanResponsePerConUnex3 mean(responseUnex3)];
            end
            responseUnex3=[responseUnex3 contrastJudgement3(i,4)];


        elseif  contrastJudgement3(i-1,1)~=contrastJudgement3(i,1)
            meanResponsePerConUnex3=[meanResponsePerConUnex3 mean(responseUnex3)];
            responseUnex3=[];
            responseUnex3=[contrastJudgement3(i,4)];
            if i==length(responses3)
               meanResponsePerConUnex3=[meanResponsePerConUnex3 mean(responseUnex3)];
            end
        end
    end
end

contrastJudgementEx1=contrastJudgement1(contrastJudgement1(:,2)==1);
contrastJudgementUnex1=contrastJudgement1(contrastJudgement1(:,2)==2);
contrastJudgementEx2=contrastJudgement2(contrastJudgement2(:,2)==1);
contrastJudgementUnex2=contrastJudgement2(contrastJudgement2(:,2)==2);
contrastJudgementEx3=contrastJudgement3(contrastJudgement3(:,2)==1);
contrastJudgementUnex3=contrastJudgement3(contrastJudgement3(:,2)==2);

alldiffEx1=unique(contrastJudgementEx1(:,1));
alldiffUnex1=unique(contrastJudgementUnex1(:,1));

alldiffEx2=unique(contrastJudgementEx2(:,1));
alldiffUnex2=unique(contrastJudgementUnex2(:,1));

alldiffEx3=unique(contrastJudgementEx3(:,1));
alldiffUnex3=unique(contrastJudgementUnex3(:,1));

figure();
plot(alldiffEx1,meanResponsePerConEx1,'-v');
hold on
plot(alldiffUnex1,meanResponsePerConUnex1,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p  (reported +45 with respect to difference) for expected vs unexpected, blocks 1-8 ",subjectID);
legend("expected 45","unexpected 45");
saveas(figure(f),[plotdir_sub,'/preported45_difference_split_blocks18.jpg']); 

f=f+1;
figure();
plot(alldiffEx2,meanResponsePerConEx2,'-v');
hold on
plot(alldiffUnex2,meanResponsePerConUnex2,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p  (reported +45 with respect to difference) for expected vs unexpected, blocks 9-12 ",subjectID);
legend("expected 45","unexpected 45");
saveas(figure(f),[plotdir_sub,'/preported45_difference_split_blocks912.jpg']); 

f=f+1;
figure();
plot(alldiffEx3,meanResponsePerConEx3,'-v');
hold on
plot(alldiffUnex3,meanResponsePerConUnex3,'--v');
xlabel("difference value (with respect to +45)")
ylabel("p (reported +45 with respect to difference)");
title("p  (reported +45 with respect to difference) for expected vs unexpected, blocks 13-16 ",subjectID);
legend("expected 45","unexpected 45");
saveas(figure(f),[plotdir_sub,'/preported45_difference_split_blocks1316.jpg']); 

end
