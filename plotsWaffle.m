function plotsWaffle()

clear;
close all;

dataDir = sprintf('%s/data', pwd);
subjectID="newS0130";
choice = 1; % a single session (1) or many (2)
sessionNum=1;
responseKeyIds=[1 2];

% Load data based on user choice
if choice == 1
    % Get session number for single file
    datafileDir=sprintf('%s/%s/session_%d/',dataDir,subjectID,sessionNum);
    date="240502_0931";
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

correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
correct_exp=correct_r==1 & precueValidities==1;
correct_unexp= correct_r==1 & precueValidities==2;

orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
nTrialsGrating=length(stair);
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
skipNotWaffle=isnan(waffleContrastsIdx);
waffleContrastsIdx(skipNotWaffle)=[];

waffleContrasts=[];
difference=[];
waffleOrientation=file.d.plaidOrientation; 
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
precueValidities(skipNotWaffle)=[];

response= file.d.targetResponseKey;
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


%% 


% %% plot contrast judgement task from waffles
% f=f+1;
% waffleContrastsIdx=file.d.plaidContrast;
% skipNotWaffle=isnan(waffleContrastsIdx);
% waffleContrastsIdx(skipNotWaffle)=[];
% waffleContrasts=[];
% for i=1:length(waffleContrastsIdx)
%     index=waffleContrastsIdx(i);
%     waffleContrasts(i)=file.p.plaidContrasts1(index);
% end
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
% correctResponse=[];
% for i=1:length(responses)
%     if responses(i)==responseKeyIds(1) && waffleOrientation(i)==1
%         correctResponse=[correctResponse 1];
%     elseif responses(i)==responseKeyIds(2) && waffleOrientation(i)==2
%         correctResponse=[correctResponse 1];
%     else
%         correctResponse=[correctResponse 0];
%     end
% end
% contrastJudgement= cat(4,waffleContrasts,precueValidities,waffleOrientation,correctResponse);
% nTrialsWaffle=length(waffleContrastsIdx);
% contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,4]);
% %sort based on contrast levels
% contrastJudgement=sortrows(contrastJudgement,1); 
% % sort based on precue validity
% contrastJudgement=sortrows(contrastJudgement,2);
% 
% sum(precueValidities==2)
% sum(precueValidities==1)
% uniqueContrastsEx=unique(contrastJudgement(contrastJudgement(:,2)==1))
% uniqueContrastsUnex=unique(contrastJudgement(contrastJudgement(:,2)==2))
% sum_ex=[];
% sum_ex_45=[];
% sum_ex_135=[];
% for i=1:length(uniqueContrastsEx)
%     sum_ex=[sum_ex sum(contrastJudgement(:,2)==1 & contrastJudgement(:,1)==uniqueContrastsEx(i))];
%     sum_ex_45=[sum_ex_45 sum(contrastJudgement(:,2)==1 & contrastJudgement(:,3)==1 & contrastJudgement(:,1)==uniqueContrastsEx(i))];
%     sum_ex_135=[sum_ex_135 sum(contrastJudgement(:,2)==1 & contrastJudgement(:,3)==2 & contrastJudgement(:,1)==uniqueContrastsEx(i))];
% 
% end
% sum_unex=[];
% sum_unex_45=[];
% sum_unex_135=[];
% for i=1:length(uniqueContrastsUnex)
%     sum_unex=[sum_unex sum(contrastJudgement(:,2)==2 & contrastJudgement(:,1)==uniqueContrastsUnex(i))];
%     sum_unex_45=[sum_unex_45 sum(contrastJudgement(:,2)==2 & contrastJudgement(:,3)==1 & contrastJudgement(:,1)==uniqueContrastsUnex(i))];
%     sum_unex_135=[sum_unex_135 sum(contrastJudgement(:,2)==2 & contrastJudgement(:,3)==2 & contrastJudgement(:,1)==uniqueContrastsUnex(i))];
% end
% correctsEx=[];
% correctsUnex=[];
% meanCorrectPerConEx=[];
% meanCorrectPerConUnex=[];
% for i=1:length(responses)
%     if contrastJudgement(i,2)==1%% expected
%         if i==1
%             correctsEx=[correctsEx contrastJudgement(i,4)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             correctsEx=[correctsEx contrastJudgement(i,4)];
%             if i==length(responses)
%                 correctsEx=[correctsEx contrastJudgement(i,4)];
%                 meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             end
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsEx=[];
%             correctsEx=[correctsEx contrastJudgement(i,4)];
%         end
%     elseif contrastJudgement(i,2)==2 
%         if contrastJudgement(i-1,2)==1
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsUnex=[correctsUnex contrastJudgement(i,4)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             if i==length(responses)
%                 correctsUnex=[correctsUnex contrastJudgement(i,4)];
%                 meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             end
%             correctsUnex=[correctsUnex contrastJudgement(i,4)];
% 
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             correctsUnex=[];
%             correctsUnex=[contrastJudgement(i,4)];
%         end
%     end
% end
% contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
% contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);
% 
% allconsEx=unique(contrastJudgementEx);
% allconsUnex=unique(contrastJudgementUnex);
% figure();
% plot(allconsEx,meanCorrectPerConEx,'blue-v');
% hold on
% plot(allconsUnex,meanCorrectPerConUnex,'red-v');
% xlabel("contrasts")
% ylabel("p(correct stronger)");
% title("p(correct stronger) for expected vs unexpected ",subjectID);
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrectstronger_split.jpg']);
% 
% %% p (correct tilt discrimination)
% f=f+1;
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% 
% response= file.d.targetResponseKey;
% response(skipWaffle)=[];
% 
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% 
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% correct_exp=correct_r==1 & precueValidities==1;
% correct_unexp= correct_r==1 & precueValidities==2;
% 
% 
% sum(precueValidities==2)
% sum(correct_unexp)
% mean_unexp=sum(correct_unexp)/sum(precueValidities==2)
% sum(precueValidities==1)
% sum(correct_exp)
% mean_exp=sum(correct_exp)/sum(precueValidities==1)
% 
% orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
% nTrialsGrating=length(stair);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
% %sort based on stair levels
% orientationJudgement=sortrows(orientationJudgement,1); 
% 
% 
% corrects=[];
% meanCorrectPerCon=[];
% for i=1:length(response)
%     if i==1
%         corrects=[corrects orientationJudgement(i,3)];
% 
%     elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
%         corrects=[corrects orientationJudgement(i,3)];
%         if i==length(response)
%             corrects=[corrects orientationJudgement(i,3)];
%             meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
%         end
% 
%     elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
%         meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
%         corrects=[];
%         corrects=[corrects orientationJudgement(i,3)];
%         if i==length(response)
%            meanCorrectPerCon=[meanCorrectPerCon mean(corrects)];
%         end
% 
%     end
% end
% 
% allStairs=unique(orientationJudgement(:,1));
% 
% allstairsId=[];
% for i=1:length(allStairs)
%     allstairsId=[allstairsId find(file.p.stairs==allStairs(i))];
% end
% 
% 
% figure();
% plot(allStairs,meanCorrectPerCon,'blue-v');
% xlabel("stair value")
% ylabel("p (correct tilt discrimination)");
% title("p (correct tilt discrimination) ",subjectID);
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_stairs.jpg']);
% f=f+1;
% figure();
% plot(allstairsId,meanCorrectPerCon,'blue-v');
% xlabel("stair id")
% ylabel("p (correct tilt discrimination)");
% title("p (correct tilt discrimination) ",subjectID);
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_stairsID.jpg']);
% 
% mean(correct_r)
% 
% %% p (correct tilt discrimination)
% f=f+1;
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% 
% response= file.d.targetResponseKey;
% response(skipWaffle)=[];
% 
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% 
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
% nTrialsGrating=length(stair);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
% %sort based on stair levels
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% correctsEx=[];
% correctsUnex=[];
% meanCorrectPerConEx=[];
% meanCorrectPerConUnex=[];
% for i=1:length(response)
%     if orientationJudgement(i,2)==1%% expected
%         if i==1
%             correctsEx=[correctsEx orientationJudgement(i,3)];
% 
%         elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
%             correctsEx=[correctsEx orientationJudgement(i,3)];
%             if i==length(response)
%                 correctsEx=[correctsEx orientationJudgement(i,3)];
%                 meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             end
% 
%         elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsEx=[];
%             correctsEx=[correctsEx orientationJudgement(i,3)];
%             if i==length(response)
%                meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             end
%         end
%     elseif orientationJudgement(i,2)==2 
%         if orientationJudgement(i-1,2)==1
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsUnex=[correctsUnex orientationJudgement(i,3)];
% 
%         elseif orientationJudgement(i-1,1)==orientationJudgement(i,1) 
%             if i==length(response)
%                 correctsUnex=[correctsUnex orientationJudgement(i,3)];
%                 meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             end
%             correctsUnex=[correctsUnex orientationJudgement(i,3)];
% 
% 
%         elseif  orientationJudgement(i-1,1)~=orientationJudgement(i,1)
%             meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             correctsUnex=[];
%             correctsUnex=[orientationJudgement(i,3)];
%             if i==length(response)
%                meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             end
%         end
%     end
% end
% 
% orientationJudgementEx=orientationJudgement(orientationJudgement(:,2)==1);
% orientationJudgementUnex=orientationJudgement(orientationJudgement(:,2)==2);
% 
% allstairsEx=unique(orientationJudgementEx);
% allstairsUnex=unique(orientationJudgementUnex);
% allstairsIdEx=[];
% allstairsIdUnex=[];
% for i=1:length(allstairsEx)
%     allstairsIdEx=[allstairsIdEx find(file.p.stairs==allstairsEx(i))];
% end
% 
% for i=1:length(allstairsUnex)
%     allstairsIdUnex=[allstairsIdUnex find(file.p.stairs==allstairsUnex(i))];
% end
% 
% 
% figure();
% plot(allstairsEx,meanCorrectPerConEx,'blue-v');
% hold on
% plot(allstairsUnex,meanCorrectPerConUnex,'red-v');
% xlabel("stair value")
% ylabel("p (correct tilt discrimination)");
% title("p (correct tilt discrimination) for expected vs unexpected ",subjectID);
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_split_stairs.jpg']);
% f=f+1;
% figure();
% plot(allstairsIdEx,meanCorrectPerConEx,'blue-v');
% hold on
% plot(allstairsIdUnex,meanCorrectPerConUnex,'red-v');
% xlabel("stair id")
% ylabel("p (correct tilt discrimination)");
% title("p (correct tilt discrimination) for expected vs unexpected ",subjectID);
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_split_stairsID.jpg']);
% 
% 
% %% p (correct contrast discrimination with respect to difference)
% f=f+1;
% waffleContrastsIdx=file.d.plaidContrast;
% skipNotWaffle=isnan(waffleContrastsIdx);
% waffleContrastsIdx(skipNotWaffle)=[];
% 
% waffleContrasts=[];
% difference=[];
% waffleOrientation=file.d.plaidOrientation; 
% waffleOrientation(skipNotWaffle)=[];
% 
% for i=1:length(waffleContrastsIdx)
%     index=waffleContrastsIdx(i);
%     waffleContrasts(i)=file.p.plaidContrasts1(index);
%     if waffleOrientation(i)==2 %if orientation is +45 then difference is positive
%         difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
%     elseif waffleOrientation(i)==1 %if orientation is -45 then difference is positive, 
%         difference(i)=file.p.plaidContrasts2-waffleContrasts(i); %so difference based on
%     end
% end
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipNotWaffle)=[];
% 
% response= file.d.targetResponseKey;
% response(skipNotWaffle)=[];
% 
% correct_r=(response==responseKeyIds(1) & waffleOrientation==1) | (response==responseKeyIds(2) & waffleOrientation==2)
% 
% contrastJudgement= cat(4,difference,precueValidities,correct_r,response);
% nTrialsWaffle=length(waffleContrastsIdx);
% contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,4]);
% %sort based on difference values
% contrastJudgement=sortrows(contrastJudgement,1); 
% % sort based on precue validity
% contrastJudgement=sortrows(contrastJudgement,2);
% 
% correctsEx=[];
% correctsUnex=[];
% meanCorrectPerConEx=[];
% meanCorrectPerConUnex=[];
% for i=1:length(response)
%     if contrastJudgement(i,2)==1%% expected
%         if i==1
%             correctsEx=[correctsEx contrastJudgement(i,3)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             correctsEx=[correctsEx contrastJudgement(i,3)];
%             if i==length(response)
%                 correctsEx=[correctsEx contrastJudgement(i,3)];
%                 meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             end
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsEx=[];
%             correctsEx=[correctsEx contrastJudgement(i,3)];
%             if i==length(response)
%                meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             end
%         end
%     elseif contrastJudgement(i,2)==2 
%         if contrastJudgement(i-1,2)==1
%             meanCorrectPerConEx=[meanCorrectPerConEx mean(correctsEx)];
%             correctsUnex=[correctsUnex contrastJudgement(i,3)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             if i==length(response)
%                 correctsUnex=[correctsUnex contrastJudgement(i,3)];
%                 meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             end
%             correctsUnex=[correctsUnex contrastJudgement(i,3)];
% 
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             correctsUnex=[];
%             correctsUnex=[contrastJudgement(i,3)];
%             if i==length(response)
%                meanCorrectPerConUnex=[meanCorrectPerConUnex mean(correctsUnex)];
%             end
%         end
%     end
% end
% 
% contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
% contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);
% 
% alldiffEx=unique(contrastJudgementEx);
% alldiffUnex=unique(contrastJudgementUnex);
% 
% 
% figure();
% plot(alldiffEx,meanCorrectPerConEx,'blue-v');
% hold on
% plot(alldiffUnex,meanCorrectPerConUnex,'red-v');
% xlabel("difference value (with respect to +45)")
% ylabel("p (correct contrast discrimination with respect to difference)");
% title("p (correct contrast discrimination with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrectstronger_difference_split.jpg']);
% 
% %% p (responses with +45/-45 with respect to difference)
% f=f+1;
% waffleContrastsIdx=file.d.plaidContrast;
% skipNotWaffle=isnan(waffleContrastsIdx);
% waffleContrastsIdx(skipNotWaffle)=[];
% 
% waffleContrasts=[];
% difference=[];
% waffleOrientation=file.d.plaidOrientation; 
% waffleOrientation(skipNotWaffle)=[];
% 
% for i=1:length(waffleContrastsIdx)
%     index=waffleContrastsIdx(i);
%     waffleContrasts(i)=file.p.plaidContrasts1(index);
%     if waffleOrientation(i)==2 %if orientation is +45 then difference is positive
%         difference(i)=waffleContrasts(i)-file.p.plaidContrasts2;
%     elseif waffleOrientation(i)==1 %if orientation is -45 then difference is positive, 
%         difference(i)=file.p.plaidContrasts2-waffleContrasts(i); %so difference based on
%     end
% end
% 
% precueValidities=file.d.precueValidity;
% precueValidities(skipNotWaffle)=[];
% 
% response= file.d.targetResponseKey;
% response(skipNotWaffle)=[];
% 
% response_45=(response==responseKeyIds(2));
% response_135=(response==responseKeyIds(1));
% 
% correct_r=(response==responseKeyIds(1) & waffleOrientation==1) | (response==responseKeyIds(2) & waffleOrientation==2);
% 
% contrastJudgement= cat(7,difference,precueValidities,correct_r,response_45,response_135,response,waffleContrasts);
% nTrialsWaffle=length(waffleContrastsIdx);
% contrastJudgement=reshape(contrastJudgement,[nTrialsWaffle,7]);
% %sort based on difference values
% contrastJudgement=sortrows(contrastJudgement,1); 
% % sort based on precue validity
% contrastJudgement=sortrows(contrastJudgement,2);
% 
% responseEx=[];
% responseUnex=[];
% meanResponsePerConEx=[];
% meanResponsePerConUnex=[];
% for i=1:length(response)
%     if contrastJudgement(i,2)==1%% expected
%         if i==1
%             responseEx=[responseEx contrastJudgement(i,4)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             responseEx=[responseEx contrastJudgement(i,4)];
%             if i==length(response)
%                 responseEx=[responseEx contrastJudgement(i,4)];
%                 meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
%             end
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
%             responseEx=[];
%             responseEx=[responseEx contrastJudgement(i,4)];
%             if i==length(response)
%                meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
%             end
%         end
%     elseif contrastJudgement(i,2)==2 
%         if contrastJudgement(i-1,2)==1
%             meanResponsePerConEx=[meanResponsePerConEx mean(responseEx)];
%             responsesUnex=[responseUnex contrastJudgement(i,4)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             if i==length(response)
%                 responseUnex=[responseUnex contrastJudgement(i,4)];
%                 meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
%             end
%             responseUnex=[responseUnex contrastJudgement(i,4)];
% 
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
%             responseUnex=[];
%             responseUnex=[contrastJudgement(i,4)];
%             if i==length(response)
%                meanResponsePerConUnex=[meanResponsePerConUnex mean(responseUnex)];
%             end
%         end
%     end
% end
% 
% responseEx135=[];
% responsesUnex135=[];
% meanResponsePerConEx_135=[];
% meanResponsePerConUnex_135=[];
% 
% for i=1:length(response)
%     if contrastJudgement(i,2)==1%% expected
%         if i==1
%             responseEx135=[responseEx135 contrastJudgement(i,5)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             responseEx135=[responseEx135 contrastJudgement(i,5)];
%             if i==length(response)
%                 responseEx135=[responseEx135 contrastJudgement(i,5)];
%                 meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
%             end
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
%             responseEx135=[];
%             responseEx135=[responseEx135 contrastJudgement(i,5)];
%             if i==length(response)
%                meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
%             end
%         end
%     elseif contrastJudgement(i,2)==2 
%         if contrastJudgement(i-1,2)==1
%             meanResponsePerConEx_135=[meanResponsePerConEx_135 mean(responseEx135)];
%             responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];
% 
%         elseif contrastJudgement(i-1,1)==contrastJudgement(i,1) 
%             if i==length(response)
%                 responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];
%                 meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
%             end
%             responsesUnex135=[responsesUnex135 contrastJudgement(i,5)];
% 
% 
%         elseif  contrastJudgement(i-1,1)~=contrastJudgement(i,1)
%             meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
%             responsesUnex135=[];
%             responsesUnex135=[contrastJudgement(i,5)];
%             if i==length(response)
%                meanResponsePerConUnex_135=[meanResponsePerConUnex_135 mean(responsesUnex135)];
%             end
%         end
%     end
% end
% 
% contrastJudgementEx=contrastJudgement(contrastJudgement(:,2)==1);
% contrastJudgementUnex=contrastJudgement(contrastJudgement(:,2)==2);
% 
% alldiffEx=unique(contrastJudgementEx);
% alldiffUnex=unique(contrastJudgementUnex);
% 
% % x-axis is difference between grating contrasts inside plaid with respect
% % to +45
% figure();
% yyaxis left
% plot(alldiffEx,meanResponsePerConEx,'-v');
% hold on
% plot(alldiffUnex,meanResponsePerConUnex,'--v');
% ylabel("p (reported +45 with respect to difference)");
% 
% yyaxis right
% plot(alldiffEx,meanResponsePerConEx_135,'-v');
% hold on
% plot(alldiffEx,meanResponsePerConUnex_135,'--v');
% xlabel("difference value (with respect to +45)")
% ylabel("p (reported -45 with respect to difference)");
% title("p  (reported -45/+45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected 45","unexpected 45","expected -45","unexpected -45");
% saveas(figure(f),[plotdir_sub,'/preported_difference_split.jpg']);
% f=f+1;

% figure();
% plot(alldiffEx,meanResponsePerConEx,'-v');
% hold on
% plot(alldiffUnex,meanResponsePerConUnex,'--v');
% xlabel("difference value (with respect to +45)")
% ylabel("p (reported +45 with respect to difference)");
% title("p  (reported +45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected 45","unexpected 45");
% saveas(figure(f),[plotdir_sub,'/preported45_difference_split.jpg']);
% f=f+1;
% figure();
% plot(alldiffEx,meanResponsePerConEx_135,'-v');
% hold on
% plot(alldiffEx,meanResponsePerConUnex_135,'--v');
% xlabel("difference value (with respect to -45)")
% ylabel("p (reported -45 with respect to difference)");
% title("p (reported -45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected -45","unexpected -45");
% saveas(figure(f),[plotdir_sub,'/preported135_difference_split.jpg']);
% 
% % x-axis is strength of plaid 
% f=f+1;
% allStrength=[flip(file.p.plaidContrasts1) file.p.plaidContrasts1 ];
% 
% figure();
% yyaxis left
% plot(allStrength,meanResponsePerConEx,'-v');
% hold on
% plot(allStrength,meanResponsePerConUnex,'--v');
% ylabel("p (reported +45 with respect to difference)");
% 
% yyaxis right
% plot(allStrength,meanResponsePerConEx_135,'-v');
% hold on
% plot(allStrength,meanResponsePerConUnex_135,'--v');
% xlabel("difference value (with respect to +45)")
% ylabel("p (reported -45 with respect to difference)");
% title("p (reported -45/+45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected 45","unexpected 45","expected -45","unexpected -45");
% saveas(figure(f),[plotdir_sub,'/preported_strength_split.jpg']);
% f=f+1;
% figure();
% plot(allStrength,meanResponsePerConEx,'-v');
% hold on
% plot(allStrength,meanResponsePerConUnex,'--v');
% xlabel("difference value (with respect to +45)")
% ylabel("p (reported +45 with respect to difference)");
% title("p (reported +45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected 45","unexpected 45");
% saveas(figure(f),[plotdir_sub,'/preported45_strength_split.jpg']);
% f=f+1;
% figure();
% plot(allStrength,meanResponsePerConEx_135,'-v');
% hold on
% plot(allStrength,meanResponsePerConUnex_135,'--v');
% xlabel("difference value (with respect to -45)")
% ylabel("p (reported -45 with respect to difference)");
% title("p (reported -45 with respect to difference) for expected vs unexpected ",subjectID);
% legend("expected -45","unexpected -45");
% saveas(figure(f),[plotdir_sub,'/preported135_strength_split.jpg']);
% 
% 
% %% New adjustments
% 
% % adjust for ceiling or floor values
% % nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% % nh(nh==0) = 1;
% % nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% % nfa(nfa==0) = 1;
% 
% % %signal present= +45 (1)
% % %signal absent = -45 (2)
% 
% f=f+1;
% signalPresent=1;
% signalAbsent=2;
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% % get all stair values and remove waffle trials
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% % get all validities and remove waffle trials
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% % get all responses and remove waffle trials
% response= file.d.targetResponseKey;
% response(skipWaffle(1:length(response)))=[];
% % get all difference multipliers and remove waffle trials
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% % get all static grating values and remove waffle trials
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% % a correct response occurs when the response is 9 and static is standard
% % and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% % response is 0 and static is test and multiplier is 1, response is 0 and
% % static is standard and multiplier is 2
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% response_45=(response==responseKeyIds(2)); %response was 0 ie CW
% response_135=(response==responseKeyIds(1)); %response was 9 ie CCW 
% 
% orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
% nTrialsGrating=length(gratingContrastsIdx);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
% %sort based on stair values
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% 
% % Expected
% %response present is CW
% responsePresent=response_45;
% %orientation present is CW
% orientationPresent=(static==1 & diffMult==1) | (static==2 & diffMult==-1);
% %precue validity is valid ie 1
% precueValid=(precueValidities==1);
% % get hits for expected
% hitsE=responsePresent==1 & orientationPresent==1 & precueValid==1;
% 
% %orientation absent: is CCW
% orientationAbsent=(static==2 & diffMult==1) | (static==1 & diffMult==-1);
% % get false alarms for expected
% fa_ex=responsePresent==1 & orientationAbsent==1 & precueValid==1;
% 
% % response absent is CCW
% responseAbsent=response_135;
% % misses for expected
% miss_ex1= responseAbsent==1 & orientationPresent==1 & precueValid==1;
% %correct rejection
% cr_ex1=responseAbsent==1 & orientationAbsent==1;
% 
% % Unexpected
% %precue validity is valid ie 1
% precueInvalid=(precueValidities==2);
% % get hits for unexpected
% hitsU=responsePresent==1 & orientationPresent==1 & precueInvalid==1;
% % get false alarms for unexpected
% fa_unex=responsePresent==1 & orientationAbsent==1  & precueInvalid==1;
% % get misses  for unexpected
% miss_unex1= responseAbsent==1 & orientationPresent==1 & precueInvalid==1;    
% % get correct rejection for unexpected
% cr_unex1=responseAbsent==1 & orientationAbsent==1  & precueInvalid==1;    
% 
% orientationJudgementEx=orientationJudgement(orientationJudgement(:,2)==1);
% orientationJudgementUnex=orientationJudgement(orientationJudgement(:,2)==2);
% 
% allstairsEx=unique(orientationJudgementEx);
% allstairsUnex=unique(orientationJudgementUnex);
% allstairsIdEx=[];
% allstairsIdUnex=[];
% 
% for i=1:length(allstairsEx)
%     allstairsIdEx=[allstairsIdEx find(file.p.stairs==allstairsEx(i))];
% end
% 
% for i=1:length(allstairsUnex)
%     allstairsIdUnex=[allstairsIdUnex find(file.p.stairs==allstairsUnex(i))];
% end
% 
% stairExp=[];
% stairExpId=file.d.stairIdxExp;
% skipsSavedEx=isnan(stairExpId);
% stairExpId(skipsSavedEx)=[];
% for i=1:length(stairExpId)
%     stairExp=[stairExp file.p.stairs(stairExpId(i))];
% end
% stairUnexp=[];
% stairUnexpId=file.d.stairIdxUn;
% skipsSavedUn=isnan(stairUnexpId);
% stairUnexpId(skipsSavedUn)=[];
% for i=1:length(stairUnexpId)
%     stairUnexp=[stairUnexp file.p.stairs(stairUnexpId(i))];
% end
% 
% hit_rates_ex=[];
% fa_rates_ex=[];
% cons_dprime_ex=[];
% cons_c_ex=[];
% miss_ex=[];
% cr_ex=[];
% 
% 
% for i=1:length(allstairsEx)
% 
%     % Expected
%     hitsEstair=(hitsE==1) & (stair==allstairsEx(i));
%     hit_rate_ex=mean(hitsEstair);
%     hit_rates_ex(i)=hit_rate_ex;
% 
%     faEstair=(fa_ex==1) & (stair==allstairsEx(i));
%     fa_rate_ex=mean(faEstair);
%     fa_rates_ex(i)=fa_rate_ex;
% 
%     %   get dprime and criterion  
%     nh=sum(hitsEstair);
%     nfa=sum(faEstair);
%     signal=orientationPresent==1 & (stair==allstairsEx(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsEx(i));
%     nnoise=sum(noise);
%     if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_ex = nh./nsignal;
%            fa_rate_ex = nfa./nnoise;
%     end
% 
%     cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
%     cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
% 
%     missEstair=miss_ex1==1 & (stair==allstairsEx(i));
%     miss_ex(i)=mean(missEstair);
%     crEstair=cr_ex1==1 & (stair==allstairsEx(i));
%     cr_ex(i)=mean(crEstair);
% 
% end
% 
% hit_rates_unex=[];
% fa_rates_unex=[];
% cons_dprime_unex=[];
% cons_c_unex=[];
% miss_unex=[];
% cr_unex=[];
% % Unexpected    
% for i=1:length(allstairsUnex)
% 
%     hitsU_stair= hitsU==1 & (stair==allstairsUnex(i));
%     hit_rate_unex=mean(hitsU_stair);
%     hit_rates_unex(i)=hit_rate_unex;
% 
%     faU_stair=fa_unex==1 & (stair==allstairsUnex(i));
%     fa_rate_unex=mean(faU_stair);
%     fa_rates_unex(i)=fa_rate_unex;
% 
%      %   get dprime and criterion  
%     nh=sum(hitsU_stair);
%     nfa=sum(faU_stair);
%     signal=orientationPresent==1 & (stair==allstairsUnex(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsUnex(i));
%     nnoise=sum(noise);
%     if fa_rate_unex==0 || fa_rate_unex==1  || hit_rate_unex==0  || hit_rate_unex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_unex = nh./nsignal;
%            fa_rate_unex = nfa./nnoise;
%     end
% 
%     %   get dprime and criterion
%     cons_dprime_unex(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
%     cons_c_unex(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
%     missUstair=miss_unex1==1 & (stair==allstairsUnex(i));
%     miss_unex(i)=mean(missUstair);
%     crUstair= cr_unex1==1 & (stair==allstairsUnex(i));
%     cr_unex(i)=mean(crUstair);
% end
% 
% % stair ids
% figure();
% subplot(2,1,1);
% plot(allstairsIdEx,cons_dprime_ex,'blue-v');
% hold on
% plot(allstairsIdUnex,cons_dprime_unex,'red-v');
% xlabel("stairs ID")
% ylabel("d'");
% title("d' expected vs unexpected");
% legend("expected","unexpected");
% 
% subplot(2,1,2);
% plot(allstairsIdEx,cons_c_ex,'blue-v');
% hold on
% plot(allstairsIdUnex,cons_c_unex,'red-v');
% xlabel("stairs ID")
% ylabel("criterion");
% title("criterion expected vs unexpected, ",subjectID);
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/dprime_c_split_stairsID.jpg']);
% f=f+1;
% % stair values
% 
% figure();
% subplot(2,1,1);
% plot(allstairsEx,cons_dprime_ex,'blue-v');
% hold on
% plot(allstairsUnex,cons_dprime_unex,'red-v');
% xlabel("stairs")
% ylabel("d'");
% title("d' expected vs unexpected");
% legend("expected","unexpected");
% 
% subplot(2,1,2);
% plot(allstairsEx,cons_c_ex,'blue-v');
% hold on
% plot(allstairsUnex,cons_c_unex,'red-v');
% xlabel("stairs")
% ylabel("criterion");
% title("criterion expected vs unexpected, ",subjectID);
% legend("expected","unexpected");
% 
% saveas(figure(f),[plotdir_sub,'/dprime_c_split_stairs.jpg']);
% 
% %% New adjustments + Signal is expected
% 
% % adjust for ceiling or floor values
% % nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% % nh(nh==0) = 1;
% % nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% % nfa(nfa==0) = 1;
% 
% % %signal present= expected
% % %signal absent = unexpected
% f=f+1;
% signalPresent=1;
% signalAbsent=2;
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% % get all stair values and remove waffle trials
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% % get all validities and remove waffle trials
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% % get all responses and remove waffle trials
% response= file.d.targetResponseKey;
% response(skipWaffle(1:length(response)))=[];
% % get all difference multipliers and remove waffle trials
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% % get all static grating values and remove waffle trials
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% % a correct response occurs when the response is 9 and static is standard
% % and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% % response is 0 and static is test and multiplier is 1, response is 0 and
% % static is standard and multiplier is 2
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% response_45=(response==responseKeyIds(2)); %response was 0 ie +45
% response_135=(response==responseKeyIds(1)); %response was 9 ie -45 or 135
% 
% orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
% nTrialsGrating=length(gratingContrastsIdx);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
% %sort based on stair values
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% 
% % Expected
% %response present is as expected (i.e. responsed 45 when it was 45 and -45
% %when it was -45)
% %observer responded with expected orientation 
% responsePresent=(response_45==1 & diffMult==1 & static==1 & precueValidities==1) | (response_45==1 & static==2 & diffMult==-1 & precueValidities==1) | (response_135==1 & static==1 & diffMult==-1 & precueValidities==1) | (response_135==1 & static==2 & diffMult==1 & precueValidities==1);
% % signal: expected orientation
% validPresent=precueValidities==1;
% 
% % get hits 
% hits=responsePresent==1 & validPresent==1;
% %signal absent: was unexpected
% validAbsent=precueValidities==2;
% % get false alarms for expected
% fa=responsePresent==1 & validAbsent==1;
% 
% % response absent is responsed -45 when +45 and +45 when -45
% responseAbsent=(response_135==1 & diffMult==1 & static==1 & precueValidities==1) | (response_135==1 & static==2 & diffMult==-1 & precueValidities==1) | (response_45==1 & static==1 & diffMult==-1 & precueValidities==1) | (response_45==1 & static==2 & diffMult==1 & precueValidities==1);
% % misses  
% miss_ex1= responseAbsent==1 & validPresent==1;
% %correct rejection
% cr_ex1=responseAbsent==1 & validAbsent==1;
% 
% allstairs=unique(stair);
% allstairsId=[];
% 
% for i=1:length(allstairs)
%     allstairsId=[allstairsId find(file.p.stairs==allstairs(i))];
% end
% 
% 
% hit_rates=[];
% fa_rates=[];
% cons_dprime=[];
% cons_c=[];
% miss_ex=[];
% cr_ex=[];
% 
% 
% for i=1:length(allstairs)
% 
%     % Expected
%     hitsEstair=(hits==1) & (stair==allstairs(i));
%     hit_rate=mean(hitsEstair);
%     hit_rates(i)=hit_rate;
% 
%     faEstair=(fa==1) & (stair==allstairs(i));
%     fa_rate=mean(faEstair);
%     fa_rates(i)=fa_rate;
% 
%     %   get dprime and criterion  
%     nh=sum(hitsEstair);
%     nfa=sum(faEstair);
%     signal=validPresent==1 & (stair==allstairs(i));
%     nsignal=sum(signal);
%     noise=validAbsent==1 & (stair==allstairs(i));
%     nnoise=sum(noise);
%     if fa_rate==0 || fa_rate==1  || hit_rate==0  || hit_rate==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate = nh./nsignal;
%            fa_rate = nfa./nnoise;
%     end
% 
%     cons_dprime(i)=norminv(hit_rate) - norminv(fa_rate);
%     cons_c(i)=-(norminv(hit_rate) + norminv(fa_rate))/2;
% 
%     missEstair=miss_ex1==1 & (stair==allstairs(i));
%     miss_ex(i)=mean(missEstair);
%     crEstair=cr_ex1==1 & (stair==allstairs(i));
%     cr_ex(i)=mean(crEstair);
% 
% end
% 
% % stair ids
% figure();
% subplot(2,1,1);
% plot(allstairsId,cons_dprime,'blue-v');
% xlabel("stair IDs")
% ylabel("d'");
% title("d'");
% 
% subplot(2,1,2);
% plot(allstairsId,cons_c,'blue-v');
% xlabel("stair IDs")
% ylabel("criterion");
% title("criterion, ",subjectID);
% 
% saveas(figure(f),[plotdir_sub,'/dprime_c_signalExp_stairsID.jpg']);
% 
% % stair values
% f=f+1;
% figure();
% subplot(2,1,1);
% plot(allstairs,cons_dprime,'blue-v');
% xlabel("stairs")
% ylabel("d'");
% title("d'");
% 
% subplot(2,1,2);
% plot(allstairs,cons_c,'blue-v');
% xlabel("stairs")
% ylabel("criterion");
% title("criterion, ",subjectID);
% saveas(figure(f),[plotdir_sub,'/dprime_c_signalExp_stairs.jpg']);
% 
% 
% %% New adjustments + expected/unexpected by +45 or -45
% 
% % adjust for ceiling or floor values
% % nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% % nh(nh==0) = 1;
% % nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% % nfa(nfa==0) = 1;
% 
% % %signal present= +45 (1)
% % %signal absent = -45 (2)
% f=f+1;
% signalPresent=1;
% signalAbsent=2;
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% % get all stair values and remove waffle trials
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% % get all validities and remove waffle trials
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% % get all responses and remove waffle trials
% response= file.d.targetResponseKey;
% response(skipWaffle(1:length(response)))=[];
% % get all difference multipliers and remove waffle trials
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% % get all static grating values and remove waffle trials
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% overallOri=file.d.orientation;
% overallOri(skipWaffle)=[];
% 
% % a correct response occurs when the response is 9 and static is standard
% % and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% % response is 0 and static is test and multiplier is 1, response is 0 and
% % static is standard and multiplier is 2
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% response_45=(response==responseKeyIds(2)); %response was 0 ie CW
% response_135=(response==responseKeyIds(1)); %response was 9 ie CCW 
% 
% orientationJudgement= cat(9,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult,overallOri);
% nTrialsGrating=length(gratingContrastsIdx);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,9]);
% %sort based on stair values
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% 
% % Expected
% %response present is CW
% responsePresent=response_45;
% %orientation present is CW
% orientationPresent=(static==1 & diffMult==1) | (static==2 & diffMult==-1);
% %precue validity is valid ie 1
% precueValid=(precueValidities==1);
% % get hits for expected +45
% hitsE45=responsePresent==1 & orientationPresent==1 & precueValid==1 & overallOri==2;
% 
% % get hits for expected -45/135
% hitsE135=responsePresent==1 & orientationPresent==1 & precueValid==1 & overallOri==1;
% 
% %orientation absent: is CCW
% orientationAbsent=(static==2 & diffMult==1) | (static==1 & diffMult==-1);
% % get false alarms for expected +45
% fa_ex45=responsePresent==1 & orientationAbsent==1 & precueValid==1 & overallOri==2;
% 
% % get false alarms for expected -45/135
% fa_ex135=responsePresent==1 & orientationAbsent==1 & precueValid==1 & overallOri==1;
% 
% 
% % response absent is CCW
% responseAbsent=response_135;
% % misses for expected 45
% miss_ex1_45= responseAbsent==1 & orientationPresent==1 & precueValid==1 & overallOri==2;
% % misses for expected 135
% miss_ex1_135= responseAbsent==1 & orientationPresent==1 & precueValid==1 & overallOri==1;
% 
% %correct rejection for expected 45
% cr_ex1_45=responseAbsent==1 & orientationAbsent==1 & precueValid==1 & overallOri==2;
% %correct rejection for expected 135 
% cr_ex1_135=responseAbsent==1 & orientationAbsent==1 & precueValid==1 & overallOri==1;
% 
% % Unexpected
% %precue validity is valid ie 1
% precueInvalid=(precueValidities==2);
% % get hits for unexpected 45
% hitsU45=responsePresent==1 & orientationPresent==1 & precueInvalid==1 & overallOri==2;
% % get hits for unexpected 135
% hitsU135=responsePresent==1 & orientationPresent==1 & precueInvalid==1 & overallOri==1;
% % get false alarms for unexpected 45
% fa_unex45=responsePresent==1 & orientationAbsent==1  & precueInvalid==1 & overallOri==2;
% % get false alarms for unexpected 135
% fa_unex135=responsePresent==1 & orientationAbsent==1  & precueInvalid==1 & overallOri==1;
% % get misses  for unexpected 45
% miss_unex1_45= responseAbsent==1 & orientationPresent==1 & precueInvalid==1 & overallOri==2;   
% % get misses  for unexpected 135
% miss_unex1_135= responseAbsent==1 & orientationPresent==1 & precueInvalid==1 & overallOri==1;  
% % get correct rejection for unexpected 45
% cr_unex1_45=responseAbsent==1 & orientationAbsent==1  & precueInvalid==1 & overallOri==2;    
% % get correct rejection for unexpected 135
% cr_unex1_135=responseAbsent==1 & orientationAbsent==1  & precueInvalid==1 & overallOri==1; 
% 
% orientationJudgementEx45=orientationJudgement(orientationJudgement(:,2)==1 & orientationJudgement(:,9)==2);
% orientationJudgementEx135=orientationJudgement(orientationJudgement(:,2)==1 & orientationJudgement(:,9)==1);
% 
% orientationJudgementUnex45=orientationJudgement(orientationJudgement(:,2)==2 & orientationJudgement(:,9)==2);
% orientationJudgementUnex135=orientationJudgement(orientationJudgement(:,2)==2 & orientationJudgement(:,9)==1);
% 
% allstairsEx45=unique(orientationJudgementEx45);
% allstairsEx135=unique(orientationJudgementEx135);
% 
% allstairsUnex45=unique(orientationJudgementUnex45);
% allstairsUnex135=unique(orientationJudgementUnex135);
% 
% allstairsIdEx45=[];
% allstairsIdEx135=[];
% allstairsIdUnex45=[];
% allstairsIdUnex135=[];
% 
% for i=1:length(allstairsEx45)
%     allstairsIdEx45=[allstairsIdEx45 find(file.p.stairs==allstairsEx45(i))];
% end
% for i=1:length(allstairsEx135)
%     allstairsIdEx135=[allstairsIdEx135 find(file.p.stairs==allstairsEx135(i))];
% end
% 
% for i=1:length(allstairsUnex45)
%     allstairsIdUnex45=[allstairsIdUnex45 find(file.p.stairs==allstairsUnex45(i))];
% end
% for i=1:length(allstairsUnex135)
%     allstairsIdUnex135=[allstairsIdUnex135 find(file.p.stairs==allstairsUnex135(i))];
% end
% 
% stairExp=[];
% stairExpId=file.d.stairIdxExp;
% skipsSavedEx=isnan(stairExpId);
% stairExpId(skipsSavedEx)=[];
% for i=1:length(stairExpId)
%     stairExp=[stairExp file.p.stairs(stairExpId(i))];
% end
% stairUnexp=[];
% stairUnexpId=file.d.stairIdxUn;
% skipsSavedUn=isnan(stairUnexpId);
% stairUnexpId(skipsSavedUn)=[];
% for i=1:length(stairUnexpId)
%     stairUnexp=[stairUnexp file.p.stairs(stairUnexpId(i))];
% end
% 
% hit_rates_ex_45=[];
% fa_rates_ex_45=[];
% cons_dprime_ex_45=[];
% cons_c_ex_45=[];
% miss_ex_45=[];
% cr_ex_45=[];
% 
% 
% for i=1:length(allstairsEx45)
% 
%     % Expected
%     hitsEstair=(hitsE45==1) & (stair==allstairsEx45(i));
%     hit_rate_ex=mean(hitsEstair);
%     hit_rates_ex_45(i)=hit_rate_ex;
% 
%     faEstair=(fa_ex45==1) & (stair==allstairsEx45(i));
%     fa_rate_ex=mean(faEstair);
%     fa_rates_ex_45(i)=fa_rate_ex;
% 
%     %   get dprime and criterion  
%     nh=sum(hitsEstair);
%     nfa=sum(faEstair);
%     signal=orientationPresent==1 & (stair==allstairsEx45(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsEx45(i));
%     nnoise=sum(noise);
%     if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_ex = nh./nsignal;
%            fa_rate_ex = nfa./nnoise;
%     end
% 
%     cons_dprime_ex45(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
%     cons_c_ex45(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
% 
%     missEstair=miss_ex1_45==1 & (stair==allstairsEx45(i));
%     miss_ex_45(i)=mean(missEstair);
%     crEstair=cr_ex1_45==1 & (stair==allstairsEx45(i));
%     cr_ex_45(i)=mean(crEstair);
% 
% end
% 
% hit_rates_ex_135=[];
% fa_rates_ex_135=[];
% cons_dprime_ex_135=[];
% cons_c_ex_135=[];
% miss_ex_135=[];
% cr_ex_135=[];
% 
% 
% for i=1:length(allstairsEx135)
% 
%     % Expected
%     hitsEstair=(hitsE135==1) & (stair==allstairsEx135(i));
%     hit_rate_ex=mean(hitsEstair);
%     hit_rates_ex_135(i)=hit_rate_ex;
% 
%     faEstair=(fa_ex135==1) & (stair==allstairsEx135(i));
%     fa_rate_ex=mean(faEstair);
%     fa_rates_ex_135(i)=fa_rate_ex;
% 
%     %   get dprime and criterion  
%     nh=sum(hitsEstair);
%     nfa=sum(faEstair);
%     signal=orientationPresent==1 & (stair==allstairsEx135(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsEx135(i));
%     nnoise=sum(noise);
%     if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_ex = nh./nsignal;
%            fa_rate_ex = nfa./nnoise;
%     end
% 
%     cons_dprime_ex135(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
%     cons_c_ex135(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
% 
%     missEstair=miss_ex1_135==1 & (stair==allstairsEx135(i));
%     miss_ex_135(i)=mean(missEstair);
%     crEstair=cr_ex1_135==1 & (stair==allstairsEx135(i));
%     cr_ex_135(i)=mean(crEstair);
% 
% end
% 
% hit_rates_unex45=[];
% fa_rates_unex45=[];
% cons_dprime_unex45=[];
% cons_c_unex45=[];
% miss_unex45=[];
% cr_unex45=[];
% % Unexpected    
% for i=1:length(allstairsUnex45)
% 
%     hitsU_stair= hitsU45==1 & (stair==allstairsUnex45(i));
%     hit_rate_unex=mean(hitsU_stair);
%     hit_rates_unex45(i)=hit_rate_unex;
% 
%     faU_stair=fa_unex45==1 & (stair==allstairsUnex45(i));
%     fa_rate_unex=mean(faU_stair);
%     fa_rates_unex45(i)=fa_rate_unex;
% 
%      %   get dprime and criterion  
%     nh=sum(hitsU_stair);
%     nfa=sum(faU_stair);
%     signal=orientationPresent==1 & (stair==allstairsUnex45(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsUnex45(i));
%     nnoise=sum(noise);
%     if fa_rate_unex==0 || fa_rate_unex==1  || hit_rate_unex==0  || hit_rate_unex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_unex = nh./nsignal;
%            fa_rate_unex = nfa./nnoise;
%     end
% 
%     %   get dprime and criterion
%     cons_dprime_unex45(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
%     cons_c_unex45(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
%     missUstair=miss_unex1_45==1 & (stair==allstairsUnex45(i));
%     miss_unex45(i)=mean(missUstair);
%     crUstair= cr_unex1_45==1 & (stair==allstairsUnex45(i));
%     cr_unex45(i)=mean(crUstair);
% end
% 
% hit_rates_unex135=[];
% fa_rates_unex135=[];
% cons_dprime_unex135=[];
% cons_c_unex135=[];
% miss_unex135=[];
% cr_unex135=[];
% % Unexpected    
% for i=1:length(allstairsUnex135)
% 
%     hitsU_stair= hitsU135==1 & (stair==allstairsUnex135(i));
%     hit_rate_unex=mean(hitsU_stair);
%     hit_rates_unex135(i)=hit_rate_unex;
% 
%     faU_stair=fa_unex135==1 & (stair==allstairsUnex135(i));
%     fa_rate_unex=mean(faU_stair);
%     fa_rates_unex135(i)=fa_rate_unex;
% 
%      %   get dprime and criterion  
%     nh=sum(hitsU_stair);
%     nfa=sum(faU_stair);
%     signal=orientationPresent==1 & (stair==allstairsUnex135(i));
%     nsignal=sum(signal);
%     noise=orientationAbsent==1 & (stair==allstairsUnex135(i));
%     nnoise=sum(noise);
%     if fa_rate_unex==0 || fa_rate_unex==1  || hit_rate_unex==0  || hit_rate_unex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_unex = nh./nsignal;
%            fa_rate_unex = nfa./nnoise;
%     end
% 
%     %   get dprime and criterion
%     cons_dprime_unex135(i)=norminv(hit_rate_unex) - norminv(fa_rate_unex);
%     cons_c_unex135(i)=-(norminv(hit_rate_unex) + norminv(fa_rate_unex))/2;
%     missUstair=miss_unex1_135==1 & (stair==allstairsUnex135(i));
%     miss_unex135(i)=mean(missUstair);
%     crUstair= cr_unex1_135==1 & (stair==allstairsUnex135(i));
%     cr_unex135(i)=mean(crUstair);
% end
% 
% 
% % stair ids
% figure();
% subplot(2,1,1);
% plot(allstairsIdEx45,cons_dprime_ex45,'-v');
% hold on
% plot(allstairsIdEx135,cons_dprime_ex135,'-v');
% hold on
% plot(allstairsIdUnex45,cons_dprime_unex45,'-v');
% hold on
% plot(allstairsIdUnex135,cons_dprime_unex135,'-v');
% xlabel("stairs ID")
% ylabel("d'");
% title("d' expected vs unexpected");
% legend("expected 45","expected 135","unexpected 45","unexpected 135");
% 
% subplot(2,1,2);
% plot(allstairsIdEx45,cons_c_ex45,'-v');
% hold on
% plot(allstairsIdEx135,cons_c_ex135,'-v');
% hold on
% plot(allstairsIdUnex45,cons_c_unex45,'-v');
% hold on
% plot(allstairsIdUnex135,cons_c_unex135,'-v');
% xlabel("stairs ID")
% ylabel("criterion");
% title("criterion expected vs unexpected, ",subjectID);
% legend("expected 45","expected 135","unexpected 45","unexpected 135");
% saveas(figure(f),[plotdir_sub,'/dprime_c_split45-135_stairsID.jpg']);
% 
% % stair values
% f=f+1;
% figure();
% subplot(2,1,1);
% plot(allstairsEx45,cons_dprime_ex45,'-v');
% hold on
% plot(allstairsEx135,cons_dprime_ex135,'-v');
% hold on
% plot(allstairsUnex45,cons_dprime_unex45,'-v');
% hold on
% plot(allstairsUnex135,cons_dprime_unex135,'-v');
% xlabel("stairs")
% ylabel("d'");
% title("d' expected vs unexpected");
% legend("expected 45","expected 135","unexpected 45","unexpected 135");
% 
% subplot(2,1,2);
% plot(allstairsEx45,cons_c_ex45,'-v');
% hold on
% plot(allstairsEx135,cons_c_ex135,'-v');
% hold on
% plot(allstairsUnex45,cons_c_unex45,'-v');
% hold on
% plot(allstairsUnex135,cons_c_unex135,'-v');
% xlabel("stairs")
% ylabel("criterion");
% title("criterion expected vs unexpected, ",subjectID);
% legend("expected 45","expected 135","unexpected 45","unexpected 135");
% 
% saveas(figure(f),[plotdir_sub,'/dprime_c_split45-135_stairs.jpg']);
% 
% 
% 
% %% p (correct tilt discrimination) - split by 16
% f=f+1;
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% stair=file.d.stair;
% %stair(skipWaffle(1:length(stair)))=[];
% 
% precueValidities=file.d.precueValidity;
% %precueValidities(skipWaffle)=[];
% 
% response= file.d.targetResponseKey;
% %response(skipWaffle)=[];
% 
% diffMult=file.d.differenceMultiplier;
% %diffMult(skipWaffle(1:length(diffMult)))=[];
% 
% static=file.d.staticGrating;
% %static(skipWaffle(1:length(static)))=[];
% 
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% orientationJudgement= cat(4,stair,precueValidities,correct_r,response);
% nTrialsGrating=length(stair);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,4]);
% orientationJudgement_blocks=[];
% for i=1:16
%     orientationJudgement_blocks{i}=orientationJudgement(i*46-46+1:i*46,:);
%     orientationJudgement_blocks{i}(skipWaffle(i*46-46+1:i*46),:)=[];
%     %sort based on stair levels
%     orientationJudgement_blocks{i}=sortrows(orientationJudgement_blocks{i},1); 
%     % sort based on precue validity
%     orientationJudgement_blocks{i}=sortrows(orientationJudgement_blocks{i},2);
% end
% orientationJudgement(skipWaffle(1:length(orientationJudgement)),:)=[];
% 
% %sort based on stair levels
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% correctsEx=[];
% correctsUnex=[];
% meanCorrectPerConEx=[];
% meanCorrectPerConUnex=[];
% 
% for j=1:length(orientationJudgement_blocks)
%     correctsEx{j}=[];
%     correctsUnex{j}=[];
%     meanCorrectPerConEx{j}=[];
%     meanCorrectPerConUnex{j}=[];
%     for i=1:length(orientationJudgement_blocks{j})
% 
%         if orientationJudgement_blocks{j}(i,2)==1%% expected
%             if i==1
%                 correctsEx{j} =[correctsEx{j} orientationJudgement_blocks{j}(i,3)];
% 
%             elseif orientationJudgement_blocks{j}(i-1,1)==orientationJudgement_blocks{j}(i,1) 
%                 correctsEx{j}=[correctsEx{j} orientationJudgement_blocks{j}(i,3)];
%                 if i==length(orientationJudgement_blocks{j})
%                     meanCorrectPerConEx{j}=[meanCorrectPerConEx{j} mean(correctsEx{j})];
%                 end
% 
%             elseif  orientationJudgement_blocks{j}(i-1,1)~=orientationJudgement_blocks{j}(i,1)
%                 meanCorrectPerConEx{j}=[meanCorrectPerConEx{j} mean(correctsEx{j})];
%                 correctsEx{j}=[];
%                 correctsEx{j}=[correctsEx{j} orientationJudgement_blocks{j}(i,3)];
%                 if i==orientationJudgement_blocks{j}
%                    meanCorrectPerConEx{j}=[meanCorrectPerConEx{j} mean(correctsEx{j})];
%                 end
%             end
%         elseif orientationJudgement_blocks{j}(i,2)==2 
%             if orientationJudgement_blocks{j}(i-1,2)==1
%                 meanCorrectPerConEx{j}=[meanCorrectPerConEx{j} mean(correctsEx{j})];
%                 correctsUnex{j}=[ correctsUnex{j} orientationJudgement_blocks{j}(i,3)];
% 
%             elseif orientationJudgement_blocks{j}(i-1,1)==orientationJudgement_blocks{j}(i,1) 
%                 if i==length(orientationJudgement_blocks{j})
%                      correctsUnex{j}=[ correctsUnex{j} orientationJudgement_blocks{j}(i,3)];
%                     meanCorrectPerConUnex{j}=[meanCorrectPerConUnex{j} mean(correctsUnex{j})];
%                 end
%                  correctsUnex{j}=[ correctsUnex{j} orientationJudgement_blocks{j}(i,3)];
% 
% 
%             elseif  orientationJudgement_blocks{j}(i-1,1)~=orientationJudgement_blocks{j}(i,1)
%                 meanCorrectPerConUnex{j}=[meanCorrectPerConUnex{j} mean( correctsUnex{j})];
%                  correctsUnex{j}=[];
%                  correctsUnex{j}=[orientationJudgement_blocks{j}(i,3)];
%                 if i==length(orientationJudgement_blocks{j})
%                    meanCorrectPerConUnex{j}=[meanCorrectPerConUnex{j} mean( correctsUnex{j})];
%                 end
%             end
%         end
%     end
% end
% 
% orientationJudgementEx=[];
% orientationJudgementUnex=[];
% for j=1:length(orientationJudgement_blocks)
%     orientationJudgementEx{j}=orientationJudgement_blocks{j}(orientationJudgement_blocks{j}(:,2)==1);
%     orientationJudgementUnex{j}=orientationJudgement_blocks{j}(orientationJudgement_blocks{j}(:,2)==2);
% end
% allstairsEx=[];
% allstairsUnex=[];
% for j=1:length(orientationJudgement_blocks)
%     allstairsEx{j}=unique(orientationJudgementEx{j});
%     allstairsUnex{j}=unique(orientationJudgementUnex{j});
% end
% 
% allstairsIdEx=[];
% allstairsIdUnex=[];
% for j=1:length(orientationJudgement_blocks)
%     allstairsIdEx{j}=[];
%     allstairsIdUnex{j}=[];
%     for i=1:length(allstairsEx{j})
%         allstairsIdEx{j}=[allstairsIdEx{j} find(file.p.stairs==allstairsEx{j}(i))];
%     end
%     for i=1:length(allstairsUnex{j})
%         allstairsIdUnex{j}=[allstairsIdUnex{j} find(file.p.stairs==allstairsUnex{j}(i))];
%     end
% end
% 
% figure(); 
% tiledlayout(2,8)
% for j=1:length(orientationJudgement_blocks)
%     nexttile(j)
%     plot(allstairsEx{j},meanCorrectPerConEx{j},'blue-v');
%     hold on
%     plot(allstairsUnex{j},meanCorrectPerConUnex{j},'red-v');
%     xlabel("stair value")
%     ylabel("p(correct tilt discrimination)");
%     title("block " + j);
% end
% sgtitle("p (correct tilt discrimination) for expected vs unexpected in the 16 blocks");
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_blocks_stairs.jpg']);
% f=f+1;
% figure();
% tiledlayout(2,8)
% for j=1:length(orientationJudgement_blocks)
%     nexttile(j)
%     plot(allstairsIdEx{j},meanCorrectPerConEx{j},'blue-v');
%     hold on
%     plot(allstairsIdUnex{j},meanCorrectPerConUnex{j},'red-v');
%     xlabel("stair ID")
%     ylabel("p(correct tilt discrimination)");
%     title("block " + j);
% end
% sgtitle("p (correct tilt discrimination) for expected vs unexpected in the 16 blocks");
% legend("expected","unexpected");
% saveas(figure(f),[plotdir_sub,'/pcorrecttilt_blocks_stairsID.jpg']);
% 

% %% d' criterion signal: expected, split by orientation (needs to be completed)
% 
% % %signal present= expected
% % %signal absent = unexpected
% signalPresent=1;
% signalAbsent=2;
% 
% gratingContrastsIdx=file.d.gratingContrast;
% skipWaffle=isnan(gratingContrastsIdx);
% gratingContrastsIdx(skipWaffle(1:length(gratingContrastsIdx)))=[];
% 
% % get all stair values and remove waffle trials
% stair=file.d.stair;
% stair(skipWaffle(1:length(stair)))=[];
% % get all validities and remove waffle trials
% precueValidities=file.d.precueValidity;
% precueValidities(skipWaffle)=[];
% % get all responses and remove waffle trials
% response= file.d.targetResponseKey;
% response(skipWaffle(1:length(response)))=[];
% % get all difference multipliers and remove waffle trials
% diffMult=file.d.differenceMultiplier;
% diffMult(skipWaffle(1:length(diffMult)))=[];
% % get all static grating values and remove waffle trials
% static=file.d.staticGrating;
% static(skipWaffle(1:length(static)))=[];
% 
% % a correct response occurs when the response is 9 and static is standard
% % and multiplier is -1, response is 9 and static is test and multiplier is 1, 
% % response is 0 and static is test and multiplier is 1, response is 0 and
% % static is standard and multiplier is 2
% correct_r=(response==responseKeyIds(1) & static==1 & diffMult==-1) | (response==responseKeyIds(1) & static==2 & diffMult==1) | (response==responseKeyIds(2) & static==1 & diffMult==1) | (response==responseKeyIds(2) & static==2 & diffMult==-1);
% 
% response_45=(response==responseKeyIds(2)); %response was 0 ie +45
% response_135=(response==responseKeyIds(1)); %response was 9 ie -45 or 135
% 
% orientationJudgement= cat(8,stair,precueValidities,correct_r,response_45,response_135,response,static,diffMult);
% nTrialsGrating=length(gratingContrastsIdx);
% orientationJudgement=reshape(orientationJudgement,[nTrialsGrating,8]);
% %sort based on stair values
% orientationJudgement=sortrows(orientationJudgement,1); 
% % sort based on precue validity
% orientationJudgement=sortrows(orientationJudgement,2);
% 
% 
% % Expected
% %response present is as expected (i.e. responsed 45 when it was 45 and -45
% %when it was -45)
% %observer responded with expected orientation 
% responsePresent=(response_45==1 & diffMult==1 & static==1) | (response_45==1 & static==2 & diffMult==-1) | (response_135==1 & static==1 & diffMult==-1) | (response_135==1 & static==2 & diffMult==1);
% % signal: expected orientation
% validPresent=precueValidities==1;
% 
% % get hits 
% hits=responsePresent==1 & validPresent==1;
% %signal absent: was unexpected
% validAbsent=precueValidities==2;
% % get false alarms for expected
% fa_ex=responsePresent==1 & validAbsent==1;
% 
% % response absent is responsed -45 when +45 and +45 when -45
% responseAbsent=(response_135==1 & diffMult==1 & static==1) | (response_135==1 & static==2 & diffMult==-1) | (response_45==1 & static==1 & diffMult==-1) | (response_45==1 & static==2 & diffMult==1);
% % misses  
% miss_ex1= responseAbsent==1 & validPresent==1;
% %correct rejection
% cr_ex1=responseAbsent==1 & validAbsent==1;
% 
% allstairs=unique(stair);
% allstairsId=[];
% 
% for i=1:length(allstairs)
%     allstairsId=[allstairsId find(file.p.stairs==allstairs(i))];
% end
% 
% 
% hit_rates_ex=[];
% fa_rates_ex=[];
% cons_dprime_ex=[];
% cons_c_ex=[];
% miss_ex=[];
% cr_ex=[];
% 
% 
% for i=1:length(allstairs)
% 
%     % Expected
%     hitsEstair=(hits==1) & (stair==allstairs(i));
%     hit_rate_ex=mean(hitsEstair);
%     hit_rates_ex(i)=hit_rate_ex;
% 
%     faEstair=(fa_ex==1) & (stair==allstairs(i));
%     fa_rate_ex=mean(faEstair);
%     fa_rates_ex(i)=fa_rate_ex;
% 
%     %   get dprime and criterion  
%     nh=sum(hitsEstair);
%     nfa=sum(faEstair);
%     signal=validPresent==1 & (stair==allstairs(i));
%     nsignal=sum(signal);
%     noise=validAbsent==1 & (stair==allstairs(i));
%     nnoise=sum(noise);
%     if fa_rate_ex==0 || fa_rate_ex==1  || hit_rate_ex==0  || hit_rate_ex==1 % % loglinear adjustment (Stanislaw & Todorov 1999) 
%            nh = nh + 0.5;
%            nfa = nfa + 0.5;
%            nnoise = nnoise + 1;
%            nsignal = nsignal + 1;
%            % proportions
%            hit_rate_ex = nh./nsignal;
%            fa_rate_ex = nfa./nnoise;
%     end
% 
%     cons_dprime_ex(i)=norminv(hit_rate_ex) - norminv(fa_rate_ex);
%     cons_c_ex(i)=-(norminv(hit_rate_ex) + norminv(fa_rate_ex))/2;
% 
%     missEstair=miss_ex1==1 & (stair==allstairs(i));
%     miss_ex(i)=mean(missEstair);
%     crEstair=cr_ex1==1 & (stair==allstairs(i));
%     cr_ex(i)=mean(crEstair);
% 
% end
% 
% 
% % stair ids
% figure();
% subplot(2,1,1);
% plot(allstairsId,cons_dprime_ex,'blue-v');
% xlabel("stair IDs")
% ylabel("d'");
% title("d'");
% 
% subplot(2,1,2);
% plot(allstairsId,cons_c_ex,'blue-v');
% xlabel("stair IDs")
% ylabel("criterion");
% title("criterion, ",subjectID);
% 
% % stair values
% 
% figure();
% subplot(2,1,1);
% plot(allstairs,cons_dprime_ex,'blue-v');
% xlabel("stairs")
% ylabel("d'");
% title("d'");
% 
% subplot(2,1,2);
% plot(allstairs,cons_c_ex,'blue-v');
% xlabel("stairs")
% ylabel("criterion");
% title("criterion, ",subjectID);

end
