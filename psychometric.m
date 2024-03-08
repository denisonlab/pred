function psychometric(file)

subjectID = input('Enter subject ID:  ','s');
sessionNum = input('Enter session number (1,2...) ');
dataDir= sprintf('%s/data',pwd);
datafile=sprintf('%s/%s_s%d_pred.mat',dataDir,subjectID,sessionNum);
file=load(datafile);


%contrasts
allcons=file.p.gratingContrasts(1:7);

%% P(test stronger) for unexpected vs expected

phigher= file.d.targetResponseKey == (1 || 4);
phigher= cat(3,file.d.testContrast,phigher,file.d.precueValidity);
nTrials=448;
phigher=reshape(phigher,[nTrials,3]);
phigher=sortrows(phigher,1);
phigher=sortrows(phigher,3);
prop_exp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==2));
cons_highex=[];
cons_highunex=[];
for i = 1:length(file.p.gratingContrasts)
    cons_highex(i)=mean(phigher(prop_exp*i-prop_exp+1:prop_exp*i,2)==1);
    start=prop_exp*length(file.p.gratingContrasts)+1;
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
nTrials=448;
corrects=reshape(corrects,[nTrials,3]);
corrects=sortrows(corrects,1);
corrects=sortrows(corrects,3);
prop_exp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==1));
prop_unexp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==2));
cons_ex=[];
cons_unex=[];

for i = 1:length(file.p.gratingContrasts)
    cons_ex(i)=mean(corrects(prop_exp*i-prop_exp+1:prop_exp*i,2)==1);
    start=prop_exp*length(file.p.gratingContrasts)+1;
    cons_unex(i)=mean(corrects((start+prop_unexp*i)-prop_unexp:start+prop_unexp*i-1,2)==1);
end

figure();
plot(allcons,cons_ex,'blue');
hold on
plot(allcons,cons_unex,'red');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct) for expected/valid vs unexpected/invalid");
legend("expected","unexpected");

%% P(correct)

corrects=cat(2,file.d.testContrast,file.d.correct);
nTrials=448;
corrects=reshape(corrects,[nTrials,2]);
corrects=sortrows(corrects,1);
prop_con=nTrials/length(file.p.gratingContrasts);
%prop_exp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==1));
%prop_unexp=(nTrials/length(file.p.gratingContrasts))*(mean(file.p.precueValidities==2));
cons_cor=[];

for i = 1:length(file.p.gratingContrasts)
    cons_cor(i)=mean(corrects(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons,cons_cor,'blue');
xlabel("contrasts")
ylabel("percent correct");
title("p(correct)");

%% P(test stronger)

phigher= file.d.targetResponseKey == (1 || 4);
phigher= cat(2,file.d.testContrast,phigher);
nTrials=448;
phigher=reshape(phigher,[nTrials,2]);
phigher=sortrows(phigher,1);
prop_con=nTrials/length(file.p.gratingContrasts);
cons_high=[];
for i = 1:length(file.p.gratingContrasts)
    cons_high(i)=mean(phigher(prop_con*i-prop_con+1:prop_con*i,2)==1);
end

figure();
plot(allcons,cons_high,'blue');
xlabel("contrasts")
ylabel("percent that answered 0 or 1 aka perceived test as higher contrast");
title("p(test stronger)");
end
