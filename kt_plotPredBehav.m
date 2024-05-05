function kt_plot(d,p)

figure
set(gcf,'Position',[100 100 1000 600])
t = 1:numel(d.stairIdx); % number of trials
varOIs = {'stairIdx','stair','correct'};

for iV = 1:numel(varOIs)
    subplot (numel(varOIs),1,iV) % plot each of the variables of interest 
    figureStyle
    for iE = 1:2 % two expectation conditions (1 = expected, 2 = unexpected)
        idx = d.precueValidity==iE & ~d.stopThisTrial; 
        var = d.(varOIs{iV}); 
        pl = plot(t(idx), var(idx),'LineWidth',1,'Marker','o','MarkerSize',5);
        pl.Color(4) = 0.3; % alpha
    end
    legend({'Expected','Unexpected'})
    xlabel('Trial')
    ylabel(varOIs{iV})

    % xlim([1 80]) % restrict x lims to look more closely
end

%% Check presentation order of plaids/koks and valid/invalid 
figure
subplot 211
imagesc(d.plaidStatus)

subplot 212
imagesc(d.precueValidity)

%% Check precue distributions 
for iStim = 1:2
    for iE = 1:2
        val = numel(d.precueValidity(d.precueValidity==iE & d.plaidStatus==iStim & ~d.stopThisTrial))/numel(d.precueValidity(d.plaidStatus==iStim & ~d.stopThisTrial));
        pE(iStim,iE) = val; % percentage of expected trials per stimulus type
    end
end

%% Check mean stair idx by valid/invalid 
figure
figureStyle
for iE = 1:2 % two expectation conditions (1 = expected, 2 = unexpected)
    blockIdx = d.block~=9 & d.block~=10 & d.block~=11 & d.block~=12; 
    idx = d.precueValidity==iE & ~d.stopThisTrial & d.plaidStatus==1 & blockIdx; % plaid status 1 is Kok
    var = d.stairIdx(idx);
    bar(iE, mean(var));
end
ylabel('Mean stair index')
xticks([1 2])
xticklabels({'Valid','Invalid'})

%% Check mean stair idx by valid/invalid (just block 9-12)
figure
figureStyle
for iE = 1:2 % two expectation conditions (1 = expected, 2 = unexpected)
    blockIdx = d.block==9 | d.block==10 | d.block==11 | d.block==12; 
    idx = d.precueValidity==iE & ~d.stopThisTrial & d.plaidStatus==1 & blockIdx; % plaid status 1 is Kok
    var = d.stairIdx(idx);
    bar(iE, mean(var));
end
ylabel('Mean stair index')
xticks([1 2])
xticklabels({'Valid','Invalid'})

%% Plaids 
% Check histogram of contrasts 
figure
figureStyle 
edges = 1:8; 
histogram(d.plaidContrast,edges)

%% Plaid plot subjective curve 
% Get average visibility 
for iE = 1:2
    for iC = 1:numel(p.plaidContrasts1) % 8 contrasts
        blockIdx = d.block==9 | d.block==10 | d.block==11 | d.block==12; % bad blocks
        idx = d.plaidContrast==iC & d.precueValidity==iE & ~d.stopThisTrial & d.plaidStatus==2 & ~blockIdx & d.plaidOrientation==2;
        var = d.targetResponseKey(idx);
        valSub{iE,iC} = var; % numel(var==1)/sum(idx);
        valSubMean(iE,iC) = mean(valSub{iE,iC})-1; 
    end
end
% Plot 
figure
figureStyle 
for iE = 1:2
    plot(valSubMean(iE,:))
end
legend({'Expected','Unexpected'})
xlabel('Plaid contrast idx')
ylabel('p(report orientation 2')



