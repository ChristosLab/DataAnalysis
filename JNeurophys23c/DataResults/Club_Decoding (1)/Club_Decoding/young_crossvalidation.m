%% Cross-Validation for Young Neurons

%% Load the Young Data

young_files = fullfile('Neuron_Data_young', '*_1_*.mat');
youngdirectory = dir(young_files);

young_neurons = cell(1,length(youngdirectory));

for i = 1:length(youngdirectory)
    matFileName = fullfile('Neuron_Data_young', youngdirectory(i).name);
    young_neurons{i} = load(matFileName);
end

%% Obtain neuron cue rates from the dataset
% Switch adult and young in variable names to obtain both data sets
% Set parameters
num_neurons = length(young_neurons);
num_classes = 8;
num_trials = 8;
cue_rate_matrix = NaN(num_neurons,num_classes, num_trials);

for i = 1:num_neurons
    neuron = young_neurons{1,i};
    
    % Check to ensure that neuron has 8 classes
    if length(neuron.MatData.class) == 8
        
        for j = 1:length(neuron.MatData.class)
            neuron_class = neuron.MatData.class(j).ntr;
            num_cues = length(neuron_class);
            
            % Check to ensure that each class has 8 cues
            if num_cues >= 8
                
                for k = 1:8
                    cuerate = neuron_class(k).cuerate;
                    if isempty(cuerate)
                        cuerate = NaN;
                    end
                    cue_rate_matrix(i,j,k) = cuerate;
                    
                end
            end
        end
    end
end

%% Reshape the cuerate matrix

names = cell(1,num_neurons);
all_cues = NaN(num_classes*num_trials,num_neurons);

% Obtain neuron names
for ii = 1:num_neurons
    names{ii} = youngdirectory(ii).name;
end

% Reshape cue rate values
for cidx = 0:7
    mult = cidx*num_trials;
    for ii = (1 + mult):(8 + mult)
        for jj = 1:num_neurons
            all_cues(ii,jj) = cue_rate_matrix(jj,1+cidx,ii - mult);
        end
    end
end

% Remove all nan
all_cues = rmmissing(all_cues,2);

%% Random Permutation
% % This method is commented out for initial testing of the decoding method,
% % but is used for second trial when adult and young # of neurons is held
% % constant for testing

num_select = 200;
cue_selection =randperm(length(all_cues), num_select);

% change all_cues to be a constant # of randomly selected group of neurons
all_cues = all_cues(:,cue_selection);

%% Initialize a label vector

all_labels = [];
for ii = 1:8
    mult = ones(8,1);
    all_labels = [all_labels;mult*ii];
end

%% Initialize fitcecoc

fitted = fitcecoc(all_cues, all_labels);

%% Implement Kai's code for cross-validation
% Testing the KFold parameter, which increases the number of folds of the
% cross-validation model. The folds = the number of subsets the dataset is
% subdivided into, where 1 subset is then used as the testing data and k-1
% subsets are used as the training data.

% NOTE: at kfold > 8, some subsets do not have all class values represented
kfolds = 2:20;
meankfloss = NaN(1,length(kfolds));

for ii = 1:length(kfolds)
    tic;
    kfloss = 0;
    kfold = kfolds(ii)
    for jj = 1:10
        CV = crossval(fitted,'KFold',kfold);
        kfloss = kfloss + kfoldLoss(CV)*100;
    end
    meankfloss(ii) = kfloss/jj;
    toc;
end

kfcorrect = 100 - meankfloss;
%% Holdout Parameter
% The holdout is the percentage of data points that are "held out" of the
% training set, to be used in the testing set. Thus, holdout*100 is the
% percentage of data points that are reserved for testing, while 1 -
% holdout*100 is used as the training set for the model.

holdouts = 0.02:0.02:1;
meanholdoutloss = NaN(1,length(holdouts));

for ii = 1:length(holdouts)
    tic;
    holdoutloss = 0;
    holdout = holdouts(ii)
    for jj = 1:10
        CV = crossval(fitted,'Holdout',holdout);
        holdoutloss = holdoutloss + kfoldLoss(CV)*100;
    end
    meanholdoutloss(ii) = holdoutloss/jj;
    toc;
end

holdoutcorrect = 100 - meanholdoutloss;
%% Leave-One-Out 

leaveoutloss = NaN(1, 10);
for ii = 1:10
    tic;
    CV = crossval(fitted, 'Leaveout', 'on');
    leaveoutloss(ii) = kfoldLoss(CV)*100;
    toc;
end
meanleaveoutloss = mean(leaveoutloss);
    
leaveoutcorrect = 100 - leaveoutloss;

%% Next Steps - plotting for all three categories

figure(1)
plot(kfolds, kfcorrect)
title("K-Fold Parameter")
xlabel 'Number of Subsets Created from Dataset'
ylabel 'Average Percent Correct'

figure(2)
plot(holdouts*100, holdoutcorrect)
title ("Holdout Parameter")
xlabel 'Percent Data Held Out for Testing'
ylabel 'Average Percent Correct'

figure(3)
plot(1:10, leaveoutcorrect)
title("Leave-One-Out")
xlabel 'Trial Number'
ylabel 'Average Percent Correct'

%% Identify best parameters

[kfoldmax,kfoldidx] = max(kfcorrect);
[holdoutmax,holdoutidx] = max(holdoutcorrect);

bestkfold = kfolds(kfoldidx)
bestholdout = holdouts(holdoutidx)

%% Save variables
% 
% parameters = matfile('young_crossvalparams','Writable',true);
% 
% parameters.yfitted = fitted;
% 
% parameters.ykfolds = kfolds;
% parameters.ykfcorrect = kfcorrect;
% parameters.ybestkfold = bestkfold;
% 
% parameters.yholdouts = holdouts;
% parameters.yholdoutcorrect = holdoutcorrect;
% parameters.ybestholdout = bestholdout;
% 
% parameters.yleaveoutcorrect = leaveoutcorrect;

% Save variables for method with random permutation

parameters = matfile('yperm_crossvalparams','Writable',true);

parameters.ypermfitted = fitted;

parameters.ypermkfolds = kfolds;
parameters.ypermkfcorrect = kfcorrect;
parameters.ypermbestkfold = bestkfold;

parameters.ypermholdouts = holdouts;
parameters.ypermholdoutcorrect = holdoutcorrect;
parameters.ypermbestholdout = bestholdout;

parameters.ypermleaveoutcorrect = leaveoutcorrect;
