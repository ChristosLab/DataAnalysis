%% Generate Figures for Decoder Analysis and Comparison

%% Compare parameter performance for young and adult - no permutation
clear; clc; close all
load young_crossvalparams.mat
load adult_crossvalparams.mat

figure(1)
plot(akfolds,akfcorrect,'r')
hold on
plot(ykfolds, ykfcorrect, 'b')
title("K-Fold Performance")
xlabel 'Number of K-Folds'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

figure(2)
plot(aholdouts,aholdoutcorrect,'r')
hold on
plot(yholdouts, yholdoutcorrect, 'b')
title("Holdout Performance")
xlabel 'Percent Held Out'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

figure(3)
plot(1:10,aleaveoutcorrect,'r')
hold on
plot(1:10, yleaveoutcorrect, 'b')
title("Leave-One-Out Performance")
xlabel 'Trial Number'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

%% Compare parameter performance for young and adult - with permutation

clear; clc; close all
load yperm_crossvalparams.mat
load aperm_crossvalparams.mat

figure(1)
plot(apermkfolds,apermkfcorrect,'r')
hold on
plot(ypermkfolds, ypermkfcorrect, 'b')
title("K-Fold Performance")
xlabel 'Number of K-Folds'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

figure(2)
plot(apermholdouts,apermholdoutcorrect,'r')
hold on
plot(ypermholdouts, ypermholdoutcorrect, 'b')
title("Holdout Performance")
xlabel 'Percent Held Out'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

figure(3)
plot(1:10,apermleaveoutcorrect,'r')
hold on
plot(1:10, ypermleaveoutcorrect, 'b')
title("Leave-One-Out Performance")
xlabel 'Trial Number'
ylabel 'Percent Correct Identification'
legend("Adult", "Young")

%% Compare parameter performance within adult group

clear; clc; close all
load adult_crossvalparams.mat
load aperm_crossvalparams.mat

figure(1)
plot(akfolds,akfcorrect,'r')
hold on
plot(apermkfolds, apermkfcorrect, 'b')
title("K-Fold Performance")
xlabel 'Number of K-Folds'
ylabel 'Percent Correct Identification'
legend("All Adult Neurons", "200 Adult Neurons")

figure(2)
plot(aholdouts,aholdoutcorrect,'r')
hold on
plot(apermholdouts, apermholdoutcorrect, 'b')
title("Holdout Performance")
xlabel 'Percent Held Out'
ylabel 'Percent Correct Identification'
legend("All Adult Neurons", "200 Adult Neurons")

figure(3)
plot(1:10,aleaveoutcorrect,'r')
hold on
plot(1:10, apermleaveoutcorrect, 'b')
title("Leave-One-Out Performance")
xlabel 'Trial Number'
ylabel 'Percent Correct Identification'
legend("All Adult Neurons", "200 Adult Neurons")

%% Compare parameter performance for young and adult - with permutation

clear; clc; close all
load young_crossvalparams.mat
load yperm_crossvalparams.mat

figure(1)
plot(ykfolds,ykfcorrect,'r')
hold on
plot(ypermkfolds, ypermkfcorrect, 'b')
title("K-Fold Performance")
xlabel 'Number of K-Folds'
ylabel 'Percent Correct Identification'
legend("All Young Neurons", "200 Young Neurons")

figure(2)
plot(yholdouts,yholdoutcorrect,'r')
hold on
plot(ypermholdouts, ypermholdoutcorrect, 'b')
title("Holdout Performance")
xlabel 'Percent Held Out'
ylabel 'Percent Correct Identification'
legend("All Young Neurons", "200 Young Neurons")

figure(3)
plot(1:10,yleaveoutcorrect,'r')
hold on
plot(1:10, ypermleaveoutcorrect, 'b')
title("Leave-One-Out Performance")
xlabel 'Trial Number'
ylabel 'Percent Correct Identification'
legend("All Young Neurons", "200 Young Neurons")

