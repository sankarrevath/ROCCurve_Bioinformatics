% Get Data from CSV Compilation
entryData = csvread('ROC_Curve_PersonalPrediction.csv', 1);
prediction = entryData(:, 2);%Uses the second column as the Prediction values
truth = entryData(:, 1);%Uses the first column as the Truth values
random = entryData(:,3);

%Vars Instantiation
bestThresholdValue = 1;
bestThresholdValue_rand = 1;

%Instantiation of Threshold Matrices - meant to help with the creation of
%the best matrix
threshSpecs = zeros(1, 6);
threshSens = zeros(1, 6);
threshSpecs_random = zeros(1,6);
threshSens_random = zeros(1,6);


% Iterate through each threshold and calculate variables - PREDICTION CURVE
for threshold = (1:6)
    cp = 0; %Correct Positive Prediction
    cn = 0; %Correct Negative Prediction
    inc_p = 0; %Incorrect Pos
    inc_n = 0; %Incorrect Neg
    
    %iterate and assign each to a particular square
    for j = (1 : length(prediction))
        if truth(j) == 1
            if prediction(j) >= threshold
                cp = cp + 1;
            else
                inc_n = inc_n + 1;
            end
        else
            if prediction(j) >= threshold
                inc_p = inc_p + 1;
            else
                cn = cn + 1;
            end
        end
    end
    %assemble final matrix
    conMatrix_random = [cp,inc_p; inc_n,cn];
    threshSens(threshold) = cp / (cp + inc_n);
    threshSpecs(threshold) = 1 - (cn / (cn + inc_p));
end

%Iterate to get Random Curve Sens and Specs
for threshold = (1:6)
    cp_rand = 0; %Correct Positive Prediction
    cn_rand = 0; %Correct Negative Prediction
    inc_p_rand = 0; %Incorrect Pos
    inc_n_rand = 0; %Incorrect Neg
    
    for j = (1 : length(random_data))
        if truth(j) == 1
            if random_data(j) >= threshold
                cp_rand = cp_rand + 1;
            else
                inc_n_rand = inc_n_rand + 1;
            end
        else
            if random_data(j) >= threshold
                inc_p_rand = inc_p_rand + 1;
            else
                cn_rand = cn_rand + 1;
            end
        end
    end
    %assemble final random matrix
    conMatrix_random = [cp_rand,inc_p_rand; inc_n_rand,cn_rand];
    %array of sensitivity and specificities in each instance of the matrix
    threshSens_random(threshold) = cp_rand / (cp_rand + inc_n_rand);
    threshSpecs_random(threshold) = 1 - (cn_rand / (cn_rand + inc_p_rand));
end



% Find best threshold level, discussed in class
% use the x^2 + y^2 = dist^2 [pythagorean theorem]
shortestDist = sqrt((threshSpecs(1)) ^ 2 + (threshSens(1) - 1) ^ 2);
for i = (2 : length(threshSens))
    currentDist = sqrt((threshSpecs(i)) ^ 2 + (threshSens(i) - 1) ^ 2);
    if currentDist < shortestDist
        shortestDist = currentDist;
        bestThresholdValue = i;
    end
end
bestThresholdValue%predictive value 

%Find best threshold for random data - same method
shortestDist_rand = sqrt((threshSpecs_random(1)) ^ 2 + (threshSens_random(1) - 1) ^ 2);
for i = (2 : length(threshSens_random))
    currentDist_rand = sqrt((threshSpecs_random(i)) ^ 2 + (threshSens_random(i) - 1) ^ 2);
    if currentDist_rand < shortestDist_rand
        shortestDist_rand = currentDist_rand;
        bestThresholdValue_rand = i;
    end
end
bestThresholdValue_rand%random value


%Calculating Necessary Values at Best Threshold Value

%same iterator complex
threshold = bestThresholdValue;
cp = 0;
cn = 0;
inc_p = 0;
inc_n = 0;
for j = (1 : length(prediction))
    if truth(j) == 1
        if prediction(j) >= threshold
            cp = cp + 1;
        else
            inc_n = inc_n + 1;
        end
    else
        if prediction(j) >= threshold
            inc_p = inc_p + 1;
        else
            cn = cn + 1;
        end
    end
end

%calculate and iterate the Random curve

threshold_r = bestThresholdValue_rand;
cp_r = 0;
cn_r = 0;
inc_p_r = 0;
inc_n_r = 0;
for j = (1 : length(random))
    if truth(j) == 1
        if random(j) >= threshold
            cp_r = cp_r + 1;
        else
            inc_n_r = inc_n_r + 1;
        end
    else
        if random(j) >= threshold
            inc_p_r = inc_p_r + 1;
        else
            cn_r = cn_r + 1;
        end
    end
end



%make the best confusion matrix
bestThresholdConfusionMatrix = [cp,inc_p; inc_n, cn]; %first is the correct prediction, then incorrect positive then incorrect negative and finally correct negative
bestThresholdConfusionMatrix

bestRandomMatrix = [cp_r,inc_p_r; inc_n_r,cn_r];
bestRandomMatrix

%ALL CALCS OF STATISTICS FOR BOTH PREDICTION AND RANDOM _r indicates random
%statistic
prevalance = (cp + inc_n) / (cp + inc_p + inc_n + cn);
prevalance
prevalence_random = (cp_r + inc_n_r)/(cp_r + inc_p_r + inc_n_r + cn_r)
prevalence_random
%tried to convert to percentage - did not work

accuracy = (cp + cn) / (cp + inc_p + inc_n + cn);
accuracy
accuracy_random = (cp_r + cn_r)/(cp_r + inc_n_r + inc_p_r + cn_r);
accuracy_random


sensitivity = cp / (cp + inc_n);
sensitivity

sensitivity_random = cp_r/(cp_r + inc_n_r);
sensitivity_random

specificity = cn / (cn + inc_p);
specificity

specificity_random = cn_r/ (cn_r + inc_p_r);
specificity_random

ppv = cp / (cp + inc_p);
ppv

ppv_r = cp_r/(cp_r+inc_p_r);
ppv_r

npv = cn / (cn + inc_n);
npv

npv_r = cn_r/(cn_r + inc_n_r);
npv_r

% Graph Specificities vs Sensitivities spec = tpf, 1-sens = fpf
figure();


hold off;
plot(threshSpecs, threshSens);
hold on;
plot(threshSpecs_random, threshSens_random,'Marker','x');


title('ROC Curve')
xlabel('1-Specificity or False Positive Fraction, FPF')
ylabel('Sensitivity or True Positive Fraction, TPF')
legend({'ROC Curve', 'Closest Point to (0,1)'})
axis([0 1 0 1])
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])

areaUnderROCCurve = -1*trapz(threshSpecs, threshSens) %trapeziodal sum not sure if this is it, double check there
areaUnderROCCurve_random = -1*trapz(threshSpecs_random,threshSens_random)