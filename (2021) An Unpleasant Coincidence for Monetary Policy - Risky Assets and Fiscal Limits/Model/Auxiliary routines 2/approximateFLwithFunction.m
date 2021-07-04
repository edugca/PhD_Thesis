%% Fit the Fiscal Limit Distribution to a logistic and a linear function

% Set-up simulation
nFLSims       = 1000;
randList    = randn(nFLSims);

% Get data
stateComb   = muComb;
muFL        = fVals_mu;
stdFL       = fVals_std;
nFL         = length(muFL);

% Build sample
sampleFL = [];
for iComb = 1:nFL
    distFL    = muFL(iComb) + stdFL(iComb)*randn(nFLSims, 1);
    for iSorted = 1:nFLSims
        probFL = normcdf(distFL, muFL(iComb), stdFL(iComb));
    end
    
    sampleFL = [sampleFL; [probFL, distFL, repmat(muFL(iComb), nFLSims, 1), repmat(stateComb(iComb,:), nFLSims, 1)]];
end

if strcmp(approximateDefaultProb, 'logistic')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOGISTIC
    % Fit non-linear model
    y = sampleFL(:,1);
    X = sampleFL(:,2:end);
    testfun = @(b,x) 1./(1 + exp( b(1) + b(2)*(x(:,1) - x(:,2)) + b(3)*x(:,3) + b(4)*x(:,4) ));
    beta0 = [0.5 0.5 0.5 0.5];
    mdlFL_logistic = fitnlm(X,y,testfun,beta0);
    %mdlFL_logistic.disp

    % Compare values
    comparisonValues = [testfun(mdlFL_logistic.Coefficients.Estimate,X), fFiscalLimit_defProb(X(:,1), X(:,3),0,0,X(:,4),0)];
    diff = comparisonValues(:,1)-comparisonValues(:,2);
    disp(['Mean difference of the fiscal limit approximation with a logistic function:' num2str(mean(diff))]);

    % Save
    save([pathSaved filesep 'fiscalLimit_approxLogistic.mat'], 'mdlFL_logistic');
    
elseif strcmp(approximateDefaultProb, 'linear')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR
    % Fit non-linear model
    y = sampleFL(:,1);
    X = sampleFL(:,2:end);
    testfun = @(b,x) b(1) + b(2)*(x(:,1) - x(:,2)) + b(3)*x(:,3) + b(4)*x(:,4) ;
    beta0 = [0.5 0.5 0.5 0.5];
    mdlFL_linear = fitlm([X(:,1) - X(:,2), X(:,3), X(:,4)],y);
    %mdlFL_logistic.disp

    % Compare values
    comparisonValues = [testfun(mdlFL_linear.Coefficients.Estimate,X), fFiscalLimit_defProb(X(:,1), X(:,3),0,0,X(:,4),0)];
    diff = comparisonValues(:,1)-comparisonValues(:,2);
    disp(['Mean difference of the fiscal limit approximation with a logistic function:' num2str(mean(diff))]);

    % Save
    save([pathSaved filesep 'fiscalLimit_approxLinear.mat'], 'mdlFL_linear');
    
end