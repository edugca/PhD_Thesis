%% Calibrate debt level to a certain default probability

% Set-up simulation
nFLSims       = 1000;
randList    = randn(nFLSims);

% Get data
stateComb   = muComb;
muFL        = fVals_mu;
stdFL       = fVals_std;
nFL         = length(muFL);

% Build sample
%bCandidate = simStruct.params.bYBar;
bCandidate = debtLevelGuess;
defMean = 1.01;
tol = 0.001;
while abs(defMean - defPolTarget) > tol
    disp(['Testing debt to default probability: ' num2str(bCandidate)]);
    
    sampleDef = [];
    for iComb = 1:nFL
        distFL    = (bCandidate*simStruct.params.yBar) + stdFL(iComb)*randn(nFLSims, 1);
        probFL    = normcdf(distFL, muFL(iComb), stdFL(iComb));

        sampleDef = [sampleDef; probFL];
    end
    defMean = mean(sampleDef);
    
    if abs(defMean - defPolTarget) > tol && defMean - defPolTarget > tol
       bCandidate = bCandidate - 0.001;
    elseif abs(defMean - defPolTarget) > tol && defPolTarget - defMean > tol
        bCandidate = bCandidate + 0.001;
    end
end

debtTargetDefProb = bCandidate;
gammaTauTargetDefProb = fzero(@(x) log(debtTargetDefProb/simStruct.params.yBar)/log(simStruct.params.bYBar) - x/simStruct.params.gammaTau, simStruct.params.gammaTau);
