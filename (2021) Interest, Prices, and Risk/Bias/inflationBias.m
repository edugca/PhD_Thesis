%% INFLATION BIAS

%% Housekeeping
clear;
clc;
close all;

pathMain    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) Monetary Policy with Risky Assets/Bias', 'C:\');
pathTables  = [pathMain filesep 'Tables'];

%% Price-level Bias

rN      = 0.04;
ddelta  = 0.60;

fPLBias = @(def, phiPii, rN, ddelta) (def.*ddelta)./((1-def.*ddelta).*phiPii) .* (1 + rN) ;

defVec      = linspace(0, 0.1, 5);
ddeltaVec   = ones(size(defVec)) .* ddelta;
phiPiiVec   = [0.2 0.5 1.0 1.5 2.0];
fMat        = NaN(length(defVec), length(phiPiiVec));

for iDef = 1:length(defVec)
    for iPhiPii = 1:length(phiPiiVec)
        
        % Check determinacy
        if phiPiiVec(iPhiPii) > 0
            fMat(iDef, iPhiPii) = fPLBias(defVec(iDef), phiPiiVec(iPhiPii), rN, ddelta) ;
        else
            fMat(iDef, iPhiPii) = NaN;
        end
    end
end

t = array2table(arrayfun(@(x) num2str(round(x .* 100, 2), '%9.1f'), fMat, 'UniformOutput', false),...
    'VariableNames', arrayfun(@(x) ['$\phi=$ ' num2str(x, '%9.1f')], phiPiiVec, 'UniformOutput', false), ...
    'RowNames', arrayfun(@(x) ['$\mathcal{D}=$ ' num2str(x .* 100, '%9.1f'), '\%'], defVec, 'UniformOutput', false));
disp(t);

colAlignment = '|ccccccccc';
tabWidth = '0.8\\textwidth';
colNames = strjoin(t.Properties.VariableNames,' & ');
edu_Table2Latex(t, [pathTables filesep 't_priceLevelBias_delta_' num2str(ddelta) '.tex'], 'colNames', colNames, 'colAlignment', colAlignment, 'tabWidth', tabWidth);

%% Inflation Bias

rN      = 0.04;
ddelta  = 0.60;

fInfBias = @(def, phiPii, rN, ddelta) (def.*ddelta)./((1-def.*ddelta).*phiPii - 1) .* (1 + rN) ;

defVec      = linspace(0, 0.1, 5);
ddeltaVec   = ones(size(defVec)) .* ddelta;
phiPiiVec   = [1.2 1.5 2.0 2.5 3.0];
fMat        = NaN(length(defVec), length(phiPiiVec));

for iDef = 1:length(defVec)
    for iPhiPii = 1:length(phiPiiVec)
        
        % Check determinacy
        if phiPiiVec(iPhiPii) > 1/(1-defVec(iDef).*ddeltaVec(iDef))
            fMat(iDef, iPhiPii) = fInfBias(defVec(iDef), phiPiiVec(iPhiPii), rN, ddelta) ;
        else
            fMat(iDef, iPhiPii) = NaN;
        end
    end
end

t = array2table(arrayfun(@(x) num2str(round(x .* 100, 2), '%9.1f'), fMat, 'UniformOutput', false),...
    'VariableNames', arrayfun(@(x) ['$\phi^{\pi}=$ ' num2str(x, '%9.1f')], phiPiiVec, 'UniformOutput', false), ...
    'RowNames', arrayfun(@(x) ['$\mathcal{D}=$ ' num2str(x .* 100, '%9.1f'), '\%'], defVec, 'UniformOutput', false));
disp(t);

colAlignment = '|ccccccccc';
tabWidth = '0.8\\textwidth';
colNames = strjoin(t.Properties.VariableNames,' & ');
edu_Table2Latex(t, [pathTables filesep 't_inflationBias_delta_' num2str(ddelta) '.tex'], 'colNames', colNames, 'colAlignment', colAlignment, 'tabWidth', tabWidth);