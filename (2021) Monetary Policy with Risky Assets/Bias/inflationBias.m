%% INFLATION BIAS

%% Housekeeping
clear;
clc;
close all;

pathMain    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) Monetary Policy with Risky Assets/Bias', 'C:\');
pathTables  = [pathMain filesep 'Tables'];

%% Inflation Bias

rN      = 0.04;
iBad    = 0.4 * rN;

fInfBias = @(def, phiPii, rN, iBad) def./((1-def).*phiPii - 1) .* (rN -iBad) ;

defVec      = linspace(0, 0.1, 5);
phiPiiVec   = [1.2 1.5 2.0 2.5 3.0];
fMat        = NaN(length(defVec), length(phiPiiVec));

for iDef = 1:length(defVec)
    for iPhiPii = 1:length(phiPiiVec)
        
        % Check determinacy
        if phiPiiVec(iPhiPii) > 1/(1-defVec(iDef))
            fMat(iDef, iPhiPii) = fInfBias(defVec(iDef), phiPiiVec(iPhiPii), rN, iBad) ;
        else
            fMat(iDef, iPhiPii) = NaN;
        end
    end
end

t = array2table(arrayfun(@(x) num2str(round(x .* 100, 2), 2), fMat, 'UniformOutput', false),...
    'VariableNames', arrayfun(@(x) ['$\phi^{\pi}=$ ' num2str(x, '%9.2f')], phiPiiVec, 'UniformOutput', false), ...
    'RowNames', arrayfun(@(x) ['$\delta=$ ' num2str(x .* 100, '%9.2f'), '\%'], defVec, 'UniformOutput', false));

colAlignment = '|ccccccccc';
tabWidth = '0.8\\textwidth';
colNames = strjoin(t.Properties.VariableNames,' & ');
edu_Table2Latex(t, [pathTables filesep 't_inflationBias.tex'], 'colNames', colNames, 'colAlignment', colAlignment, 'tabWidth', tabWidth);