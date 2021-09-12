%%%%%%%%% Amaral(2021): Monetary Policy with Risky Assets

% This script plots figures 7, 8, 9, and 10

%% Set-up

isFinal = true;

if isFinal
    imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Partial Equilibrium/';
else
    imagesFolder = 'Power/Images/';
end

%% Parameterization

phi     = 1.0;

nSumPeriods = 10;

defLow  = 0.02;
defMed  = 0.08;
defHigh = 0.20;
rNRef   = 0.04;
rNLow   = 0.03;
rNHigh  = 0.05;
iotaRef = rNRef;

%% Monetary policy power reduction

iotaWithLowDefault  = NaN(nSumPeriods, 1);
iotaWithMedDefault  = NaN(nSumPeriods, 1);
iotaWithHighDefault = NaN(nSumPeriods, 1);
iotaWithoutDefault  = NaN(nSumPeriods, 1);

rNWithLowDefault    = NaN(nSumPeriods, 1);
rNWithMedDefault    = NaN(nSumPeriods, 1);
rNWithHighDefault   = NaN(nSumPeriods, 1);
rNWithoutDefault    = NaN(nSumPeriods, 1);

rNLowWithLowDefault    = NaN(nSumPeriods, 1);
rNLowWithMedDefault    = NaN(nSumPeriods, 1);
rNLowWithHighDefault   = NaN(nSumPeriods, 1);
rNLowWithoutDefault    = NaN(nSumPeriods, 1);

rNHighWithLowDefault    = NaN(nSumPeriods, 1);
rNHighWithMedDefault    = NaN(nSumPeriods, 1);
rNHighWithHighDefault   = NaN(nSumPeriods, 1);
rNHighWithoutDefault    = NaN(nSumPeriods, 1);

% Economy is compared to no-default steady state
for ii = 1:nSumPeriods
    denProdNo   = 1;
    denProdLow  = 1;
    denProdMed  = 1;
    denProdHigh = 1;
    for jj = 1:ii
        denProdNo    = denProdNo    * (1 + phi);
        denProdLow   = denProdLow   * (1 + (1-defLow)*phi);
        denProdMed   = denProdMed   * (1 + (1-defMed)*phi);
        denProdHigh  = denProdHigh  * (1 + (1-defHigh)*phi);
    end
    
    iotaWithoutDefault(ii)      = iotaRef * 1 / denProdNo;
    iotaWithLowDefault(ii)      = iotaRef * (1-defLow) / denProdLow;
    iotaWithMedDefault(ii)      = iotaRef * (1-defMed) / denProdMed;
    iotaWithHighDefault(ii)     = iotaRef * (1-defHigh) / denProdHigh;
    
    rNWithoutDefault(ii)      = rNRef * 1 / denProdNo;
    rNWithLowDefault(ii)      = rNRef * 1 / denProdLow;
    rNWithMedDefault(ii)      = rNRef * 1 / denProdMed;
    rNWithHighDefault(ii)     = rNRef * 1 / denProdHigh;
    
    % Simulating changes in rN
    rNLowWithoutDefault(ii)      = rNLow * 1 / denProdNo;
    rNLowWithLowDefault(ii)      = rNLow * 1 / denProdLow;
    rNLowWithMedDefault(ii)      = rNLow * 1 / denProdMed;
    rNLowWithHighDefault(ii)     = rNLow * 1 / denProdHigh;
    
    % Simulating changes in rN
    rNHighWithoutDefault(ii)      = rNHigh * 1 / denProdNo;
    rNHighWithLowDefault(ii)      = rNHigh * 1 / denProdLow;
    rNHighWithMedDefault(ii)      = rNHigh * 1 / denProdMed;
    rNHighWithHighDefault(ii)     = rNHigh * 1 / denProdHigh;
end

%% Effect on the steady-state price level
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [...
     rNWithoutDefault - iotaWithoutDefault, ...
     rNWithLowDefault - iotaWithLowDefault, ...
     rNWithMedDefault - iotaWithMedDefault, ...
     rNWithHighDefault - iotaWithHighDefault]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Effect on the steady-state price level per period ($\phi$ = ' num2str(phi,'%.1f') ')'], 'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['effectOnSteadyStatePriceLevel_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Cummulative effect on the steady-state price level
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [ ...
     cumsum(rNWithoutDefault) - cumsum(iotaWithoutDefault), ...
     cumsum(rNWithLowDefault) - cumsum(iotaWithLowDefault), ...
     cumsum(rNWithMedDefault) - cumsum(iotaWithMedDefault), ...
     cumsum(rNWithHighDefault) - cumsum(iotaWithHighDefault) ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Cummulative effect on the steady-state price level ($\phi$ = ' num2str(phi,'%.1f') ')'], ...
    'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['cummEffectOnSteadyStatePriceLevel_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Effect on the steady-state price level (rNLow)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [...
     rNLowWithoutDefault - iotaWithoutDefault, ...
     rNLowWithLowDefault - iotaWithLowDefault, ...
     rNLowWithMedDefault - iotaWithMedDefault, ...
     rNLowWithHighDefault - iotaWithHighDefault]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Effect on the steady-state price level per period ($\phi$ = ' num2str(phi,'%.1f') ')'], 'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['effectOnSteadyStatePriceLevel_rNLow_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Cummulative effect on the steady-state price level (rNLow)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [ ...
     cumsum(rNLowWithoutDefault) - cumsum(iotaWithoutDefault), ...
     cumsum(rNLowWithLowDefault) - cumsum(iotaWithLowDefault), ...
     cumsum(rNLowWithMedDefault) - cumsum(iotaWithMedDefault), ...
     cumsum(rNLowWithHighDefault) - cumsum(iotaWithHighDefault) ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Cummulative effect on the steady-state price level ($\phi$ = ' num2str(phi,'%.1f') ')'], ...
    'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['cummEffectOnSteadyStatePriceLevel_rNLow_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Effect on the steady-state price level (rNHigh)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [...
     rNHighWithoutDefault - iotaWithoutDefault, ...
     rNHighWithLowDefault - iotaWithLowDefault, ...
     rNHighWithMedDefault - iotaWithMedDefault, ...
     rNHighWithHighDefault - iotaWithHighDefault]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Effect on the steady-state price level per period ($\phi$ = ' num2str(phi,'%.1f') ')'], 'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['effectOnSteadyStatePriceLevel_rNHigh_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Cummulative effect on the steady-state price level (rNHigh)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [ ...
     cumsum(rNHighWithoutDefault) - cumsum(iotaWithoutDefault), ...
     cumsum(rNHighWithLowDefault) - cumsum(iotaWithLowDefault), ...
     cumsum(rNHighWithMedDefault) - cumsum(iotaWithMedDefault), ...
     cumsum(rNHighWithHighDefault) - cumsum(iotaWithHighDefault) ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['Cummulative effect on the steady-state price level ($\phi$ = ' num2str(phi,'%.1f') ')'], ...
    'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['cummEffectOnSteadyStatePriceLevel_rNHigh_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Effect on the steady-state price level (rN falls)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [...
     rNLowWithoutDefault - iotaWithoutDefault, ...
     rNLowWithLowDefault - iotaWithLowDefault, ...
     rNLowWithMedDefault - iotaWithMedDefault, ...
     rNLowWithHighDefault - iotaWithHighDefault, ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['$r^n$ falls by 1 p.p.: effect on the steady-state price level per period ($\phi$ = ' num2str(phi,'%.1f') ')'], 'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default $( r^n_t = \overline{\iota}_t )$ price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['effectOnSteadyStatePriceLevel_rNFalls_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Cummulative effect on the steady-state price level (rN falls)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [ ...
     cumsum(rNLowWithoutDefault) - cumsum(iotaWithoutDefault), ...
     cumsum(rNLowWithLowDefault) - cumsum(iotaWithLowDefault), ...
     cumsum(rNLowWithMedDefault) - cumsum(iotaWithMedDefault), ...
     cumsum(rNLowWithHighDefault) - cumsum(iotaWithHighDefault) ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['$r^n$ falls by 1 p.p.: cummulative effect on the steady-state price level ($\phi$ = ' num2str(phi,'%.1f') ')'], ...
    'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default $( r^n_t = \overline{\iota}_t )$ price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['cummEffectOnSteadyStatePriceLevel_rNFalls_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Effect on the steady-state price level (rN rises)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [...
     rNHighWithoutDefault - iotaWithoutDefault, ...
     rNHighWithLowDefault - iotaWithLowDefault, ...
     rNHighWithMedDefault - iotaWithMedDefault, ...
     rNHighWithHighDefault - iotaWithHighDefault, ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['$r^n$ rises by 1 p.p.: effect on the steady-state price level per period ($\phi$ = ' num2str(phi,'%.1f') ')'], 'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default $( r^n_t = \overline{\iota}_t )$ price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['effectOnSteadyStatePriceLevel_rNRises_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Cummulative effect on the steady-state price level (rN rises)
f = figure;

previousInterpreter = edu_GraphSetInterpreter('Latex');

p = plot(1:nSumPeriods, ...
    100 * [ ...
     cumsum(rNHighWithoutDefault) - cumsum(iotaWithoutDefault), ...
     cumsum(rNHighWithLowDefault) - cumsum(iotaWithLowDefault), ...
     cumsum(rNHighWithMedDefault) - cumsum(iotaWithMedDefault), ...
     cumsum(rNHighWithHighDefault) - cumsum(iotaWithHighDefault) ...
     ]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '-.';
p(4).LineStyle = '--';
set(gca,'FontSize', 12);
title(['$r^n$ rises by 1 p.p.: cummulative effect on the steady-state price level ($\phi$ = ' num2str(phi,'%.1f') ')'], ...
    'interpreter', 'latex');
xlabel('term period', 'interpreter', 'latex');
ylabel('\% dev. from no-default $( r^n_t = \overline{\iota}_t )$ price level', 'interpreter', 'latex');
legend('Without default risk', ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defLow)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defMed)], ...
    ['With constant $\mathcal{D}^{Policy}$ = ', num2str(defHigh)], ...
    'location', 'best', 'interpreter', 'latex');
xlim([0 nSumPeriods]);

edu_GraphSetInterpreter(previousInterpreter);

simGraphName = ['cummEffectOnSteadyStatePriceLevel_rNRises_phi', replace(num2str(phi),'.',''), '.png'];
set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);
