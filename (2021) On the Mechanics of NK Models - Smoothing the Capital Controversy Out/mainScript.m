%%%%%%%%% Amaral(2021): On the Mechanics of New-Keynesian Models: Smoothing the Capital Controversy Out
%%%%%%%%% Replication of Rupert and Sustek (2019) with the addition of Interest-rate-smoothing

clear all; clc;

pathData = 'C:\Users\Eduardo\OneDrive\MATLAB\My Thesis\Paper 7\Data';
pathTables = 'C:\Users\Eduardo\OneDrive\MATLAB\My Thesis\Paper 7\Tables';
pathImages = 'C:\Users\Eduardo\OneDrive\MATLAB\My Thesis\Paper 7\Images';

dynare smoothingCapital;

%% Sweep parameters

%%%%%%
% The steady-state capital and depreciation must be changed in the script .mod
%%%%%%

params_names = {'rrhoM', 'rrhoNr', 'kkSS'};
params_tex_names = {'$\rho^m$', '$\rho^i$', '$\overline{K}$'};
params_1 = 0:0.1:0.9;
params_2 = 0:0.1:0.9; %0:0.1:0.9;
params_3 = [5.5 15.4]; %15.4

matSign_r = NaN(length(params_1), length(params_2), length(params_3));
matSign_nr = NaN(length(params_1), length(params_2), length(params_3));
first_time = 1;
for iParam = 1:length(params_1)
    for jParam = 1:length(params_2)
        for kParam = 1:length(params_3)
       
            if first_time
                dynare smoothingCapital noclearall;
                first_time = 0;

                kkSS = get_param_by_name('kkSS');
                ddelta = get_param_by_name('ddelta');
                kkappa = get_param_by_name('kkappa');
            end

            set_param_value(params_names{1}, params_1(iParam));
            set_param_value(params_names{2}, params_2(jParam));
            set_param_value(params_names{3}, params_3(kParam));

            % Check Blanchard and Khan
            [eigVals, BK] = check(M_,options_,oo_);
            
            if BK
                disp(['Checking rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ]);
                %dynare smoothingCapital;
                [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);
                
                if info
                    disp(['Computation fails for rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ]);
                else
                    matSign_r(iParam, jParam, kParam) = sign(oo_.irfs.r_eepsM(1));
                    matSign_nr(iParam, jParam, kParam) = sign(oo_.irfs.nr_eepsM(1));
                end
            else
                disp(['Blanchard-Khan fails for rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ]);
                % 2 indicates that the BK consitions are not satisfied
                matSign_r(iParam, jParam, kParam) = 2;
                matSign_nr(iParam, jParam, kParam) = 2;
            end
            
        end
    end
end

%% Build table of parameter sweep

idx = 1;
kkSS = params_3(idx);

arrTable = matSign_r(:,:,idx);
arrTable = string(arrTable);
arrTable = strrep(arrTable, "-1", "-");
arrTable = strrep(arrTable, "1", "+");
arrTable = strrep(arrTable, "2", "");


t = array2table(arrTable, ...
    'VariableNames', string(strcat([params_tex_names{2} ' = '], num2str(params_2'))), ...
    'RowNames', string(strcat([params_tex_names{1} ' = '], num2str(params_1')))) ;

disp(t);

n_col = size(t,2);
tabWidth     = '1.15\\textwidth';
colAlignment = repmat('c',1,n_col);
edu_Table2Latex(t, [pathTables filesep 'ParameterSweep_capital_smoothing_' num2str(kkSS) '_' num2str(ddelta) '_' num2str(kkappa) '.tex'], ...
    'tabWidth', tabWidth, 'colAlignment', colAlignment);

%% Plot IRFs

varList         = {'nr', 'r', 'ppi', 'y', 'c', 'kk_lag'};
varTexList      = {'$\hat{i}_t$', '$\hat{r}_t$', '$\pi_t$', '$\hat{y}_t$', '$\hat{c}_t$', '$\hat{k}_t$'};

params_names = {'rrhoM', 'rrhoNr', 'kkSS', 'kkappa'};
params_tex_names = {'$\rho^m$', '$\rho^i$', '$\overline{K}$', '$\kappa$'};
params_1 = [0, 0.1, 0.5, 0.9];
params_2 = [0, 0.1, 0.5, 0.9];
params_3 = [5.5 15.4]; %15.4
params_4 = [0 0.1 0.2 0.5 0.9];

for iParam = 1:length(params_1)
    for jParam = 1:length(params_2)
        for kParam = 1:length(params_3)
            for mParam = 1:length(params_4)
            
                set_param_value(params_names{1}, params_1(iParam));
                set_param_value(params_names{2}, params_2(jParam));
                set_param_value(params_names{3}, params_3(kParam));
                set_param_value(params_names{4}, params_4(mParam));

                % Check Blanchard and Khan
                [eigVals, BK] = check(M_,options_,oo_);

                if BK
                    disp(['Checking rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ]);
                    %dynare smoothingCapital;
                    [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);

                    if info
                        disp(['Computation fails for rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ]);
                    else

                        f = figure;
                        edu_GraphSetInterpreter('latex');
                        nLins = 2;
                        nCols = ceil(length(varList)/nLins);
                        tl = tiledlayout(nLins, nCols);
                        tl.TileSpacing  = 'compact';
                        tl.Padding      = 'compact';

                        for iVar = 1:length(varList)
                            nexttile();
                            p = plot(oo_.irfs.([ varList{iVar} '_eepsM']));
                            p.LineWidth = 1.5;
                            edu_GraphDrawZeroAxis(p);
                            title(varTexList{iVar});
                            xlabel('');
                            set(gca, 'FontSize', 12);
                            
                            if mod(iVar, nCols) == 1
                               ylabel('\%'); 
                            end
                        end

                        graphName = ['IRF_' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ' and ' num2str(params_4(mParam)) '.png'];
                        set(f, 'Position',  [100, 0, 800, 500]); % resize figure
                        exportgraphics(f, [pathImages filesep 'IRFs' filesep graphName]);

                    end
                else
                    disp(['Blanchard-Khan fails for rho = ' num2str(params_1(iParam)) ' and ' num2str(params_2(jParam)) ' and ' num2str(params_3(kParam)) ' and ' num2str(params_4(mParam)) ]);
                    % 2 indicates that the BK consitions are not satisfied

                end
            end
        end
    end
end
            
%% Load empirical data

dataset = readtable([pathData filesep 'pwt100.xlsx'], 'Sheet', 'Data');

%% Read empirical data

y = dataset.rgdpna;
k = dataset.rnna;

kY = k ./ (y/4); % Capital to quarterly GDP

% Find specific countries
idx = strcmp(dataset.countrycode, 'USA');
kY_USA = k(idx) ./ (y(idx)/4) ;
idx = strcmp(dataset.countrycode, 'BRA');
kY_BRA = k(idx) ./ (y(idx)/4) ;
idx = strcmp(dataset.countrycode, 'DEU');
kY_DEU = k(idx) ./ (y(idx)/4) ;
idx = strcmp(dataset.countrycode, 'JPN');
kY_JPN = k(idx) ./ (y(idx)/4) ;

% Remove NaN
kY = rmmissing(kY);
kY_USA = rmmissing(kY_USA);
kY_BRA = rmmissing(kY_BRA);
kY_DEU = rmmissing(kY_DEU);
kY_JPN = rmmissing(kY_JPN);

f = figure;
binLimits = [prctile(kY, 5) prctile(kY, 95)];
p = histogram(kY, 'BinLimits', binLimits);
p.FaceAlpha = 0.4;
edu_GraphSetInterpreter('latex');
title('Histogram of capital-to-quarterly-output ratio');
xlabel('$\frac{K_t}{Y_t}$');
ylabel('frequency');
set(gca, 'FontSize', 14);

xl = xline(median(kY, 'omitnan'), '--r');
xl.LineWidth = 2;

set(f, 'Position',  [100, 0, 600, 400]); % resize figure
exportgraphics(f, [pathImages filesep 'Histogram.png']);

% Build table

tKY = table(dataset.year, dataset.country, dataset.rgdpna, dataset.rnna, dataset.rnna./(dataset.rgdpna./4), ...
    'VariableNames', {'Year', 'Country', 'Annual GDP', 'Capital Stock', '$\frac{\overline{K}}{\overline{Y}}$'});

tKY.("Annual GDP") = arrayfun(@(c) ThousandSep(c, 'US\\$ %.2f', ','),tKY.("Annual GDP"),'UniformOutput',false);
tKY.("Capital Stock") = arrayfun(@(c) ThousandSep(c, 'US\\$ %.2f', ','),tKY.("Capital Stock"),'UniformOutput',false);
tKY.("$\frac{\overline{K}}{\overline{Y}}$") = arrayfun(@(c) ThousandSep(round(c, 1), '%.1f', ','), tKY.("$\frac{\overline{K}}{\overline{Y}}$"), 'UniformOutput',false);

idx = ( ismember(tKY.('Year'), [1960, 2019] ) ) ;
idx = idx .* ( ismember(tKY.('Country'), {'Brazil', 'United States', 'Canada', 'United Kingdom', 'Germany', 'France', 'Japan', 'Mexico', 'South Africa', 'China'})) ;

tSelected = tKY(find(idx), :);

n_col = size(tSelected,2);
tabWidth     = '0.8\\textwidth';
colAlignment = 'clrrr';
edu_Table2Latex(tSelected, [pathTables filesep 'CapitalToOutput_SelectedCountries' '.tex'], ...
    'tabWidth', tabWidth, 'colAlignment', colAlignment);

%% Undetermined coefficients

% Return parameter values to their benchmark values
dynare smoothingCapital;

syms a0 a1 a2 b0 b1 b2 d0 d1 d2 f0 f1 f2;
syms nnu aalpha eeta bbeta ddelta ttheta rrhoM rrhoNr ppsi kkappa kkappaLine;
syms ySS cSS kkSS rSS cchiSS;

eqns_manual = [
                -a0 == -a0*f0 + nnu*b0*(1-rrhoNr)*(1 - a2 - b2) - b0*f0 ;
                -a1 == -a0*f1 - a1*rrhoM + b1*(1-rrhoNr)*nnu*(1 - a2 - b2) - a2 - b0*f1 - b1*rrhoM - b2 + 1  ;
                -a2 == -a0*f2 - a2*rrhoNr - a2*(1-rrhoNr)*nnu*b2 + rrhoNr + (1-rrhoNr)*nnu*b2 - b0*f2 - b2*rrhoNr - b2*(1-rrhoNr)*nnu*b2 ;
                
                0 == - a0 + (1-rSS)*a0*f0 - kkappaLine - kkappaLine*f0^2 + 2*kkappaLine*f0 - rSS*(1 + eeta)/(1 - aalpha)*d0*f0 + rSS*(1+aalpha*eeta)/(1-aalpha)*f0 - kkappaLine*f2*(1-rrhoNr)*nnu*b0 - rSS*(1+eeta)/(1-aalpha)*(1-rrhoNr)*nnu*b0;
                0 == - a1 + (1-rSS)*(a0*f1 + a1*rrhoM) - rSS*(1+eeta)/(1-aalpha)*(d0*f1 + d1*rrhoM + d2*(1-rrhoNr)*nnu*b1 + d2) - kkappaLine*(f0*f1 + rrhoM*f1 - 2*f1 + f2*(1-rrhoNr)*nnu*b1 + f2 ) + rSS*(1+aalpha*eeta)/(1-aalpha)*f1 ;
                0 == - a2 + (1-rSS)*(a0*f2 + a2*(rrhoNr + (1-rrhoNr)*nnu*b2)) - rSS*(1+eeta)/(1-aalpha)*(d0*f2 + d2*(rrhoNr + (1-rrhoNr)*nnu*b2)) + rSS*(1+aalpha*eeta)/(1-aalpha)*f2 - kkappaLine*(f0*f2 + f2*(rrhoNr + (1-rrhoNr)*nnu*b2) - 2*f2) ;
                
                b0 == ppsi*(eeta + aalpha)/(1-aalpha)*d0 - ppsi*(aalpha*eeta+aalpha)/(1-aalpha) + ppsi*a0 + bbeta*b0*f0 + bbeta*b2*(1-rrhoNr)*nnu*b0;
                b1 == ppsi*(eeta + aalpha)/(1-aalpha)*d1 + ppsi*a1 + bbeta*b0*f1 + bbeta*b1*rrhoM + bbeta*b2*((1-rrhoNr)*nnu*b1 + 1) ;
                b2 == ppsi*( (eeta+aalpha)/(1-aalpha)*d2 + a2 ) + bbeta*b0*f2 + bbeta*b2*(rrhoNr + (1-rrhoNr)*nnu*b2 ) ;
                
                d0 == cSS/ySS*a0 + kkSS/ySS*f0 - (1 - ddelta)*kkSS/ySS ;
                d1 == cSS/ySS*a1 + kkSS/ySS*f1 ;              
                d2 == cSS/ySS*a2 + kkSS/ySS*f2 ;
];

% ELIMINATE STATE VARIABLE i_{t-1}
%eqns = subs(eqns, [a2 b2 d2 f2], [0 0 0 0]);

eqns_manual_0 = eqns_manual([1,4,7,10]);
eqns_manual_1 = eqns_manual([2,5,8,11]);
eqns_manual_2 = eqns_manual([3,6,9,12]);

% LOAD EQUATIONS
eqns_stepByStep;

% COMPARE MANUAL WITH PROG
comp = [
    isequal(expand(eqns_manual(1)), expand(eqns_0(1)));
    isequal(expand(eqns_manual(2)), expand(eqns_1(1)));
    isequal(expand(eqns_manual(3)), expand(eqns_2(1)));

    isequal(expand(eqns_manual(4)) + expand(eqns_0(2)), expand(eqns_0(2)) + expand(eqns_manual(4)));
    isequal(expand(eqns_manual(5)) + expand(eqns_1(2)), expand(eqns_1(2)) + expand(eqns_manual(5)));
    isequal(expand(eqns_manual(6)) + expand(eqns_2(2)), expand(eqns_2(2)) + expand(eqns_manual(6)));

    isequal(expand(eqns_manual(7)), expand(eqns_0(3)));
    isequal(expand(eqns_manual(8)), expand(eqns_1(3)));
    isequal(expand(eqns_manual(9)), expand(eqns_2(3)));

    isequal(expand(eqns_manual(10)), expand(eqns_0(4)));
    isequal(expand(eqns_manual(11)), expand(eqns_1(4)));
    isequal(expand(eqns_manual(12)), expand(eqns_2(4)))
]

% % Solve by order
% eqns_soln = eqns_0;
% eqns_soln(1) = isolate(eqns_0(1), a0); % isolate a0 in (1)
% eqns_soln(2) = isolate(eqns_0(2), d0); % isolate d0 in (1)
% eqns_soln(3) = isolate(eqns_0(3), b0); % isolate a0 in (1)
% eqns_soln(4) = isolate(eqns_0(4), f0); % isolate f0 in (4)
% 
% eqns_soln(1) = subs(eqns_soln(1), f0, rhs(eqns_soln(4))); % replace f0 in (1)
% eqns_soln(2) = subs(eqns_soln(2), f0, rhs(eqns_soln(4))); % replace f0 in (2)
% eqns_soln(3) = subs(eqns_soln(3), f0, rhs(eqns_soln(4))); % replace f0 in (3)
% 
% eqns_soln(1) = subs(eqns_soln(1), b0, rhs(eqns_soln(3))); % replace b0 in (1)
% eqns_soln(2) = subs(eqns_soln(2), b0, rhs(eqns_soln(3))); % replace b0 in (2)
% 
% % Isolate d0
% eqns_soln(2) = isolate(eqns_soln(2), d0); % isolate d0 in (2)
% 
% % % Find a0 (only numerical solution)
% eqns_soln(1) = subs(eqns_soln(1), d0, rhs(eqns_soln(2)));  % replace d0 in (1)


% Numerical
aalpha_val = get_param_by_name('aalpha');
bbeta_val = get_param_by_name('bbeta');
ddelta_val = get_param_by_name('ddelta');
eeta_val = get_param_by_name('eeta');
ttheta_val = get_param_by_name('ttheta');
rrhoM_val = get_param_by_name('rrhoM');
rrhoNr_val = get_param_by_name('rrhoNr');
nnu_val = get_param_by_name('nnu');
kkappa_val = get_param_by_name('kkappa');

ySS_val     = get_param_by_name('ySS');
kkSS_val     = get_param_by_name('kkSS');


% assume(a0,'clear')
% assume(f1,'clear')

% assume(bbeta>0 & bbeta<1);
% assumeAlso(ddelta>0 & ddelta<1);
% assumeAlso(aalpha>0 & aalpha<1);
% assumeAlso(eeta>0);
% assumeAlso(ppsi>0);
% assumeAlso(rrhoM>=0 & rrhoM<1);
% assumeAlso(nnu>1);
% assumeAlso(ySS>0 & cSS>0 & kkSS>0 & rSS>0);

%%%%%%%%%%%%%%%%%%%%%%% NO ADJUSTMENT COSTS and NO SMOOTHING

kkappa_val = 0;
rrhoNr_val = 0;

kkappaLine_val      = kkappa_val*kkSS_val;
rSS_val             = 1/bbeta_val - 1 + ddelta_val;
cSS_val             = ySS_val - ddelta_val*kkSS_val;
wSS_val             = cSS_val*(ySS_val/(kkSS_val^aalpha_val))^(eeta_val/(1-aalpha_val));
cchiSS_val          = wSS_val/(1-aalpha_val)*(ySS_val/kkSS_val)^(aalpha_val/(1-aalpha_val));
ppsi_val            = cchiSS_val*(1-ttheta_val)*(1-ttheta_val*bbeta_val)/ttheta_val;
varVals = { 
            aalpha      aalpha_val
            bbeta       bbeta_val
            ddelta      ddelta_val
            eeta        eeta_val
            ttheta      ttheta_val
            ppsi        ppsi_val
            
            nnu         nnu_val
            kkSS        kkSS_val
            rSS         rSS_val
            
            kkappaLine  kkappaLine_val
            cchiSS      cchiSS_val
            ySS         ySS_val
            cSS         cSS_val
            
            %rrhoM       rrhoM_val
            rrhoNr       rrhoNr_val
            };

solveUC_plotGraph;


%%%%%%%%%%%%%%%%%%%%%%% WITH ADJUSTMENT COSTS and NO SMOOTHING

kkappa_val = 0.1;
rrhoNr_val = 0;

kkappaLine_val      = kkappa_val*kkSS_val;
rSS_val             = 1/bbeta_val - 1 + ddelta_val;
cSS_val             = ySS_val - ddelta_val*kkSS_val;
wSS_val             = cSS_val*(ySS_val/(kkSS_val^aalpha_val))^(eeta_val/(1-aalpha_val));
cchiSS_val          = wSS_val/(1-aalpha_val)*(ySS_val/kkSS_val)^(aalpha_val/(1-aalpha_val));
ppsi_val            = cchiSS_val*(1-ttheta_val)*(1-ttheta_val*bbeta_val)/ttheta_val;
varVals = { 
            aalpha      aalpha_val
            bbeta       bbeta_val
            ddelta      ddelta_val
            eeta        eeta_val
            ttheta      ttheta_val
            ppsi        ppsi_val
            
            nnu         nnu_val
            kkSS        kkSS_val
            rSS         rSS_val
            
            kkappaLine  kkappaLine_val
            cchiSS      cchiSS_val
            ySS         ySS_val
            cSS         cSS_val
            
            %rrhoM       rrhoM_val
            rrhoNr       rrhoNr_val
            };

solveUC_plotGraph;


%%%%%%%%%%%%%%%%%%%%%%% WITH NO ADJUSTMENT COSTS and WITH SMOOTHING

kkappa_val = 0.0;
rrhoNr_val = 0.5;
%rrhoM_val  = 0.0;
%kkSS_val   = 5.5;

kkappaLine_val      = kkappa_val*kkSS_val;
rSS_val             = 1/bbeta_val - 1 + ddelta_val;
cSS_val             = ySS_val - ddelta_val*kkSS_val;
wSS_val             = cSS_val*(ySS_val/(kkSS_val^aalpha_val))^(eeta_val/(1-aalpha_val));
cchiSS_val          = wSS_val/(1-aalpha_val)*(ySS_val/kkSS_val)^(aalpha_val/(1-aalpha_val));
ppsi_val            = cchiSS_val*(1-ttheta_val)*(1-ttheta_val*bbeta_val)/ttheta_val;
varVals = { 
            aalpha      aalpha_val
            bbeta       bbeta_val
            ddelta      ddelta_val
            eeta        eeta_val
            ttheta      ttheta_val
            ppsi        ppsi_val
            
            nnu         nnu_val
            kkSS        kkSS_val
            rSS         rSS_val
            
            kkappaLine  kkappaLine_val
            cchiSS      cchiSS_val
            ySS         ySS_val
            cSS         cSS_val
            
            %rrhoM       rrhoM_val
            rrhoNr       rrhoNr_val
            };

solveUC_plotGraph;



%%%%%%%%%%%%%%%%%%%%%%% WITH SOME ADJUSTMENT COSTS and WITH SMOOTHING

kkappa_val = 0.1;
rrhoNr_val = 0.5;
%rrhoM_val  = 0.0;
%kkSS_val   = 5.5;

kkappaLine_val      = kkappa_val*kkSS_val;
rSS_val             = 1/bbeta_val - 1 + ddelta_val;
cSS_val             = ySS_val - ddelta_val*kkSS_val;
wSS_val             = cSS_val*(ySS_val/(kkSS_val^aalpha_val))^(eeta_val/(1-aalpha_val));
cchiSS_val          = wSS_val/(1-aalpha_val)*(ySS_val/kkSS_val)^(aalpha_val/(1-aalpha_val));
ppsi_val            = cchiSS_val*(1-ttheta_val)*(1-ttheta_val*bbeta_val)/ttheta_val;
varVals = { 
            aalpha      aalpha_val
            bbeta       bbeta_val
            ddelta      ddelta_val
            eeta        eeta_val
            ttheta      ttheta_val
            ppsi        ppsi_val
            
            nnu         nnu_val
            kkSS        kkSS_val
            rSS         rSS_val
            
            kkappaLine  kkappaLine_val
            cchiSS      cchiSS_val
            ySS         ySS_val
            cSS         cSS_val
            
            %rrhoM       rrhoM_val
            rrhoNr       rrhoNr_val
            };

solveUC_plotGraph;
