%% ESTIMATE MODEL

% Model
listMdls = [1];
estMdl      = mdlVector(listMdls);
modelName   = {riseModelNames{listMdls}};
setUpName   = {riseSetUpNames{listMdls}};
policyName  = {riseConfigVector(listMdls).conf_policyRule};
defProcess  = {riseConfigVector(listMdls).conf_defRProcessType};

% Estimated parameters
estimParamsNames = cell(length(listMdls), 1);
for ii = 1:length(listMdls)
    estimParamsNames{ii} = {estMdl(ii).estimation.priors.name};
end

%% Plot priors

plot_priors(estMdl);

%% Compute the bounds from the mean and standard deviations from dynare
prob=0.95;
Params={
	'sigmaShocks'  ,  'inv_gamma'     ,    0.01, 0.1	  
	'rhoShocks'  ,  'beta'     ,    0.5, 0.15
	};
npars=size(Params,1);
bounds=nan(npars,2);
for ii=1:npars
    bounds(ii,:)=distributions.find_bounds(Params{ii,2},Params{ii,3},Params{ii,4},prob);
end
disp(bounds)

%% Estimating the model

% Set initial values
% if exist('modeResults', 'var')
%    estMdl = set(estMdl, 'estim_start_vals', modeResults);
% end

% One thing that can often occur, especially in regime switching models,
% is that the initial covariance matrix at which to initialize the filter
% does not exist. In this case, one thing that can help is to use a
% diffuse prior, which is implemented in RISE under option kf_init_variance.
% You have to set this option to some arbitrary number.
% To not use diffuse prior, set kf_init_variance = [];
kf_init_variance        = 5;
kf_presample            = 16;
estim_start_date        = obs2date(firstFullDate_date,1); % retrieve the date of the 1st observation
estim_start_from_mode   = true;
solve_check_stability   = true;

optimizerName = 'fmincon'; %fmincon, fminsearch, bee_gate, particleswarm
opt_useParallel = false;
StepTolerance = 1e-16;
opt_maxFunEvals = 50000;

if ismember(optimizerName, {'fmincon', 'fminsearch'})
    opt = optimset(optimizerName);
    opt.Algorithm = 'interior-point'; %'interior-point'
    opt.UseParallel = opt_useParallel;
    opt.MaxFunEvals = opt_maxFunEvals;
    opt.TolFun      = StepTolerance;
    opt.TolX        = StepTolerance;
    opt.Diagnostics = 'off';
elseif ismember(optimizerName, {'particleswarm'})
    opt = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon);
    opt.UseParallel = opt_useParallel;
    opt.MaxIterations = opt_maxFunEvals;
elseif ismember(optimizerName, {'surrogateopt'})
    opt = optimoptions(optimizerName);
    opt.UseParallel = opt_useParallel;
    opt.MaxFunctionEvaluations = opt_maxFunEvals;
else
    opt = optimoptions(optimizerName);
    opt.UseParallel = opt_useParallel;
end

%global fitModel_mu
%fitModel_mu = fitModel_mu;
%addpath('C:\Users\Eduardo\OneDrive\MATLAB\Resources\RISE_toolbox-master\classes\models\@dsge\');
[estMdl, filtration] = estimate(estMdl, ...
                            'optimizer', optimizerName, ...
                            'kf_presample', kf_presample, ...
                            'kf_init_variance', kf_init_variance, ...
                            'estim_start_date', estim_start_date, ...
                            'estim_start_from_mode', estim_start_from_mode, ...
                            'solve_check_stability', solve_check_stability, ...
                            'optimset', opt);

create_containers = @(n)arrayfun(@(x)containers.Map(), 1:n, 'UniformOutput', false);
mapEstimParams = create_containers(length(listMdls));
                    
for iMap = 1:length(listMdls)
    mapEstimParams{iMap} = containers.Map(estimParamsNames{iMap}, ...
        estMdl(iMap).estimation.posterior_maximization.mode);
end
                       
                        
% <--- mdlDSGE=mdlDSGE.estimate;
%[mdlDSGE, filtration] = estimate(mdlDSGE, 'optimizer', 'fmincon', 'optimset', opt, 'estim_start_from_mode', true); % <--- mdlDSGE=mdlDSGE.estimate;

% m_new=estimate(m_old,'optimizer','fmincon','estim_start_from_mode', true)

%% Model comparison

mdlComparison = estMdl;

modelComparison;

%% Update priors
%%%% UNDER CONSTRUCTION

%mdlPriors = estMdl.estimation.priors;

%estMdl = set(estMdl,'estim_priors', mdlPriors);

%% Print and save results

% Print results
estMdl.print_estimation_results;

modeResults = cell2struct(...
    num2cell(estMdl.estimation.posterior_maximization.mode), ...
    {estMdl.estimation.priors.name}', 1);

iEstMdl = iEstMdl + 1;
estMdl_Final{iEstMdl} = estMdl;

%% Do posterior simulation

% Initiate randomizer
rng(1900)

myParPool = gcp;
for iMdl = 1:length(listMdls)
    [objective,lb,ub,x0,SIG] = pull_objective(estMdl(iMdl));

    SIG = utils.cov.nearest(SIG);

    nChains_mcmc = 8;
    draws_mcmc = 200000; % number of effective parameter draws through MCMC.
    ndraws_burnin = floor(0.1*draws_mcmc);

    mcmc_options=struct('burnin',ndraws_burnin,'N',draws_mcmc,'thin',1,...
        'nchain', nChains_mcmc);
    postSimResults.(['Model_' num2str(iMdl)]) = mh_sampler(objective,lb,ub,mcmc_options,x0,SIG);

    for iChain = 1:nChains_mcmc
        postSimResults.(['Model_' num2str(iMdl)]){iChain}.stats
    end
end
    
%% Plot priors, posteriors, priors and posteriors

%plot_priors(estMdl);

%plot_posteriors(estMdl, Results);
%plot_priors_and_posteriors(estMdl, postSimResults.Model_1,['alphaG'],true)
% [draw,m1]=draw_parameter(estMdl,postSimResults.Model_1{1}.pop)

tEstimation = struct();

paramsOrder = {
    '$\sigma$',
    ...%'$\chi$',
    '$\alpha_G$',
    '$\gamma_{G\Psi}$',
    '$\phi$',
    '$\phi^Y$',
    %'$\phi^{dY}$',
    '$\phi^i$',
    %'$\overline{\mathcal{D}^{\ast}} - \overline{\mathcal{D}}$',
    '$\rho^A$',
    '$\rho^\beta$',
    '$\rho^{GY}$',
    '$\rho^{GG}$',
    '$\rho^{\mathcal{M}}$',
    %'$\rho^{\mathcal{D}}$',
    '$\sigma^A$',
    '$\sigma^\beta$',
    '$\sigma^{G}$',
    '$\sigma^M$',
    %'$\sigma^{\mathcal{D}}$',
    '$\sigma^{me,Y}$'
    %'$\sigma^{me,\mathcal{D}}$',
};

for iMdl = 1:length(listMdls)
    estimPostSimData = plot_priors_and_posteriors(estMdl(iMdl), postSimResults.(['Model_' num2str(iMdl)]));

    nEstimParams = length(estimParamsNames{iMdl});
    nCols = 4;
    nLins = min(4, ceil(nEstimParams / nCols));
    
    mdlPriorsTexName    = {estMdl(iMdl).estimation.priors.tex_name};
    mdlPriorsDist       = {estMdl(iMdl).estimation.priors.prior_distrib};
    mdlPriorsLQ         = {estMdl(iMdl).estimation.priors.lower_quantile};
    mdlPriorsUQ         = {estMdl(iMdl).estimation.priors.upper_quantile};
    mdlPriorsProb         = {estMdl(iMdl).estimation.priors.prior_prob};
    
    tEstimation.(['Model_' num2str(iMdl)]) = cell2table(cell(nEstimParams, 8));
    tEstimation.(['Model_' num2str(iMdl)]).Properties.VariableNames = {'Parameter', 'Prior Dist.', 'CI Min', 'CI Max', '\% CI', 'Post. Mode', 'Post. Mean', 'Post. Std.'};
    tEstimation.(['Model_' num2str(iMdl)]).Properties.RowNames      = mdlPriorsTexName';
    
    tEstimation.(['Model_' num2str(iMdl)]){:,'Parameter'}           = mdlPriorsTexName';
    tEstimation.(['Model_' num2str(iMdl)]){:,'Prior Dist.'}         = strrep(mdlPriorsDist', '_', '. ');
    tEstimation.(['Model_' num2str(iMdl)]){:,'CI Min'}              = mdlPriorsLQ';
    tEstimation.(['Model_' num2str(iMdl)]){:,'CI Max'}              = mdlPriorsUQ';
    tEstimation.(['Model_' num2str(iMdl)]){:,'\% CI'}               = mdlPriorsProb';
    
    nFigures = 0;
    for iParam=1:nEstimParams
        
        % Create new figure
        if mod(iParam, nLins*nCols) == 1
           
            nFigures = nFigures + 1;
            f = figure;
            previousInterpreter = pub_GraphSetInterpreter('latex');
            
            if nEstimParams - iParam + 1 >= nLins*nCols
                nEffectiveLins = nLins;
            else
                nEffectiveLins = ceil(mod(nEstimParams - iParam + 1, nLins*nCols) / nCols);
            end
            
            % Set layout and reduce empty spaces in the graph
            tl = tiledlayout(nEffectiveLins, nCols);
            tl.TileSpacing  = 'compact';
            tl.Padding      = 'compact';
            
        end
        
        tl_h = nexttile();

        pStruct = estimPostSimData.(estimParamsNames{iMdl}{iParam});
        x = [pStruct.x_prior pStruct.x_kdens'];
        y = [pStruct.f_prior pStruct.f_kdens'];
        p = plot(x, y);
        p(1).LineWidth = 2;
        p(2).LineWidth = 2;
        p(2).LineStyle = ':';
        
        lPostModeSim = xline(pStruct.post_mode_sim, '--g');
        lPostModeSim.LineWidth = 2;
        lMeanSim = xline(pStruct.mean_sim, '--k');
        lMeanSim.LineWidth = 2;
        
        title(pStruct.tex_name{1});
        
        %legend(pStruct.tex_name, 'Location', 'northoutside');
        
        
        % Fill in table
        tEstimation.(['Model_' num2str(iMdl)]){iParam, 'Post. Mode'}   = {pStruct.post_mode_sim};
        tEstimation.(['Model_' num2str(iMdl)]){iParam, 'Post. Mean'}   = {pStruct.mean_sim};
        tEstimation.(['Model_' num2str(iMdl)]){iParam, 'Post. Std.'}   = {std(pStruct.x_kdens, pStruct.f_kdens)};
        
        % Format and save graph
        if mod(iParam, nLins*nCols) == 0 || iParam == nEstimParams
           
            pub_GraphSetInterpreter(previousInterpreter);

            % Custom font
            set(f.Children.Children, 'FontSize', 12);
            set(f.Children.Children, 'FontWeight', 'bold');

            % Resize title
            for iGraph = 1:length(tl.Children)
                tl.Children(iGraph).Title.FontSize = 18;
            end

            %set(gca, 'FontSize', 12);
            set(f, 'Position',  [100, 100, 800, 200 * nEffectiveLins]); % resize figure
            exportgraphics(f, ...
                strjoin({pathImages 'Estimation' ['Graph - Prior and Posterior - ' riseFullSetUpNames{iMdl} ' Part ' num2str(nFigures) '.png']}, filesep));
            pub_GraphSetInterpreter(previousInterpreter);
            
        end
    end
    
    % Format cells
    tEstimation.(['Model_' num2str(iMdl)]){:,'CI Min'} = cellfun(@(x) num2str(x,'%.3f'), tEstimation.(['Model_' num2str(iMdl)]){:,'CI Min'}, 'UniformOutput', false);
    tEstimation.(['Model_' num2str(iMdl)]){:,'CI Max'} = cellfun(@(x) num2str(x,'%.3f'), tEstimation.(['Model_' num2str(iMdl)]){:,'CI Max'}, 'UniformOutput', false);
    tEstimation.(['Model_' num2str(iMdl)]){:,'\% CI'} = cellfun(@(x) num2str(round(x,3)*100,'%.1f\\%%'), tEstimation.(['Model_' num2str(iMdl)]){:,'\% CI'}, 'UniformOutput', false);
    tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Mode'} = cellfun(@(x) num2str(x,'%.3f'), tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Mode'}, 'UniformOutput', false);
    tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Mean'} = cellfun(@(x) num2str(x,'%.3f'), tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Mean'}, 'UniformOutput', false);
    tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Std.'} = cellfun(@(x) num2str(x,'%.3f'), tEstimation.(['Model_' num2str(iMdl)]){:,'Post. Std.'}, 'UniformOutput', false);
    
    % Sort rows
    tEstimation.(['Model_' num2str(iMdl)]) = tEstimation.(['Model_' num2str(iMdl)])(paramsOrder, :);
    
    % Remove row names
    tEstimation.(['Model_' num2str(iMdl)]).Properties.RowNames = {};
    
    % Column names
    colNames = tEstimation.(['Model_' num2str(iMdl)]).Properties.VariableNames;
    tabWidth = '1.0\\textwidth';
    pub_Table2Latex(tEstimation.(['Model_' num2str(iMdl)]), ...
        [pathTables filesep 'Estimation' filesep 'modelEstimation_Parameters_' num2str(iMdl) '.tex'], ...
        'colAlignment', 'llccccccc', 'colNames', colNames, 'tabWidth', tabWidth);
    
    % Display estimation table
    disp(tEstimation.(['Model_' num2str(iMdl)]));
end

% Save results from estimation
save('estimPosteriorResults', 'estimPostSimData', 'tEstimation');

%% Check curvature at the mode

for iMdl = 1:length(listMdls)
    curvData = mode_curvature(estMdl(iMdl));

    f = figure;
    previousInterpreter = pub_GraphSetInterpreter('latex');

    nEstimParams = length(estimParamsNames{iMdl});
    nCols = 4;
    nLins = ceil(nEstimParams / nCols);

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';

    for iParam=1:nEstimParams
        tl_h = nexttile();

        pStruct = curvData.(estimParamsNames{iMdl}{iParam});
        x = [pStruct.x' pStruct.x'];
        y = [pStruct.log_lik' pStruct.log_post'];
        p = plot(x, y);
        p(1).LineWidth = 2;
        p(2).LineWidth = 2;
        p(2).LineStyle = ':';
        title(pStruct.tex_name);
    end
    pub_GraphSetInterpreter(previousInterpreter);

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    set(gca, 'FontSize', 12);
    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Estimation' ['Graph - Curvature - ' riseFullSetUpNames{listMdls(iMdl)} '.png']}, filesep));
    pub_GraphSetInterpreter(previousInterpreter);

end

%% Estimation diagnostics

% preparing the variables for posterior sampling [objective,lb,ub,mu,SIG]=pull_objective(m);
% running the Metropolis-Hastings algorithm results=mh_sampler(objective,lb,ub,options,mu);
% preparing the environment for convergence diagnostics obj=mcmc(results,pnames), where pnames is the list of names of the estimated parameters
% then if you want to do the trace plots traceplot(obj,name,chain_id) where 'name' is the name of one particular parameter and 'chain_id' is the ID of the chain you want to do the traceplot from
% if you want to plot the density densplot(obj,name,chain_id,N); where N is the number of points to consider on the x-axis. The default is 250
% if you want to plot the recursive means meanplot(obj,name,chain_id), where chain_id can be empty, in which case the recursive means of all the chains are plotted.
% If you want to plot the autocorrelation of a parameter autocorrplot(obj,name,chain_id,order) where 'order' is the maximum order of the autocorrelation to consider. It can be left empty in which case the default value is used.
% if you want to plot the Potential Scale Reduction Factor (PSRF), psrf_plot(obj,name)
% if you want to do a scatter plot for two parameters, scatterplot(obj,name1,name2,chain_id)

fontSize = 12;
fontSizeOutside = 14;

analysisName = 'psrf_plot'; %'traceplot'; 'densplot'; 'meanplot'; 'autocorrplot'; 'psrf_plot'; 'scatterplot'

for iMdl = 1:length(listMdls)
    estimResults = postSimResults.(['Model_' num2str(iMdl)]);
    obj = mcmc(estimResults, estimParamsNames{iMdl});
    
    f = figure;
    oldInterpreter = pub_GraphSetInterpreter('latex');
    nCols = 4;
    nLins = ceil(nParams / nCols);
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    nPoints = 250;
    maxOrder = 10000;
    iterationStart = 1; %iteration at which to start the plot of the PSRF. If a function handle is used then it should take as input the total number of observations and return the point at which to start. e.g. @(x)round(0.5*x)
    refParam = 1;
    pName1 = estimParamsNames{iMdl}{1};
    for iParam = 1:length(estimParamsNames{iMdl})
        pName = estimParamsNames{iMdl}{iParam};
        pName2 = estimParamsNames{iMdl}{refParam};
        nexttile();
        for iChain = 1:obj.nchains
            hold on;
            %traceplot(obj, pName, iChain);
            %densplot(obj, pName, iChain, nPoints);
            %meanplot(obj, pName, iChain);
            %autocorrplot(obj, pName, iChain, maxOrder);
            if strcmp(analysisName, 'psrf_plot')
                h = psrf_plot(obj, pName, iterationStart);
                h.LineWidth = 1.5;
                yCurrentLim = h.Parent.YLim;
                nDigits = 2;
                h.Parent.YTick = linspace(1, round(10^nDigits*yCurrentLim(2))/10^nDigits, 5);
                h.Parent.YTickLabel = linspace(1, round(10^nDigits*yCurrentLim(2))/10^nDigits, 5);
                h.Parent.YLim = [1 max(h.Parent.YTick)];
            end
            %scatterplot(obj, pName1, pName2, iChain);
            hold off;
        end
        
        title(estMdl.estimation.priors(iParam).tex_name);
        set(gca, 'FontSize', fontSize);
       
    end
end

set(f, 'Position',  [100, 100, 800, 800]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Estimation' ['Graph - ' analysisName ' - ' riseFullSetUpNames{listMdls(iMdl)} '.png']}, filesep));
pub_GraphSetInterpreter(previousInterpreter);

%% My estimation diagnostics: Acceptance ratio

acceptanceRatios = NaN(nChains_mcmc,1);
for iChain = 1:nChains_mcmc
    acceptanceRatios(iChain) = postSimResults.Model_1{iChain}.stats.accept_ratio;
end

% Display
acceptanceRatios
mean(acceptanceRatios)

%% My estimation diagnostics: Autocorrelation of the MC (finds minimum size of the chain)

% What is smallest lag to give an ρlag ≈ 0?
% One of several methods for estimating how many iterations of
% Markov chain are needed for effectively independent samples

fontSize = 12;
fontSizeOutside = 14;

for iMdl = 1:length(listMdls)
    estimResults = postSimResults.(['Model_' num2str(iMdl)]);
    obj = mcmc(estimResults, estimParamsNames{iMdl});
    
    f = figure;
    oldInterpreter = pub_GraphSetInterpreter('latex');
    nParams = length(estimParamsNames{iMdl});
    nCols = 4;
    nLins = ceil(nParams / nCols);
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';
    nPoints = 250;
    maxOrder = 39000;
    pName1 = estimParamsNames{iMdl}{1};
    for iParam = 1:length(estimParamsNames{iMdl})
        pName = estimParamsNames{iMdl}{iParam};
        pName2 = estimParamsNames{iMdl}{iParam};
        tile = nexttile();
        for iChain = 1:obj.nchains
            hold on;
            [acf,lags,bounds] = autocorr(obj.draws(iParam, 1:end, iChain), 'NumLags', size(obj.draws,2)-1);
            p = plot(acf);
            hold off;
        end
        title(estMdl.estimation.priors(iParam).tex_name);
        pub_GraphDrawZeroAxis(p);
        
        set(gca, 'FontSize', fontSize);
    end
end

%title(tl,'Autocorrelation of Markov Chains')
xlabel(tl,'Periods', 'FontSize', fontSizeOutside);
ylabel(tl,'Autocorrelation', 'FontSize', fontSizeOutside);

set(f, 'Position',  [100, 100, 800, 800]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Estimation' ['Graph - Autocorrelation MC - ' riseFullSetUpNames{listMdls(iMdl)} '.png']}, filesep));
pub_GraphSetInterpreter(previousInterpreter);

%% Save estimation results

%load([pathSaved filesep 'preliminary - estimationResults - ' riseFullSetUpNames{listMdls} '.mat']);
save([pathSaved filesep 'Final - estimationResults - ' riseFullSetUpNames{listMdls} '.mat']);

%load([pathSaved filesep 'estimationResults - ' riseFullSetUpNames{iMdl} '.mat']);
%save([pathSaved filesep 'estimationResults - ' riseFullSetUpNames{listMdls} '.mat']);
%load([pathSaved filesep 'estimationResults.mat']);
%save([pathSaved filesep 'estimationResults.mat']);

%% Plot filtered and smoothed series

% 1- filtered_variables refer to one-step ahead forecasts: a_{t|t-1}
% 2- updated_variables refer to the updates a_{t|t}
% 3- smoothed_variables refer to the final estimates conditional on all
% available information a_t{t|n}
% 
% As for the expected counterparts, they are averages across all regimes.
% This means that in a regime switching model, you will have for every
% variable, as many series as the number of regimes. The weights used to
% compute those averages are the probabilities.

close all;

plotObsAndSmoothed = true;

for iMdl = 1:length(listMdls)

    % Select variables
    obsVars             = {'obs_y', 'obs_c', 'obs_g', 'obs_pii', 'obs_swap_PreDI_3m'};
    filterSmoothVars    = {'y', 'c', 'g', 'Pii', 'nrPolicy'};
    vlocs               = locate_variables(filterSmoothVars, estMdl(iMdl).endogenous.name);
    vtexNames           = estMdl(iMdl).endogenous.tex_name(vlocs);
    plotVars            = filterSmoothVars;
    plotVarNames        = vtexNames;

    % Filter and smooth series
    myFiltSmooth = filter(estMdl(iMdl),'data',myData);

    nCols = 2;
    nLins = ceil(length(plotVars) / nCols);

    f = figure;
    previousInterpreter = pub_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';


    for iVar = 1:numel(plotVars)

        tl_h = nexttile();
        %tl_h = [tl_h, nexttile()];

        obsData       = myData(obsVars{iVar}).values;
        filtData      = myFiltSmooth.filtered_variables.(filterSmoothVars{iVar}).values;
        smoothData    = myFiltSmooth.smoothed_variables.(filterSmoothVars{iVar}).values;


        seriesRef = myFiltSmooth.filtered_variables.(filterSmoothVars{iVar});
        dates = (datetime(seriesRef.start, 'InputFormat', 'yyyyQQQ'):calquarters(1):datetime(seriesRef.finish, 'InputFormat', 'yyyyQQQ'))';
        [obsData, filtData, smoothData] = transformVariables(filterSmoothVars{iVar}, mapEstimParams{iMdl}, obsData, filtData, smoothData);

        nReg = size(filtData, 2);
        plotData = [];
        legStr = {'observed'};
        for iReg = 1:nReg
            plotData = [plotData, real([ ...
                            [obsData; NaN], ...
                            filtData(:,iReg), ...
                            [smoothData(:,iReg); NaN] ])
                         ];
            legStr = [legStr {'filtered'} {'smoothed'}];
        end
        graph = plot(dates, plotData);
        graph(1).LineWidth = 1.5;
        graph(2).LineWidth = 1.5;
        graph(2).LineStyle = ':';
        graph(2).Color = 'green';
        graph(3).LineWidth = 1.5;
        graph(3).LineStyle = '-.';
        graph(3).Color = 'red';
        title(plotVarNames{iVar});
        legend(legStr, 'Location', 'best', 'FontSize', 14);
        ylabel('\%');

        pub_GraphDrawZeroAxis(graph);
        
        if plotObsAndSmoothed
            for iGraph = 2:2:(nMdl+1)
                delete(graph(iGraph));
            end
        end

        % Row label
        %text(graph(1), mean(xlim(graph(1))), mean(ylim(graph(1))), ...
        %    shocksToPlotNames{iShock}, ...
        %    'Interpreter', 'latex', ...
        %    'FontSize', 20);

        % First row
        %if iVar == 1
        %    graph(1).Legend.Location = 'North';
        %end

        % Hide after the first graph
        if iVar > 1
        %    title(graph(1:end), '');
            legend('hide');
        end

        % Between lines 1 and 3
        %if iShock >= 1 && iShock <= 3
        %    set(graph(3:5), 'Color', [0.9 0.9 0.9]);
        %end

        % Hide after column 2
        %ylabel(graph(3:end), '');

        % Hide after column 2
        %yticks(graph(3:end), []);

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    ax = nexttile(1);
    ax.Legend.FontSize = 14;

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Estimation' 'Graph - Estimation - Filtered and Smoothed Shocks.png'}, filesep));
    pub_GraphSetInterpreter(previousInterpreter);
    
end

%% Plot smoothed shocks
close all;

plotShocks = {'epsA', 'epsBeta', 'epsG', 'epsM'};

for iMdl = 1:length(listMdls)
    
    % Filter and smooth series
    myFiltSmooth = filter(estMdl(iMdl),'data',myData);
    mdlPlotShocks = shock_names{iMdl};
    [~, idxShocks] = ismember(mdlPlotShocks, plotShocks);
    idxShocks = idxShocks ~= 0;
    mdlPlotShocks = mdlPlotShocks(idxShocks);
    mdlPlotTexNames = shock_texnames{iMdl}(idxShocks);
    
    nCols = 2;
    nLins = ceil(length(mdlPlotShocks) / nCols);

    f = figure;
    previousInterpreter = pub_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';

    % Smoothed shocks
    smoothedShocks      = myFiltSmooth.smoothed_shocks;

    for iShock = 1:numel(mdlPlotShocks)

        tl_h = nexttile();
        %tl_h = [tl_h, nexttile()];

        smoothedShocksData      = myFiltSmooth.smoothed_shocks.(mdlPlotShocks{iShock}).values;

        seriesRef = myFiltSmooth.smoothed_shocks.(mdlPlotShocks{iShock});
        dates = (datetime(seriesRef.start, 'InputFormat', 'yyyyQQQ'):calquarters(1):datetime(seriesRef.finish, 'InputFormat', 'yyyyQQQ'))';

        plotData = seriesRef.values;
        graph = bar(dates, plotData);
        title(mdlPlotTexNames{iShock});
        ylabel('std.');

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    set(f, 'Position',  [100, 100, 800, 400]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Estimation' 'Graph - Individual Smoothed Shocks.png'}, filesep));
    pub_GraphSetInterpreter(previousInterpreter);

end
    
%% Historical decomposition

close all;

plotConfig_HistDec;

plotVars = {'y', 'c', 'g', 'Pii', 'nrPolicy'};
plotShocks = {'epsA', 'epsBeta', 'epsG', 'epsM', 'init'};

nCols = 2;
nLins = ceil((length(plotVars) + 1) / nCols);

% Plot the historical decomposition
histdec = struct();
for listMdls = 1:length(estMdl)
    
    modelIdxName = ['Model_' char(string(listMdls))];
    histdec.(modelIdxName) = historical_decomposition(estMdl(listMdls));
    
    f = figure('name','Historical decomposition of shocks and initial conditions');
    previousInterpreter = pub_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';
    legendInFirstSubPlot = true;

    if legendInFirstSubPlot
       subGraph_legend = nexttile(); 
    end

    nBars = size(histdec.(modelIdxName).( plotVars{1}), 2);
    for iVar = 1:numel(plotVars)
        
        vname           = plotVars{iVar};
        [~, idxShockSelected] = ismember([histdec.(modelIdxName).(vname).varnames]', plotShocks);
        idxShockSelected = idxShockSelected ~=0;
        [~, idxShock] = ismember([histdec.(modelIdxName).(vname).varnames]', shockFullName(:,1));
        legendStr       = shockFullName(idxShock(idxShockSelected),2);
        
        
        ssValue = endo_ss{listMdls}(find(ismember(endo_names{iMdl}, vname)));
        
        nexttile();
        %smoothedVar = myFiltSmooth.smoothed_variables.(vname);
        p = bar(histdec.(modelIdxName).(vname)(plotShocks), 'stacked','BarWidth', 1);
        for iBar = 1:length(p)
            p(iBar).FaceAlpha = 1.0;
        end
        set(p, 'LineStyle', 'none');
%         p(1).FaceColor = 'blue';
%         p(2).FaceColor = 'red';
%         p(3).FaceColor = 'green';
%         p(4).FaceColor = 'blue';
%         p(5).FaceColor = 'red';
        
        p(1).YData = p(1).YData - ssValue; % Subtract steady state
        yAxisLims_left = ylim();
        
%         % Detrend
%         yTrendEstimated = estMdl(iMdl).estimation.posterior_maximization.mode( ...
%              find(strcmp({estMdl(iMdl).estimation.priors.name},'ytrend')));
%         if ismember(vname, {'y', 'c'})
%         p(1).YData = p(1).YData - yTrendEstimated.*(0:(length(p(1).YData)-1)); % Subtract trend
%         end
        
        % Plot smoothed series
        plotShocks_woutInit = plotShocks(~strcmp(plotShocks, 'init'));
        initSeries = ssValue - histdec.(modelIdxName).(vname)('init');
        ySmoothed = sum(histdec.(modelIdxName).(vname)(plotShocks),2) - ssValue;
        %[~, ~, ySmoothed] = transformVariables(vname, mapEstimParams{iMdl}, [], [], ySmoothed);
        hold on;
        %yyaxis right;
        pSmoothed = plot(p(1).XData, ySmoothed);
        pSmoothed.Color = 'black';
        pSmoothed.LineWidth = 2;
        pSmoothed.LineStyle = '--';
        hold off;
        
        %%align axes for left and right
%         yyaxis left; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
%         if ratio ~= 0
%             yyaxis right; yliml = get(gca,'Ylim');
%             if yliml(2)*ratio<yliml(1)
%                 set(gca,'Ylim', ssValue + [yliml(2)*ratio yliml(2)])
%             else
%                 set(gca,'Ylim', ssValue + [yliml(1) yliml(1)/ratio])
%             end
%         end
        
        %hold on;
        %p2 = plot(smoothedVar);
        %hold off;
        
        title(endo_texnames{iMdl}(find(ismember(endo_names{iMdl}, vname))));
        
        if iVar == 1 && legendInFirstSubPlot
            %legendAxes = [p2; flip(p)];
            legendAxes = [ flip(p)];
            pub_GraphPutLegendInFirstSubplot(...
                subGraph_legend, legendAxes, legendStr, ...
                'Interpreter', 'latex', ...
                'Background', 'boxoff');
        end

        pub_GraphDrawZeroAxis(pSmoothed);
        
    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    ax = nexttile(1);
    ax.Legend.FontSize = 14;
    ax.Legend.NumColumns = ceil(length(legendStr) / 4);
    
    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'IRFs' modelName{listMdls} 'DefR' defProcess{iMdl} ['Graph - Historical Decomposition - ' riseConfigVector(listMdls).conf_policyRule '.png']}, filesep));
    pub_GraphSetInterpreter(previousInterpreter);
    
end

%% Variance decomposition

close all;

refMdl = 1;
plotConfig_HistDec;

plotVars = {'y', 'c', 'g', 'Pii', 'nrPolicy'};
plotShocks = {'epsA', 'epsBeta', 'epsG', 'epsM'};

plotPeriods = 1:40;

nCols = 2;
nLins = ceil((length(plotVars) + 1) / nCols);


% Plot the historical decomposition
vardec = struct();
for listMdls = 1:length(estMdl)
    
    modelIdxName = ['Model_' char(string(listMdls))];
    
    %vardec_shocks = ; % List of shocks
    vardec_theoretical  = true;
    vardec_ergodic      = true;
    [vardec.(modelIdxName), obj] = variance_decomposition(estMdl(listMdls), ...
        'vardec_theoretical',   vardec_theoretical, ...
        'vardec_ergodic',       vardec_ergodic);
    
    %theoAutoCov     = theoretical_autocovariances(estMdl(iMdl));
    %theoAutoCorrel  = theoretical_autocorrelations(estMdl(iMdl));
    
    f = figure('name','Variance decomposition of shocks');
    previousInterpreter = pub_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';
    legendInFirstSubPlot = true;

    if legendInFirstSubPlot
       subGraph_legend = nexttile(); 
    end

    for iVar = 1:numel(plotVars)

        vname           = plotVars{iVar};
        [~, idxShockSelected] = ismember(vardec.(modelIdxName).conditional.(vname).varnames', plotShocks);
        idxShockSelected = idxShockSelected ~=0;
        [~, idxShock] = ismember(vardec.(modelIdxName).conditional.(vname).varnames', shockFullName(:,1));
        legendStr       = shockFullName(idxShock(idxShockSelected),2);
        
        nexttile();
        barWidth = 1;
        plotData = 100 * vardec.(modelIdxName).conditional.(vname)(plotPeriods,:).values;
        plotData = plotData(:,idxShockSelected);
        plotData = [plotData; NaN(2,size(plotData,2))];
        plotDataInf = vardec.(modelIdxName).infinity.(vname)(1,:).values;
        plotDataInf = plotDataInf(idxShockSelected);
        plotData = [plotData; 100 * plotDataInf];
        p = bar(plotData, barWidth, 'stacked');
        
        title(endo_texnames{iMdl}(find(ismember(endo_names{iMdl}, vname))));
        xlabel('periods');
        tickPeriods = [1, plotPeriods(mod(plotPeriods,5)==0), 43];
        xticks(tickPeriods);
        xticklabels([string(tickPeriods(1:end-1)) '$\infty$']);
        xlim([1 45]);
        ylabel('\%')
        
        if iVar == 1 && legendInFirstSubPlot
            legendAxes = p;
            pub_GraphPutLegendInFirstSubplot(...
                subGraph_legend, legendAxes, legendStr, ...
                'Interpreter', 'latex', ...
                'Background', 'boxoff');
        end

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    ax = nexttile(1);
    ax.Legend.FontSize = 14;
    ax.Legend.NumColumns = ceil(length(legendStr) / 4);

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Estimation' ['Graph - Variance Decomposition - ' riseConfigVector(listMdls).conf_policyRule '.png']}, filesep));
    pub_GraphSetInterpreter(previousInterpreter);
    
end

%% Estimation sampling

estimationSampling;

%% Forecasting over all sampled parameters

myfkst=struct();
% date_start=[];
date_start='2017Q1';
nsteps=12;
shock_uncertainty=false;
Rfunc=[];

conditions=struct();
conditions.rtwi={'2017Q1','2019Q1'}; % range over which we want to condition

myfkst.ve=forecast(models.ve,db,date_start,params.ve,nsteps,...
    shock_uncertainty,Rfunc,conditions);

myfkst.ve_lr=forecast(models.ve_lr,db,date_start,params.ve_lr,nsteps,...
    shock_uncertainty,Rfunc,conditions);

% Set environment
ci=[30,50,68,90];

% Forecast plots
modelnames=fieldnames(models);

for jj=1:numel(modelnames)
    
    modname=modelnames{jj};
    
    figure('name',['model (',modname,') Forecasts of Norwegian Data']);
    
    for ii=1:numel(endog)
        
        subplot(3,2,ii)
        
        d=myfkst.(modname).(endog{ii});
        
        out=fanchart(d,ci);
        
        plot_fanchart(out,[244, 122, 66]/255)%'c'
        
        %     hold on
        %
        %     plot('1995:2004',db.(endog{ii}))
        
        title(tex.(endog{ii}))
        
        axis tight
        
    end
    
    xrotate(45)
    
end

%% Auxiliary functions

function [obsData, filtData, smoothData] = transformVariables(varName, mapEstimParams, obsData, filtData, smoothData)

if ismember(varName, {'c', 'y', 'g'})
    obsData     = obsData * 100;
    filtData    = (log(filtData./lagmatrix(filtData,1)) + pub_IIF(mapEstimParams.isKey('ytrend'), @() mapEstimParams('ytrend'), 0)) * 100;
    smoothData  = (log(smoothData./lagmatrix(smoothData,1)) + pub_IIF(mapEstimParams.isKey('ytrend'), @() mapEstimParams('ytrend'), 0)) * 100 ;
elseif ismember(varName, {'n'})
    obsData     = obsData * 100;
    filtData    = log(filtData./lagmatrix(filtData,1)) * 100;
    smoothData  = log(smoothData./lagmatrix(smoothData,1)) * 100;
elseif ismember(varName, {'Pii'})
    obsData     = obsData * 100;
    filtData    = log(filtData) * 100;
    smoothData  = log(smoothData) * 100;
elseif ismember(varName, {'nrPolicy'})
    obsData     = obsData * 100;
    filtData    = ((1+filtData).^4 - 1) .* 100;
    smoothData  = ((1+smoothData).^4 - 1) .* 100;
elseif ismember(varName, {'polDef'})
    obsData     = obsData * 100;
    filtData    = ( (filtData - lagmatrix(filtData,1)) + pub_IIF(mapEstimParams.isKey('polDefWedge'), @() mapEstimParams('polDefWedge'), 0) ) * 100;
    smoothData  = ( (smoothData - lagmatrix(smoothData,1)) + pub_IIF(mapEstimParams.isKey('polDefWedge'), @() mapEstimParams('polDefWedge'), 0) ) * 100;
end

end