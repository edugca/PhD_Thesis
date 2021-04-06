%%%%%%%%% Amaral(2021): New-Keynesian Models with Defaultables Policy Assets in Zero-Net-Supply

clear all; clc;

pathImages = 'C:\Users\Eduardo\OneDrive\MATLAB\My Thesis\Paper 7\Images';

%% Run dynare

flexiblePrices = false; % You need to change in the Dynare mod file as well

cd 'Closed Economy'
dynare risky_ClosedEconomy
cd ../

%% Collect impulse response functions

if flexiblePrices
    priceRegime = 'flexible';
else
    priceRegime = 'sticky';
end

varNames = {
                '$\xi^m_t$'
                '$\xi^d_t$'
                '$\xi^\omega_t$'
                '$\hat{i}_t$'
                '$\pi_t$'
                '$\hat{r}^{RF}_t$'
                '$\hat{y}_t$'
                '$E_t \hat{\mathcal{D}}_{t+1}$'
                '$E_t \hat{\omega}_{t+1}$'
            };

irfs_epsM = [   oo_.irfs.xxiM_eepsM', ...
                oo_.irfs.xxiD_eepsM', ...
                oo_.irfs.xxiR_eepsM', ...
                oo_.irfs.nr_eepsM', ...
                oo_.irfs.ppi_eepsM', ...
                oo_.irfs.rRF_eepsM', ...
                oo_.irfs.y_eepsM', ...
                oo_.irfs.defExp_eepsM', ...
                oo_.irfs.recExp_eepsM'
             ];

irfs_epsD = [   oo_.irfs.xxiM_eepsD', ...
                oo_.irfs.xxiD_eepsD', ...
                oo_.irfs.xxiR_eepsD', ...
                oo_.irfs.nr_eepsD', ...
                oo_.irfs.ppi_eepsD', ...
                oo_.irfs.rRF_eepsD', ...
                oo_.irfs.y_eepsD', ...
                oo_.irfs.defExp_eepsD', ...
                oo_.irfs.recExp_eepsD'
             ];

irfs_epsR = [   oo_.irfs.xxiM_eepsR', ...
                oo_.irfs.xxiD_eepsR', ...
                oo_.irfs.xxiR_eepsR', ...
                oo_.irfs.nr_eepsR', ...
                oo_.irfs.ppi_eepsR', ...
                oo_.irfs.rRF_eepsR', ...
                oo_.irfs.y_eepsR', ...
                oo_.irfs.defExp_eepsR', ...
                oo_.irfs.recExp_eepsR'
             ];
         
% Numerical accuracy
minTol = 1e-10;
irfs_epsM(abs(irfs_epsM) < minTol) = 0;
irfs_epsD(abs(irfs_epsD) < minTol) = 0;
irfs_epsR(abs(irfs_epsR) < minTol) = 0;

if flexiblePrices
    flexible_irfs_epsM = irfs_epsM;
    flexible_irfs_epsD = irfs_epsD;
    flexible_irfs_epsR = irfs_epsR;
else
    sticky_irfs_epsM = irfs_epsM;
    sticky_irfs_epsD = irfs_epsD;
    sticky_irfs_epsR = irfs_epsR;
end

%% Plot IRFs

%%% epsM
f = figure;
edu_GraphSetInterpreter('latex');
tl = tiledlayout(3,3);
tl.TileSpacing      = 'compact';
tl.Padding          = 'compact';

for iCol = 1:length(varNames)
    nexttile();
    p1 = plot(flexible_irfs_epsM(:,iCol), 'DisplayName', 'Flexible');
    hold on;
    p2 = plot(sticky_irfs_epsM(:,iCol), 'DisplayName', 'Sticky');
    p1.LineWidth = 2;
    p2.LineWidth = 1.5;
    p1.LineStyle = ':';
    p2.LineStyle = '-';
    %grid('on');
    title(varNames{iCol});
    set(gca, 'FontSize', 14);
    
    if mod(iCol, 3) == 1
       ylabel('\%') ; 
    end
end

lg  = legend([p1 p2], 'Orientation', 'horizontal');
lg.Layout.Tile = 'North'; % <-- place legend east of tiles

set(f, 'Position',  [100, 0, 800, 600]); % resize figure
exportgraphics(f, [pathImages filesep 'IRFs_epsM.png']);

%%% epsD
f = figure;
edu_GraphSetInterpreter('latex');
tl = tiledlayout(3,3);
tl.TileSpacing      = 'compact';
tl.Padding          = 'compact';

for iCol = 1:length(varNames)
    nexttile();
    p1 = plot(flexible_irfs_epsD(:,iCol), 'DisplayName', 'Flexible');
    hold on;
    p2 = plot(sticky_irfs_epsD(:,iCol), 'DisplayName', 'Sticky');
    p1.LineWidth = 2;
    p2.LineWidth = 1.5;
    p1.LineStyle = ':';
    p2.LineStyle = '-';
    %grid('on');
    title(varNames{iCol});
    set(gca, 'FontSize', 14);
    
    if mod(iCol, 3) == 1
       ylabel('\%') ; 
    end
end

lg  = legend([p1 p2], 'Orientation', 'horizontal');
lg.Layout.Tile = 'North'; % <-- place legend east of tiles

set(f, 'Position',  [100, 0, 800, 600]); % resize figure
exportgraphics(f, [pathImages filesep 'IRFs_epsD.png']);

%%% epsR
f = figure;
edu_GraphSetInterpreter('latex');
tl = tiledlayout(3,3);
tl.TileSpacing      = 'compact';
tl.Padding          = 'compact';

for iCol = 1:length(varNames)
    nexttile();
    p1 = plot(flexible_irfs_epsR(:,iCol), 'DisplayName', 'Flexible');
    hold on;
    p2 = plot(sticky_irfs_epsR(:,iCol), 'DisplayName', 'Sticky');
    p1.LineWidth = 2;
    p2.LineWidth = 1.5;
    p1.LineStyle = ':';
    p2.LineStyle = '-';
    %grid('on');
    title(varNames{iCol});
    set(gca, 'FontSize', 14);
    
    if mod(iCol, 3) == 1
       ylabel('\%') ; 
    end
end

lg  = legend([p1 p2], 'Orientation', 'horizontal');
lg.Layout.Tile = 'North'; % <-- place legend east of tiles

set(f, 'Position',  [100, 0, 800, 600]); % resize figure
exportgraphics(f, [pathImages filesep 'IRFs_epsR.png']);

%% Calculate coefficients

%%%%%%%%%%% DEEP PARAMS
bbeta = get_param_by_name('bbeta');
eeta = get_param_by_name('eeta');
ttheta = get_param_by_name('ttheta');
nnu = get_param_by_name('nnu');

rrhoM = get_param_by_name('rrhoM');
rrhoD = get_param_by_name('rrhoD');
rrhoR = get_param_by_name('rrhoR');

defSS = get_param_by_name('defSS');
recSS = get_param_by_name('recSS');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% VECTOR
defSS_vals = 0:0.01:0.99;
recSS_vals = 0:0.01:0.99;

defSS_Bench = find(defSS_vals == defSS);
recSS_Bench = find(recSS_vals == recSS);

cM_vec = NaN(length(defSS_vals),length(recSS_vals));
cD_vec = NaN(length(defSS_vals),length(recSS_vals));
cR_vec = NaN(length(defSS_vals),length(recSS_vals));

dM_vec = NaN(length(defSS_vals),length(recSS_vals));
dD_vec = NaN(length(defSS_vals),length(recSS_vals));
dR_vec = NaN(length(defSS_vals),length(recSS_vals));

for iParam = 1:length(defSS_vals)
    for jParam = 1:length(recSS_vals)
    
        defSS = defSS_vals(iParam);
        recSS = recSS_vals(jParam);

        oomega = (1+eeta)*(1-ttheta)*(1-ttheta*bbeta) / ttheta ;
        pphii = 1 - bbeta*defSS*recSS ;
        pphid = (-defSS + bbeta*defSS*recSS) / (1 - defSS) ;
        pphir = bbeta * defSS * recSS ;

        cM = -pphii*(1-bbeta*rrhoM) / (oomega*(pphii*nnu - rrhoM) + (1-bbeta*rrhoM)*(1-rrhoM) ) ;
        cD = -pphid*(1-bbeta*rrhoD) / (oomega*(pphii*nnu - rrhoD) + (1-bbeta*rrhoD)*(1-rrhoD) ) ;
        cR = -pphir*(1-bbeta*rrhoR) / (oomega*(pphii*nnu - rrhoR) + (1-bbeta*rrhoR)*(1-rrhoR) ) ;

        dM = -pphii*oomega / (oomega*(pphii*nnu - rrhoM) + (1-bbeta*rrhoM)*(1-rrhoM) ) ;
        dD = -pphid*oomega / (oomega*(pphii*nnu - rrhoD) + (1-bbeta*rrhoD)*(1-rrhoD) ) ;
        dR = -pphir*oomega / (oomega*(pphii*nnu - rrhoR) + (1-bbeta*rrhoR)*(1-rrhoR) ) ;

        cM_vec(iParam, jParam) = cM;
        cD_vec(iParam, jParam) = cD;
        cR_vec(iParam, jParam) = cR;

        dM_vec(iParam, jParam) = dM;
        dD_vec(iParam, jParam) = dD;
        dR_vec(iParam, jParam) = dR;
    end
end

%% GRAPHS defSS

f = figure;
tl = tiledlayout(2,3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

%%%%%%%%%%%%%% cM
nexttile();
p = plot(defSS_vals, cM_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), cM_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_m$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% cD
nexttile();
p = plot(defSS_vals, cD_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), cD_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_d$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% cR
nexttile();
p = plot(defSS_vals, cR_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), cR_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_\omega$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%% GRAPHS

%%%%%%%%%%%%%% dM
nexttile();
p = plot(defSS_vals, dM_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), dM_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_m$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% dD
nexttile();
p = plot(defSS_vals, dD_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), dD_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_d$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% dR
nexttile();
p = plot(defSS_vals, dR_vec(:, recSS_Bench) );
p.LineWidth = 1.5;
hold on;
p2 = plot(defSS_vals(1), dR_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_\omega$');
xlabel('$\overline{\mathcal{D}}$');
ylabel('value');
set(gca, 'FontSize', 14);

set(f, 'Position',  [100, 0, 800, 500]); % resize figure
exportgraphics(f, [pathImages filesep 'defaultableNK_coeffs_defSS.png']);

%% GRAPHS recSS

f = figure;
tl = tiledlayout(2,3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

%%%%%%%%%%%%%% cM
nexttile();
p = plot(recSS_vals, cM_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), cM_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_m$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% cD
nexttile();
p = plot(recSS_vals, cD_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), cD_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_d$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% cR
nexttile();
p = plot(recSS_vals, cR_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), cR_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$c_\omega$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% dM
nexttile();
p = plot(recSS_vals, dM_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), dM_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_m$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% dD
nexttile();
p = plot(recSS_vals, dD_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), dD_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_d$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

%%%%%%%%%%%%%% dR
nexttile();
p = plot(recSS_vals, dR_vec(defSS_Bench, :)' );
p.LineWidth = 1.5;
hold on;
p2 = plot(recSS_vals(1), dR_vec(1, 1),'r*');
p2.MarkerSize = 10;
grid('on');
edu_GraphSetInterpreter('latex');
title('$d_\omega$');
xlabel('$\overline{\omega}$');
ylabel('value');
set(gca, 'FontSize', 14);

set(f, 'Position',  [100, 0, 800, 500]); % resize figure
exportgraphics(f, [pathImages filesep 'defaultableNK_coeffs_recSS.png']);