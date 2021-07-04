%% Parameterize model
function [paramsStruct, priorsStruct] = createParameters(mdl, riseConfig, parameterization, setInitialPriorValuesToParams, varargin)

if length(varargin) == 1
   simStruct = varargin{1};
end

% Specify alternative calibrations
paramsStruct = struct();
priorsStruct = struct();
    
paramNames  = fieldnames(simStruct.params);
ssNames     = fieldnames(simStruct.ss);

iParam = 0;
for strParam = mdl.parameters.name
    iParam = iParam + 1;
    if ismember(strParam, paramNames)
        paramsStruct.(strParam{:}) = simStruct.params.(strParam{:});
    elseif ismember(strParam, ssNames)
        paramsStruct.(strParam{:}) = simStruct.ss.(strParam{:});
    elseif mdl.parameters.is_auxiliary(iParam) == false
        if mdl.parameters.is_switching(iParam) == false
            disp(['Parameter ' strParam{:} ' is empty!']);
            paramsStruct.(strParam{:}) = [];
        end
    end
end

if strcmp(riseConfig.conf_prodProcessType, 'Single_LogNormal')
    paramsStruct.defSS         = 0.0;
end

% Occasionally binding constraints
if strcmp(riseConfig.conf_occBinConstraint, 'FiscalLimit')
    
    % Conditional parameters
    paramsStruct.bindFis_r1fisLim_1   = 0;
    paramsStruct.bindFis_r1fisLim_2   = 1;
    
    paramsStruct.bindTax_r2taxLim_1   = 0;
    paramsStruct.bindTax_r2taxLim_2   = 1;
    
    paramsStruct.r1fisLim_tp_1_2      = 0;
    paramsStruct.r1fisLim_tp_2_1      = 0;
    
    paramsStruct.r2taxLim_tp_1_2      = 0;
    paramsStruct.r2taxLim_tp_2_1      = 0;
    
    %paramsStruct                    = rmfield(paramsStruct, 'taxLim_tp_1_2');
    %paramsStruct                    = rmfield(paramsStruct, 'taxLim_tp_2_1');
    
    % Tax lim
    %paramsStruct.sigmaTau_taxLim_1   = paramsStruct.sigmaTau;
    %paramsStruct.sigmaTau_taxLim_2   = 0;
    %paramsStruct                    = rmfield(paramsStruct, 'sigmaTau');
    
    % Conditional parameters
    
    
    %paramsStruct.phi_fisLim_1       = paramsStruct.phi;
    %paramsStruct.phi_fisLim_2       = paramsStruct.phi;
    %paramsStruct                    = rmfield(paramsStruct, 'phi');
    
    %paramsStruct.gammaTau_fisLim_1  = paramsStruct.gammaTau;
    %paramsStruct.gammaTau_fisLim_2  = paramsStruct.gammaTau;
    %paramsStruct                    = rmfield(paramsStruct, 'gammaTau');
    
    %paramsStruct                    = rmfield(paramsStruct, 'fisLim_tp_1_2');
    %paramsStruct                    = rmfield(paramsStruct, 'fisLim_tp_2_1');
    
else
    paramsStruct.bindTax   = 0;
    paramsStruct.bindFis   = 0;
end

if riseConfig.conf_effectiveLowerBound
    paramsStruct.bindELB_r3elbLim_1   = 0;
    paramsStruct.bindELB_r3elbLim_2   = 1;
else
    paramsStruct.bindELB   = 0;
end

%% CALCULATE STEADY STATE polDef
aux_bindFis = 0; % Value in the Regime for which to calculate the steady state
syms sym_Psi sym_tax sym_polDef sym_b sym_nr sym_r

if strcmp(riseConfig.conf_policyRule, "polDefAdjusted_priceLevel") || strcmp(riseConfig.conf_policyRule, "polDefAdjusted_inflation")
    eqns = [
        sym_Psi             == (1-paramsStruct.defSS)*paramsStruct.aTildeBar + paramsStruct.defSS*paramsStruct.recTildeSS ;
        sym_tax             == paramsStruct.tauBar * paramsStruct.kBar * ((1-paramsStruct.tauBar) * paramsStruct.kBar/paramsStruct.eta * sym_Psi)^(1/paramsStruct.chi) * sym_Psi ;
        sym_polDef          == fFiscalLimit_defProb(sym_b, paramsStruct.aTildeBar, paramsStruct.recTildeSS, paramsStruct.defSS, paramsStruct.gYSS*paramsStruct.yBar, paramsStruct.bbetaSS) ;
        sym_b               == (1 + sym_nr) * (paramsStruct.gYSS*paramsStruct.yBar + paramsStruct.zYSS*paramsStruct.yBar - sym_tax) / (1  -(1 + sym_nr)*(1-aux_bindFis*paramsStruct.deltaBar)/(1+paramsStruct.piiBar)) ;
        sym_nr              == -1 + ( sym_r + sym_polDef ) / (1 - aux_bindFis*paramsStruct.deltaBar) ;
        sym_r               == -1 + 1/(paramsStruct.beta * (1 - aux_bindFis*paramsStruct.deltaBar)) ;
    ];
else
    eqns = [
        sym_Psi             == (1-paramsStruct.defSS)*paramsStruct.aTildeBar + paramsStruct.defSS*paramsStruct.recTildeSS ;
        sym_tax             == paramsStruct.tauBar * paramsStruct.kBar * ((1-paramsStruct.tauBar) * paramsStruct.kBar/paramsStruct.eta * sym_Psi)^(1/paramsStruct.chi) * sym_Psi ;
        sym_polDef          == fFiscalLimit_defProb(sym_b, paramsStruct.aTildeBar, paramsStruct.recTildeSS, paramsStruct.defSS, paramsStruct.gYSS*paramsStruct.yBar, paramsStruct.bbetaSS) ;
        sym_b               == (1 + sym_nr) * (paramsStruct.gYSS*paramsStruct.yBar + paramsStruct.zYSS*paramsStruct.yBar - sym_tax) / (1  -(1 + sym_nr)*(1-aux_bindFis*paramsStruct.deltaBar)/(1+paramsStruct.piiBar)) ;
        sym_nr              == - 1 + (1 + sym_r) * (1+paramsStruct.piiBar) ;
        sym_r               == -1 + 1/(paramsStruct.beta * (1 - aux_bindFis*paramsStruct.deltaBar)) ;
    ];
end
soln = vpasolve(eqns, symvar(eqns));
paramsStruct.polDefSS = double(soln.sym_polDef);

%% Add the initial conditions of the priors to the parameters
%------------------------------------------------------------
if setInitialPriorValuesToParams
    fields = fieldnames(priorsStruct);
    for ip = 1:numel(fields)
        name = fields{ip};
        paramsStruct.(name) = priorsStruct.(name){1};
    end
end

end