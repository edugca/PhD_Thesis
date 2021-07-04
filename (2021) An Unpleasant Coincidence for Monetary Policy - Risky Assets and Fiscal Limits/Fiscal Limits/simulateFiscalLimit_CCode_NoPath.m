function [fiscalLimit_mu, fiscalLimit_std, fiscalLimit_dbl] = simulateFiscalLimit_CCode_NoPath(...
    nSimulations, nStates, nFiscalLimit, ...
    sim_epsA, sim_epsR, sim_epsD, sim_epsG, sim_epsBeta, ...
    rng_epsA, rng_epsR, rng_epsD, rng_epsG, rng_epsBeta, ...
    aTildeSS, recTildeSS, gSS, defSS, ySS, zSS, tLSSS, yMaxSS, tauMaxSS, ...
    beta, chi, eta, kkSS, sigma, alphaG, gammaGPSI, fracNR, ...
    rhoBeta, rhoA, rhoR, rhoD, rhoGG, rhoGY, ...
    sigmaBeta, sigmaA, sigmaR, sigmaD, sigmaG, ...
    psiShockObservedWithDelay, ...
    govCeiling) %#codegen

% Asserts to avoid dynamic memory allocation
%assert(nStates <= 21);
%assert(nSimulations <= 200000);

% Build arrays
fiscalLimit_mu  = NaN(nStates,nStates,nStates,nStates,nStates);
fiscalLimit_std = NaN(nStates,nStates,nStates,nStates,nStates);
fiscalLimit_dbl = NaN(nSimulations, 1);

% Define variables
base_sim_uCMax = NaN(nSimulations, 1);

% Prepare steady states
sim_bbeta      = beta .* complex(ones(nSimulations,2));
sim_aTilde      = aTildeSS .* complex(ones(nSimulations,2));
%sim_aTildeExp   = aTildeSS .* complex(ones(nSimulations,2));
sim_recTilde    = recTildeSS .* complex(ones(nSimulations,2));
%sim_recTildeExp = recTildeSS .* complex(ones(nSimulations,2));
sim_g           = gSS .* complex(ones(nSimulations,2));
%sim_gExp        = gSS .* complex(ones(nSimulations,2));
sim_def         = defSS .* complex(ones(nSimulations,2));
%sim_defExp      = defSS .* complex(ones(nSimulations,2));
sim_z           = zSS .* complex(ones(nSimulations,2));
sim_tLS         = tLSSS .* complex(ones(nSimulations,2));
%sim_psi         = complex(NaN(nSimulations,2));
%sim_psiExp      = complex(NaN(nSimulations,2));
%sim_nMax        = complex(NaN(nSimulations,2));
sim_yMax        = yMaxSS .* complex(ones(nSimulations,2));
sim_tauMax      = tauMaxSS .* complex(ones(nSimulations,2));
%sim_taxMax      = complex(NaN(nSimulations,2));
%sim_cMax        = complex(NaN(nSimulations,2));
%sim_uCMax       = complex(NaN(nSimulations,2));

nStates_epsA    = length(rng_epsA);
nStates_epsR    = length(rng_epsR);
nStates_epsD    = length(rng_epsD);
nStates_epsG    = length(rng_epsG);
nStates_epsBeta = length(rng_epsBeta);

% Iterate the state space
for iA = 1:nStates_epsA
    
    sim_epsA(:,2)   = rng_epsA(iA);
    
    for iR = 1:nStates_epsR
        
        sim_epsR(:,2)   = rng_epsR(iR);
        for iD = 1:nStates_epsD
            
            sim_epsD(:,2)   = rng_epsD(iD);
            
            for iG = 1:nStates_epsG
                
                disp(['Iterating (aTilde, recTilde, def, g, bbeta): ', sprintf('%d',iA), ' ', sprintf('%d',iR), ' ', sprintf('%d',iD), ' ', sprintf('%d',iG), ' ', sprintf('%d',1)]);
                tic
                
                sim_epsG(:,2)   = rng_epsG(iG);
                
                for iBeta = 1:nStates_epsBeta
            
                    sim_epsBeta(:,2)   = rng_epsBeta(iBeta);
                
                    %sim_pii(:,2)   = piiSS;

                    %%% Fiscal limit
                    t = int64(2);
                    j = 0;
                    fiscalLimit = complex(zeros(nSimulations,1));

                    % Start values
                    sim_bbeta_old         = sim_bbeta(:,1);
                    sim_aTilde_old         = sim_aTilde(:,1);
                    %sim_aTildeExp_old      = sim_aTildeExp(:,1);
                    sim_g_old              = sim_g(:,1);
                    %sim_gExp_old           = sim_gExp(:,1);
                    sim_def_old            = sim_def(:,1);
                    %sim_defExp_old         = sim_defExp(:,1);
                    sim_recTilde_old       = sim_recTilde(:,1);
                    %sim_recTildeExp_old    = sim_recTildeExp(:,1);
                    %sim_psi_old            = sim_psi(:,1);
                    %sim_psiExp_old         = sim_psiExp(:,1);
                    sim_tauMax_old         = sim_tauMax(:,1);
                    sim_z_old              = sim_z(:,1);
                    sim_tLS_old            = sim_tLS(:,1);
                    %sim_nMax_old           = sim_nMax(:,1);
                    sim_yMax_old           = sim_yMax(:,1);
                    %sim_taxMax_old         = sim_taxMax(:,1);
                    %sim_cMax_old           = sim_cMax(:,1);
                    %sim_uCMax_old          = sim_uCMax(:,1);

                    % Fiscal limit variables
                    for tt=2:(nFiscalLimit+1)

                        sim_aTilde_new        = exp( (1 - rhoA).* log(aTildeSS) + rhoA.*log(sim_aTilde_old) + sigmaA.*sim_epsA(:,tt) );
                        sim_aTildeExp_new     = exp( (1 - rhoA).* log(aTildeSS) + rhoA.*log(sim_aTilde_old) );

                        sim_g_new             = exp( (1 - rhoGG).* log(gSS) + rhoGG.*log(sim_g_old) + rhoGY.*log(sim_yMax_old./ySS) + sigmaG.*sim_epsG(:,tt) );
                        
                        % Impose government ceiling
                        sim_g_new             = min(sim_g_new, govCeiling);
                        
                        %sim_gExp_new          = exp( (1 - rhoGG).* log(gSS) + rhoGG.*log(sim_g_old) + rhoGY.*log(sim_yMax_old./ySS) );

                        sim_kk_new            = kkSS .* (sim_g_new./gSS).^gammaGPSI;

                        sim_bbeta_new        = exp( (1 - rhoBeta).* log(beta) + rhoBeta.*log(sim_bbeta_old) + sigmaBeta.*sim_epsBeta(:,tt) );
                        
                        sim_def_new           = exp( (1 - rhoD).* log(defSS) + rhoD.*log(sim_def_old) + sigmaD.*sim_epsD(:,tt) );
                        sim_defExp_new        = exp( (1 - rhoD).* log(defSS) + rhoD.*log(sim_def_old) );

                        sim_recTilde_new      = exp( (1 - rhoR).* log(recTildeSS) + rhoR.*log(sim_recTilde_old) + sigmaR.*sim_epsR(:,tt) );
                        sim_recTildeExp_new   = exp( (1 - rhoR).* log(recTildeSS) + rhoR.*log(sim_recTilde_old) );

                        sim_psi_new       = (1 - sim_def_new).*sim_aTilde_new + sim_def_new.*sim_recTilde_new;
                        if psiShockObservedWithDelay
                             sim_psiExp_new    = (1 - sim_defExp_new).*sim_aTildeExp_new + sim_defExp_new.*sim_recTildeExp_new;  
                        else
                            sim_psiExp_new    = sim_psi_new; 
                        end

                        sim_tauMax_new    = sim_tauMax_old;
                        sim_z_new         = sim_z_old;
                        sim_tLS_new       = sim_tLS_old;

                        sim_nMax_new      = ( (1-sim_tauMax_new) .* sim_kk_new./eta .* sim_psiExp_new ).^(1/chi) ;
                        sim_wpMax_new     = (1./(1-sim_tauMax_new)) .* eta .* sim_nMax_new.^chi ;

                        sim_yMax_new      = sim_kk_new .* ( (1-sim_tauMax_new) .* sim_kk_new./eta .* sim_psiExp_new ).^(1./chi) .* sim_psi_new;                    
                        sim_taxMax_new    = sim_tauMax_new .* sim_kk_new .* ((1-sim_tauMax_new) .* sim_kk_new./eta .* sim_psiExp_new).^(1/chi) .* sim_psi_new;

                        sim_cMax_new      = sim_yMax_new - sim_g_new;
                        sim_cNRMax_new    = (1-sim_tauMax_new) .* sim_wpMax_new .* sim_nMax_new + sim_z_new;
                        sim_cRMax_new     = (1/(1-fracNR)) .* ( sim_cMax_new - fracNR .* sim_cNRMax_new );
                        sim_uCMax_new     = ( sim_cRMax_new + alphaG.*sim_g_new - eta .* sim_nMax_new.^(1+chi)./(1+chi) ).^(-sigma) ;

                        j = j + 1;
                        % Get marginal utility at the present
                        if t == tt
                            base_sim_uCMax = sim_uCMax_new;
                        end
                        fiscalLimit = fiscalLimit + (sim_bbeta_new.^double((j-1)) .* sim_uCMax_new./base_sim_uCMax .* (sim_taxMax_new - sim_g_new - sim_z_new - sim_tLS_new));

                        % Move forward
                        sim_bbeta_old         = sim_bbeta_new;
                        sim_aTilde_old         = sim_aTilde_new;
                        %sim_aTildeExp_old      = sim_aTildeExp_new;
                        sim_g_old              = sim_g_new;
                        %sim_gExp_old           = sim_gExp_new;
                        sim_def_old            = sim_def_new;
                        %sim_defExp_old         = sim_defExp_new;
                        sim_recTilde_old       = sim_recTilde_new;
                        %sim_recTildeExp_old    = sim_recTildeExp_new;
                        %sim_psi_old            = sim_psi_new;
                        %sim_psiExp_old         = sim_psiExp_new;
                        %sim_tauMax_old         = sim_tauMax_new;
                        sim_z_old              = sim_z_new;
                        sim_tLS_old            = sim_tLS_new;
                        %sim_nMax_old           = sim_nMax_new;
                        sim_yMax_old           = sim_yMax_new;
                        %sim_taxMax_old         = sim_taxMax_new;
                        %sim_cMax_old           = sim_cMax_new;
                        %sim_uCMax_old          = sim_uCMax_new;
                    end


    %                 for j=1:nFiscalLimit
    %                     fiscalLimit = fiscalLimit + (beta^double((j-1)) * sim_uCMax(:,t+j-1)./sim_uCMax(:,t) .* (sim_taxMax(:,t+j-1) - sim_g(:,t+j-1) - sim_z(:,t+j-1) - sim_tLS(:,t+j-1)));
    %                 end
                    fiscalLimit_dbl = real(fiscalLimit);
                    %%%

                    % Fit normal distribution to simulation data
                    %pdf_delta                       = fitdist(fiscalLimit,'Normal');
                    %fiscalLimit_mu(iA,iR,iD,iG)     = pdf_delta.mu;
                    %fiscalLimit_std(iA,iR,iD,iG)    = pdf_delta.sigma;
                    % Faster method
                    [fiscalLimit_mu(iA,iR,iD,iG,iBeta), fiscalLimit_std(iA,iR,iD,iG,iBeta)] = normfit(fiscalLimit_dbl);
                end
                toc
            end
        end
    end
end