function pdf = pickFiscalLimit(simStruct, fiscalLimit_mu, fiscalLimit_std, ...
    v_epsA, v_epsR, v_epsD, v_epsG, v_epsBeta)

interpValues = true;

nStates = 11;
rng_epsA        = pub_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaA,simStruct.params.rhoA,"Tauchen1986");
rng_epsR        = pub_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaR,simStruct.params.rhoR,"Tauchen1986");
rng_epsD        = pub_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaD,simStruct.params.rhoD,"Tauchen1986");
rng_epsG        = pub_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaG,simStruct.params.rhoGG,"Tauchen1986");
rng_epsBeta     = pub_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaBeta,simStruct.params.rhoBeta,"Tauchen1986");

[~, idxA]       = min(abs(rng_epsA - v_epsA));
[~, idxR]       = min(abs(rng_epsR - v_epsR));
[~, idxD]       = min(abs(rng_epsD - v_epsD));
[~, idxG]       = min(abs(rng_epsG - v_epsG));
[~, idxBeta]    = min(abs(rng_epsBeta - v_epsBeta));

if interpValues
    
    if any(rng_epsA)
        interpA = interp1(rng_epsA - v_epsA, 1:nStates, v_epsA, 'linear', 'extrap');
    else
        interpA = ceil(length(rng_epsA)/2);
    end
    
    if any(rng_epsR)
        interpR     = interp1(rng_epsR - v_epsR, 1:nStates, v_epsR, 'linear', 'extrap');
    else
        interpR = ceil(length(rng_epsR)/2);
    end
    
    if any(rng_epsD)
        interpD     = interp1(rng_epsD - v_epsD, 1:nStates, v_epsD, 'linear', 'extrap');
    else
        interpD = ceil(length(rng_epsD)/2);
    end
    
    if any(rng_epsG)
        interpG     = interp1(rng_epsG - v_epsG, 1:nStates, v_epsG, 'linear', 'extrap');
    else
        interpG = ceil(length(rng_epsG)/2);
    end
    
    if any(rng_epsBeta)
        interpBeta  = interp1(rng_epsBeta - v_epsBeta, 1:nStates, v_epsBeta, 'linear', 'extrap');
    else
        interpBeta = ceil(length(rng_epsBeta)/2);
    end
end

% Verify the size of the vector
if length(fiscalLimit_mu) > 1
    
    if interpValues
        mu_A = interp1(1:nStates, vec(fiscalLimit_mu(:,idxR,idxD,idxG,idxBeta)), interpA, 'spline', 'extrap');
        mu_R = interp1(1:nStates, vec(fiscalLimit_mu(idxA,:,idxD,idxG,idxBeta)), interpR, 'spline', 'extrap');
        mu_D = interp1(1:nStates, vec(fiscalLimit_mu(idxA,idxR,:,idxG,idxBeta)), interpD, 'spline', 'extrap');
        mu_G = interp1(1:nStates, vec(fiscalLimit_mu(idxA,idxR,idxD,:,idxBeta)), interpG, 'spline', 'extrap');
        mu_Beta = interp1(1:nStates, vec(fiscalLimit_mu(idxA,idxR,idxD,idxG, :)), interpBeta, 'spline', 'extrap');
        
        std_A = interp1(1:nStates, vec(fiscalLimit_std(:,idxR,idxD,idxG,idxBeta)), interpA, 'spline', 'extrap');
        std_R = interp1(1:nStates, vec(fiscalLimit_std(idxA,:,idxD,idxG,idxBeta)), interpR, 'spline', 'extrap');
        std_D = interp1(1:nStates, vec(fiscalLimit_std(idxA,idxR,:,idxG,idxBeta)), interpD, 'spline', 'extrap');
        std_G = interp1(1:nStates, vec(fiscalLimit_std(idxA,idxR,idxD,:,idxBeta)), interpG, 'spline', 'extrap');
        std_Beta = interp1(1:nStates, vec(fiscalLimit_std(idxA,idxR,idxD,idxG,:)), interpBeta, 'spline', 'extrap');
        
        pdf = makedist('Normal','mu', mean([mu_A,mu_R,mu_D,mu_G,mu_Beta]), ...
                        'sigma', mean([std_A,std_R,std_D,std_G,std_Beta]));
    else
        
        pdf = makedist('Normal','mu', fiscalLimit_mu(idxA,idxR,idxD,idxG,idxBeta), ...
                        'sigma',fiscalLimit_std(idxA,idxR,idxD,idxG,idxBeta));
                
    end
    
else
    pdf = makedist('Normal','mu', fiscalLimit_mu, ...
                    'sigma',fiscalLimit_std);
end

end