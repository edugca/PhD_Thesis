%% rN and def are stochastic; Mon pol does not track the natural interest rate

disp([num2str(iCorrel), ' ', num2str(iPhi)]);
disp('rN and def are stochastic; Mon pol does not track the natural interest rate');
iExp = iExp + 1;
expNames = [expNames {'Intercept is fixed at $\overline{r^{n}}$'}];
for iSimul = 1:nSimulations
    shocks_1 = randn(size(rN,1),1);
    shocks_2 = correl_rN_defPol*shocks_1 + (1-correl_rN_defPol)*randn(size(rN,1),1);
    rN(:, iSimul)        = real(log(1 + 0.02*shocks_1));
    iota(:, iSimul)      = 0;
    defPol(:, iSimul)    = max(0,0.03 + 0.02*shocks_2);
    nrBad(:, iSimul)     = -0.40;
    ddelta(:, iSimul)     = 0.40;

    nrGap(:, iSimul) = 0;
    for ii=1:nPeriods
        for j=ii:(ii+nFuture)
            % Shock rule
            %

            if strcmp(policyTargetRule, 'priceLevel')
                gamma = 1;
                
                for k=1:(j-ii+1)
                        gamma = gamma * 1/(1+(1-defPol(ii+k-1, iSimul)*ddelta(ii+k-1, iSimul))*phi);
                end
                
                nrGap(j, iSimul)   = gamma * (rN(j, iSimul) - (1-defPol(j, iSimul)*ddelta(j, iSimul))*iota(j, iSimul) + defPol(j, iSimul)*ddelta(j, iSimul));
            elseif strcmp(policyTargetRule, 'inflation')
                gammaPii = 1;
                
                for k=1:(j-ii+1)
                        gammaPii = gammaPii * 1/((1-defPol(ii+k-1, iSimul)*ddelta(ii+k-1, iSimul))*phiPii);
                end
                
                nrGap(j, iSimul)   = gammaPii * (rN(j, iSimul) - (1-defPol(j, iSimul)*ddelta(j, iSimul))*iota(j, iSimul) + defPol(j, iSimul)*ddelta(j, iSimul));
            end
            
        end
       
        if strcmp(policyTargetRule, 'priceLevel')
            p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                pii(ii, iSimul) =  p(ii, iSimul) - 0;
            else
                pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
            end
        elseif strcmp(policyTargetRule, 'inflation')
            pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                p(ii, iSimul) =  pii(ii, iSimul) - 0;
            else
                p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
            end
        end
       
    end
end

p_Vec(num2str(thisPhi))       = [p_Vec(num2str(thisPhi)), p(:)];
pii_Vec(num2str(thisPhi))     = [pii_Vec(num2str(thisPhi)), pii(:)];
rN_Vec(num2str(thisPhi))      = [rN_Vec(num2str(thisPhi)),  reshape(rN(1:nPeriods, :), [], 1)];
defPol_Vec(num2str(thisPhi))  = [defPol_Vec(num2str(thisPhi)), reshape(defPol(1:nPeriods, :), [], 1)];

%% rN and def are stochastic; Mon pol tracks the natural interest rate

disp('rN and def are stochastic; Mon pol tracks the natural interest rate');
iExp = iExp + 1;
expNames = [expNames {'Intercept tracks $r^n_t$'}];
for iSimul = 1:nSimulations
    rN(:, iSimul)         = rN(:, iSimul); % use the same sequence of shocks
    iota(:, iSimul)       = rN(:, iSimul);
    defPol(:, iSimul)     = defPol(:, iSimul); % use the same sequence of shocks
    nrBad(:, iSimul)      = -0.40;
    ddelta(:, iSimul)     = 0.40;
    
    nrGap(:, iSimul) = 0;
    for ii=1:nPeriods
        for j=ii:(ii+nFuture)
            % Shock rule
            %
            
            if strcmp(policyTargetRule, 'priceLevel')
                gamma = 1;
                
                for k=1:(j-ii+1)
                        gamma = gamma * 1/(1+(1-defPol(ii+k-1, iSimul)*ddelta(ii+k-1, iSimul))*phi);
                end
                
                nrGap(j, iSimul)   = gamma * (rN(j, iSimul) - (1-defPol(j, iSimul)*ddelta(j, iSimul))*iota(j, iSimul) + defPol(j, iSimul)*ddelta(j, iSimul));
            elseif strcmp(policyTargetRule, 'inflation')
                gammaPii = 1;
                
                for k=1:(j-ii+1)
                        gammaPii = gammaPii * 1/((1-defPol(ii+k-1, iSimul)*ddelta(ii+k-1, iSimul))*phiPii);
                end
                
                nrGap(j, iSimul)   = gammaPii * (rN(j, iSimul) - (1-defPol(j, iSimul)*ddelta(j, iSimul))*iota(j, iSimul) + defPol(j, iSimul)*ddelta(j, iSimul));
            end
            
        end
       
       if strcmp(policyTargetRule, 'priceLevel')
            p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                pii(ii, iSimul) =  p(ii, iSimul) - 0;
            else
                pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
            end
        elseif strcmp(policyTargetRule, 'inflation')
            pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                p(ii, iSimul) =  pii(ii, iSimul) - 0;
            else
                p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
            end
        end
       
    end
end


p_Vec(num2str(thisPhi))       = [p_Vec(num2str(thisPhi)), p(:)];
pii_Vec(num2str(thisPhi))     = [pii_Vec(num2str(thisPhi)), pii(:)];
rN_Vec(num2str(thisPhi))      = [rN_Vec(num2str(thisPhi)),  reshape(rN(1:nPeriods, :), [], 1)];
defPol_Vec(num2str(thisPhi))  = [defPol_Vec(num2str(thisPhi)), reshape(defPol(1:nPeriods, :), [], 1)];

%% rN and def are stochastic; Mon pol tracks rN adjusted to underlying risk

disp('rN and def are stochastic; Mon pol tracks rN adjusted to underlying default risk');
iExp = iExp + 1;
expNames = [expNames {'Optimal intercept'}];
for iSimul = 1:nSimulations

    rN(:, iSimul)        = rN(:, iSimul); % use the same sequence of shocks
    defPol(:, iSimul)    = defPol(:, iSimul); % use the same sequence of shocks
    nrBad(:, iSimul)      = -0.40;
    ddelta(:, iSimul)     = 0.40;
    iota(:, iSimul)      = (rN(:,iSimul) - defPol(:,iSimul).*nrBad(:,iSimul)) ./ (1-defPol(:,iSimul));

    nrGap(:, iSimul) = 0;
    for ii=1:nPeriods
        for j=ii:(ii+nFuture)
            % Shock rule
            %

            if strcmp(policyTargetRule, 'priceLevel')
                gamma = 1;
                
                for k=1:(j-ii+1)
                        gamma = gamma * 1/(1+(1-defPol(ii+k-1, iSimul))*phi);
                end
                
                nrGap(j, iSimul)   = gamma * (rN(j, iSimul) - (1-defPol(j, iSimul))*iota(j, iSimul) - defPol(j, iSimul)*nrBad(j, iSimul));
            elseif strcmp(policyTargetRule, 'inflation')
                gammaPii = 1;
                
                for k=1:(j-ii+1)
                        gammaPii = gammaPii * 1/((1-defPol(ii+k-1, iSimul))*phiPii);
                end
                
                nrGap(j, iSimul)   = gammaPii * (rN(j, iSimul) - (1-defPol(j, iSimul))*iota(j, iSimul) - defPol(j, iSimul)*nrBad(j, iSimul));
            end
            
        end
       
       if strcmp(policyTargetRule, 'priceLevel')
            p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                pii(ii, iSimul) =  p(ii, iSimul) - 0;
            else
                pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
            end
        elseif strcmp(policyTargetRule, 'inflation')
            pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture, iSimul));
            if ii == 1
                p(ii, iSimul) =  pii(ii, iSimul) - 0;
            else
                p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
            end
        end
       
    end
end

p_Vec(num2str(thisPhi))       = [p_Vec(num2str(thisPhi)), p(:)];
pii_Vec(num2str(thisPhi))     = [pii_Vec(num2str(thisPhi)), pii(:)];
rN_Vec(num2str(thisPhi))      = [rN_Vec(num2str(thisPhi)),  reshape(rN(1:nPeriods, :),[],1)];
defPol_Vec(num2str(thisPhi))  = [defPol_Vec(num2str(thisPhi)), reshape(defPol(1:nPeriods, :),[],1)];