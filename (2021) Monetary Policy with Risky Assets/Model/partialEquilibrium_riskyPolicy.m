%% Wicksellian model augmented with policy-asset default

isFinal = false;

if isFinal
    imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Partial Equilibrium/';
else
    imagesFolder = 'Model/Images/';
end

%%% To replicate the Paper's figures, change phi (or phiPii) and set the respective policyTarget rule
phi         = 0.15; % price-level targeting: monetary policy parameter > 0
phiPii      = 1.15; % inflation targeting:monetary policy parameter > 1

policyTargetRule = 'inflation'; % priceLevel, inflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha       = 1/(1+phi);
alphaPii    = 1/phiPii;

nPeriods        = 100;
nFuture         = 100;
nSimulations    = 4;
p               = nan(nPeriods, nSimulations);
pii             = nan(nPeriods, nSimulations);
rN              = nan(nPeriods + nFuture,1);
iota            = nan(nPeriods + nFuture,1);
nrGap           = nan(nPeriods + nFuture,1);
defPol          = nan(nPeriods + nFuture,1);
nrBad           = nan(nPeriods + nFuture,1);

fontSize        = 12;
fontSizeLegend  = 12;
titleFontSize   = 14;

%% Partial equilibrium with policy default probability

%%% Steady state
iSimul        = 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
    
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
end
%%%%%%%%

%%% Steady state (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end

        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
        
    end

    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% rN

%%% Jump in rN
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            rN(j)        = 0.05;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
        
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rN (MIT shocks)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            rN(ii)        = 0.05;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rN (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.05;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            rN(j)        = 0.05;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rN (MIT shocks) (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.05;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            rN(ii)        = 0.05;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% iota

%%% Jump in iota
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            iota(j)        = 0.06;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in iota (MIT shocks)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            iota(ii)        = 0.06;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in iota (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            iota(j)        = 0.06;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in iota (MIT shocks) (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            iota(ii)        = 0.06;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% def

%%% Jump in def
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            defPol(j)        = 0.03;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in def (MIT shocks)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            defPol(ii)        = 0.03;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in def (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            defPol(j)        = 0.08;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in def (MIT shocks) (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            defPol(ii)        = 0.08;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% nrBad

%%% Jump in nrBad
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.40;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            nrBad(j)        = -0.50;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in nrBad (MIT shocks)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.00;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            nrBad(ii)        = 0.50;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in nrBad (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            nrBad(j)        = -0.50;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in nrBad (MIT shocks) (def>)
iSimul        = iSimul + 1;
rN(:)         = 0.04;
iota(:)       = 0.04;
defPol(:)     = 0.03;
nrBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            nrBad(ii)        = -0.50;
        end
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% Stochastic rN

%%% rN and def are stochastic; Mon pol does not track the natural interest rate
iSimul       = iSimul + 1;
rN(:)        = real(log(1 + 0.01*randn(length(rN),1)));
iota(:)      = 0;
defPol(:)    = max(0,0.03 + 0.01*randn(length(rN),1));
nrBad(:)      = -0.40;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        %
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% def is stochastic; Mon pol tracks the natural interest rate
iSimul       = iSimul + 1;
rN(:)        = rN(:); % use the same sequence of shocks
iota(:)      = rN(:);
defPol(:)    = defPol(:); % use the same sequence of shocks
nrBad(:)      = -0.40;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        %
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% rN stochastic; Mon pol tracks the risky natural interest rate
iSimul       = iSimul + 1;
rN(:)        = rN(:); % use the same sequence of shocks
defPol(:)    = defPol(:); % use the same sequence of shocks
nrBad(:)      = -0.40;
iota(:)      = (rN - defPol.*nrBad) ./ (1-defPol);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        %
        
        gamma = 1;
        gammaPii = 1;
        for k=1:(j-ii+1)
            gamma = gamma * 1/(1+(1-defPol(ii+k-1))*phi);
            gammaPii = gammaPii * 1/((1-defPol(ii+k-1))*phiPii);
        end
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = gamma * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = gammaPii * (rN(j) - (1-defPol(j))*iota(j) - defPol(j)*nrBad(j));
        end 
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%% Plot graphs

nPlotPeriods = 20;

f = figure;

previousInterpreter = pub_GraphSetInterpreter('Latex');

if strcmpi(policyTargetRule, 'priceLevel')
    yData = p;
    ruleLabel = '$p_t$   (\%)';
    simGraphName = ['graph_WicksellianTheoryDefault_', policyTargetRule, '_phi', replace(num2str(phi), '.', ''), '.png'];
elseif strcmpi(policyTargetRule, 'inflation')
    yData = pii;
    ruleLabel = '$\pi_t$   (\%)';
    simGraphName = ['graph_WicksellianTheoryDefault_', policyTargetRule, '_phi', replace(num2str(phiPii), '.', ''), '.png'];
end

tl = tiledlayout(3,2, 'tilespacing', 'compact', 'padding', 'compact');

iPlot = 1;
ax1 = nexttile();
graph = plot(100*[yData(:,1), yData(:,2)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("Steady state", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-10,50]);
pub_GraphDrawZeroAxis(graph);
legend('$\mathcal{D}^{Policy}_0$=0\%', '$\mathcal{D}^{Policy}_0$=3\%', ...
    'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

iPlot = 2;
ax2 = nexttile();
graph = plot(100*[yData(:,3), yData(:,4), yData(:,5), yData(:,6)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
graph(3).LineStyle = '-.';
graph(4).LineStyle = '--';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^n$ goes up by 5 p.p. from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-10,10]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated ($\mathcal{D}^{Policy}_0$=0\%)', 'MIT ($\mathcal{D}^{Policy}_0$=0\%)', ...
    'anticipated ($\mathcal{D}^{Policy}_0$=3\%)', 'MIT ($\mathcal{D}^{Policy}_0$=3\%)', ...
     'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

iPlot = 3;
ax3 = nexttile();
graph = plot(100*[yData(:,7), yData(:,8), yData(:,9), yData(:,10)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
graph(3).LineStyle = '-.';
graph(4).LineStyle = '--';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$\bar{\iota}$ goes up by 2 p.p. from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-10,50]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated ($\mathcal{D}^{Policy}_0$=0\%)', 'MIT ($\mathcal{D}^{Policy}_0$=0\%)', ...
    'anticipated ($\mathcal{D}^{Policy}_0$=3\%)', 'MIT ($\mathcal{D}^{Policy}_0$=3\%)', ...
     'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

iPlot = 4;
ax4 = nexttile();
graph = plot(100*[yData(:,11), yData(:,12), yData(:,13), yData(:,14)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
graph(3).LineStyle = '-.';
graph(4).LineStyle = '--';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$\mathcal{D}^{Policy}$ goes up by 3 p.p. from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('');
xlim([0,nPlotPeriods]);
ylim([-10,70]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated ($\mathcal{D}^{Policy}_0$=0\%)', 'MIT ($\mathcal{D}^{Policy}_0$=0\%)', ...
    'anticipated ($\mathcal{D}^{Policy}_0$=3\%)', 'MIT ($\mathcal{D}^{Policy}_0$=3\%)', ...
     'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

iPlot = 5;
ax5 = nexttile();
graph = plot(100*[yData(:,15), yData(:,16), yData(:,17), yData(:,18)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
graph(3).LineStyle = '-.';
graph(4).LineStyle = '--';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$i^{Bad}$ goes up by 10 p.p. from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-10,50]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated ($\mathcal{D}^{Policy}_0$=0\%)', 'MIT ($\mathcal{D}^{Policy}_0$=0\%)', ...
    'anticipated ($\mathcal{D}^{Policy}_0$=3\%)', 'MIT ($\mathcal{D}^{Policy}_0$=3\%)', ...
     'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

iPlot = 6;
ax6 = nexttile();
graph = plot(100*[yData(:,19), yData(:,20), yData(:,21)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
graph(3).LineStyle = '-.';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^n$ and $\mathcal{D}^{Policy}$ are independent stoch. processes", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-10,30]);
pub_GraphDrawZeroAxis(graph);
legend('fixed intercept (ignoring risk)', '$r^n$ in the intercept', 'optimal intercept', ...
     'Location','southoutside', 'FontSize', fontSizeLegend, 'Interpreter', 'latex', 'NumColumns', 2);
grid 'on';

%linkaxes([ax1, ax2, ax3, ax4, ax5, ax6]);

%edu_Suptitle(['Neo-Wicksellian theory of the price level with policy asset default ($\phi$ = ', num2str(phi), ')'], ...
%    'FontSize', 20, 'Interpreter', 'latex');

pub_GraphSetInterpreter(previousInterpreter);

set(gcf, 'Position',  [100, 100, 1000, 800]); % resize figure
saveas(f,[imagesFolder, simGraphName]);
