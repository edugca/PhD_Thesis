function [grid, transitionMatrix] = pub_arToMarkov(points, nStdDev, distributionName, mean, valStdDev, arCoeff, method)

% I use the method of Tauchen (1986) to covert the AR process into a
% finite discrete Markov Chain
%
% points            = number of points in the grid
% nStdDev           = number of standard deviations to each side for the grid
% distributionName  = name of the distribution of the process
% mean              = mean value of the process
% valStdDev         = value of standard deviation
% arCoeff           = AR(1) coefficient
% method            = default is Tauchen 1986 ("Tauchen1986", "TauchenHussey1991")

if isempty(method) || method == "Tauchen1986"
    [grid, transitionMatrix] = pub_arToMarkov_Tauchen1986(points, nStdDev, distributionName, mean, valStdDev, arCoeff);
elseif isempty(method) || method == "Tauchen1986_SPLANAR"
    [grid, transitionMatrix] = pub_arToMarkov_Tauchen1986_SPLANAR(points, nStdDev, distributionName, mean, valStdDev, arCoeff);
elseif method == "TauchenHussey1991"
    [grid, transitionMatrix] = pub_arToMarkov_TauchenHussey1991(points, nStdDev, distributionName, mean, valStdDev, arCoeff);
end

end


function [grid, transitionMatrix] = pub_arToMarkov_Tauchen1986(points, nStdDev, distributionName, mean, valStdDev, arCoeff)

% I use the method of Tauchen (1986) to covert the AR process into a
% finite discrete Markov Chain
% z = arCoeff*z(-1) + eps ~ D(mean, valStdDev)
%
% points            = number of points in the grid
% nStdDev           = number of standard deviations to each side for the grid
% distributionName  = name of the distribution of the process
% mean              = mean value of the process
% valStdDev         = value of standard deviation
% arCoeff           = AR(1) coefficient


% Create the grids
zstar = mean/(1-arCoeff); % expected value of z
sigmaz = valStdDev/sqrt(1-arCoeff^2); %stddev of z

minLim = -nStdDev*sigmaz;
maxLim = nStdDev*sigmaz;
interval = (maxLim-minLim)/(points-1); % Note that all the points are equidistant by construction.

if interval ~= 0
    grid = zstar*ones(points,1) + ( minLim:interval:maxLim )';
else
    grid = zstar*ones(points,1) ;
end

% Create the transition matrix
transitionMatrix = NaN(points,points);
for i=1:points
	for j=1:points
		if j==1		
			transitionMatrix(i,j) = ...
					cdf(distributionName, grid(j)-arCoeff*grid(i) + interval/2, mean, valStdDev);
        elseif j==points
			transitionMatrix(i,j) = ...
					1 - cdf(distributionName, grid(j)-arCoeff*grid(i) - interval/2, mean, valStdDev);
        else
			transitionMatrix(i,j) = ...
					cdf(distributionName, grid(j)-arCoeff*grid(i) + interval/2, mean, valStdDev) ...
					- cdf(distributionName, grid(j)-arCoeff*grid(i) - interval/2, mean, valStdDev);
        end
    end
end

end

function [grid, transitionMatrix] = pub_arToMarkov_Tauchen1986_SPLANAR(points, nStdDev, distributionName, mean, valStdDev, arCoeff)

% SPLANAR VERSION
% I use the method of Tauchen (1986) to covert the AR process into a
% finite discrete Markov Chain
% z = arCoeff*z(-1) + eps ~ D(mean, valStdDev)
%
% points            = number of points in the grid
% nStdDev           = number of standard deviations to each side for the grid
% distributionName  = name of the distribution of the process
% mean              = mean value of the process
% valStdDev         = value of standard deviation
% arCoeff           = AR(1) coefficient

% Create the grids
zstar = mean/(1-arCoeff); % expected value of z
sigmaz = valStdDev/sqrt(1-arCoeff^2); %stddev of z

minLim = -nStdDev*sigmaz;
maxLim = nStdDev*sigmaz;
interval = (maxLim-minLim)/(points-1); % Note that all the points are equidistant by construction.

grid = zstar + minLim;
for ii = 2:points
    grid = [grid; zstar + minLim + interval*(ii-1)];
end

% Create the transition matrix
transitionMatrix(points,points) = splanar;
for i=1:points
	for j=1:points
		if j==1		
			transitionMatrix(i,j) = ...
					cdf(distributionName, grid(j)-arCoeff*grid(i) + interval/2, mean, valStdDev);
        elseif j==points
			transitionMatrix(i,j) = ...
					1 - cdf(distributionName, grid(j)-arCoeff*grid(i) - interval/2, mean, valStdDev);
        else
			transitionMatrix(i,j) = ...
					cdf(distributionName, grid(j)-arCoeff*grid(i) + interval/2, mean, valStdDev) ...
					- cdf(distributionName, grid(j)-arCoeff*grid(i) - interval/2, mean, valStdDev);
        end
    end
end

end


function [grid, transitionMatrix] = pub_arToMarkov_TauchenHussey1991(points, nStdDev, distributionName, mean, valStdDev, arCoeff)

%%% INCOMPLETE!!!!!

% I use the method of Tauchen and Hussey (1991) to covert the AR process into a
% finite discrete Markov Chain using Gauss-Hermite quadrature
%
% points            = number of points in the grid
% nStdDev           = number of standard deviations to each side for the grid
% distributionName  = name of the distribution of the process
% arCoeff           = AR(1) coefficient
% mean              = mean value of the process
% valStdDev         = value of standard deviation
% method            = default is Tauchen 1986 ("Tauchen1986", "TauchenHussey1991")

% Size of each grid interval
%interval = (nStdDev*valStdDev-(-nStdDev*valStdDev))/(points-1);

% Hermitean roots
syms x;
gridHermitean = sort(real(eval(roots(fliplr(coeffs(hermiteH(points,x),'All'))))));

% Create the grids
grid = mean + sqrt(2)*valStdDev*gridHermitean;

% Create the transition matrix
transitionMatrix = NaN(points,points);
for i=1:points
	for j=1:points
		if j==1		
			transitionMatrix(i,j) = ...
					cdf(distributionName, grid(j)-arCoeff*grid(i) + interval/2, mean, valStdDev);
        elseif j==points
			transitionMatrix(i,j) = ...
					1 - cdf(distributionName, grid(j)-arCoeff*grid(i) - interval/2, mean, valStdDev);
        else
			transitionMatrix(i,j) = ...
					1/sqrt(pi) * ...
                    gridHermitean(j) * ...
                    pdf(distributionName, grid(j)-arCoeff*grid(i), mean, valStdDev) / ...
                    pdf(distributionName, grid(j)-arCoeff*mean, mean, valStdDev);
        end
    end
end

% Normalize the probabilities
transitionMatrix = transitionMatrix ./ repmat(sum(transitionMatrix, 2), 1, points);

end