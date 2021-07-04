% Check steady state
%switchingParamsNames = {mdlVector(ii).parameters.name{mdlVector(ii).parameters.is_switching}};
%edu_RISE_checkSteadyState(ssFileAddress, paramsStruct, ...
%    riseConfigVector(ii), mdlVector(ii).markov_chains.state_names, ...
%    switchingParamsNames);


% "RISE is allergic to cases where the variables entering the
%  transition probabilities do not have a unique steady state. 
%  If the variables entering the transition probabilities
%  have multiple steady states, the transition probabilities will
%  not be unique. What that means is that we will not know how to
%  form expectations because for a given state of the economy, we
%  have different probabilities for a given event to occur.
%  When you impose a steady state you do resolve the confusion
% above but at the same time you create a new set of issues.
% Imposing a steady state implies that the residuals evaluated in
% a given regime might be nonzero. But RISE still needs to find a
% solution to the problem such that for the dynamic model we still
% have E[f(x_{t+1},x_t,x_{t-1},e_t)]=0 even if for the static model
% f(x_t,x_t,x_t,0) is different from zero.
% Loosely speaking, RISE approximates the system around a point that
% is not the steady state. RISE can only achieve this through solving
% for a constant that compensates for the fact that the
% regime-specific steady state residuals are different from zero:
% This is the role of @sig, which then appears already at a
% first-order approximation because in that case certainty equivalence
% does not hold. You can easily check this by running exogenous
% probabilities without imposing the steady state.
% Summing up, if you want to use endogenous probabilities, you will
% have to ensure the variables entering the time-varying transition
% probabilities have a unique steady state. But this is mostly a
% weakness of the solution approach rather than a fundamental property
% of the problem: if we knew how to solve the problem exactly, we would
% not need to approximate it and hence we would not need to refer to
% the steady state.

% (a) imposed(default=false): RISE computes the solution at the
%                   specified point without checking that the point solves for the
%                   steady state
% (b) unique (default=false): RISE computes the steady state at
%                   the ergodic distribution of the parameters. In case the
%                   probabilities are endogenous, the ergodic distribution of the
%                   parameters is itself a function of the steady state of the
%                   variables.
% (c) loop(default=false): RISE considers the equations
%                   calculating the steady state as true and just solves for the
%                   missing variables by looping over the steady state program. The
%                   user can then use the values pushed into the steady state
%                   program to calculate the steady state for the included
%                   variables.