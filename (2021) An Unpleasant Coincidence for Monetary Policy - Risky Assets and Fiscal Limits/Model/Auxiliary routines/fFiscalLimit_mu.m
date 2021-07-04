function mu = fFiscalLimit_mu(arg1, arg2, arg3, arg4, arg5, varargin)

global fitModel_mu

%diffOrder = 0;

%f_mu  = evalin('base', 'fitModel_mu.PolynomialExpression');
f_mu  = fitModel_mu(arg1,arg2,arg3,arg4,arg5);

%if isa(arg1, 'double')
    mu = f_mu;
%else
%     if ~isempty(varargin)
%         for ii=1:lentgh(varargin)
%             if strcmp(varargin(ii), 'diff')
%                 diffOrder = varargin(ii+1);
%             end
%         end
%     end
%     
%     fSymbolic = sym(f_mu);
%     f_mu_sym = symfun(fSymbolic, symvar(fSymbolic));
%     
%     if diffOrder == 0
%         mu = f_mu_sym(arg1,arg2,arg3,arg4);
%     else
%         mu = f_mu_sym(arg1,arg2,arg3,arg4);
%     end
% end

end