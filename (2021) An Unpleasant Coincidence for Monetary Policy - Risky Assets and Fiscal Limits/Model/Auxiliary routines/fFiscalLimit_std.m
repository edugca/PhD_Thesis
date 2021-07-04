function std = fFiscalLimit_std(x1,x2,x3,x4,x5, varargin)

global fitModel_std

%f_std  = evalin('base', 'fitModel_std.PolynomialExpression');
f_std  = fitModel_std(x1,x2,x3,x4,x5);

%if isa(x1, 'double')
    std = f_std;
% else
%    f_std_sym = sym(f_std);
%    std = f_std_sym(x1,x2,x3,x4);
% 
% end

end