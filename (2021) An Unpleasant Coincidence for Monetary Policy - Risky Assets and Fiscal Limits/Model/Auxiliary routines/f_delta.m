function result = f_delta(x, epsA, epsR, epsD, epsG, epsBeta)

dist = pdf_delta(epsA, epsR, epsD, epsG, epsBeta);
            
result = normcdf(x, ...
    dist.mu, ...
    dist.std);

end