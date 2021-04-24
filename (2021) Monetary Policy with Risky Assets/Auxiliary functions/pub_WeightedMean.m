function result = pub_WeightedMean(values, weights, dim)

wMean = NaN(size(values,1), 1);

if dim == 2
    if size(weights,2) == 1; weights = weights'; end
        
    for iLin = 1:size(values,1)
        validWeights = find(~isnan(values(iLin,:)));
        wMean(iLin) = sum(values(iLin,:) .* weights ./ sum(weights(validWeights), 'omitnan'), 'omitnan');
    end
end

result = wMean;

end