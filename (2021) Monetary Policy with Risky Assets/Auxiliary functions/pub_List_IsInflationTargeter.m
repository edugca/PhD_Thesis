% Source: Wikipedia: https://en.wikipedia.org/wiki/Inflation_targeting
function response = pub_List_IsInflationTargeter(country, date)

tb_infTarget = pub_List_InflationTargeters();

idxCountry = find(strcmp(tb_infTarget.Country, country));
if isempty(idxCountry)
    response = zeros(size(date));
else
    response = ...
        date >= ...
        repmat(tb_infTarget.Adoption(idxCountry), length(date), 1);
end

end