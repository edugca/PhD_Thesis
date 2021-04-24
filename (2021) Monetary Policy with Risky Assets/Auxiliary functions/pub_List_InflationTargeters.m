% Source: Wikipedia: https://en.wikipedia.org/wiki/Inflation_targeting
function tb = pub_List_InflationTargeters()

% Name, Date of Adoption of Inflation Targeting
cellArray = ...
    {
    'New Zealand',      datetime(1989,12,1);
    'Canada',           datetime(1991,2,1);
    'United Kingdom',	datetime(1992,10,1);
    'Sweden',           datetime(1993,1,1);
    'Australia',        datetime(1993,6,1);
    'Israel',           datetime(1997,6,1);
    'Czech Republic',	datetime(1997,12,1);
    'Poland',           datetime(1998,1,1); % Just the year
    'Brazil',           datetime(1999,6,1);
    'Chile',            datetime(1999,9,1);
    'Colombia',         datetime(1999,10,1);
    'South Africa',     datetime(2000,2,1);
    'Thailand',         datetime(2000,5,1);
    'Mexico',           datetime(2001,1,1); % Just the year
    'Norway',           datetime(2001,3,1);
    'Iceland',          datetime(2001,3,1);
    'Peru',             datetime(2002,1,1);
    'Philippines',      datetime(2002,1,1);
    'Guatemala',        datetime(2005,1,1);
    'Indonesia',        datetime(2005,7,1);
    'Romania',          datetime(2005,8,1);
    'Armenia',          datetime(2006,1,1);
    'Turkey',           datetime(2006,1,1);
    'Ghana',            datetime(2007,5,1);
    'Georgia',          datetime(2009,1,1);
    'Serbia',                   datetime(2009,1,1);
    'United States',            datetime(2012,1,1);
    'Japan',                    datetime(2013,1,1);
    'Russian Federation',       datetime(2014,1,1);
    'Republic of Kazakhstan',	datetime(2015,8,1);
    'India',                    datetime(2016,8,1);
    'Argentina',                datetime(2016,9,1);
    };

tb = cell2table(cellArray, 'VariableNames',{'Country', 'Adoption'});

end