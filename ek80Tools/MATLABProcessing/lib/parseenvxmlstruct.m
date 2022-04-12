% Parse Environment XML data to an Environment Matlab structure
% Simrad, Lars Nonboe Andersen, 10/10-13

function envdata = parseenvxmlstruct(xmldata)

nattributes = length(xmldata.Attributes);

for i = 1:nattributes,
    if ~strcmp(xmldata.Attributes(i).Name, 'Copyright')
        envdata.(xmldata.Attributes(i).Name) = str2double(xmldata.Attributes(i).Value);
    else
        envdata.(xmldata.Attributes(i).Name) = xmldata.Attributes(i).Value;
    end
end

% Make salinity have units of PSU.
for j = 1:length(envdata)
    % If the salinity value looks like 'per 1', convert back to
    % parts per thousand (or perhaps PSU).
    if envdata(j).Salinity < 0.1 && envdata(j).Salinity > 0
        envdata(j).Salinity = envdata(j).Salinity * 1000; %#ok<AGROW>
    end
end
