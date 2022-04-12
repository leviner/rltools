% Parse Parameter XML data to a Parameter Matlab structure

function paramdata = parseparamxmlstruct(xmldata)

nchildren = length(xmldata.Children);

for i = 1:nchildren
    nattributes = length(xmldata.Children(i).Attributes);
    for j = 1:nattributes
        name = xmldata.Children(i).Attributes(j).Name;
        paramdata(i).(name) = xmldata.Children(i).Attributes(j).Value;
        if ~strcmp(xmldata.Children(i).Attributes(j).Name, 'ChannelID')
            paramdata(i).(name) = str2double(paramdata(i).(name));
        end
    end
end
