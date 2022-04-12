% Parse Sensor XML data to a Parameter Matlab structure
%
%  New 7/2020
%
% xmldata
%  
%           Name: 'Sensor'
%     Attributes: [1×3 struct]
%           Data: ''
%       Children: []
%       
%  xmldata.Attributes
%     Name
%     Value  
%  
% xmldata.Attributes.Name
%     'IsManual'
%     'ManualValue'
%     'Type'
%     
    


function paramdata = parseSensorxmlstruct(xmldata)

    
    nattributes = length(xmldata.Attributes);
    
    
    for j = 1:nattributes
        name = xmldata.Attributes(j).Name;
        paramdata.(name) = xmldata.Attributes(j).Value;
        
        % change these to numbers
        if strcmp(xmldata.Attributes(j).Name, 'IsManual')
            paramdata.(name) = str2double(paramdata.(name));
        end
        if strcmp(xmldata.Attributes(j).Name, 'ManualValue')
            paramdata.(name) = str2double(paramdata.(name));
        end      
        
    end
    
end   