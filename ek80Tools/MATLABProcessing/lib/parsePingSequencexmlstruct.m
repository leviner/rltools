% Parse PingSequence XML data to a Parameter Matlab structure
%
%  New 7/2020
%
%  xmldata  
%  
%           Name: 'PingSequence'
%     Attributes: []
%           Data: ''
%       Children: [1×1 struct] 
%       
%       
% xmldata.Children
%           Name: 'Ping'
%     Attributes: [1×1 struct]
%           Data: ''
%       Children: [] 
%       
% xmldata.Children.Attributes
%     Name: 'ChannelID'
%     Value: 'WBT Tube 247715-15 ES70-7CD_ES'
%       


function paramdata = parsePingSquencexmlstruct(xmldata)

    nattributes = length(xmldata.Children.Attributes);
    
    for j = 1:nattributes
        name = xmldata.Children.Attributes(j).Name;
        paramdata.(name) = xmldata.Children.Attributes(j).Value;      
    end
    
end   