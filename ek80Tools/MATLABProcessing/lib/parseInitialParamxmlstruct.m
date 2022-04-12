% Parse InitialParameter XML data to a Parameter Matlab structure
%
%  New 7/2020
%
% xmldata = 
%           Name: 'InitialParameter'
%     Attributes: []
%           Data: ''
%       Children: [1×1 struct]
%       
% xmldata.Children
%           Name: 'Channels'
%     Attributes: []
%           Data: ''
%       Children: [1×1 struct] 
%       
% xmldata.Children.Children
%           Name: 'Channel'
%     Attributes: [1×11 struct]
%           Data: ''
%       Children: []
%
% xmldata.Children.Children.Attributes
%     Name
%     Value
%     
% xmldata.Children.Children.Attributes.Name
%     'ChannelID'
%     'ChannelMode'
%     'FrequencyEnd'
%     'FrequencyStart'
%     'PingId' 
%     'PulseDuration' 
%     'PulseForm' 
%     'SampleInterval' 
%     'Slope' 
%     'SoundVelocity' 
%     'TransmitPower'      
%       
%       
% 
%
%



function paramdata = parseInitialParamxmlstruct(xmldata)


if xmldata.Children.Children.Name == 'Channel'
    
    nattributes = length(xmldata.Children.Children.Attributes);
    
    
    for j = 1:nattributes
        name = xmldata.Children.Children.Attributes(j).Name;
        paramdata.(name) = xmldata.Children.Children.Attributes(j).Value;
        if ~strcmp(xmldata.Children.Children.Attributes(j).Name, 'ChannelID')
            paramdata.(name) = str2double(paramdata.(name));
        end
        
    end
    
end   % if