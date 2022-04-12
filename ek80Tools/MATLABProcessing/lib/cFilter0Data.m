classdef cFilter0Data
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        channel
        channelid
        filterdata
    end
    
    methods
        function obj = read(obj,fid,length)
            obj.channel     = fread(fid,1,'int16');
            obj.channelid   = char(fread(fid,128,'char')');
            obj.filterdata  = char(fread(fid,length-128-2,'char')');
        end
    end
    
end

