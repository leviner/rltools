classdef cHeader
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        datetime
    end
    
    methods
        function obj = read(obj,fid)
            obj.type = char(fread(fid,4,'char')');
            
            lowdatetime = fread(fid,1,'uint32');
            highdatetime = fread(fid,1,'uint32');
            
            obj.datetime = NTTime2Mlab(highdatetime*2^32 + lowdatetime);
        end
    end
    
end

