classdef cNMEA
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        text
        timestamp
    end
    
    methods
        function obj = read(obj,fid,length)
            obj.text = char(fread(fid,length,'char')');
        end
    end
    
end

