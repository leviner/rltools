classdef cMRUData
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        heave
        roll
        pitch
        heading
        timestamp
    end
    
    methods
        function x = asStruct(obj)
            x = struct(...
                'heave', obj.heave, ...
                'roll', obj.roll, ...
                'pitch', obj.pitch, ...
                'heading', obj.heading);
        end
        
        function obj = read(obj,fid)
            obj.heave = fread(fid,1,'float32');
            obj.roll = fread(fid,1,'float32');
            obj.pitch = fread(fid,1,'float32');
            obj.heading = fread(fid,1,'float32');
        end
    end
    
end

