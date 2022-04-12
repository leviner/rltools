classdef cConfig
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        surveyname
        transectname
        soundername
        spare
        transducercount
        
        configtransducer = cConfigTransducer;
    end
    
    methods
        function obj = read(obj,fid)
            obj.surveyname      = char(fread(fid,128,'char')');
            obj.transectname    = char(fread(fid,128,'char')');
            obj.soundername     = char(fread(fid,128,'char')');
            obj.spare           = char(fread(fid,128,'char')');
            obj.transducercount = fread(fid,1,'int32');
            
            obj.configtransducer(obj.transducercount,1) = cConfigTransducer;
            for i=1:obj.transducercount,
                obj.configtransducer(i) = cConfigTransducer;
                obj.configtransducer(i) = obj.configtransducer(i).read(fid);
            end
        end
    end
    
end

