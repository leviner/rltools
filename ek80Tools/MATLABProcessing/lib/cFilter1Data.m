classdef cFilter1Data
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Stage
        channel
        ChannelID
        Decimation
        FilterData
    end
    
    methods
        function obj = read(obj,fid)
            obj.Stage       = fread(fid,1,'int16');
            obj.channel     = fread(fid,1,'int16');
            obj.ChannelID   = deblank((fread(fid,128,'*char')'));
            ncoeff          = fread(fid,1,'int16');
            obj.Decimation  = fread(fid,1,'int16');
            % The real and imag filter coefficients are stored as:
            % real0 imag0 real1 imag1 ... realN imagN
            coeff           = fread(fid,[2, ncoeff],'float32'); 
            obj.FilterData  = coeff(1,:)+1i*coeff(2,:);
        end
    end
    
end

