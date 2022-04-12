classdef cSampleDataRAW3_DT
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        timestamp % [Matlab time number of ping]
        channelID % 
        datatype
        spare
        offset 
        count 
        
        power
        angle
        alongship
        athwartship
        
        complexsamples
        compressed
        range
        minCount
        maxCount
        minRange
        maxRange
        dR
    end
    
    methods
    
        function x = asStruct(obj)
            x = struct(...
                'timestamp', obj.timestamp, ...
                'channelID', obj.channelID, ...
                'datatype', obj.datatype, ...
                'spare', obj.spare, ...
                'offset', obj.offset, ...
                'count', obj.count, ...
                'complexsamples', obj.complexsamples, ...
                'minCount', obj.minCount, ...
                'maxCount', obj.maxCount, ...
                'minRange', obj.minRange, ...
                'maxRange', obj.maxRange, ...
                'dR', obj.dR ...
                );
        end
        
        function obj = read(obj,fid)
            obj.channelID               = deblank(fread(fid,128,'*char')');
            obj.datatype                = fread(fid,1,'int16');
            obj.spare                   = fread(fid,2,'*char');
            obj.offset                  = fread(fid,1,'int32');
            obj.count                   = fread(fid,1,'int32');
            if (bitget(obj.datatype, 1, 'int16')) % power values
                obj.power               = fread(fid,obj.count,'int16');
                obj.power               = obj.power*10*log10(2)/256;
                if (bitget(obj.datatype, 2, 'int16')) % angle values
                    obj.angle           = fread(fid,[2 obj.count],'int8');
                    obj.angle           = obj.angle(1,:) + obj.angle(2,:)*256;
                    obj.alongship       = obj.angle(2,:)';
                    obj.athwartship     = obj.angle(1,:)';
                end
            elseif (bitget(obj.datatype, 4, 'int16')) % complex float 32 bit
                ncomplex     = bitshift(obj.datatype, -8, 'int16');
                samples      = fread(fid,2*ncomplex*obj.count,'float32=>single');
                samples      = reshape(samples,[2 ncomplex obj.count]);
                obj.complexsamples  = (squeeze(complex(samples(1,:,:),samples(2,:,:)))).';
            else
                error('Unknown sample mode');
            end
        end
        
    end
    
end

