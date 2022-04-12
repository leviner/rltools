classdef cSampleData
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        timestamp % [Matlab time number of ping]
        channel % [channel number]
        mode_low % [1 = Power, 2 = angle, 3 = Both, XXXX]
        mode_high 
        transducerdepth % [m]
        fstart % [Hz]
        transmitpower % [W]
        pulselength % [s]
        bandwidth % [Hz]
        sampleinterval % [s]
        soundvelocity % [m/s]
        absorptioncoefficient % [dB/m]
        heave % [m]
        roll % [degrees]
        pitch % [degrees]
        temperature % [deg C]
        heading % [degrees]
        transmitmode % 0 = Active, 1 = Passive, 2 = Test, -1 = Unknown
        spare1 
        sweep % [Hz/s]
        offset 
        count 
        slope % [%]
        
        power
        angle
        alongship
        athwartship
        
        complexsamples
        compressed
        range
    end
    
    methods
    
        function x = asStruct(obj)
            x = struct(...
                'timestamp', obj.timestamp, ...
                'channel', obj.channel, ...
                'transducerdepth', obj.transducerdepth, ...
                'fstart', obj.fstart, ...
                'transmitpower', obj.transmitpower, ...
                'pulselength', obj.pulselength, ...
                'bandwidth', obj.bandwidth, ...
                'sampleinterval', obj.sampleinterval, ...
                'soundvelocity', obj.soundvelocity, ...
                'absorptioncoefficient', obj.absorptioncoefficient, ...
                'heave', obj.heave, ...
                'roll', obj.roll, ...
                'pitch', obj.pitch, ...
                'temperature', obj.temperature, ...
                'heading', obj.heading, ...
                'transmitmode', obj.transmitmode, ...
                'sweep', obj.sweep, ...
                'offset', obj.offset, ...
                'count', obj.count, ...
                'slope', obj.slope, ...
                'complexsamples', obj.complexsamples ...
                );
        end
        
        function obj = read(obj,fid)
            obj.channel                 = fread(fid,1,'int16');
            obj.mode_low                = fread(fid,1,'int8');
            obj.mode_high               = fread(fid,1,'int8');
            obj.transducerdepth         = fread(fid,1,'float32');
            obj.fstart                  = fread(fid,1,'float32');
            obj.transmitpower           = fread(fid,1,'float32');
            obj.pulselength             = fread(fid,1,'float32');
            obj.bandwidth               = fread(fid,1,'float32');
            obj.sampleinterval          = fread(fid,1,'float32');
            obj.soundvelocity           = fread(fid,1,'float32');
            obj.absorptioncoefficient   = fread(fid,1,'float32');
            obj.heave                   = fread(fid,1,'float32');
            obj.roll                    = fread(fid,1,'float32');
            obj.pitch                   = fread(fid,1,'float32');
            obj.temperature             = fread(fid,1,'float32');
            obj.heading                 = fread(fid,1,'float32');
            obj.transmitmode            = fread(fid,1,'int16');
            obj.spare1                  = char(fread(fid,2,'char')');
            obj.sweep                   = fread(fid,1,'float32');
            obj.offset                  = fread(fid,1,'int32');
            obj.count                   = fread(fid,1,'int32');
            if (obj.mode_low<4)
                obj.power               = fread(fid,obj.count,'int16');
                obj.power               = obj.power*10*log10(2)/256;
                if (obj.mode_low==3)
                    obj.angle           = fread(fid,[2 obj.count],'int8');
                    obj.angle           = obj.angle(1,:) + obj.angle(2,:)*256;
                    obj.alongship       = obj.angle(2,:)';
                    obj.athwartship     = obj.angle(1,:)';
                end
            elseif (obj.mode_low==8)
                ncomplex            = obj.mode_high;
                samples      = fread(fid,2*ncomplex*obj.count,'float32=>single');
                samples      = reshape(samples,[2 ncomplex obj.count]);
                obj.complexsamples  = (squeeze(complex(samples(1,:,:),samples(2,:,:)))).';
            else
                error('Unknown sample mode');
            end
                
            % slope is not currently in this telegram (it should be), so
            % hard-code it in at the most common setting...
            obj.slope = 0.1;
        end
        
    end
    
end

