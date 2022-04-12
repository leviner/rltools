classdef cConfigTransducer
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        channelid
        beamtype
        frequency
        gain
        equivalentbeamangle
        beamwidthalongship
        beamwidthathwartship
        anglesensitivityalongship
        anglesensitivityathwartship
        angleoffsetalongship
        angleoffsetathwartship
        posx
        posy
        posz
        dirx
        diry
        dirz
        pulselengthtable
        spare2
        gaintable
        spare3
        sacorrectiontable
        slope
        directivitydrop
        transceiverswver
        spare4
    end
    
    methods
        function obj = read(obj,fid)
            obj.channelid                   = char(fread(fid,128,'char')');
            obj.beamtype                    = fread(fid,1,'int32');
            obj.frequency                   = fread(fid,1,'float32');
            obj.gain                        = fread(fid,1,'float32');
            obj.equivalentbeamangle         = fread(fid,1,'float32');
            obj.beamwidthalongship          = fread(fid,1,'float32');
            obj.beamwidthathwartship        = fread(fid,1,'float32');
            obj.anglesensitivityalongship   = fread(fid,1,'float32');
            obj.anglesensitivityathwartship = fread(fid,1,'float32');
            obj.angleoffsetalongship        = fread(fid,1,'float32');
            obj.angleoffsetathwartship      = fread(fid,1,'float32');
            obj.posx                        = fread(fid,1,'float32');
            obj.posy                        = fread(fid,1,'float32');
            obj.posz                        = fread(fid,1,'float32');
            obj.dirx                        = fread(fid,1,'float32');
            obj.diry                        = fread(fid,1,'float32');
            obj.dirz                        = fread(fid,1,'float32');
            obj.pulselengthtable            = fread(fid,5,'float32');
            obj.spare2                      = char(fread(fid,8,'char')');
            obj.gaintable                   = fread(fid,5,'float32');
            obj.spare3                      = char(fread(fid,8,'char')');
            obj.sacorrectiontable           = fread(fid,5,'float32');
            obj.slope                       = fread(fid,1,'float32');
            obj.directivitydrop             = fread(fid,1,'float32');
            obj.transceiverswver            = char(fread(fid,16,'char')');
            obj.spare4                      = char(fread(fid,28,'char')');
        end
    end
    
end

