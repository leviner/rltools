classdef cRawFile
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        fid
        
        headerlength = 12;
        dglength
        
        header      = cHeader;
        config      = cConfig;
        filter0data  = cFilter0Data;
        filter1data  = cFilter1Data;
        nmea        = cNMEA;
        annotation  = cAnnotation;
        sampledata  = cSampleData;
        sampledata3 = cSampleDataRAW3;
        xml         = cXML;
        mru         = cMRUData;
    end
    
    methods
        
        function obj = cRawFile(name)
            obj.name = name;
        end
        
        function obj = open(obj)
            obj.fid = fopen(obj.name,'r');
        end
        
        function obj = readdg(obj)
            % Read configuration datagram
            obj.dglength = fread(obj.fid,1,'int32');

            if (feof(obj.fid))
                return
            end
            obj.header = obj.header.read(obj.fid);
            
            %disp(['Read datagram of type ' obj.header.type])

            switch (obj.header.type)

                case 'CON0' % Configuration datagram
                    obj.config = obj.config.read(obj.fid);

                case 'FIL0' % Filter datagram
                    obj.filter0data = obj.filter0data.read(obj.fid,obj.dglength-obj.headerlength);

                case 'FIL1' % Filter datagram
                    obj.filter1data = obj.filter1data.read(obj.fid);
                    
                case 'NME0' % NMEA datagram
                    obj.nmea = obj.nmea.read(obj.fid,obj.dglength-obj.headerlength);
                    
                case 'TAG0' % Annotation datagram
                    obj.annotation = obj.annotation.read(obj.fid,obj.dglength-obj.headerlength);

                case 'RAW0' % Sample datagram
                    obj.sampledata = obj.sampledata.read(obj.fid);
                    obj.sampledata.timestamp = obj.header.datetime;
                    
                case 'XML0' % XML datagram
                    obj.xml = obj.xml.read(obj.fid,obj.dglength-obj.headerlength);
                    
                case 'MRU0' % MRU datagram
                    obj.mru = obj.mru.read(obj.fid);
                    
                case 'RAW3'
                    obj.sampledata = obj.sampledata3.read(obj.fid);
                    obj.sampledata.timestamp = obj.header.datetime;
                    
                otherwise
                    warning(strcat('Unknown datagram ''',obj.header.type,''' in file'));
                    obj.header.type = 'unknown';

            end
            
            obj.dglength = fread(obj.fid,1,'int32');

        end
        
    end
    
end
