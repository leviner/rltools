function data = EK80readRawV3(fname)

% Reads an entire .raw file and returns the data in various variables
% in a structure.
%
% data=[];
%
%  Revised Perkins
%
%   Revised to work with partial data files
%     Newhall 2019

pingno = 0;
envCount = 1;
paramCount = 1;
nmeaCount = 1;
mruCount = 1;
filterCount = 1;
headerlength = 12;

 % initialize hash tables
channels = containers.Map;
channelsInverse = containers.Map('KeyType', 'int32', 'ValueType', 'char');

channelNum = 1;
minRange=[];maxRange=[];
fid = fopen(fname,'r');

if (fid == -1)
    error(['Could not open file: ' fname]);
else
    while(1)
        dglength = fread(fid,1,'int32');
        
        if feof(fid)
            break
        end
        
        header = cHeader;                % make a new header class
        header = header.read(fid);       % read method
        
        switch(header.type)
            case 'XML0'
                xml = cXML;
                xml = xml.read(fid, dglength-headerlength);
                xmldata = xmlreadstring(deblank(xml.text));
                xmldata = parseXML(xmldata);
                %                     disp(xmldata.Name)
                switch (xmldata.Name)
                    case 'Configuration'
                        configData = parseconfxmlstruct(xmldata);
                    case 'Environment'
                        e = parseenvxmlstruct(xmldata);
                        e.timestamp = header.datetime;
                        envData(envCount) = e;
                        envCount = envCount + 1;
                    case 'Parameter'
                        p = parseparamxmlstruct(xmldata);
                        % Assumes that this datagram always comes
                        % before its associated RAW3 datagram
                        if ~channels.isKey(p.ChannelID)
                            channels(p.ChannelID) = channelNum;
                            
                               % debugging....
                            fprintf('p.ChannelID (hdr) %s\n',p.ChannelID);
                            
                            channelsInverse(channelNum) = p.ChannelID;
                            channelNum = channelNum + 1;
                        end
                        channel = channels(p.ChannelID);  % hash to get channel #
                        if channel == 1
                            pingno = pingno + 1;
                            no_reads = 1;
                        end
                        
                        p.timestamp = header.datetime;
                        
                        %%%%% ------------------------------------------
                        % Perkins edit to load new EK80 files from AR10 or
                        % later
                        if isfield(p,'Frequency')
                            p.FrequencyEnd=p.Frequency;
                            p.FrequencyStart=p.Frequency;
                            p=rmfield(p,'Frequency');
                            p=orderfields(p,[1;2;9;10;3;4;5;6;7;8]);
                        end
                        %%%%%% ----------------------------------------
                        
                        paramData(channel, pingno) = p; %#ok<*AGROW>
                        paramCount = paramCount + 1;
                    otherwise
                        disp(['Unknown XML datagram with toplevel element name of ' ...
                            xmldata.Name])
                end
            case 'FIL1'
                filter = cFilter1Data;
                filterData(filterCount) = filter.read(fid);
                filterCount = filterCount + 1;
            case 'MRU0'
                mru = cMRUData;
                mru = mru.read(fid);
                mru.timestamp = header.datetime;
                mruData(mruCount) = mru.asStruct();
                mruCount = mruCount + 1;
            case 'RAW3'
                s = cSampleDataRAW3_AR13;     % New class for AR13 v3
                s = s.read(fid);
                s.timestamp = header.datetime;
                
                s.dR=envData(end).SoundSpeed * paramData(channel, end).SampleInterval/2;
                s.minRange=0;
                s.maxRange=size(s.complexsamples, 1)*s.dR;
                s.minCount = 1;
                s.maxCount = round(s.maxRange./s.dR);
                s.maxCount=min(s.maxCount,size(s.complexsamples, 1));
                s.count = s.maxCount - s.minCount + 1;
                
                sampleData(channel, pingno) = s.asStruct();
                sampleData(channel, pingno).timestamp = header.datetime;
                ps(1:s.count,pingno)=sum(s.complexsamples,2);
            case 'NME0'
                nmea = cNMEA;           % new nmea class
                nmea = nmea.read(fid, dglength-headerlength);
                nmea.timestamp = header.datetime;
                nmeaData(nmeaCount) = nmea;
                nmeaCount = nmeaCount + 1;
            otherwise
                disp(['Unsupported datagram of type: ' header.type])
                fread(fid, dglength-headerlength);
        end
        fread(fid, 1, 'int32'); % the trailing datagram marker
    end
end

% Rearrange the filter structures, based on the channel number and the
% filter stage.
for i = 1:length(filterData)
    % debugging
    % fprintf('filterData(%d).ChannelID: %s\n',i,filterData(i).ChannelID);
    
     % fixing partial files 2019
    if ~isKey(channels,filterData(i).ChannelID)
        data = [];
        fprintf('\t Missing Key for channels.. Partial file?\n');
        return
    end
    
   
    
    channel = channels(filterData(i).ChannelID);
    f(channel, filterData(i).Stage) = filterData(i);
end
filterData = f;

if ~exist('nmeaData', 'var')
    nmeaData = struct([]);
end

data = struct('echodata', sampleData, 'filters', filterData, ...
    'config', configData, 'environ', envData, ...
    'nmea', nmeaData, 'param', paramData, 'channelMap', channels, ...
    'mru', mruData, 'channelIDs', (channelsInverse));

