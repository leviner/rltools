function data_all = EK80readRawFiles(rawfiles, version)

% Keeping version option for future revisions
%
% Set version to 'V2' for data recorded before the Oct 2014 survey
% Set version to 'V3' for data record during and after the Oct 2014 survey.
% 'V3' files have the new XML0/parameter datagram, and the RAW3 
% datagram.full_rawfiles=char(rawfiles);
%
%  revised Perkins
%
%  revised 2019 Newhall to work with empty data files

full_rawfiles=char(rawfiles);
disp(['Loading ' full_rawfiles])

if strcmp(version, 'V2')
    data = EK80readRaw(full_rawfiles, maxRange);   % in EK80 folder
    data.version = version;
elseif strcmp(version, 'V3')
    data = EK80readRawV3(full_rawfiles);
elseif strcmp(version, 'V4')
    data = EK80readRawV4(full_rawfiles);
    
    % fixing partial files  2019
    if isempty(data)
        data_all = [];
        return
    end
    
    data.version = version;
 
elseif strcmp(version, 'V4')
    data = EK80readRawV4(full_rawfiles);
    
    % fixing partial files  2019
    if isempty(data)
        data_all = [];
        return
    end
    
    data.version = version;
   
    
else
    disp(['Unknown file version of ' version])
    data.version = 'unknown';
end
fprintf(sprintf('Finished reading the file.\n'));

% Some hard-coded EK80 parameters
parameters.fs = 1.5e6; % sampling frequency [Hz]
parameters.Rwbtrx = 1000; % Wideband transceiver impedance [Ohms]
parameters.Ztrd = 75; % Transducer quadrant nominal impedance [Ohms]
parameters.rawfilename = full_rawfiles;
[~, name, ext] = fileparts(full_rawfiles);
parameters.filename = [name ext];
data.parameters = parameters;
data_all=data;