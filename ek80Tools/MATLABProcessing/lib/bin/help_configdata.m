%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% help configdata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Much of the information in here is not specific to the settings but more
% general configuration data, e.g. the max and min frequency values
% allowed, max power allowed, etc. 
%
% Beam parameters and specific transducer
% serial numbers are in here.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configdata = 
%          header: [1x1 struct]
%    transceivers: [1x3 struct] 
%    (3 transducers and 3 WBTs in this example)
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configdata.header = 
%    ApplicationName: 'EK80'
%         GPSOffsetX: '0' NOT USED
%         GPSOffsetY: '0' NOT USED
%         GPSOffsetZ: '0' NOT USED
%          MRUAlphaX: '0' NOT USED
%          MRUAlphaY: '0' NOT USED
%          MRUAlphaZ: '0' NOT USED
%         MRUOffsetX: '0' NOT USED
%         MRUOffsetY: '0' NOT USED
%         MRUOffsetZ: '0' NOT USED
%       Multiplexing: '0' NOT USED - MULTIPLE TRANSDUCERS TO ONE TRANCEIVER
%         SurveyName: 'TEST_WC381_200kHz_TARGET_DETECTION'
%           TimeBias: '0' LOCAL TIME OFFSET
%            Version: '1.0.5030.13621'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configdata.transceivers = 
    % 1x3 struct array with fields:
    % id
    % EthernetAddress
    % IPAddress
    % MarketSegment
    % SerialNumber
    % TransceiverNumber
    % TransceiverSoftwareVersion
    % TransceiverType
    % Version
    % channels
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configdata.transceivers.id AUMATICALLY GENERATED
    % WBT__nbsp--00907207adc4
    % WBT__nbsp--00907207adc7
    % WBT__nbsp--00907207adca
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDER SET IN SOFTWARE INSTALLATION OF TRANSDUCERS
% (in this example 120 kHz, 200 kHz, 70 kHz
%
% configdata.transceivers.channels
%
%                     Name: 'WBT__nbsp--1-1__nbsp--192.168.1.10__nbsp--ES120-7C'
%                BandWidth: '33300'
%            ChannelIdLong: 'WBT 1-1 192.168.1.10 ES120-7C'
%    MaxTxPowerTransceiver: '2000'
%              PulseLength: '6.4E-05;0.000128;0.000256;0.000512;0.001024'
%           SampleInterval: '1.6E-05;3.2E-05;6.4E-05;0.000128;0.000256'
%               transducer: [1x1 struct]
%
%                     Name: 'WBT__nbsp--2-1__nbsp--192.168.1.11__nbsp--ES200-7CD'
%                BandWidth: '33300'
%            ChannelIdLong: 'WBT 2-1 192.168.1.11 ES200-7CD'
%    MaxTxPowerTransceiver: '2000'
%              PulseLength: '6.4E-05;0.000128;0.000256;0.000512;0.001024'
%           SampleInterval: '1.6E-05;3.2E-05;6.4E-05;0.000128;0.000256'
%               transducer: [1x1 struct]
%
%                    Name: 'WBT__nbsp--3-1__nbsp--192.168.1.12__nbsp--ES70-7C'
%                BandWidth: '33300'
%            ChannelIdLong: 'WBT 371790 B 3913-1 192.168.1.10'
%    MaxTxPowerTransceiver: '2000'
%              PulseLength: '0.000128;0.000256;0.000512;0.001024;0.002048'
%           SampleInterval: '3.2E-05;6.4E-05;0.000128;0.000256;0.000512'
%               transducer: [1x1 struct]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configdata.transceivers(1).channels.transducer
%           AngleOffsetAlongship: '0'
%          AngleOffsetAthwartship: '0'
%       AngleSensitivityAlongship: '23'
%     AngleSensitivityAthwartship: '23'
%                        BeamType: '1'
%              BeamWidthAlongship: '7'
%            BeamWidthAthwartship: '7'
%    DirectivityDropAt2XBeamWidth: '0'
%             EquivalentBeamAngle: '-21'
%                       Frequency: '120000'
%                FrequencyMaximum: '160000'
%                FrequencyMinimum: '95000'
%                            Gain: '25.5;26.8;27;27;27'
%            MaxTxPowerTransducer: '250'
%                    SaCorrection: '0;0;0;0;0'
%                  TransducerName: 'ES120-7C'