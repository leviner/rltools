function [Lat, Lon, GPStime] =  EK80_get_nmea_time_position(nmea_data)
% ---------------------------------------------------------
%function [Lat, Lon, GPStime] =  EK80_get_nmea_time_position(nmea_data)
%
% make GPS vectors from NMEA strings
%
%   Lat, Lon, GPStime
%   if not complete, use NaNs
%   Order for CTRiverDT
%       GPGGA, GPGSA, GPRMC
%
%   SHOULD ONLY USE GPRMC,  has all the info...
% ------------------------------------------------------------
    
    GPStime=[]; Lat=[]; Lon=[];
    for k=1:length(nmea_data)
        
        if strcmp(nmea_data(k).text(1:6),'$GPRMC')
            [lat,lon,YYYY,MM,DD,GPShh,GPSmm,GPSss] = EK80_parse_nmea_GPRMC(nmea_data(k).text); 
            GPStime(end+1)=datenum([YYYY' MM' DD' GPShh' GPSmm' GPSss'])';
            Lat(end+1)=lat;
            Lon(end+1)=lon;
        end         
   
    end  % end of parsing nmea strings    
        
    