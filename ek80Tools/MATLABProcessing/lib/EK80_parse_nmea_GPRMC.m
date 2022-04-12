function [Lat,Lon, YY,MM,DD,hh,mm,ss] = EK80_parse_nmea_GPRMC(nmea)
% ----------------------------------------------------------------
%     [Latd,Latm,Lond,Lonm, YY,MM,DD,hh,mm,ss] = EK80_parse_nmea_GPRMC(nmea)
%
%   function to parse GPRMC nmea string
%
% 
% 
% ----------------------------------------------------------------

% should not happen but just checking
if nmea(1:6) ~= '$GPRMV'
     error('Not correct nmea string');
end
    
latdir = 1; londir = 1;

itm = regexp(nmea,',','split');   % regular expression to split on comma

if (itm{5} =='S'), latdir = -1; end

Latd = str2num(itm{4}(1:2));
Latm = str2num(itm{4}(3:end));
Lat = (Latd + Latm/60) * latdir;


if (itm{7} =='W'), londir = -1; end
 

Lond = str2num(itm{6}(1:3));
Lonm = str2num(itm{6}(4:end));
Lon = (Lond + Lonm/60) * londir;

YY   = str2num(itm{10}(5:6))+ 2000;
MM   = str2num(itm{10}(3:4));
DD   = str2num(itm{10}(1:2));
hh   = str2num(itm{2}(1:2));
mm   = str2num(itm{2}(3:4));
ss   = str2num(itm{2}(5:end));






