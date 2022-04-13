%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% help filterdatavec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There is information for two filterm stage 1 and stage 2. 
% These are complex coefficients.
%
% The first filter, stage 1, is actually applied digitally to the data, 
% so to accurately get the matched filter signal this filter must be applied. 
% Similarly, there is a first stage of decimation that is applied to the 
% data and shoudl be applied to the generation of the MF signal.
%
% The second filter, stage 2, is applied to tha data in the software
% display, but not to the saved "raw" data. So this does not need to be
% applied.But Lars advocates that it should!
%
% The dicimations are not the same for all the channels and the decimations
% are set by the signals chosen, that is, the bandwidth.
%
% The maximum bandwidth supported by the WBTs is 200 kHz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filterdatavec = 
% 3x2 struct array with fields:
%    stage
%    channelno
%    channelid
%    ncoeff
%    dec
%    coeff
% (There are three transducers in this example, each split beam connected
% to a separate WBT, so it is a 3x2 array, as there are two stages of
% filters.)
%
% Each field has 6 elements, and the order is 
% ch1 stage 1, ch2 stage 1, ch3 stage 1, ch1 stage 2, ch2 stage 2, ch3
% stage 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filterdatavec.ncoeff
%  These are the number of complex filter coefficients in each 
%  filter in the ordeer stated above.
%    47
%    59
%   109
%    25
%    59
%    43
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For example the coeffiecients for ch3 stage 2 are given by:
% filterdatavec(6).coeff = 
%   0.0001 - 0.0002i
%    0.0007 - 0.0005i
%   -0.0036 + 0.0004i
%    0.0042 + 0.0017i
%   -0.0001 - 0.0002i
%   -0.0011 - 0.0057i
%   -0.0013 + 0.0042i
%   -0.0043 + 0.0040i
%    0.0104 - 0.0027i
%    0.0014 + 0.0004i
%   -0.0123 - 0.0116i
%    0.0028 + 0.0087i
%   -0.0035 + 0.0184i
%    0.0164 - 0.0199i
%    0.0101 - 0.0040i
%   -0.0457 - 0.0058i
%    0.0112 + 0.0081i
%    0.0280 + 0.0595i
%    0.0045 - 0.0721i
%    0.0429 - 0.0676i
%   -0.2711 + 0.1491i
%    0.4204 + 0.0000i
%   -0.2711 - 0.1491i
%    0.0429 + 0.0676i
%    0.0045 + 0.0721i
%    0.0280 - 0.0595i
%    0.0112 - 0.0081i
%   -0.0457 + 0.0058i
%    0.0101 + 0.0040i
%    0.0164 + 0.0199i
%   -0.0035 - 0.0184i
%    0.0028 - 0.0087i
%   -0.0123 + 0.0116i
%    0.0014 - 0.0004i
%    0.0104 + 0.0027i
%   -0.0043 - 0.0040i
%   -0.0013 - 0.0042i
%   -0.0011 + 0.0057i
%   -0.0001 + 0.0002i
%    0.0042 - 0.0017i
%   -0.0036 - 0.0004i
%    0.0007 + 0.0005i
%    0.0001 + 0.0002i
