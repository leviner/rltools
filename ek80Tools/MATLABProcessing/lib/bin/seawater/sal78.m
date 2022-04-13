function [SAL78]=sal78(CND,T,P);
%  SAL78  salinity from conductivity, S/ma, temp C and p(dbar)
%
%    [sal78]=sal78(CND,T,P);
%
%     conductivity ratio is calculated right away, then used.
%
%     THE CONDUCTIVITY RATIO (CND) = 1.0000000 FOR SALINITY = 35 PSS-78
%     TEMPERATURE = 15.0 DEG. CELSIUS , AND ATMOSPHERIC PRESSURE.
%             
% FUNCTION TO CONVERT CONDUCTIVITY RATIO TO SALINITY (M = 0)
%
%   REFERENCES:   ALSO LOCATED IN UNESCO REPORT # 37 1981
%  PRACTICAL SALINITY SCALE 1978: E.L. LEWIS IEEE OCEAN ENG. JAN. 1980
% 
% UNITS:      
%       PRESSURE        P        DECIBARS
%       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
%       CONDUCTIVITY    CND      RATIO     (M=0)
%       SALINITY        SAL78    (PSS-78)  (M=0)
%  CHECKVALUES:
%      SAL78=40.00000 :CND=1.888091*4.2914,T=40 DEG C,P=10000 DECIBARS:   
 
% SAL78: RETURNS ZERO FOR SALINITY:  < 0.02    
%
%  PRACTICAL SALINITY SCALE 1978 DEFINITION WITH TEMPERATURE CORRECTION
%  XT=T-15.0 : XR=SQRT(RT)
%     IMPLICIT REAL*8 (A-H,O-Z)
%     SAL(XR,XT) =((((2.7081*XR-7.0261)*XR+14.0941)*XR+25.3851)*XR
%    X-0.1692)*  XR+0.0080
%    X  +(XT/(1.0+0.0162*XT))*(((((-0.0144*XR+
%    X   0.0636)*XR-0.0375)*XR-0.0066)*XR-0.0056)*XR+0.0005)

%  DSAL(XR,XT)  FUNCTION FOR DERIVATIVE OF SAL(XR,XT) WITH XR.
%     DSAL(XR,XT) =((((13.5405*XR-28.1044)*XR+42.2823)*XR+50.7702)*XR  ...
%       -0.1692)+(XT/(1.0+0.0162*XT))*((((-0.0720*XR+0.2544)*XR ...
%       -0.1125)*XR-0.0132)*XR-0.0056) ;
%  FUNCTION RT35 :  C(35,T,0)/C(35,15,0) VARIATION WITH TEMPERATURE
%  WITH TEMPERATURE.
%     RT35(XT) = (((1.0031E-9*XT-6.9698E-7)*XT+1.104259E-4)*XT ...
%       + 2.00564E-2)*XT + 0.6766097  ;
% POLYNOMIALS OF RP: C(S,T,P)/C(S,T,0) VARIATION WITH PRESSURE 
%  C(XP) POLYNOMIAL CORRESPONDS TO A1-A3 CONSTANTS: LEWIS 1980
%     C(XP) = ((3.989E-15*XP-6.370E-10)*XP+2.070E-5)*XP ;
%     B(XT) = (4.464E-4*XT+3.426E-2)*XT + 1.0  ;
%  A(XT) POLYNOMIAL CORRESPONDS TO B3 AND B4 CONSTANTS: LEWIS 1980
%     A(XT) = -3.107E-3*XT + 0.4215  ;

%*******************************************************************
      XT = T - 15.0 ;
% ************************************************
      R = CND./4.2914;
   %   RT = R/(RT35(T)*(1.0 + C(P)/(B(T) + A(T)*R))) ;  repalce this!
      C = ((3.989E-15.*P-6.370E-10).*P+2.070E-5).*P ;
      B = (4.464E-4.*T+3.426E-2).*T + 1.0  ;
      A = -3.107E-3.*T + 0.4215  ;

      RT35=   (((1.0031E-9.*T-6.9698E-7).*T+1.104259E-4).*T ...
        + 2.00564E-2).*T + 0.6766097  ;
      RT = R./(RT35.*(1.0 + C./(B + A.*R))) ;

      XR = sqrt(abs(RT)) ;

      SAL78 =((((2.7081 .*XR-7.0261).*XR+14.0941).*XR+25.3851).*XR ...
       -0.1692) .*  XR+0.0080  ...
       +(XT./(1.0+0.0162 .*XT)).*(((((-0.0144.*XR+  ...
        0.0636).*XR-0.0375).*XR-0.0066).*XR-0.0056).*XR+0.0005);
