function cmap = CMap_ek80
% --------------------------------------------
%   cmap = CMap_ek80
%
%   close to original colormap on ek80 software
%      taken from screendump using gimp
% 
% from gimp and ek80 dashboard
% ----------------------------------------------

 

% % % % another attempt --------------------------------
% Pos = [1    2   3   4    5   6   7   8   9   10  11  12   13   17   23   25   27  28  30   31  33  35    37  39    41  43   45    47  49   51 53   55   57  61   63   65   66];
% R = [ 156  141 126 112  97  82  68  53   39  24   9  02    9   36   37   24   10  21  114  161 255 253  252  252  252  252  252  253  254 255 226 198  176  153  151 153  155]';
% G = [ 138  125 113 100  88  76  76  83   90  96  103 102  84   12  200  171  141 139 185   208 255 204  153  116  110   99   85   61   36  37  37  51   57   44   31  23    7]';
% B = [ 168  150 132 114  96  78  94  129 163  197 232 249  234  174 122  111   99  92  72    62  42  44   46   63   85  130  160  118   75  51  51  49   49   58   45   33   11]';
% R = R/255;
% G = G/255;
% B = B/255;
 
% N = 64;
% R2 = linterp(Pos, R, 1: (N-1)/N : Pos(end));
% G2 = linterp(Pos, G, 1 :(N-1)/N : Pos(end));
% B2 = linterp(Pos, B, 1: (N-1)/N : Pos(end));
% cmap2 = [R2' G2' B2'];
% cmap2 = [ 1 1 1; cmap2];
% cmap = cmap2;
% % % 
% ----------------------------------------------------------
 
% taken from screendump using gimp
cmap = [     
156	138	168; ...
141	125	150; ...
126	113	132; ...
112	100	114; ...
97	88	96; ...
82	76	78; ...
68	76	94; ...
53	83	129; ...
39	90	163; ...
24	96	197; ...
9	103	232; ...
2	102	249; ...
9	84	234; ...
15	66	219; ...
22	48	204; ...
29	30	189; ...
36	12	174; ...
37	49	165; ...
38	86	156; ...
39	123	147; ...
40	160	138; ...
41	197	129; ...
37	200	122; ...
30	185	116; ...
24	171	111; ...
17	156	105; ...
10	141	99; ...
21	139	92; ...
68	162	82; ...
114	185	72; ...
161	208	62; ...
208	231	52; ...
255	255	42; ...
254	229	43; ...
253	204	44; ...
253	179	45; ...
252	153	46; ...
252	128	47; ...
252	116	63; ...
252	110	85; ...
252	105	108; ...
252	99	130; ...
252	93	153; ...
252	85	160; ...
252	73	139; ...
253	61	118; ...
253	48	96; ...
254	36	75; ...
255	24	54; ...
240	30	52; ...
226	37	51; ...
212	44	50; ...
198	51	49; ...
184	57	48; ...
176	57	49; ...
170	54	51; ...
165	51	54; ...
159	47	56; ...
153	44	58; ...
150	39	56; ...
151	31	45; ...
153	23	33; ...
154	15	22; ...
155	7	11 ];

cmap = cmap/255;
cmap = [ 1 1 1; cmap];


if 0
    figure; clf 
    
   % make up some data
   Z = linspace(0,1,100)'*ones(1,100);
   surf(Z,'edgealpha',0);
   
   colormap(cmap);
  
end
