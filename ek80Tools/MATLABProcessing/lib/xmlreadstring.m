%% 
function [parseResult,p] = xmlreadstring(stringToParse,varargin)
%XMLREADSTRING Modified XMLREAD function to read XML data from a string.
% Author: Luis Cantero.
% The MathWorks.

p = locGetParser(varargin);
locSetEntityResolver(p,varargin);
locSetErrorHandler(p,varargin);

% remove hex newlines... '&#x0A;','&#x0D;' and others   Newhall 2020
 
stringToParse = regexprep(stringToParse,'&#x..;',' ');  
stringToParse = regexprep(stringToParse,'ÖD1ÿéC',' '); % added by RK 7/2020
stringToParse = regexprep(stringToParse,'@{ ®Gáz @ôqÄ§',' '); % added by RK 7/2020 
stringToParse = regexprep(stringToParse,'Ã¿Ã¿â‚¬',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã–D1ÃŸÃ©C',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã—DÃ€\?Ã§C',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'m:Â J9Å¡  Â¸;Ã®â€º9Â¤Ã¶VÃˆÂ³&lt;',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Â§9Â±Ã‰â€ž9Ã¡Ã¦ 9',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Â¥9AÂ³29Å½Ãž 9Å’kÃ–8;6(5',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Â¥9AÂ³29Å½Ãž 9Å’kÃ–8  X4',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã–D1Ã¿Ã©C',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'â‚¬DÂ±u@B',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'DÃ’Â¯&amp;B Ã¡AI ÃÂ¼&gt;ÂºGÃ\\Â»',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã–DÃ®uÃªC',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã–DnÃ¬Ã©C',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Â¥9AÂ³29Å½Ãž 9Å’kÃ–8 aâ€”â€œ',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Â¥9AÂ³29Å½Ãž 9Å’kÃ–8naâ€°â€œb',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'dÃ‚`w',' '); % added by RK 6/2021
stringToParse = regexprep(stringToParse,'Ã¬Ã‚Ã  ',' '); % added by RK 6/2021


% Parse and return.
parseStringBuffer = java.io.StringBufferInputStream(stringToParse);
parseResult = p.parse(parseStringBuffer);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = locGetParser(args)

p = [];
%%
for i=1:length(args)
    if isa(args{i},'javax.xml.parsers.DocumentBuilderFactory')
        javaMethod('setValidating',args{i},locIsValidating(args));
        p = javaMethod('newDocumentBuilder',args{i});
        break;
    elseif isa(args{i},'javax.xml.parsers.DocumentBuilder')
        p = args{i};
        break;
    end
end

if isempty(p)
    parserFactory = javaMethod('newInstance',...
        'javax.xml.parsers.DocumentBuilderFactory');
        
    javaMethod('setValidating',parserFactory,locIsValidating(args));
    %javaMethod('setIgnoringElementContentWhitespace',parserFactory,1);
    %ignorable whitespace requires a validating parser and a content model
    p = javaMethod('newDocumentBuilder',parserFactory);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf=locIsValidating(args)

tf=any(strcmp(args,'-validating'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locSetEntityResolver(p,args)

for i=1:length(args)
    if isa(args{i},'org.xml.sax.EntityResolver')
        p.setEntityResolver(args{i});
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locSetErrorHandler(p,args)

for i=1:length(args)
    if isa(args{i},'org.xml.sax.ErrorHandler')
        p.setErrorHandler(args{i});
        break;
    end
end
