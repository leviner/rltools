function cal = parseCals(calfilename, calfilepath);
for j=1:length(calfilename)
    cal(j) = readstruct(fullfile(calfilepath,string(calfilename(j))));
end