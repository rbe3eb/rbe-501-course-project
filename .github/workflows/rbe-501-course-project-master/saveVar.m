function [] = saveVar(filename,var)
    fid = fopen(filename, 'w');
    fprintf(fid, '[\n');
    arrayfun(@(ROWIDX) fprintf(fid,'%s,',var(ROWIDX,1:end-1)) + fprintf(fid,'%s;\n', var(ROWIDX,end)), (1:size(var,1)).');
    fprintf(fid,']\n');
    fclose(fid);
end