function save_lines_povray(fname, s, A, B, diam, colors)

assert(size(A,2)==3)
assert(size(B,2)==3)
assert(size(A,1)==size(B,1))

fid = fopen(fname,'w');

fprintf(fid, 'union {\n');

for k=1:size(A,1)
    
    p = A(k,:);
    q = B(k,:);
    
    sz = num2str(diam);
    
    if p==q % degenerate cylinder
        continue
    end
    
    if nargin<6
        script = sprintf('cylinder\n{\n\t<%.4f,%.4f,%.4f>*%f, <%.4f,%.4f,%.4f>*%f, \t%s\n}',...
            p(1), p(2), p(3), s, q(1), q(2), q(3), s, sz);
    else
        script = sprintf('cylinder\n{\n\t<%.4f,%.4f,%.4f>*%f, <%.4f,%.4f,%.4f>*%f, \t%s\n\tpigment {color rgb<%.3f,%.3f,%.3f>}\n\tfinish { ambient 0.9 }\n}',...
            p(1), p(2), p(3), s, q(1), q(2), q(3), s, sz, colors(k,1), colors(k,2), colors(k,3));
    end
    
    fprintf(fid, '%s\n', script);
    
end

fprintf(fid, '}\n');

fclose(fid);

end
