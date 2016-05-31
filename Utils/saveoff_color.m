function saveoff_color(fname,v,f,vc)
%
% saveoff(fname,v,vc,f)
%
% save a surface mesh or pointcloud to Object File Format (OFF)
%
% this is a faster implementation and color extention to a 'saveoff' code published by Qianqian Fang
% (fangq<at> nmr.mgh.harvard.edu) 
% written by Or Litany (orlitany <at> gmail <dot> com )
% 
% input:
% fname: output file name (including the 'off' suffix)%      
% v: input, surface node list, dimension (nn,3)
% vc: input, vertices color
% f: input, surface face element list, dimension (be,3)%


fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

if exist('vc')
    fprintf(fid,'COFF\n');
else
    fprintf(fid,'OFF\n');
end

fprintf(fid,'%d %d %d\n',length(v),length(f),0);

if exist('vc')
    fprintf(fid,'%f %f %f %d %d %d %d\n',[v';vc';255*ones(1,size(vc,1))]);
    

else
    fprintf(fid,'%f %f %f\n',v');
end

fprintf(fid,'%d %d %d %d \n',[ones(1,size(f,1))*size(f,2); f'-1]);
fclose(fid);

