function[tod]=read_tod(dirfile)
% sample: /cita/d/raid-nolta/act/data/season1/1197424216.1197424241

rows=repmat(0:31,[32 1]);
cols=rows';
rowvec=reshape(rows,[numel(rows) 1]);
colvec=reshape(cols,[numel(rows) 1]);
if ~exist('dirfile')
    dirfile='/cita/d/raid-nolta/act/data/season1/1197424216.1197424241';
end

dirfile=[dirfile '/'];  %make sure we have a trailing slash.
ndet=length(rowvec);
for j=1:ndet,
    tag=sprintf('tesdatar%02dc%02d',rows(j),cols(j));
    fid=fopen([dirfile tag]);
    vec=fread(fid,inf,'int');
    fclose(fid);
        
    if (j==1),
        tod.data(length(vec),ndet)=0;
    end
    
    tod.data(:,j)=vec;
end
tod.rows=rowvec;
tod.cols=colvec;

 




