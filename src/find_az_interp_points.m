function[ind]=find_az_interp_points(az,tol)
if ~exist('tol')
  tol=1.0;
end

keep_ind(1:length(az))=false;
keep_ind(1)=true;
keep_ind(end)=true;
i1=1;
i2=2;
n=length(az);
while (1) 
  aa=interp1([i1 i2],az([i1 i2]),i1:i2,'linear');
  maxerr=max(abs(aa'-az(i1:i2)));
  if (maxerr>tol)
    i1=i2-1;
    keep_ind(i1)=true;
  else
    i2=i2+1;
    if (i2==n)
      keep_ind(i2)=true;
      break;    
    end
  end
end

ii=1:length(az);
ind=ii(keep_ind);
