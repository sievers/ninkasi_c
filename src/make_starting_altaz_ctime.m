crap=load('for_sievers.txt');
alt=crap(:,3);
az=crap(:,4);
ctime=crap(:,1);
len=crap(:,2);

plot(alt,az,'.');
is_strip=((az>100)&(az<250));
scan_len=410;  %everybody is going to be this long

fid_equitorial=fopen('altaz_ctime_equitorial.txt','w');
fid_strip=fopen('altaz_ctime_strip.txt','w');
for j=1:length(ctime),
    if is_strip(j),
        fid=fid_strip;
    else
        fid=fid_equitorial;
    end
    n=ceil(len(j)/scan_len);
    for k=0:n-1,
        fprintf(fid,'%12.4f %12.4f %14.2f\n',alt(j),az(j),ctime(j)+k*scan_len);
    end
end
fclose(fid_equitorial);
fclose(fid_strip);


        