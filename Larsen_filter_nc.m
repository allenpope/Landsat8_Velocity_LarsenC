%%
%Filtering based on stack
cd ~/Desktop/pycorr/Out/LarsenC/Crop/

time = ncread('LarsenC_vx.nc','time');
x = ncread('LarsenC_vx.nc','x');
y = ncread('LarsenC_vx.nc','y');
time = ncread('LarsenC_vx.nc','time');
duration = ncread('LarsenC_vx.nc','duration');
year_start = ncread('LarsenC_vx.nc','year_start');
doy_start = ncread('LarsenC_vx.nc','doy_start');
year_end = ncread('LarsenC_vx.nc','year_end');
doy_end = ncread('LarsenC_vx.nc','doy_end');
path = ncread('LarsenC_vx.nc','path');
row = ncread('LarsenC_vx.nc','row');

today = date;

%%
%create new vx netcdf files to keep reproducible
vx_file_out = strcat('LarsenC_vx_',today,'.nc');

nccreate(vx_file_out,'time','Dimensions',{'time' Inf},'Datatype','double'); 
ncwriteatt(vx_file_out,'time','long_name','decimal date at middle of velocity image pair');
ncwriteatt(vx_file_out,'time','units','decimal years');
ncwrite(vx_file_out,'time',time);

nccreate(vx_file_out,'duration','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'duration','long_name','Duration between velocity image pairs')
ncwriteatt(vx_file_out,'duration','units','days');
ncwrite(vx_file_out,'duration',duration);

nccreate(vx_file_out,'year_start','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'year_start','long_name','Year of first velocity image')
ncwriteatt(vx_file_out,'year_start','units','year');
ncwrite(vx_file_out,'year_start',year_start);

nccreate(vx_file_out,'doy_start','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'doy_start','long_name','Day of year of first velocity image')
ncwriteatt(vx_file_out,'doy_start','units','day of year');
ncwrite(vx_file_out,'doy_start',doy_start);

nccreate(vx_file_out,'year_end','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'year_end','long_name','Year of second velocity image')
ncwriteatt(vx_file_out,'year_end','units','year');
ncwrite(vx_file_out,'year_end',year_end);

nccreate(vx_file_out,'doy_end','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'doy_end','long_name','Day of year of second velocity image')
ncwriteatt(vx_file_out,'doy_end','units','day of year');
ncwrite(vx_file_out,'doy_end',doy_end);

nccreate(vx_file_out,'path','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'path','long_name','WRS2 path of velocity image pairs')
ncwriteatt(vx_file_out,'path','units','path');
ncwrite(vx_file_out,'path',path);

nccreate(vx_file_out,'row','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vx_file_out,'row','long_name','WRS2 row of velocity image pairs')
ncwriteatt(vx_file_out,'row','units','row');
ncwrite(vx_file_out,'row',row);

nccreate(vx_file_out, 'x','Dimensions',{'x' size(x,1)},'Datatype','double');
ncwriteatt(vx_file_out,'x','long_name','Easting');
ncwriteatt(vx_file_out,'x','units','meters easting, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vx_file_out,'x',x);

nccreate(vx_file_out, 'y','Dimensions',{'y' size(y,1)},'Datatype','double');
ncwriteatt(vx_file_out,'y','long_name','Northing');
ncwriteatt(vx_file_out,'y','units','meters northing, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vx_file_out,'y',y);

nccreate(vx_file_out,'vx','Dimensions',{'x' 'y' 'time'},'Datatype','double');
ncwriteatt(vx_file_out,'vx','long_name','velocity in the x direction')
ncwriteatt(vx_file_out,'vx','units','meters per day');
ncwriteatt(vx_file_out,'vx','code','PyCorr v1.11');

%%
vy_file_out = strcat('LarsenC_vy_',today,'.nc');

nccreate(vy_file_out,'time','Dimensions',{'time' Inf},'Datatype','double'); 
ncwriteatt(vy_file_out,'time','long_name','decimal date at middle of velocity image pair');
ncwriteatt(vy_file_out,'time','units','decimal years');
ncwrite(vy_file_out,'time',time);

nccreate(vy_file_out,'duration','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'duration','long_name','Duration between velocity image pairs')
ncwriteatt(vy_file_out,'duration','units','days');
ncwrite(vy_file_out,'duration',duration);

nccreate(vy_file_out,'year_start','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'year_start','long_name','Year of first velocity image')
ncwriteatt(vy_file_out,'year_start','units','year');
ncwrite(vy_file_out,'year_start',year_start);

nccreate(vy_file_out,'doy_start','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'doy_start','long_name','Day of year of first velocity image')
ncwriteatt(vy_file_out,'doy_start','units','day of year');
ncwrite(vy_file_out,'doy_start',doy_start);

nccreate(vy_file_out,'year_end','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'year_end','long_name','Year of second velocity image')
ncwriteatt(vy_file_out,'year_end','units','year');
ncwrite(vy_file_out,'year_end',year_end);

nccreate(vy_file_out,'doy_end','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'doy_end','long_name','Day of year of second velocity image')
ncwriteatt(vy_file_out,'doy_end','units','day of year');
ncwrite(vy_file_out,'doy_end',doy_end);

nccreate(vy_file_out,'path','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'path','long_name','WRS2 path of velocity image pairs')
ncwriteatt(vy_file_out,'path','units','path');
ncwrite(vy_file_out,'path',path);

nccreate(vy_file_out,'row','Dimensions',{'time'},'Datatype','double');
ncwriteatt(vy_file_out,'row','long_name','WRS2 row of velocity image pairs')
ncwriteatt(vy_file_out,'row','units','row');
ncwrite(vy_file_out,'row',row);

nccreate(vy_file_out, 'x','Dimensions',{'x' size(x,1)},'Datatype','double');
ncwriteatt(vy_file_out,'x','long_name','Easting');
ncwriteatt(vy_file_out,'x','units','meters easting, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vy_file_out,'x',x);

nccreate(vy_file_out, 'y','Dimensions',{'y' size(y,1)},'Datatype','double');
ncwriteatt(vy_file_out,'y','long_name','Northing');
ncwriteatt(vy_file_out,'y','units','meters northing, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vy_file_out,'y',y);

nccreate(vy_file_out,'vy','Dimensions',{'x' 'y' 'time'},'Datatype','double');
ncwriteatt(vy_file_out,'vy','long_name','velocity in the y direction')
ncwriteatt(vy_file_out,'vy','units','meters per day');
ncwriteatt(vy_file_out,'vy','code','PyCorr v1.11');

%%
%
info_file_out = 'LarsenC_stack_info.nc';

nccreate(info_file_out, 'x','Dimensions',{'x' size(x,1)},'Datatype','double');
ncwriteatt(info_file_out,'x','long_name','Easting');
ncwriteatt(info_file_out,'x','units','meters easting, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(info_file_out,'x',x);

nccreate(info_file_out, 'y','Dimensions',{'y' size(y,1)},'Datatype','double');
ncwriteatt(info_file_out,'y','long_name','Northing');
ncwriteatt(info_file_out,'y','units','meters northing, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(info_file_out,'y',y);

nccreate(info_file_out,'count','Dimensions',{'x' 'y'},'Datatype','double');
ncwriteatt(info_file_out,'count','long_name','Number of data points in stack for this region')
ncwriteatt(info_file_out,'count','units','datapoints');

nccreate(info_file_out,'mean','Dimensions',{'x' 'y'},'Datatype','double');
ncwriteatt(info_file_out,'mean','long_name','Average velocity of data points for this region')
ncwriteatt(info_file_out,'mean','units','meters per day');

nccreate(info_file_out,'std','Dimensions',{'x' 'y'},'Datatype','double');
ncwriteatt(info_file_out,'std','long_name','Standard deviation of velocities of data points for this region')
ncwriteatt(info_file_out,'std','units','meters per day');


%%
%Filtering outside 3 std devs, over 5 m/day, too small sample size
%3x3 km box, aka 5 pixels
%could change to step of 2*pix_1 so no overlap of boxes
pix = 5;

parr = floor((size(y,1)-pix-pix)/11); %so that the parfor loop can be consecutive integers but still step 11 in real life

for par_r = 1:parr
	r = par_r*11+pix; %this should actually be r = 1+pix+(2*pix+1)*(par_r-1)
    a = strcat('r is _',num2str(par_r),'_of_',num2str(parr)); %num2str(size(y,1))
 	disp(a);
     
    for c = (pix+1):(2*pix+1):(size(x,1)-pix)
        a = strcat('c is _',num2str(c),'_of_',num2str(size(x,1)));
        disp(a);
    
        temp_x = ncread('LarsenC_vx.nc','vx',[c-pix r-pix 1],[11 11 size(time,1)]); %dims x, y, time
        temp_x(find(temp_x==0)) = nan;
        
        temp_y = ncread('LarsenC_vy.nc','vy',[c-pix r-pix 1],[11 11 size(time,1)]); %dims x, y, time
        temp_y(find(temp_y==0)) = nan;
        
        temp = (temp_x.^2 + temp_y.^2).^0.5;
        
        %going WAY too fast
        temp(find(temp > 4)) = nan;
        
        %first pass
        m = nanmean(temp(find(isfinite(temp))));
        s = nanstd(temp(find(isfinite(temp))));
        
        %big enough sample size
        if length(find(isfinite(temp))) < 150
            temp(find(isfinite(temp))) = nan;
        end
        
        %outside of 3 standard deviations
        temp(find(temp > m+3*s)) = nan;
        temp(find(temp < m-3*s)) = nan;

        
        %iterate a second time
        m = nanmean(temp(find(isfinite(temp))));
        s = nanstd(temp(find(isfinite(temp))));
        
        %big enough sample size
        if length(find(isfinite(temp))) < 150
            temp(find(isfinite(temp))) = nan;
        end
        
        %outside of 3 standard deviations
        temp(find(temp > m+3*s)) = nan;
        temp(find(temp < m-3*s)) = nan;

        
        %iterate a third time
        m = nanmean(temp(find(isfinite(temp))));
        s = nanstd(temp(find(isfinite(temp))));
        
        %big enough sample size
        if length(find(isfinite(temp))) < 150
            temp(find(isfinite(temp))) = nan;
        end
        
        %outside of 3 standard deviations
        temp(find(temp > m+3*s)) = nan;
        temp(find(temp < m-3*s)) = nan;
        
        
        %rewrite temps
        temp_x2 = zeros(size(temp_x));
        temp_x2(find(isfinite(temp))) = temp_x(find(isfinite(temp)));
        temp_y2 = zeros(size(temp_y));
        temp_y2(find(isfinite(temp))) = temp_y(find(isfinite(temp)));
        
        %creat count, std, and mean for writing
        m_block = m*ones(2*pix+1, 2*pix+1);
        s_block = s*ones(2*pix+1, 2*pix+1);
        c_block = length(find(isfinite(temp)))*ones(2*pix+1, 2*pix+1);
        
        %write back out to new files
        ncwrite(vx_file_out,'vx',temp_x2,[c-pix r-pix 1]);
        ncwrite(vy_file_out,'vy',temp_y2,[c-pix r-pix 1]);
        
        ncwrite(info_file_out,'count',c_block,[c-pix r-pix]);
        ncwrite(info_file_out,'mean',m_block,[c-pix r-pix]);
        ncwrite(info_file_out,'std',s_block,[c-pix r-pix]);
    end
end
%clear temp m s r c a

%%
clear all
