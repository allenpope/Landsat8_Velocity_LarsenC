%%
%Building a mosiac (based on times?)
%Directory with files
%edit commented/uncommented below to change which mosaic is made (could do this as function/parameter...)
cd ~/Desktop/pycorr/Out/LarsenC/Crop/

vx_file_in = 'LarsenC_vx_03-Nov-2016.nc';
vy_file_in = 'LarsenC_vy_03-Nov-2016.nc';

%%
%
time = ncread(vx_file_in,'time');
x = ncread(vx_file_in,'x');
y = ncread(vx_file_in,'y');

duration = ncread(vx_file_in,'duration');

%2013-2014 mosaic; 6
%index_time = find(time<2014.5);

%2014-2015 mosaic; 78
%index_time = find(time>2014.5 & time<2015.5); 

%2015-2016 mosaic; 209
%index_time = find(time>2015.5);% & time<2016.5); 

%entire mosaic; 293
index_time = find(time);


LarsenC_vx = zeros(length(y),length(x));
LarsenC_vy = zeros(length(y),length(x));
weight = zeros(length(y),length(x));

for i = 1:length(index_time)
    disp(strcat('processing_',num2str(i),'_of_',num2str(length(index_time))));
    
    vx_filtered = ncread(vx_file_in,'vx',[1 1 index_time(i)],[length(x) length(y) 1])'; %dims x, y, time
    vy_filtered = ncread(vy_file_in,'vy',[1 1 index_time(i)],[length(x) length(y) 1])'; %dims x, y, time

        if duration(index_time(i)) == 16
           LarsenC_vx = LarsenC_vx + vx_filtered*0.3;
           LarsenC_vy = LarsenC_vy + vy_filtered*0.3;
           weight(find(vx_filtered)) = weight(find(vx_filtered)) + 0.3;
        end
        if duration(index_time(i)) == 32
           LarsenC_vx = LarsenC_vx + vx_filtered*0.6;
           LarsenC_vy = LarsenC_vy + vy_filtered*0.6;
           weight(find(vx_filtered)) = weight(find(vx_filtered)) + 0.6;
        end
        if duration(index_time(i)) == 48
           LarsenC_vx = LarsenC_vx + vx_filtered*0.9;
           LarsenC_vy = LarsenC_vy + vy_filtered*0.9;
           weight(find(vx_filtered)) = weight(find(vx_filtered)) + 0.9;
        end
        if duration(index_time(i)) > 48
           LarsenC_vx = LarsenC_vx + vx_filtered;
           LarsenC_vy = LarsenC_vy + vy_filtered;
           weight(find(vx_filtered)) = weight(find(vx_filtered)) + 1;
        end 
end

LarsenC_vx = LarsenC_vx ./ weight;
LarsenC_vy = LarsenC_vy ./ weight;
LarsenC_vv = (LarsenC_vx.*LarsenC_vx+LarsenC_vy.*LarsenC_vy).^0.5;

clear i vy_filtered vx_filtered
clear time x y duration ans
clear weight index_time

%%
%Remove non-data points

LarsenC_vv(find(LarsenC_vv > 5)) = 0;
LarsenC_vv(find(isnan(LarsenC_vv))) = 0;

LarsenC_vx(find(LarsenC_vx > 5)) = 0;
LarsenC_vx(find(isnan(LarsenC_vx))) = 0;

LarsenC_vy(find(LarsenC_vy > 5)) = 0;
LarsenC_vy(find(isnan(LarsenC_vy))) = 0;

%%
%Export Geotiffs

vx_files = dir('*vx.tif');
[vx, vx_info] = geotiffread(vx_files(1).name);
tiffinfo= geotiffinfo(vx_files(1).name);

%geotiffwrite('LarsenC_vx_1314.tif', LarsenC_vx, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vy_1314.tif', LarsenC_vy, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vv_1314.tif', LarsenC_vv, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);

%geotiffwrite('LarsenC_vx_1415.tif', LarsenC_vx, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vy_1415.tif', LarsenC_vy, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vv_1415.tif', LarsenC_vv, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);

%geotiffwrite('LarsenC_vx_1516.tif', LarsenC_vx, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vy_1516.tif', LarsenC_vy, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%geotiffwrite('LarsenC_vv_1516.tif', LarsenC_vv, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);

geotiffwrite('LarsenC_vx.tif', LarsenC_vx, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LarsenC_vy.tif', LarsenC_vy, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('LarsenC_vv.tif', LarsenC_vv, vx_info, 'GeoKeyDirectoryTag',tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);

%%
%

clear