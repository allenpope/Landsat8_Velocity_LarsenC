%%
%Directory with files
cd ~/Desktop/pycorr/Out/LarsenC/Crop/

vx_file_out = 'LarsenC_vx.nc';
vy_file_out = 'LarsenC_vy.nc';

%%
%List of files
vx_files = dir('*VX.tif');
vy_files = dir('*VY.tif');
delcorr_files = dir('*delcorr.tif');

year_start = zeros(size(vx_files));
doy_start = zeros(size(vx_files));
year_end = zeros(size(vx_files));
doy_end = zeros(size(vx_files));
path = zeros(size(vx_files));
row = zeros(size(vx_files));
duration = zeros(size(vx_files));
time = zeros(size(vx_files));

for file = 1:size(vx_files);
    filename = vx_files(file).name;
    
    %Create year_start
    year_start(file) = str2num(filename(27:30));
    %create doy_start
    doy_start(file) = str2num(filename(32:34));
    %Create year_end
    year_end(file) = str2num(filename(36:39));
    %create doy_end
    doy_end(file) = str2num(filename(41:43));
    %create path
    path(file) = str2num(filename(15:17));
    %create row
    row(file) = str2num(filename(19:21));
    %create duration
    duration(file) = str2num(filename(23:25));
end

%create time
time = (year_start+doy_start/365+year_end+doy_end/365)/2;

clear filename file

%%
%Create grid
[vx, vx_info] = geotiffread(vx_files(1).name);

x = [vx_info.XWorldLimits(1)+vx_info.CellExtentInWorldX/2:vx_info.CellExtentInWorldX:vx_info.XWorldLimits(2)-vx_info.CellExtentInWorldX/2];
y = [vx_info.YWorldLimits(2)-vx_info.CellExtentInWorldY/2:-vx_info.CellExtentInWorldY:vx_info.YWorldLimits(1)+vx_info.CellExtentInWorldY/2];

clear vx vx_info


%%
%Start writing vx nc file

nccreate(vx_file_out,'time','Dimensions',{'time' Inf},'Datatype','double'); %100 as placeholder...
ncwriteatt(vx_file_out,'time','long_name','decimal date at middle of velocity image pair');
ncwriteatt(vx_file_out,'time','units','decimal years');

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

nccreate(vx_file_out, 'x','Dimensions',{'x' size(x,2)},'Datatype','double');
ncwriteatt(vx_file_out,'x','long_name','Easting');
ncwriteatt(vx_file_out,'x','units','meters easting, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vx_file_out,'x',x);

nccreate(vx_file_out, 'y','Dimensions',{'y' size(y,2)},'Datatype','double');
ncwriteatt(vx_file_out,'y','long_name','Northing');
ncwriteatt(vx_file_out,'y','units','meters northing, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vx_file_out,'y',y);

nccreate(vx_file_out,'vx','Dimensions',{'x' 'y' 'time'},'Datatype','double');
ncwriteatt(vx_file_out,'vx','long_name','velocity in the x direction')
ncwriteatt(vx_file_out,'vx','units','meters per day');
ncwriteatt(vx_file_out,'vx','code','PyCorr v1.11');

%%
%Start writing vy nc file

nccreate(vy_file_out,'time','Dimensions',{'time' Inf},'Datatype','double'); %100 as placeholder...
ncwriteatt(vy_file_out,'time','long_name','decimal date at middle of velocity image pair');
ncwriteatt(vy_file_out,'time','units','decimal years');

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

nccreate(vy_file_out, 'x','Dimensions',{'x' size(x,2)},'Datatype','double');
ncwriteatt(vy_file_out,'x','long_name','Easting');
ncwriteatt(vy_file_out,'x','units','meters easting, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vy_file_out,'x',x);

nccreate(vy_file_out, 'y','Dimensions',{'y' size(y,2)},'Datatype','double');
ncwriteatt(vy_file_out,'y','long_name','Northing');
ncwriteatt(vy_file_out,'y','units','meters northing, Antarctic Polar Stereographic, EPSG: 3031');
ncwrite(vy_file_out,'y',y);

nccreate(vy_file_out,'vy','Dimensions',{'x' 'y' 'time'},'Datatype','double');
ncwriteatt(vy_file_out,'vy','long_name','velocity in the y direction')
ncwriteatt(vy_file_out,'vy','units','meters per day');
ncwriteatt(vy_file_out,'vy','code','PyCorr v1.11');

%%

for file = 1:size(vx_files);
    a = strcat('Processing file_',num2str(file),'_of_',num2str(size(vx_files,1)));
    disp(a);
    
    %read in vx
    expression = strcat('[vx, vx_info] = geotiffread(''',vx_files(file).name,''');');
    eval(expression);
    
    %read in vy
    expression = strcat('[vy, vy_info] = geotiffread(''',vy_files(file).name,''');');
    eval(expression);
    
    %read in delcorr
    expression = strcat('[delcorr, delcorr_info] = geotiffread(''',delcorr_files(file).name,''');');
    eval(expression);

    %filter with a delcor > 0.15
    delcor_threshold = 0.15;
    expression = strcat('index = find(delcorr > delcor_threshold);');
    eval(expression);
    
    vx_filtered = zeros(size(vx));
    vx_filtered(index) = vx(index);

    vy_filtered = zeros(size(vy));
    vy_filtered(index) = vy(index);
    
    %write out to appropriate place
    ncwrite(vx_file_out,'time',time(file),[file]);
    ncwrite(vx_file_out,'vx',vx_filtered',[1 1 file]);
    ncwrite(vy_file_out,'time',time(file),[file]);
    ncwrite(vy_file_out,'vy',vy_filtered',[1 1 file]);
end

clear ans expression file index a
clear expression delcor_threshold delcorr delcorr_info
clear vx vx_filtered vx_info
clear vy vy_filtered vy_info

clear time
clear x* y*

clear vx_file_out vy_file_out
clear delcorr_files vx_files vy_files
clear doy_end clear doy_start duration path row