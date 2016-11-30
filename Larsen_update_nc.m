%%
%Directory with files
cd ~/Desktop/pycorr/Out/LarsenC/Crop/

%netCDF files
vx_file_out = 'LarsenC_vx_03-Nov-2016.nc';
vy_file_out = 'LarsenC_vy_03-Nov-2016.nc';

%%
%List of all files
vx_files = dir('L8*VX.tif');
vy_files = dir('L8*VY.tif');
delcorr_files = dir('L8*delcorr.tif');

%%
%Identifying new files
expression = strcat('doy_start = ncread(''',vx_file_out,''', ''doy_start'');');
eval(expression);

expression = strcat('doy_end = ncread(''',vx_file_out,''', ''doy_end'');');
eval(expression);

expression = strcat('year_start = ncread(''',vx_file_out,''', ''year_start'');');
eval(expression);

expression = strcat('year_end = ncread(''',vx_file_out,''', ''year_end'');');
eval(expression);

expression = strcat('duration = ncread(''',vx_file_out,''', ''duration'');');
eval(expression);

expression = strcat('path = ncread(''',vx_file_out,''', ''path'');');
eval(expression);

expression = strcat('row = ncread(''',vx_file_out,''', ''row'');');
eval(expression);

expression = strcat('time = ncread(''',vx_file_out,''', ''time'');');
eval(expression);

old_files = [];
for file = 1:length(path);
    if duration(file) < 100
        duration_s = strcat('0',num2str(duration(file)));
    else
        duration_s = num2str(duration(file));
    end
    
    if path(file)<10
        path_s = strcat('00',num2str(path(file)));
    elseif path(file) < 100
        path_s = strcat('0',num2str(path(file)));
    else
        path_s = num2str(path(file));
    end
    
    if row(file) < 10
        row_s = strcat('00',num2str(row(file)));
    elseif row(file) < 100
        row_s = strcat('0',num2str(row(file)));
    else
        row_s = num2str(row(file));
    end
    
    if doy_start(file) < 10
        doy_start_s = strcat('00',num2str(doy_start(file)));
    elseif doy_start(file) < 100
        doy_start_s = strcat('0',num2str(doy_start(file)));
    else
        doy_start_s = num2str(doy_start(file));
    end
    
    if doy_end(file) < 10
        doy_end_s = strcat('00',num2str(doy_end(file)));
    elseif doy_end(file) < 100
        doy_end_s = strcat('0',num2str(doy_end(file)));
    else
        doy_end_s = num2str(doy_end(file));
	end    
        
    old_files = [old_files; strcat('L8_Ant_v00_S8_',path_s,'_',row_s,'_',duration_s,'_',num2str(year_start(file)),'_',doy_start_s,'_',num2str(year_end(file)),'_',doy_end_s,'_hp_filt_3.0_vx.tif')];
end
clear *_s file expression

new_files = [];
for file = 1:length(vx_files)
    new_files = [new_files; vx_files(file).name];
end
clear file

new_files = setdiff(new_files, old_files, 'rows');
new_files = cellstr(new_files);

%%
%write new vy and delcorr file lists
new_files_vx = new_files;
new_files_vy = new_files;
new_files_delcorr = new_files;

for file = 1:size(new_files_vx,1)
    new_files_vy(file) = strrep(new_files_vy(file), 'vx', 'vy');
end

for file = 1:size(new_files_vx,1)
    new_files_delcorr(file) = strrep(new_files_delcorr(file, :), 'vx', 'delcorr');
end

clear file

%%
%background data

year_start_new = zeros(size(new_files));
doy_start_new = zeros(size(new_files));
year_end_new = zeros(size(new_files));
doy_end_new = zeros(size(new_files));
path_new = zeros(size(new_files));
row_new = zeros(size(new_files));
duration_new = zeros(size(new_files));
time_new = zeros(size(new_files));

for file = 1:size(new_files);
    filename = char(new_files(file));
    
    %Create year_start
    year_start_new(file) = str2num(filename(27:30));
    %create doy_start
    doy_start_new(file) = str2num(filename(32:34));
    %Create year_end
    year_end_new(file) = str2num(filename(36:39));
    %create doy_end
    doy_end_new(file) = str2num(filename(41:43));
    %create path
    path_new(file) = str2num(filename(15:17));
    %create row
    row_new(file) = str2num(filename(19:21));
    %create duration
    duration_new(file) = str2num(filename(23:25));
end

%create time
time_new = (year_start_new+doy_start_new/365+year_end_new+doy_end_new/365)/2;

%%
%delcor filter
%standard deviation filter
%add to nc file

stack_count = ncread('LarsenC_stack_info.nc','count')';
stack_mean = ncread('LarsenC_stack_info.nc','mean')';
stack_std = ncread('LarsenC_stack_info.nc','std')';
x = ncread('LarsenC_stack_info.nc','x')';
y = ncread('LarsenC_stack_info.nc','y')';

for file = 1:size(new_files_vx,1);
    a = strcat('Processing file_',num2str(file),'_of_',num2str(size(new_files_vx,1)));
    disp(a);
    
    %read in vx
    expression = strcat('[vx, vx_info] = geotiffread(''',new_files_vx(file),''');');
    eval(expression{:});
    
    %read in vy
    expression = strcat('[vy, vy_info] = geotiffread(''',new_files_vy(file),''');');
    eval(expression{:});
    
    %read in delcorr
    expression = strcat('[delcorr, delcorr_info] = geotiffread(''',new_files_delcorr(file),''');');
    eval(expression{:});

    %filter with a delcor > 0.15
    delcor_threshold = 0.15;
    expression = strcat('index = find(delcorr > delcor_threshold);');
    eval(expression);
    
    vx_filtered = zeros(size(vx));
    vx_filtered(index) = vx(index);

    vy_filtered = zeros(size(vy));
    vy_filtered(index) = vy(index);
    
    %mean/std filter
    vv_filtered = (vx.^2 + vy.^2).^0.5;
    vv_filtered(find(vv_filtered > 4)) = 0; %going WAY too fast
        
    index1 = find(vv_filtered < (stack_mean+stack_std*3));
    index2 = find(vv_filtered > (stack_mean-stack_std*3));
    index3 = find(vv_filtered);
    index = intersect(index1,index2);
    index = intersect(index,index3);
        
    vx = zeros(size(vx_filtered));
    vx(index) = vx_filtered(index);
    
    vy = zeros(size(vy_filtered));
    vy(index) = vy_filtered(index);
    
    vx_filtered = vx;
    vy_filtered = vy;
    
    %re-calculate count and mean
    pix = 5;
    parr = floor((size(y,2)-pix-pix)/11); %so that the parfor loop can be consecutive integers but still step 11 in real life

    for par_r = 1:parr
        r = par_r*11+pix; %this should actually be r = 1+pix+(2*pix+1)*(par_r-1)
     
        for c = (pix+1):(2*pix+1):(size(x,2)-pix)
            %vx_temp = vx_filtered((r-1-pix):(r-1+pix),(c+1-pix):(c+1+pix)); %because nc is zero indexed and matlab is 1-indexed
            vx_temp = vx_filtered((r-pix):(r+pix),(c-pix):(c+pix));
            vy_temp = vy_filtered((r-pix):(r+pix),(c-pix):(c+pix));
            vv_temp = (vx_temp.^2+vy_temp.^2).^0.5;
            
            %correct the mean
            if length(find(vx_temp)) > 0 & length(find(vv_temp)) > 0;
                stack_mean((r-pix):(r+pix),(c-pix):(c+pix)) = [(stack_mean((r-pix):(r+pix),(c-pix):(c+pix)).*stack_count((r-pix):(r+pix),(c-pix):(c+pix))+length(find(vx_temp))*mean(vv_temp(find(vv_temp))))/(stack_count(r,c)+length(find(vx_temp)))];
            end
            
            %correct the count
            stack_count((r-pix):(r+pix),(c-pix):(c+pix)) = stack_count((r-pix):(r+pix),(c-pix):(c+pix)) + ones(size(vx_temp))*length(find(vx_temp));
            
            %can't correct std
            %if updating many times, should re-run everything
        end
    end
    
    %write out to appropriate place
    ncwrite(vx_file_out,'time',time_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'duration',duration_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'year_start',year_start_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'doy_start',doy_start_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'year_end',year_end_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'doy_end',doy_end_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'path',path_new(file),[(size(old_files,1)+file)]);
    ncwrite(vx_file_out,'row',row_new(file),[(size(old_files,1)+file)]);

    ncwrite(vy_file_out,'time',time_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'duration',duration_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'year_start',year_start_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'doy_start',doy_start_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'year_end',year_end_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'doy_end',doy_end_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'path',path_new(file),[(size(old_files,1)+file)]);
    ncwrite(vy_file_out,'row',row_new(file),[(size(old_files,1)+file)]);

    ncwrite(vx_file_out,'vx',vx_filtered',[1 1 (size(old_files,1)+file)]);
    ncwrite(vy_file_out,'vy',vy_filtered',[1 1 (size(old_files,1)+file)]);
end

    
%read out mean/count stuff...
ncwrite('LarsenC_stack_info.nc','count',stack_count');
ncwrite('LarsenC_stack_info.nc','mean',stack_mean');

%%
%rename files to tell something about file history

%rename vx
movefile(vx_file_out,strcat(vx_file_out(1:(end-3)),'__',date,'.nc'));

%rename vy
movefile(vy_file_out,strcat(vy_file_out(1:(end-3)),'__',date,'.nc'));

clear