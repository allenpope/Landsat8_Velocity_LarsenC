%%
%
vx_file = 'LarsenC_vx_03-Nov-2016.nc';
vy_file = 'LarsenC_vy_03-Nov-2016.nc';

%%
%Time series analaysis of a chunk

time = ncread(vx_file,'time');

x = ncread(vx_file,'x');
y = ncread(vy_file,'y');
duration = ncread(vx_file,'duration');


vx_temp = ncread(vx_file,'vx',[975 400 1],[60 60 size(time,1)]); %dims x, y, time
%vx_temp = ncread(vx_file,'vx',[400 220 1],[350 20 size(time,1)]); %dims x, y, time
vx_temp(find(vx_temp == 0)) = nan;

vy_temp = ncread(vy_file,'vy',[975 400 1],[60 60 size(time,1)]); %dims x, y, time
%vy_temp = ncread(vy_file,'vy',[400 220 1],[350 20 size(time,1)]); %dims x, y, time
vy_temp(find(vy_temp == 0)) = nan;


vv_temp = (vx_temp.*vx_temp + vy_temp.*vy_temp).^0.5;

vv_mean = zeros(size(time,1),1);
%vv_count = zeros(293,1);
for t = 1:size(time,1)
    temp = vv_temp(:,:,t);
    %vv_count(t) = sum(isfinite(temp(:)));
    vv_mean(t) = nanmean(temp(:));
end

start_time = time-duration/365/2;
end_time = time+duration/365/2;

figure()
hold on
for t = 1:size(time,1)
    %The below would work in new matlab, but not 2016a
    %errorbar(time(t),vv_mean(t),1.5/duration(t),1.5/duration(t),duration(t)/365/2,duration(t)/365/2,'o'); %0.1 pixels (1.5 m) uncertainty divided by days
    
    %for 2016a, plot with a workaround:
    errorbar([start_time(t), end_time(t)],[vv_mean(t), vv_mean(t)],[1.5/duration(t), 1.5/duration(t)]);
end
hold off
set(gca,'Fontsize',24);
xlabel('Date [Year]');
ylabel('Velocity [meters per day]');

%this then goes into Illustrator to make it pretty...