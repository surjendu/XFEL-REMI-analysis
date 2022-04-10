tic;
% parpool('local',8);

clear all;
clc;
close all;
addpath('r0064')
data_file='CORR-R0064-REMI01-S0000';

nfiles=1;

t_bin = [0 5000]; 

% extract all hits multiple files

t0=  -108.032990;
k0=9.607321;
shift = 0;

tic
t=[];x=[];y=[];m=[];
t_i=[];x_i=[];y_i=[];m_i=[];

 parfor i=1:nfiles 
    
filename=num2str(i-1);
files= [data_file,filename,'.h5'];

% dat_dir = '/gpfs/exfel/exp/SQS/202101/p002448/proc/'+run+'/'
% datf = h5py.File(dat_dir+'CORR-R0018-REMI01-S00000.h5','r')
% 
% datf = h5py.File('CORR-R0018-REMI01-S00000.h5','r')

data = h5read(files,'/INSTRUMENT/SQS_REMI_DLD6/DET/TOP:output/rec/hits');
t_1= data.t';
x_1= data.x';
y_1 = data.y';
m_1= data.m';

t_i=[t_i;t_1]; x_i=[x_i;x_1]; y_i=[y_i;y_1]; m_i=[m_i;m_1];

 end
 
t=[t_i(:,1);t_i(:,2);t_i(:,3);t_i(:,4);t_i(:,5);t_i(:,6);t_i(:,7);t_i(:,8);t_i(:,9);t_i(:,10)];
x=[x_i(:,1);x_i(:,2);x_i(:,3);x_i(:,4);x_i(:,5);x_i(:,6);x_i(:,7);x_i(:,8);x_i(:,9);x_i(:,10)];
y=[y_i(:,1);y_i(:,2);y_i(:,3);y_i(:,4);y_i(:,5);y_i(:,6);y_i(:,7);y_i(:,8);y_i(:,9);y_i(:,10)];
m=[m_i(:,1);m_i(:,2);m_i(:,3);m_i(:,4);m_i(:,5);m_i(:,6);m_i(:,7);m_i(:,8);m_i(:,9);m_i(:,10)];
toc
%%
close all
hit_counts = max(t_i./t_i,0);
hit_counts=sum(hit_counts,1);
plot(hit_counts, 'b', 'linewidth',2)
xlabel('Detetcted ion per shot ','FontWeight', 'normal','FontName', 'Arial');
ylabel('No. of shots', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20)
xlim([0 50]);

%% plot tof of for different hits
%plot tof in ns
binsize_t=2; %us   
% edges_t=[0:binsize_t:12000];
edges_t=[0:binsize_t:t_bin(end)];
i_t=1:length(edges_t)-1;
bincent_t=[];
bincent_t(i_t)=(edges_t(i_t)+edges_t(i_t+1))/2;
%
%plot first-five hits

t1=t_i(:,1+shift);  %j_gate = t_1 > t_bin(1)  &  t_1 < t_bin(end);  t_1=t_1(j_gate);
x1=x_i(:,1+shift);   %x_1=x_1(j_gate);
y1=y_i(:,1+shift);  % y_1=y_1(j_gate);
m1=m_i(:,1+shift);  % m_1=m_1(j_gate);
%
t2=t_i(:,2+shift);  %j_gate = t_2 > t_bin(1)  &  t_2 < t_bin(end);  t_2=t_2(j_gate);
x2=x_i(:,2+shift);  % x_2=x_2(j_gate);
y2=y_i(:,2+shift);  % y_2=y_2(j_gate);
m2=m_i(:,2+shift);   %m_2=m_2(j_gate);


%
t3=t_i(:,3+shift);  %j_gate = t_3 > t_bin(1)  &  t_3 < t_bin(end);  t_3=t_3(j_gate);
x3=x_i(:,3+shift);  % x_3=x_3(j_gate);
y3=y_i(:,3+shift);  % y_3=y_3(j_gate);
m3=m_i(:,3+shift);  % m_3=m_3(j_gate);

t4=t_i(:,4+shift);  % j_gate = t_4 > t_bin(1)  &  t_4 < t_bin(end); t_4=t_4(j_gate);
x4=x_i(:,4+shift);  % x_4=x_4(j_gate);
y4=y_i(:,4+shift);  % y_4=y_4(j_gate);
m4=m_i(:,4+shift);  % m_4=m_4(j_gate);

t5=t_i(:,5+shift);  % j_gate = t_5 > t_bin(1)  &  t_5 < t_bin(end); t_5=t_5(j_gate);
x5=x_i(:,5+shift);  % x_5=x_5(j_gate);
y5=y_i(:,5+shift);  % y_5=y_5(j_gate);
m5=m_i(:,5+shift);  % m_5=m_5(j_gate);
%

[t1_counts,t_edges]=histcounts(t1,edges_t);
[t2_counts,t_edges]=histcounts(t2,edges_t);
[t3_counts,t_edges]=histcounts(t3,edges_t);
[t4_counts,t_edges]=histcounts(t4,edges_t);
[t5_counts,t_edges]=histcounts(t5,edges_t);

%
ker_allhits_hist=[bincent_t' (t1_counts+t2_counts+t3_counts+t4_counts+t5_counts)'];
dlmwrite('ker_allhits_hist.csv',ker_allhits_hist);

figure 
plot(bincent_t,t1_counts,bincent_t,t2_counts,bincent_t,t3_counts,bincent_t,t4_counts,bincent_t,t5_counts,'LineWidth',1) 
legend({'first hit','second hit','third hit','fourth hit','fifth hit'}, 'FontSize', 25, 'FontWeight', 'normal','FontName', 'Arial')
% figure 
% plot(bincent_t,t1_counts+t2_counts+t3_counts+t4_counts+t5_counts,'LineWidth',1) 


  
% grid on
xlabel('TOF/ns','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20)
xlim([ t_bin(1)  t_bin(end)]);



%%  xtof

close all

tic;
% x=x(j_gate);
tof_range = 0;
tedges =   t_bin(1): 1 : t_bin(end);
xedges =  -140:0.5:140;
yedges =  -140:0.5:140;

%plot xtof
count_xtof_1 = histcounts2(t1,x1,tedges,xedges);
count_mod_xtof_1 = max(count_xtof_1,1); 

count_xtof_2 = histcounts2(t2,x2,tedges,xedges);
count_mod_xtof_2 = max(count_xtof_1,1); 

count_xtof_3 = histcounts2(t3,x3,tedges,xedges);
count_mod_xtof_3 = max(count_xtof_3,1); 

count_xtof_4 = histcounts2(t4,x4,tedges,xedges);
count_mod_xtof_4 = max(count_xtof_4,1); 

count_xtof_5 = histcounts2(t5,x5,tedges,xedges);
count_mod_xtof_5 = max(count_xtof_5,1); 

figure
subplot(2,1,1)
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
imagesc( tedges,xedges, (count_xtof_1'+count_xtof_2'+count_xtof_3'+count_xtof_4'+count_xtof_5'));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
% xlabel('TOF /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(X) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',20)
set(gca,'colorscale','log');
xticks([ ])


%plot xtof
count_ytof_1 = histcounts2(t1,y1,tedges,xedges);
count_mod_ytof_1 = max(count_ytof_1,1); 

count_ytof_2 = histcounts2(t2,y2,tedges,xedges);
count_mod_ytof_2 = max(count_ytof_1,1); 

count_ytof_3 = histcounts2(t3,y3,tedges,xedges);
count_mod_ytof_3 = max(count_ytof_3,1); 

count_ytof_4 = histcounts2(t4,y4,tedges,xedges);
count_mod_ytof_4 = max(count_ytof_4,1); 

count_ytof_5 = histcounts2(t5,y5,tedges,xedges);
count_mod_ytof_5 = max(count_ytof_5,1); 

subplot(2,1,2)
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
imagesc( tedges,yedges, (count_ytof_1'+count_ytof_2'+count_ytof_3'+count_ytof_4'+count_ytof_5'));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
xlabel('TOF /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',20)
set(gca,'colorscale','log');
%%
close all
%plot xy
count_xy_1 = histcounts2(x1,y1,xedges,yedges);
count_mod_xy_1 = max(count_xy_1,1); 
count_xy_2 = histcounts2(x2,y2,xedges,yedges);
count_mod_xy_2 = max(count_xy_2,1); 
count_xy_3 = histcounts2(x3,y3,xedges,yedges);
count_mod_xy_3 = max(count_xy_3,1); 
count_xy_4 = histcounts2(x4,y4,xedges,yedges);
count_mod_xy_4 = max(count_xy_4,1); 
count_xy_5 = histcounts2(x5,y5,xedges,yedges);
count_mod_xy_5 = max(count_xy_5,1); 


figure
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
imagesc( xedges,yedges, (count_xy_1'+count_xy_2'+count_xy_3'+count_xy_4'+count_xy_5'));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
xlabel('ion position(X) /mm', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',20)
set(gca,'colorscale','log');
pbaspect([1 1 1]);



%%  pipico
t1=[t_1;t_1;t_1;t_1;t_2;t_2;t_2;t_3;t_3;t_4];
t2=[t_2;t_3;t_4;t_5;t_3;t_4;t_5;t_4;t_5;t_5];
% clc
% t_1_test=t_1(1:50000,1);
% combos = combntns(t_1_test,2);
% 
% t1=combos(:,1);
% t2=combos(:,2);
%%
Xedges = min(t1):5:max(t1);
Yedges = min(t2):5:max(t2);

count = histcounts2(t1,t2,Xedges,Yedges);

%sum(sum(count))
% subplot(1,2,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (count'));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
axis equal
pbaspect([1 1 1]);
hold on;
set(gca,'FontSize',25)
set(gca,'colorscale','log');
% scatter([sct_peak_tof_1], [sct_peak_tof_2],300,'x','g','LineWidth',3)

xlabel('TOF_1 /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('TOF_2 /ns', 'FontWeight', 'normal','FontName', 'Arial');





%%
for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--y', frag_m_z_str(i)); 
end

%

%plot ytof
count_ytof = histcounts2(t_reg(:,1),y_reg(:,1),tedges_reg_gate,yedges_reg_gate);
count_mod_ytof = max(count_ytof,1); 
% count_mod = count; 
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
figure;
imagesc( tedges_reg_gate,yedges_reg_gate, ((count_mod_ytof')));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
xlabel('TOF /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');
toc;

for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--y', frag_m_z_str(i));
end


%%   ytof


%%
% Pipico and TRIPICO plot then momentum calculation; Question how to choose
% coincident files?

%%
clc; close all;

disp('Importing data now')
files= dir(['Run',num2str(fileno),'_XYTallbin1.out']);
%files= dir('XYTallbin1.out');
file = files';
fileID_real=fopen(file.name);
x=fread(fileID_real,[7 Inf],'*float');
fclose(fileID_real);
x = x';


id = find(x(:,1));
measurement.files = {files.name};
measurement.data.raw.delay_position = x(id,1);
measurement.data.raw.PD_voltage = x(id,5);
measurement.data.raw.XYT.frst = x(id, 2:4); 
measurement.data.raw.laser_shot = x(id,6);


clearvars x
disp('Done Importing');
toc
%%
%binning TOF
x = measurement.data.raw.XYT.frst(:,1); %x position in mm
y = measurement.data.raw.XYT.frst(:,2); %y position in mm
t = measurement.data.raw.XYT.frst(:,3); %time of flight in ns
% laser_shot=measurement.data.raw.laser_shot = x(:,4);
t_bin = [0 max(t)]; %tof range    
% t_bin = [0 12000]; %tof range
j_gate = t > t_bin(1)  &  t < t_bin(2);

t=t(j_gate);
x=x(j_gate);
y=y(j_gate);

% clearvars measurement;

%plot tof in ns
binsize_t=0.25; %us   
% edges_t=[0:binsize_t:12000];
edges_t=[0:binsize_t:t_bin(end)-10000];
i_t=1:length(edges_t)-1;
bincent_t=[];
bincent_t(i_t)=(edges_t(i_t)+edges_t(i_t+1))/2;
 
[t_counts,t_edges]=histcounts(t,edges_t);
t_hist=[bincent_t' t_counts'];

 dlmwrite('t_hist.csv',t_hist);


figure 
plot(bincent_t,t_counts,'LineWidth',1,'Color','k')   
  
% grid on
xlabel('TOF/ns','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20)


% %Initial velocity V0x in mm/ns
% %%
% for i=1:length(frag_m_z)
%     xline(  C*(frag_m_z(i)*mpvsme)^0.5+t0, '--r', frag_m_z_str(i))
% end
%% Import data all data
tic
clc
measurement_data_raw_delay_position = [];
measurement_data_raw_PD_voltage = [];
measurement_data_raw_XYT_frst = [];
measurement_data_raw_laser_shot = [];

disp('Importing data now')

for i=1:nfiles 
filename=num2str(i);
files= dir(['Run',num2str(fileno),'_XYTallbin',filename,'.out']);
file = files';
fileID_real=fopen(file.name);
x=fread(fileID_real,[7 Inf],'*float');
fclose(fileID_real);
x = x';

id = find(x(:,1));%id = find(x(:,1)>0);
delay_position = x(id,1);
PD_voltage = x(id,5);
XYT_frst = x(id,2:4);
laser_shot =  x(id,6); 

t = XYT_frst(:,3);
j_gate = t > t_bin(1)  &  t < t_bin(2);

delay_position=delay_position(j_gate,:);
PD_voltage=PD_voltage(j_gate,:);
XYT_frst=XYT_frst(j_gate,:);
laser_shot = laser_shot(j_gate,:);

measurement_data_raw_delay_position = [measurement_data_raw_delay_position;delay_position];
measurement_data_raw_PD_voltage = [measurement_data_raw_PD_voltage;PD_voltage];
measurement_data_raw_XYT_frst = [measurement_data_raw_XYT_frst; XYT_frst];
measurement_data_raw_laser_shot = [measurement_data_raw_laser_shot;laser_shot];

clearvars x
disp(['Done Importing',' ',filename]);
pause(0.1);

end


measurement.files = {files.name};
measurement.data.raw.delay_position =measurement_data_raw_delay_position;
measurement.data.raw.PD_voltage = measurement_data_raw_PD_voltage;
measurement.data.raw.XYT.frst = measurement_data_raw_XYT_frst;
measurement.data.raw.laser_shot = measurement_data_raw_laser_shot;

clearvars measurement_data_raw_delay_position measurement_data_raw_event_number measurement_data_raw_XYT_frst delay_position event_number;

save('measurement', 'measurement', '-v7.3'); % the measurement is the file name second is the variable itself
plot(measurement.data.raw.delay_position)
toc
%%
clc
close all
load('measurement');
figure
fck=measurement.data.raw.delay_position;
plot(fck)
length(fck)
clearvars fck

x = measurement.data.raw.XYT.frst(:,1); %x position in mm
y = measurement.data.raw.XYT.frst(:,2); %y position in mm
t = measurement.data.raw.XYT.frst(:,3); %time of flight in ns
delay_step=measurement.data.raw.delay_position;
% laser_shot=measurement.data.raw.laser_shot = x(:,4);
% t_bin = [0 max(t)]; %tof range    
% t_bin = [0 12000]; %tof range
j_gate = t > t_bin(1)  &  t < t_bin(2);

t=t(j_gate);
x=x(j_gate);
y=y(j_gate);
delay_step=delay_step(j_gate);
% clearvars measurement;

%plot tof in ns
binsize_t=1; %us   
% edges_t=[0:binsize_t:12000];
edges_t=[0:binsize_t:t_bin(end)];
i_t=1:length(edges_t)-1;
bincent_t=[];
bincent_t(i_t)=(edges_t(i_t)+edges_t(i_t+1))/2;
 
[t_counts,t_edges]=histcounts(t,edges_t);
t_hist=[bincent_t' t_counts'];

% dlmwrite('t_hist.csv',t_hist);


figure 
plot(bincent_t,t_counts,'LineWidth',1,'Color','k')   
  
% grid on
xlabel('TOF/ns','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20)


% %Initial velocity V0x in mm/ns
% %%
% for i=1:length(frag_m_z)
%     xline(  C*(frag_m_z(i)*mpvsme)^0.5+t0, '--r', frag_m_z_str(i))
% end
%%
run_name = '200 mW';

mmns_to_au_conv=0.4571; %mm/ns to a.u. conversion
% auto_eV = 27.211385;
har_to_ev=  27.211396;
mpvsme=1822.888486424682;
q=1;
H = 1.00782503223;
D = 2.01410177812;
He =  4.00260325413;
Li = 7.0160034366;
C = 12;
N = 14.00307400443;
O = 15.994914619;
F = 18.99840316273;
Ne =  19.9924401762;
S = 31.9720711744;
Cl = 34.968852682;
Ar = 39.9623831237;
Br79 = 78.9183376;
Br81 = 80.9162906;
I = 126.9044719;
C13 = 13.00335483521;



t0 =19.174779;
k0=16.513980*0.99991;

% m=(C+H+Br79*3)*mpvsme-q;
          

%  frag={{'H',1,1}  {'Br79',1,1} {'H',1,'Br79',1,1} {'Br79',2,1}  {'N',2,1} {'C',1,'H',1,'Br79',1,1} {'C',1,'Br79',2,1}   {'C',1,'H',1,'Br79',2,1}      {'C',1,'H',1,'Br79',3,2} {'C',1,'H',1,'Br79',3,1}};

frag={ {'C',1,'Br79',2,1}  };

% frag={{'C',1,'H',1,'Br79',2,1}};
 
frag_m=[]; 
frag_m_z=[];
charge_z=[];
 frag_m_z_str=string.empty;
 for i=1:length(frag)
     frag_m_z_int_sum=0;
    frag_m_z_int_str=[];
     for j=1:(length(frag{i})-1)/2
          frag_m_z_int=frag{i}{2*j}*eval(frag{i}{2*j-1});
          frag_m_z_int_sum = frag_m_z_int_sum + frag_m_z_int;
          frag_m_z_int=[(frag{i}{2*j-1}),'_',num2str(frag{i}{2*j})];
          frag_m_z_int_str=[frag_m_z_int_str,frag_m_z_int];

     end
     charge_z=[charge_z,frag{i}{2*j+1}];
     frag_m=[frag_m, frag_m_z_int_sum];
     frag_m_z=[frag_m_z, frag_m_z_int_sum/frag{i}{2*j+1}];
     frag_m_z_str(i)=[frag_m_z_int_str '^' num2str(frag{i}{2*j+1}) '^' '+' ];
     
 end
[charge_z;frag_m;frag_m_z]
 frag_m=frag_m*mpvsme-1; %q=1
 frag_m_z=frag_m_z*mpvsme-1; %q=1
 
  for i=1:length(frag_m_z)
       hold on
    xline(  k0*(frag_m_z(i))^0.5+t0, '--r', frag_m_z_str(i),'LabelVerticalAlignment' ,'bottom','fontweight','bold','fontsize',16);
     end



%% Code for determining parameters
clc 
 
%  tof_gate = [9212 9236]; % for CH79Br2
tof_gate = [9200 9210];% for CH79Br81Br

%Position gates
% 
% x_gate = [-40 , 40];
% y_gate = [-40, 40];


x_gate = [-1 , 10];
y_gate = [-1, 10];


%regional gate
reg_gate = t > tof_gate(1)  &  t < tof_gate(2) &  x > x_gate(1)  &  x < x_gate(2)&  y > y_gate(1)  &  y < y_gate(2);

t_reg=t(reg_gate);
x_reg=x(reg_gate);
y_reg=y(reg_gate);
delay_step=delay_step(reg_gate);

%%
close all

tic;
% x=x(j_gate);
tof_range = 0;
tedges_reg_gate=tof_gate(1)-tof_range:1:tof_gate(2)+tof_range;
xedges_reg_gate = x_gate(1):0.2:x_gate(2);
yedges_reg_gate = y_gate(1):0.2:y_gate(2);

%plot ytof
count_xtof = histcounts2(t_reg(:,1),x_reg(:,1),tedges_reg_gate,xedges_reg_gate);
count_mod_xtof = max(count_xtof,1); 
%count_mod = count; 

clc;close all;
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
imagesc( tedges_reg_gate,xedges_reg_gate, ((count_xtof')));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
xlabel('TOF /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(X) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');

for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--y', frag_m_z_str(i)); 
end

%

%plot ytof
count_ytof = histcounts2(t_reg(:,1),y_reg(:,1),tedges_reg_gate,yedges_reg_gate);
count_mod_ytof = max(count_ytof,1); 
% count_mod = count; 
% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 0;
figure;
imagesc( tedges_reg_gate,yedges_reg_gate, ((count_mod_ytof')));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;
xlabel('TOF /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');
toc;

for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--y', frag_m_z_str(i));
end

%%
%Calculating the momentum COLTRIMS

x0_cent = 4.7-0.45 -0.25;
V0x_cent = 0;

%Center of detector in mm 
%y0_cent=0.6859;
y0_cent=-4.7-0.5 + 0.45;
V0y_cent= 1.1e-03; 
 
%z0 in mm
z0 = -10;
V0z_cent =0;    


% t0 =19.174779;
% k0=16.513980*0.99991;

%Ion Flight Distance in mm
Lion = 237.05;
 
px = frag_m(1)*((x_reg(:,1)-x0_cent)./(t_reg(:,1)-t0)-V0x_cent)*mmns_to_au_conv;
py = frag_m(1)*((y_reg(:,1)-y0_cent)./(t_reg(:,1)-t0)-V0y_cent)*mmns_to_au_conv;
pz = frag_m(1)*((Lion-z0)./(t_reg(:,1)-t0)-0.5*(t_reg(:,1)-t0)*(q * 2 * (Lion - z0) / (frag_m_z(1)*k0^2))-V0z_cent)*mmns_to_au_conv;
% frag_m(1)*(  (Lion-z0)./(t1_gate_rot-t0)  - 0.5* (t1_gate_rot-t0)*(1 * 2 * (Lion - z0) / (frag_m_z(1)*k0^2) )-v0z)*mmns_to_au_conv;
p=sqrt(px.^2+py.^2+pz.^2);   

    
%plot tof in ns
binsize_px=1; 
binsize_py=4; 
binsize_pz=0.25;
binsize_p=0.5;

edge_max=ceil(max([max([px; py; pz]); -min([px; py; pz])]));

edges_px=[-edge_max:binsize_px:edge_max];
edges_py=[-edge_max:binsize_py:edge_max]; 
edges_pz=[-edge_max:binsize_pz:edge_max];
edges_p=[min(p):binsize_p:max(p)];

i_px=1:length(edges_px)-1;
i_py=1:length(edges_py)-1;
i_pz=1:length(edges_pz)-1;
i_p=1:length(edges_p)-1;

bincent_px=[];
bincent_py=[];
bincent_pz=[];
bincent_p=[];


bincent_px(i_px)=(edges_px(i_px)+edges_px(i_px+1))/2;
bincent_py(i_py)=(edges_py(i_py)+edges_py(i_py+1))/2;
bincent_pz(i_pz)=(edges_pz(i_pz)+edges_pz(i_pz+1))/2;
bincent_p(i_p)=(edges_p(i_p)+edges_p(i_p+1))/2;


[px_counts,px_edges]=histcounts(px,edges_px);
[py_counts,py_edges]=histcounts(py,edges_py);
[pz_counts,pz_edges]=histcounts(pz,edges_pz);
[p_counts,p_edges]=histcounts(p,edges_p);

p_allhits_hist={[bincent_px', px_counts'], [bincent_py', py_counts'], [bincent_pz' pz_counts'], [bincent_p' p_counts']};
%  dlmwrite('p_allhits_hist.csv',p_allhits_hist);

close all;
%figure 
subplot(2,2,1)
plot(bincent_px,px_counts,'r',bincent_px,fliplr(px_counts),'b',bincent_px,(px_counts+fliplr(px_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-50 50]);
ylim([0 max(px_counts)*1.1]);
xlabel('px/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'px','-px','px\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

 
subplot(2,2,2)
%figure
plot(bincent_py,py_counts,'r',bincent_py,fliplr(py_counts),'b',bincent_py,(py_counts+fliplr(py_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-100 100]);
ylim([0 max(py_counts)*1.1]);
xlabel('py/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'py','-py','py\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3) 
%figure
plot(bincent_pz,pz_counts,'r',bincent_pz,fliplr(pz_counts),'b',bincent_pz,(pz_counts+fliplr(pz_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-20 20]);
ylim([0 max(pz_counts)*1.1]);
xlabel('pz/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'pz','-pz','pz\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,4) 
%figure
plot(bincent_p,p_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-100 100]);
ylim([0 max(p_counts)*1.1]);
xlabel('\Sigmap/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmap'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

%%
%plotting kinetic energy
ker =  har_to_ev*(px.^2+py.^2+pz.^2)/(2*frag_m(1));
ker_rng = [0.0 1]; 
j_ker= ker > ker_rng(1) &  ker < ker_rng(2); 

ker=ker(j_ker);
delay_step = delay_step(j_ker);

binsize_ker=0.001; 
edges_ker=[0:binsize_ker:ker_rng(end)];
i_ker=1:length(edges_ker)-1;
bincent_ker=[];
bincent_ker(i_ker)=(edges_ker(i_ker)+edges_ker(i_ker+1))/2;
[ker_counts,ker_edges]=histcounts(ker,edges_ker);

% ker_allhits_hist=[bincent_ker' ker_counts'];
% dlmwrite('ker_allhits_hist.csv',ker_allhits_hist);

close all;
figure 
plot(bincent_ker,ker_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
grid on;
xlabel('kinetic energy/eV',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'ker'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

% ker_allhits_hist=[bincent_ker' ker_counts'];
% dlmwrite('ker_allhits_hist.csv',ker_allhits_hist);
%% delay dependent
close all
delay_range=[-200 2500 5];
binsize_delay=delay_range(3);
binsize_ke=0.01;    

% delay_step=delay_step(reg_gate);
delay=(delay_step-1)*delay_range(3)+delay_range(1);



edges_delay=[min(delay):binsize_delay:max(delay)];
edges_ke=[0:binsize_ke:ker_rng(end)];

i_delay=1:length(edges_delay)-1;
bincent_delay=[];
bincent_delay(i_delay)=(edges_delay(i_delay)+edges_delay(i_delay+1))/2;

i_ke=1:length(edges_ke)-1;
bincent_ke=[];
bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2; 


figure
count_delay = histcounts2(delay(:,1),ker(:,1),edges_delay,edges_ke);
% count_delay_mod = max(count_delay,1); 

subplot(2,2,1)
% myColorMap=flipud(hot);
% myColorMap(1,:) = 1;
imagesc( edges_delay, edges_ke,((count_delay')));
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
% pbaspect([1 1 1]);
xlabel('pump-probe delay/fs', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('KER / eV', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',16)
set(gca,'colorscale','log');


subplot(2,2,2)
plot(bincent_ke,sum(count_delay),'r','LineWidth',2);
xlabel('KER / eV','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
% set(gca, 'YScale', 'log')
set(gca,'FontSize',16)
% % xlim([-1 1]);
ylim([0 max(sum(count_delay))*1.1]);
%pbaspect([1 1 1]);

subplot(2,2,3)
plot(bincent_delay,sum(count_delay'),'r','LineWidth',2);

ker_delay_allhits_hist=[bincent_delay' (sum(count_delay')')];
dlmwrite('ker_delay_allhits_hist.csv',ker_delay_allhits_hist);

xlabel('pump-probe delay / fs','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
% set(gca, 'YScale', 'log')
set(gca,'FontSize',16)
xlim([delay_range(1) delay_range(2)]);
ylim([0 max(sum(count_delay'))*1.1]);
%pbaspect([1 1 1]);

%%  delay step selection________looks like its difficult to separate the loops in this method--so right now go the old method Apr03,2021
% close all
% delay_range=[-200 2500 5];
% binsize_delay=1;% 3cause it is step size not delay in time
% binsize_ke=0.05;    
% 
% j_loop = delay_step > 0 & delay_step <= length(delay_range(1) : delay_range(3) : delay_range(2))*23;
% delay_step=delay_step(j_loop);
% ker = ker (j_loop);
% 
% 
% delay=delay_step;
% 
% 
% edges_delay=[min(delay):binsize_delay:max(delay)];
% edges_ke=[0:binsize_ke:5];
% 
% i_delay=1:length(edges_delay)-1;
% bincent_delay=[];
% bincent_delay(i_delay)=(edges_delay(i_delay)+edges_delay(i_delay+1))/2;
% 
% i_ke=1:length(edges_ke)-1;
% bincent_ke=[];
% bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2; 
% 
% 
% figure
% 
% count_delay = histcounts2(delay(:,1),ker(:,1),edges_delay,edges_ke);
% count_delay_mod = max(count_delay,1); 
% % myColorMap=flipud(hot);
% myColorMap=flipud(jet);
% myColorMap(1,:) = 0;
% imagesc( edges_delay, edges_ke,((count_delay_mod')));
% colormap(myColorMap);
% colorbar 
% axis xy;
% % axis([-300 300 -300 300]);
% % xticks([-300 -150 0 150 300]);
% % yticks([-300 -150 0 150 300]);
% % pbaspect([1 1 1]);
% xlabel('delay step', 'FontWeight', 'normal','FontName', 'Arial');
% ylabel('KER / eV', 'FontWeight', 'normal','FontName', 'Arial');
% set(gca,'FontSize',16)
% set(gca,'colorscale','log');

%%  ________looks like its difficult to separate the loops in this method--so right now go the old method Apr03,2021
% close all 
% delay_range=[-500 4000 30];
% binsize_delay=delay_range(3);
% binsize_ke=0.05;    
% 
% % delay_step=delay_step(j_KE_all);
% 
% delay_step_new=rem(delay_step,length(delay_range(1):delay_range(3):delay_range(2)));
% 
% for i=1:length(delay_step_new)
%     if delay_step_new(i)>0
%         x_step=delay_step_new(i);
%     else x_step=length(delay_range(1):delay_range(3):delay_range(2));
%     end
%     delay_step(i)=x_step;
% end
% % delay_step=delay_step(j_KE_all);
% delay=(delay_step-1)*delay_range(3)+delay_range(1);
% 
% 
% edges_delay=[min(delay)-binsize_delay/2:binsize_delay:max(delay)+binsize_delay/2];
% edges_ke=[0:binsize_ke:5];
% 
% i_delay=1:length(edges_delay)-1;
% bincent_delay=[];
% bincent_delay(i_delay)=(edges_delay(i_delay)+edges_delay(i_delay+1))/2;
% 
% i_ke=1:length(edges_ke)-1;
% bincent_ke=[];
% bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2; 
% 
% 
% figure
% 
% count_delay = histcounts2(delay(:,1),ker(:,1),edges_delay,edges_ke);
% % save('count_delay', 'count_delay', '-v7.3'); 
% % dlmwrite('count_delay.csv',count_delay);
% % delay_ker = [delay(:,1),KE12(:,1)];
% % dlmwrite('delay_ker.csv',delay_ker);
% count_delay_mod = max(count_delay,1); 
% 
% subplot(2,2,1)
% % myColorMap=flipud(hot);
% myColorMap=flipud(jet);
% myColorMap(1,:) = 0;
% imagesc( edges_delay, edges_ke,((count_delay')));
% colormap(myColorMap);
% colorbar 
% axis xy;
% % axis([-300 300 -300 300]);
% % xticks([-300 -150 0 150 300]);
% % yticks([-300 -150 0 150 300]);
% % pbaspect([1 1 1]);
% xlabel('pump-probe delay/fs', 'FontWeight', 'normal','FontName', 'Arial');
% ylabel('KE_1+KE_2 / eV', 'FontWeight', 'normal','FontName', 'Arial');
% set(gca,'FontSize',16)
% set(gca,'colorscale','log');
% %set(gca,'colorscale','linear');
% 
% 
% subplot(2,2,2)
% plot(bincent_ke,sum(count_delay),'r','LineWidth',2);
% xlabel('KER / eV','FontWeight', 'normal','FontName', 'Arial');
% ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% %set(gca, 'XScale', 'linear')
% % set(gca, 'YScale', 'log')
% %set(gca,'FontSize',16)
% % % xlim([-1 1]);
% ylim([0 max(sum(count_delay))*1.1]);
% %pbaspect([1 1 1]);
% 
% subplot(2,2,3)
% plot(bincent_delay,sum(count_delay'),'r','LineWidth',2);
% xlabel('pump-probe delay / fs','FontWeight', 'normal','FontName', 'Arial');
% ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% set(gca, 'XScale', 'linear')
% % set(gca, 'YScale', 'log')
% set(gca,'FontSize',16)
% xlim([delay_range(1) delay_range(2)]);
% % ylim([0 max(sum(count_delay'))*1.1]);
% %pbaspect([1 1 1]);
% 
% % ker_delay=[bincent_delay',sum(count_delay')'];
% % dlmwrite('ker_delay.csv',ker_delay);
% 









%%
%plot tof in microsecond (t_mis)
t_mis=t/1000;
%t_mis=t;

binsize_t=1/1000; %us   
%binsize_t=1; %us 
edges_t=[0:binsize_t:17];
i_t=1:length(edges_t)-1;
bincent_t=[];
bincent_t(i_t)=(edges_t(i_t)+edges_t(i_t+1))/2;
 
[t_counts,t_edges]=histcounts(t_mis,edges_t);
t_hist=[bincent_t' t_counts'];
% 
% dlmwrite('t_hist.csv',t_hist);


figure 
close all;
plot(bincent_t,t_counts,'LineWidth',1,'Color','k')   
  
% grid on
xlabel('TOF/\mus','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20)

    
%%
%binning  XTOF
tic;
% x=x(j_gate);
Xedges = -50:0.2:50;

count = histcounts2(t_mis(:,1),x(:,1),edges_t,Xedges);
 count_mod = max(count,1); 
% count_mod = count; 

clc;close all;
%myColorMap=dlmread('mycolormap.txt');

% myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 1;
imagesc( edges_t, Xedges, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;

axis([0 12 -50 50]);
xtickformat('%.0f');
ytickformat('%.0f');
xticks(0 :2 :12);
yticks(-50 :25 :50);

xlabel('TOF /\mus', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(X) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');

toc;



for i=1:length(frag_m_z);
    xline(  (k0*(frag_m_z(i))^0.5+t0)/1000, '--r', frag_m_z_str(i));
end






%%
%binning  YTOF
tic

% y=y(j_gate);
Yedges = -50:0.2:50;

count = histcounts2(t_mis(:,1),y(:,1),edges_t,Yedges);
 count_mod = max(count,1); 
%count_mod = count; 

clc;close all;
%myColorMap=dlmread('mycolormap.txt');

%myColorMap=flipud(hot);
myColorMap=jet;
myColorMap(1,:) = 1;
imagesc( edges_t, Yedges, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
axis xy;

axis([0 12 -50 50]);
xtickformat('%.0f');
ytickformat('%.0f');
xticks(0 :2 :12);
yticks(-50 :25 :50);


xlabel('TOF /\mus', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
%clims=[1 7];
set(gca,'colorscale','log');


for i=1:length(frag_m_z);
    xline(  (k0*(frag_m_z(i))^0.5+t0)/1000, '--r', frag_m_z_str(i));
end


toc;