%Close all windows and clear all workspace variables.
close all; clear;

%Energy Plus Output File Path (Shades Off)
filePath = 'Combined Outputs/ShadesOff.csv';
%Open a file with comma-separated values. (date/time on the first column).
fid = fopen(filePath, 'rt');
C1 = textscan(fid, ['%s' repmat('%f',1,5)],'Delimiter',',','HeaderLines',1);
fclose(fid);
 
%Energy Plus Output File Path (Shades On)
filePath = 'Combined Outputs/ShadesOn.csv';
%Open a file with comma-separated values. (date/time on the first column).
fid = fopen(filePath, 'rt');
C2 = textscan(fid, ['%s' repmat('%f',1,5)],'Delimiter',',','HeaderLines',1);
fclose(fid);

% convert date/time to serial date number
dt = datenum(C1{1}(:,1), 'mm/dd HH:MM:SS');
 nbrDays = round(dt(end) - dt(1));
  
%Combine all in one matrix
M1 = [dt C1{2} C1{3} C1{4} C1{5}];
%Combine all in one matrix
M2 = [dt C2{2} C2{3} C2{4} C2{5}];

startInt = datenum('01/01 00:00:00','mm/dd HH:MM:SS');
endInt   = datenum('01/02 00:00:00','mm/dd HH:MM:SS');
 
%Show on the command window the staring date
datestr(startInt,'mm/dd HH:MM:SS')
 
%Show the ending data
datestr(endInt,'mm/dd HH:MM:SS')
  
Int1 = getIntervalData(startInt, endInt, M1);
Int2 = getIntervalData(startInt, endInt, M2);


%+------------------------------------------------------------------------)
%|                            Time  
%+------------------------------------------------------------------------)

%'Sampling Time (seconds)'
Ts = round((M1(2,1) - M1(1,1))*3600*24)

%Serial Date Time
sDateTime = Int1(:,1); 

%Maximum time (seconds);
tMax = size(Int1,1)*Ts;
 
%Time scale (seconds)
%ts = 0:(tMax - 1);
 
%Sampling points 
tm = 0:Ts:(tMax - 1);
 
%Outdoor Ambient Temperature
Ta  = Int1(:,2);

%Zone Temperature
Tin_EP_off = Int1(:,5);
Tin_EP_on  = Int2(:,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%|>                        Show Plots                                    <|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%+---------------------------------------+
%|>  Outdoor and Zone Temperatures      <|
%+---------------------------------------+
figure(1);
hold on;
 
%Plot outdoor temperature (hourly time scale)
plot(sDateTime,Ta,'k');
 
%Energy Plus
plot(sDateTime,Tin_EP_off,'b-');
plot(sDateTime,Tin_EP_on,'r:');

xlabel('Time (hours)');
ylabel('Temperature ^\circC');
ylim([10 15.5]); 
datetick('x','hh');
title(datestr(startInt,'mm/dd'));
legend('T_a','T_{in}\_EP (ws=0)', 'T_{in}\_EP (ws=1)');
movegui('southwest');

%Save Figure to Files (eps, pdf)
mkdir('Exp4');
makeFiles('Exp4/Exp4a');

 
