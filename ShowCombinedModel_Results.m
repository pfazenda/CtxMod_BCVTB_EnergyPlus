%Close all windows and clear all workspace variables.
close all; clear;

%Energy Plus Output File Path
filePath = 'Output/Simple.csv';

%Open a file with comma-separated values. (date/time on the first column).
fid = fopen(filePath, 'rt');
 
C = textscan(fid, ['%s' repmat('%f',1,5)],'Delimiter',',','HeaderLines',1);
fclose(fid);
 
% convert date/time to serial date number
dt = datenum(C{1}(:,1), 'mm/dd HH:MM:SS');
 
nbrDays = round(dt(end) - dt(1));
  
%Combine all in one matrix
M = [dt C{2} C{3} C{4} C{5}];

startInt = datenum('01/01 00:00:00','mm/dd HH:MM:SS');
endInt   = datenum('01/02 00:00:00','mm/dd HH:MM:SS');
 
%Show on the command window the staring date
datestr(startInt,'mm/dd HH:MM:SS')
 
%Show the ending data
datestr(endInt,'mm/dd HH:MM:SS')
  
Int = getIntervalData(startInt, endInt, M);

%+------------------------------------------------------------------------)
%|                            Time  
%+------------------------------------------------------------------------)

%'Sampling Time (seconds)'
Ts = round((M(2,1) - M(1,1))*3600*24)

%Serial Date Time
sDateTime = Int(:,1); 

%Maximum time (seconds);
tMax = size(Int,1)*Ts;
 
%Time scale (seconds)
%ts = 0:(tMax - 1);
 
%Sampling points 
tm = 0:Ts:(tMax - 1);
 
%Outdoor Ambient Temperature
Ta  = Int(:,2);

%Zone Temperature
Tin_EP = Int(:,5);

%Surface Outside Face Temperature
T1_EP  = Int(:,4);

%Surface Inside Face Temperature
T4_EP  = Int(:,3);

%Initial temperatures
Tin_0 = Tin_EP(1);
T1_0  = T1_EP(1);
T4_0  = T4_EP(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%|>             Continuous - Time 3R2C Wall Model                        <|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
% Materials:
% Stucco
l1   =  0.0254;             % (m)
rho1 =    1858;             % (kg/m3)
cp1  =     837;             % (J/kg-K)
K1   =    0.69;             % (W/m-K) 

%Brick
l2   =     0.1;             % (m)
rho2 =    1922;             % (kg/m3)
cp2  =     837;             % (J/kg-K)
K2   =    0.73;             % (W/m-K) 
   
%Plaster
l3   =  0.019;              % (m)
rho3 =   1602;              % (kg/m3)
cp3  =    837;              % (J/kg-K)
K3   =  0.726;              % (W/m-K) 

%Glass
lglass =  0.003;            % (m)
Kglass =  0.8;              % (W/m-K) 

%Area of Windows
Awind = 4*1.5*4;            % (m^2)

%Area of the B. Envelope
A = 12*3*4 + 12*12 - Awind; % (m^2) 

Cz = 525765;

hext =  10.22;
hint =  3.076;
Rext = 1/(hext*A);   
Rint = 1/(hint*A);


R = l1/(K1*A) + l2/(K2*A) + l3/(K3*A) + Rext + Rint;
C = (rho1*cp1*l1 + rho2*cp2*l2 +rho3*cp3*l3)*A;

% Calculating the 3R2C wall model parameters

alpha1 = 0.23818828858317055376954213751263;
alpha2 = 0.16209971292445897672226108436541;
alpha3 = 0.59971199849237046950819677812195;

beta1 = 0.51685417244656046687598306989775;
beta2 = 0.48314582755343953312401693010225;
 
R1 = R*alpha1 - Rext;
R2 = R*alpha2;
R3 = R*alpha3 - Rint;

C1 = C*beta1;
C2 = C*beta2;

%Calculating the window parameters
%Rw = lglass/(Kglass*Awind);

UValue = 3.071;
Rshade = 0.005/(Awind*0.05);

RwOn = 1/(Awind*UValue) + Rshade;
Rw = 1/(Awind*UValue);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Model_Off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rea = Rext + R1;
Rie = Rint + R3;

a = [-(1/(C1*Rea) + 1/(C1*R2))   1/(C1*R2)                      0;
       1/(C2*R2)             -(1/(C2*R2) + 1/(C2*Rie))      1/(C2*Rie) 
         0                       1/(Cz*Rie)     -(1/(Cz*Rie) + 1/(Rw*Cz))]; 

b = [1/(C1*Rea); 0; 1/(Rw*Cz)];

c = [Rext/(Rext + R1)     0                             0;
       1                  0                             0;
       0                  1                             0;
       0                 Rint/(R3 + Rint)  R3/(R3 + Rint);
       0                  0                            1];
 
d = [R1/(Rext + R1);0 ;0 ;0 ;0];
    
%Create the state-space model.
sys = ss(a,b,c,d);

%Simulate model (obtain the Predicted Tin)
[Out] = lsim(sys,[Ta],tm,[12.5; T4_0; Tin_0]);  
    
%Model Output Temperatures
Tin_Model_Off = Out(:,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Model_On
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rea = Rext + R1;
Rie = Rint + R3;

a = [-(1/(C1*Rea) + 1/(C1*R2))   1/(C1*R2)                      0;
       1/(C2*R2)             -(1/(C2*R2) + 1/(C2*Rie))      1/(C2*Rie) 
         0                       1/(Cz*Rie)   -(1/(Cz*Rie) + 1/(RwOn*Cz))]; 

b = [1/(C1*Rea); 0; 1/(RwOn*Cz)];

c = [Rext/(Rext + R1)     0                             0;
       1                  0                             0;
       0                  1                             0;
       0                 Rint/(R3 + Rint)  R3/(R3 + Rint);
       0                  0                            1];
 
d = [R1/(Rext + R1);0 ;0 ;0 ;0];
    
%Create the state-space model.
sys = ss(a,b,c,d);

%Simulate model (obtain the Predicted Tin)
[Out] = lsim(sys,[Ta],tm,[12.5; T4_0; Tin_0]);  
    
%Model Output Temperatures
Tin_Model_On = Out(:,5);



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
plot(sDateTime,Tin_Model_Off,'b-');
plot(sDateTime,Tin_Model_On,'r:');

xlabel('Time (hours)');
ylabel('Temperature ^\circC');
ylim([10 15.5]); 
datetick('x','hh');
title(datestr(startInt,'mm/dd'));
legend('T_a','T_{in}\_Model (ws=0)', 'T_{in}\_Model (ws=1)');
movegui('southwest');

%Save Figure to Files (eps, pdf)
%mkdir('Exp4');
makeFiles('Exp4/Exp4b');

 