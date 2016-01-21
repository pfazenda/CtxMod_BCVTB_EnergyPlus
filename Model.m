%+-----------------------------------------------------------------------
%|  _____  _     _____    _______ _               _     
%| |  __ \| |   |  __ \  |__   __| |             (_)    
%| | |__) | |__ | |  | |    | |  | |__   ___  ___ _ ___ 
%| |  ___/| '_ \| |  | |    | |  | '_ \ / _ \/ __| / __|
%| | |    | | | | |__| |    | |  | | | |  __/\__ \ \__ \
%| |_|    |_| |_|_____/     |_|  |_| |_|\___||___/_|___/
%|                                                                         
%|  __  __ _____ _______         _____ _______ 
%| |  \/  |_   _|__   __|       |  __ \__   __|
%| | \  / | | |    | |  ______  | |__) | | |   
%| | |\/| | | |    | | |______| |  ___/  | |   
%| | |  | |_| |_   | |          | |      | |   
%| |_|  |_|_____|  |_|          |_|      |_|   
%|
%|
%|                      Building Simulation
%+-----------------------------------------------------------------------
%|                 Sustainable Energy Systems 
%|
%|          Energy Efficiency Monitoring and Management 
%|               to Promote Sustainable Behaviors
%+-----------------------------------------------------------------------
%| Instituto Superior Técnico  
%| (http://www.ist.utl.pt/)
%|
%| Universidade Técnica de Lisboa  
%| (http://www.utl.pt/)
%| 
%| Institute for Systems and Robotics
%| (http://welcome.isr.ist.utl.pt/home/)
%|
%| MIT-Portugal Program
%| (http://www.mitportugal.org) 
%|
%| The ALFA group: Anyscale Learning for All (ALFA)
%| Computer Science and Artificial Intelligence Laboratory (CSAIL)
%| Massachusetts Institute of Technology (MIT)
%+-----------------------------------------------------------------------
%|          Solar Gains With Windows HSyst (Example 2)
%+-----------------------------------------------------------------------

%Close all windows and clear all workspace variables.
close all; clear;

%Energy Plus Output File Path
filePath = 'Output/Simple.csv';

%Open a file with comma-separated values. (date/time on the first column).
fid = fopen(filePath, 'rt');
 
C= textscan(fid, ['%s' repmat('%f',1,14)],'Delimiter',',','HeaderLines',1);
fclose(fid);
 
% convert date/time to serial date number
dt = datenum(C{1}(:,1), 'mm/dd HH:MM:SS');
 
nbrDays = round(dt(end) - dt(1));
  
%Combine all in one matrix
M = [dt C{2} C{3}  C{4}  C{5}  C{6}  C{7}  C{8}... 
        C{9} C{10} C{11} C{12} C{13} C{14} C{15}];

startInt = datenum('01/01 00:00:00','mm/dd HH:MM:SS');
endInt   = datenum('01/03 00:00:00','mm/dd HH:MM:SS');
 
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

%Maximum time (seconds);099
tMax = size(Int,1)*Ts;
 
%Time scale (seconds)
%ts = 0:(tMax - 1);
 
%Sampling points 
tm = 0:Ts:(tMax - 1);
 
%Outdoor Ambient Temperature [C]
Ta  = Int(:,2);

%Zone Windows Total Transmitted Solar Radiation Rate [W] 
Wind_SolarHGain = Int(:,3);

%Walls - Inside Face Temperature [C]1.0
T4_EP = Int(:,4);

%Walls - Outside Face Temperature [C]
T1_EP = Int(:,5);

%Surfaces Outside Face Solar Heat Gain per Area[W/m2]
S1_SolarHGain  = Int(:,6);
S2_SolarHGain  = Int(:,7);
S3_SolarHGain  = Int(:,8);
S4_SolarHGain  = Int(:,9);
S6_SolarHGain  = Int(:,11);

%Surface5 (Ground) Outside Face Temperature [C]
SOut_EP  = Int(:,10);

%Zone Temperature [C]
Tin_EP = Int(:,12);

%AFN Linkage Node 1 to Node 2 Volume Flow Rate [kg/s]
AirFlow = Int(:,13);

%Unit Heater Outlet Node Temperature [C]
TunitHeaterOut = Int(:,14);

%Unit Heater Outlet Mass Flow Rate [kg/s]
TunitHeaterOutFRate = Int(:,15); 

%Initial temperatures [C]
Tin_0 = Tin_EP(1);
T1_0  = T1_EP(1);
T4_0  = T4_EP(2);

nbrSamples = size(Ta,1);

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
   
%Plaster*2200O
l3   =  0.019;              % (m)
rho3 =   1602;              % (kg/m3)
cp3  =    837;              % (J/kg-K)
K3   =  0.726;              % (W/m-K) 

%Concrete
lfloor   =     0.01;        % (m)
rhofloor = 2242.600;        % (kg/m3)
cpfloor  =      837;        % (J/kg-K)
Kfloor   =    1.723;        % (W/m-K) 

%Glass
lglass =  0.003;            % (m)
Kglass =  0.8;              % (W/m-K) 

%Area of Windows
Awind = 4*1.5*4;            % (m^2)
%Area of the B. Envelope
A = 12*3*4 + 12*12 - Awind; % (m^2) 

Cair = 1005;                %J/kgK
densAir = 1.205;            %kg/m³ 
Cz   = 432*densAir*Cair;

%+------------------------------------------------------------------------)
%|                             Walls and Roof
%+------------------------------------------------------------------------)
hext = 18.14;               % (W/m²K) 
hint =  3.076;              % (W/m²K) 

Rext = 1/(hext*A);   
Rint = 1/(hint*A);

R = l1/(K1*A) + l2/(K2*A) + l3/(K3*A) + Rext + Rint;
C = (rho1*cp1*l1 + rho2*cp2*l2 +rho3*cp3*l3)*A;

%+-----------------------------------------------)
%| Calculating the 3R2C wall model parameters    ) 
%+-----------------------------------------------)

alpha1 = 0.23818828858317055376954213751263;
alpha2 = 0.16209971292445897672226108436541;
alpha3 = 0.59971199849237046950819677812195;

beta1 = 0.51685417244656046687598306989775;
beta2 = 0.48314582755343953312401693010225;

R1 = R*alpha1 - Rext;
R2 = R*alpha2;
R3 = R*alpha3 - Rint;

a1 = 0.72;          
a2 = 0.67;

C1 = C*beta1*a1;
C2 = C*beta1*(1-a1);
C3 = C*beta2*(1-a2);
C4 = C*beta2*a2;

Rea = Rext + R1;
Rie = Rint + R3;

%+-----------------------------------------------)
%|               Windows and Shades              ) 
%+-----------------------------------------------)
UValue = 5.88;
ws = 0; %Shades Off
ws_history = zeros(nbrSamples,1);

Rshade = 0.02*ws/(Awind*0.05);
Rw = 1/(Awind*UValue);
RW = Rw + Rshade;
%+----------------------)
%| Floor
%+----------------------)
hgint = 1.00;

Afloor = 12*12;
Rg = lfloor/(Kfloor*Afloor);
Cg = (rhofloor*cpfloor*lfloor)*Afloor;
Rgint = 1/(hgint*Afloor);

Rg1 = Rg/2;
Rg2 = Rg/2;
Rgie = Rg2 + Rgint;

%+------------------------------------------------------------------------)
%|                           Internal Gains
%+------------------------------------------------------------------------)
%Heater
Rih = 1/(1.1942055522671*Cair);
Qh  = 8000 + 75;    %Heater + Fan

%+----------------------)
%|  Calculate Ch        )  
%-----------------------)
%Heating Ramp
ind = getTimeIndex(startInt, '01/02 09:30:00', Ts);
tempDif = TunitHeaterOut(ind) - Tin_EP(ind);
varTin = Tin_EP(ind + 1) - Tin_EP(ind);
aux = Qh - (tempDif/Rih);
Ch  = (aux/varTin)*Ts;

sgh = 0;
h   = 1;

%+-----------------)
%|   Occupancy     ) 
%+-----------------)
Qi = 132;
o = 0;

%+------------------------------------------------------------------------)
%|                          External Gains 
%+------------------------------------------------------------------------)

%Solar Gains (Building Envelope)
AvgSolarHGain = (S1_SolarHGain + ...
                 S2_SolarHGain + ...
                 S3_SolarHGain + ...
                 S4_SolarHGain + ...
                 S6_SolarHGain)./5; 

             
%Solar Gains (Through Windows)
[b,a] = butter(6,0.0099); 
fWind_SolarHGain = filter(b,a,Wind_SolarHGain);
rwww = zeros(nbrSamples,1);

%+------------------------------------------------------------------------)
%|                              Model 
%+------------------------------------------------------------------------)

%Starte Vars: [T1 T2 T3 T4 T5 Th Tin]

T1_Model  = zeros(nbrSamples,1);
T2_Model  = zeros(nbrSamples,1);
T3_Model  = zeros(nbrSamples,1);
T4_Model  = zeros(nbrSamples,1);
T5_Model  = zeros(nbrSamples,1);
Th_Model  = zeros(nbrSamples,1);
Tin_Model = zeros(nbrSamples,1);
 
%Initial Values
T1_Model(1)  = 13.012;
T2_Model(1)  = 14.6;
T3_Model(1)  = 15.32;
T4_Model(1)  = 15.32;
T5_Model(1)  = 15.41;
Th_Model(1)  = 15.00;
Tin_Model(1) = 15.3620;

%Occupancy Schedule
ti1 = getTimeIndex(startInt, '01/01 08:00:00', Ts);
ti2 = getTimeIndex(startInt, '01/01 13:00:00', Ts);
ti3 = getTimeIndex(startInt, '01/02 08:00:00', Ts);
ti4 = getTimeIndex(startInt, '01/02 13:00:00', Ts);

t01 = getTimeIndex(startInt, '01/01 12:00:00', Ts);
t02 = getTimeIndex(startInt, '01/01 18:00:00', Ts);
t03 = getTimeIndex(startInt, '01/02 12:00:00', Ts);
t04 = getTimeIndex(startInt, '01/02 18:00:00', Ts);

%Heater Schedule
tHeaterOn  = getTimeIndex(startInt, '01/02 9:00:00', Ts);
tHeaterOff = getTimeIndex(startInt, '01/02 18:00:00', Ts);

%Shades Schedule
sOpen   = getTimeIndex(startInt, '01/01 00:00:00', Ts);
sClosed = getTimeIndex(startInt, '01/02 23:59:59', Ts);

%Window Schedule
wOpen   = getTimeIndex(startInt, '01/01 15:00:00', Ts); 
wClosed = getTimeIndex(startInt, '01/01 17:00:00', Ts);


%Discrete States 
vacant       = 0;
occupied     = 1;
shadesOpened = 3;
shadesClosed = 2;
preHeat      = 4;
cTemp        = 5;
nVent        = 6;

%Current State
state = vacant;

for i=2:nbrSamples,

   Qsw = fWind_SolarHGain(i);

   %+----------------------------)
   %| Current State
   %+----------------------------)
   T1  = T1_Model(i - 1); 
   T2  = T2_Model(i - 1); 
   T3  = T3_Model(i - 1);
   T4  = T4_Model(i - 1);
   T5  = T5_Model(i - 1);
   Th  = Th_Model(i - 1);
   Tin = Tin_Model(i - 1);
   
   %+---------------------------------------------------------------------)
   %| Occupancy Out -> Working(8-12) -> Lunch -> Working(13-18) -> Out
   %+---------------------------------------------------------------------)
   if (i > ti1 && i < t01) || (i > ti2 && i < t02) ||...
      (i > ti3 && i < t03) || (i > ti4 && i < t04),
     o = 5;
   else 
     o = 0;
   end

   %+----------------------)
   %| Heater Operation
   %+----------------------)
   if (i > tHeaterOn && i <= tHeaterOff), 
     sgh = 1;
   else 
     sgh = 0;
   end
  
   %+----------------------)
   %| Shades Operation
   %+----------------------)
    if (i > sOpen && i < sClosed),     
        if (AvgSolarHGain(i) >= 150), 
            ws = 1;       %Shades On
        elseif (AvgSolarHGain(i) < 130),    
            ws = 0;       %Shades Off    
        end
    end
   ws_history(i) = ws;
   Rshade = 0.02*ws/(Awind*0.05);
   RW = Rw + Rshade;
  
   
   %+----------------------)
   %| Windows Operation
   %+----------------------)
   if (i > wOpen && i < wClosed), 
      mv =  AirFlow(i);
      
     %Smooth AirFlow and remove zeros 
     if mv == 0, mv = mean(AirFlow(find(AirFlow >1))); end
      mv = round(mv*100)/100;       %Discretize mv (100 different levels).
      if (i == wOpen + 2), mv = mv/2; end  %Smooth out transition.
      RW =  0.76/(mv*Cair);
      rwww(i) = mv;
   end 
   
   
   
   %+---------------------------------------------------------------------)
   %| State Transitions
   %+---------------------------------------------------------------------)
   switch state
     case vacant
       if (o > 0)
         i  
         state = occupied
       end
       
     case occupied
       if (AvgSolarHGain(i) >= 150), 
         i
         state = shadesClosed
       elseif (AvgSolarHGain(i) < 130),    
         i 
         state = shadesOpened
       end
       
     case shadesOpened
       if(sgh == 1),
         i
         state = preHeat
       elseif (i > wOpen && i < wClosed),
         i      
         state = nVent
       elseif (AvgSolarHGain(i) >= 150)
         i
         state = occupied
       elseif (o == 0)
         i
         state = vacant
       end
       
     case shadesClosed

       if(sgh == 1),
         i
         state = preHeat
       elseif (AvgSolarHGain(i) < 130)
         i
         state = occupied
       elseif (o == 0)
         i
         state = vacant
       end

     case preHeat     
       if Tin >= 22,
          i
          state = cTemp
       end
       
     case cTemp
       if (sgh == 0),
          i 
          if (AvgSolarHGain(i) >= 150),
            state = shadesClosed
          else
            state = shadesOpened  
          end
       end 
       
     case nVent       
         if (i <= wOpen || i >= wClosed),
            i 
            state = shadesOpened
         end
   end    
   
   
   
   %Inc
   dT1  =  -T1*(1/(C1*Rext) + 1/(C1*R1)) + T2/(C1*R1) +  Ta(i)/(C1*Rext)...
           + (AvgSolarHGain(i)*A)/C1;
   dT2  =  T1/(C2*R1) -T2*(1/(C2*R1) + 1/(C2*R2)) + T3/(C2*R2);  
   dT3  =  T2/(C3*R2) -T3*(1/(C3*R2) + 1/(C3*R3)) + T4/(C3*R3); 
   dT4  =  T3/(R3*C4) -T4*(1/(C4*R3) + 1/(C4*Rint)) + Tin/(Rint*C4); 
   dT5  =  -T5*(1/(Cg*Rg1) + 1/(Cg*Rgie)) + Tin/(Cg*Rgie)...
           + SOut_EP(i)/(Cg*Rg1);
   
   if Tin >= 22 && sgh == 1,
     dTh  =  0;
     dTin =  0;
   else
     dTh  =  -(sgh*Th)/(Ch*Rih) + (sgh*Tin)/(Ch*Rih) + (Qh*h*sgh)/Ch;                        
     dTin =  T4/(Cz*Rint) + Ta(i)/(Cz*RW) + T5/(Cz*Rgie)...
             + (Th*sgh)/(Cz*Rih)...
             -(Tin/Cz)*(1/Rgie + 1/Rint + 1/RW + sgh/Rih)...
             + (o*Qi)/Cz + (Qsw*0.18)/Cz; 
   end;
   
   %Update    
   T1_Model(i)  = T1  + dT1*Ts; 
   T2_Model(i)  = T2  + dT2*Ts; 
   T3_Model(i)  = T3  + dT3*Ts;
   T4_Model(i)  = T4  + dT4*Ts; 
   T5_Model(i)  = T5  + dT5*Ts;
   Th_Model(i)  = Th  + dTh*Ts; 
   Tin_Model(i) = Tin + dTin*Ts;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%|>                        Show Plots                                    <|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+
%+---------------------------------------+
%|>  Outdoor and Zone Temperatures      <|
%+---------------------------------------+
figure(1);
hold on;
box on;   

%Plot outdoor temperature (hourly time scale)
plot(sDateTime,Ta,'k');
    
%Energy Plus
plot(sDateTime,Tin_EP,'b:','LineWidth',1.5);
   
%Plot zone temperature (Model)
plot(Int(:,1),Tin_Model,'r--','LineWidth',1.5);
   
xlabel('Time (hours)');
ylabel('Temperature ^\circC');
datetick('x','hh');
title(datestr(startInt,'mm/dd'));
%legend('T_a','T_{in}\_EP', 'T_{in}\_Model','Location','NorthWest');
movegui('southwest');
  
%Save Figure to Files (eps, pdf)
mkdir('Example2');
makeFiles('Example2/Example2');

%+----------------------------------+
%|>  Model Error and Histogram     <|
%+----------------------------------+
figure(2);
box on;
     
error = (Tin_Model(:,1) - Tin_EP);
  
RSS = sum(abs(error))

MAE = mae(error);
  
[MAX,I] = max(abs(error))
when = I*Ts;
  
when= when/3600 
if when >= 24,
   when = when -24
end
  
hist(error)
set(get(gca,'child'),'FaceColor','none','EdgeColor','k');
xlabel('Error in ^\circC');
ylabel('Frequency');
%title(datestr(startInt,'mm/dd'));
movegui('northwest');
makeFiles('Example2/Error2');
 
% %+----------------------------------------------------+
% %|>  Surface Outside Face Solar Heat Gains           <|
% %+----------------------------------------------------+
figure(3);
hold on;
box on;

%Surfaces Outside Face Solar Heat Gain per Area[W/m2]
%plot(sDateTime,S1_SolarHGain,'b:');
%plot(sDateTime,S2_SolarHGain,'b:');
%plot(sDateTime,S3_SolarHGain,'b:');
%plot(sDateTime,S4_SolarHGain,'b:');
%plot(sDateTime,S6_SolarHGain,'b:');
 
plot(sDateTime,AvgSolarHGain,'k');
  
xlabel('Time (hours)');
ylabel('Solar Heat Gain (Wm^{-2})');
datetick('x','hh');
title(datestr(startInt,'mm/dd'));
movegui('west');
  
%Save Figure to Files (eps, pdf)
makeFiles('Example2/Solar');

%+---------------------------------------------------------+
%|>  Zone Windows Total Transmitted Solar Radiation Rate  <|
%+---------------------------------------------------------+
figure(4)
hold on;
box on;

plot(sDateTime,Wind_SolarHGain,'k');
plot(sDateTime,fWind_SolarHGain,'b--');
 
xlabel('Time (hours)');
ylabel('Windows Solar Radiation (W)');
datetick('x','hh');
%title(datestr(startInt,'mm/dd'));
%legend('Qsw','fsw');
movegui('west');

%Save Figure to Files (eps, pdf)
makeFiles('Example2/SolarWind');

 
 
% %+---------------------------------------------------------+
% %|>  Window Shades                                        <|
% %+---------------------------------------------------------+
% figure(5)
%   
% plot(sDateTime,ws_history,'b:');
% xlabel('Time (hours)');
% ylabel('ws');
% datetick('x','hh');
% title(datestr(startInt,'mm/dd'));
% movegui('west');
% makeFiles('Example2/WindShades');
% 
% 
%  
% %+---------------------------------------------------------+
% %|> AFN Linkage Node 1 to Node 2 Volume Flow Rate AirFlow <|
% %+---------------------------------------------------------+
figure(6)
box on;


plot(sDateTime,AirFlow,'k');
  
xlabel('Time (hours)');
ylabel('Windows Volume Flow Rate AirFlow (kg/s)');
datetick('x','hh');
%title(datestr(startInt,'mm/dd'));
movegui('west');
makeFiles('Example2/AirFlow');
