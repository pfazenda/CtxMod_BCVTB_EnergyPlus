function Index = getTimeIndex(dT_Start, date, Ts)

dateNum = datenum(date,'mm/dd HH:MM:SS');
Index = ceil(((dateNum - dT_Start)*86400)/Ts + 1);
