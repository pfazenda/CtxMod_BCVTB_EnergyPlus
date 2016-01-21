function M = getIntervalData(dT_Start, dT_Stop, data)

M = data(data(:,1) >= dT_Start,:,:);
M = M(M(:,1) <= dT_Stop,:,:);
 