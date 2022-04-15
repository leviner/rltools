function data = addCals(data,cal)
data.calibration = struct();
for j=1:length(data.config.transceivers)
    
    calEntry = 0;
    
    for k = 1:length(cal)
        if string([data.config.transceivers(j).channels.transducer.SerialNumber])==string(cal(k).Calibration.Common.Transducer.SerialNumber)
            calEntry = k;
        end
    end
    
    if calEntry == 0
        fprintf('Matching serial number for %s not found in calibration\n',data.config.transceivers(j).channels.transducer.TransducerName);
        calFields = fieldnames(cal(1).Calibration.CalibrationResults); % Use the first cal to enter emptys
        for i= 1:length(calFields)
            [data.calibration(j).(calFields{i})] = [];
        end
        [data.calibration(j).CalbrationEnvironment] = [];
        
    else
        fprintf('Calibration found for %s, added to data \n',data.config.transceivers(j).channels.transducer.TransducerName);
        
      
        calFields = fieldnames(cal(calEntry).Calibration.CalibrationResults);
        
        for i= 1:length(calFields)
            hold = strsplit(cal(calEntry).Calibration.CalibrationResults.(calFields{i}),';');
            for ii=1:length(hold)
                hold_vector(ii) = round(str2num(cell2mat(hold(ii))),4);
            end
            [data.calibration(j).(calFields{i})] = hold_vector;
        end
        [data.calibration(j).CalbrationEnvironment] = cal(calEntry).Calibration.Common.EnvironmentData;
    end
end
