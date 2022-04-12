function data = EK80pulseCompress(data, channel)
%
% Pulse compress the signal for FM.
% Estimate the effective pulse length for use later one
%
%   revised for CTRiver - perkins
%   revised to be more generic - Newhall

disp('Doing pulse compression...')

settingsPing = 1; % ping to get pulse parameters from...

ytx = EK80calcSentSignal(data, channel, settingsPing);

for ping = 1:size(data.echodata, 2) % loop over number of pings
    if ~isempty(data.echodata(channel,ping).channelID)
        
        % check for FM sweep
        if (data.param(channel, ping).FrequencyEnd - data.param(channel, ping).FrequencyStart) > 1
            
            
            % Do pulse compression for each quadrant
% % %             for q = 1:4
% % %                 tmp = data.echodata(channel, ping).complexsamples(:,q);
% % %                 % The PC filter and decimation (optional)
% % %                 %tmp = conv(tmp, data.filters(channel, 2).FilterData);
% % %                 %tmp = downsample(tmp, data.filters(channel, 2).Decimation);
% % %                 
% % %                 data.echodata(channel, ping).compressed(:, q) = ...
% % %                     conv(tmp, flipud(conj(ytx))) / norm(ytx)^2;
% % %                 
% % %             end


            % changed Newhall
            % check if 1 row or multiple rows, this is necessary to get correct
            %     number of channels for each transducer
            
            % only 1 channel
            if (isrow(data.echodata(channel).complexsamples))
                
               tmp = data.echodata(channel, ping).complexsamples; 
               data.echodata(channel, ping).compressed(:, 1) = ...
                    conv(tmp, flipud(conj(ytx)),'same') / norm(ytx)^2; 
                
            else   % multiple channels
                 
                for q = 1:size(data.echodata(channel, ping).complexsamples,2)
                    
                    tmp = data.echodata(channel, ping).complexsamples(:,q);
                    % The PC filter and decimation (optional)
                    %tmp = conv(tmp, data.filters(channel, 2).FilterData);
                    %tmp = downsample(tmp, data.filters(channel, 2).Decimation);
                    
                    %data.echodata(channel, ping).compressed(:, q) = ...
                    %    conv(tmp, flipud(conj(ytx))) / norm(ytx)^2;
                    
                     data.echodata(channel, ping).compressed(:, q) = ...
                        conv(tmp, flipud(conj(ytx)),'same') / norm(ytx)^2;                   
                end      
                
            end
              
            
            %             data.echodata(channel, ping).compressed_amp=sum(data.echodata(channel, ping).compressed,2);
            ptx(1:size(data.echodata(channel, ping).compressed,1), ping) = abs(sum(data.echodata(channel, ping).compressed,2));
            % Store the effective pulse length for use later on
            % ytxa = conv(ytx, conj(flipud(ytx))) / (norm(ytx) ^ 2);
            ytxa = conv(ytx, conj(flipud(ytx)),'same') / (norm(ytx) ^ 2);
            ptxa = abs(ytxa) .^ 2;
            fs_dec = 1/data.param(channel, ping).SampleInterval;
            data.param(channel, ping).tauEffective = ...
                sum(ptxa) / (max(ptxa) * fs_dec);
            sampleDataField = 'compressed';
            %             cp_ps(1:size(data.echodata(channel, ping).compressed,1),ping)=abs(sum(data.echodata(channel, ping).compressed,2));
            %             pings(ping)=ping;
            
        else     % CW
            
            % Store the effective pulse length for use later on
            ptxa = abs(ytx) .^ 2;
            fs_dec = 1/data.param(channel, ping).SampleInterval;
            data.param(channel, ping).tauEffective = ...
                sum(ptxa) / (max(ptxa) * fs_dec);
            sampleDataField = 'complexsamples';
        end
        
        % Range vector in metres, compensating for transmit pulse length
        if isfield(data.param,'PulseLength')
            data.echodata(channel, ping).range = single(data.echodata(channel, ping).minRange) + ...
                single((data.echodata(channel, ping).offset:size(data.echodata(channel, ping).(sampleDataField), 1)-1)' * ...
                data.environ.SoundSpeed * data.param(channel, ping).SampleInterval/2 - ...
                data.param(channel, ping).PulseLength * data.environ.SoundSpeed/2);
        elseif isfield(data.param,'PulseDuration')
            data.echodata(channel, ping).range = single(data.echodata(channel, ping).minRange) + ...
                single((data.echodata(channel, ping).offset:size(data.echodata(channel, ping).(sampleDataField), 1)-1)' * ...
                data.environ.SoundSpeed * data.param(channel, ping).SampleInterval/2 - ...
                data.param(channel, ping).PulseDuration * data.environ.SoundSpeed/2);
        end
        
    end   % if isempty
end  % end ping loop

 

end             % end function
