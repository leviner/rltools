classdef cTransceiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         sweeprate
%         slope
        
        trdimpedance = 75;
        trdmeasfname
        
%         filter_dir
    end
    
    properties(SetAccess='private')
        trdmeas

        filter_g;
        filter_w;
        filter_b1
        FFTb1
        fb1
        decimation;
        fsdec;
        
        txpower
        txamp
        
        sweeprate
        slope
        pulselength
        fstart
        fstop
        fc
        bw
        
        ytx
        FFTytx
        fytx
        
        ytxcomp
        FFTytxcomp
        fytxcomp
        
        ytxcompauto
        ytxcompautogate
    end

    
    methods
        function obj = init(obj,configtransducer,filterdata,sampledata)
            hv          = Constants.wbt.hvvoltage;
            transf      = Constants.wbt.transformer;
            rwbttx      = Constants.wbt.tximpedance;
            fs          = Constants.wbt.fs;

            % Transducer
            idx = strfind(configtransducer.channelid,'ES');
            trdname = deblank(configtransducer.channelid(idx:end));

            switch trdname
                case 'ES70-7C'
                    obj.trdmeasfname = 'ES70-7C_dummy.txt';
                case 'ES120-7C'
                    obj.trdmeasfname = 'ES120-7CD_106_5.txt';
                case 'ES200-7C'
                    obj.trdmeasfname = 'ES200-7CD_120_5.txt';
                case 'ES333-7C'
                    obj.trdmeasfname = 'ES333_7CD_110.txt';
            end
            
            obj.trdmeas = load(obj.trdmeasfname);
            
            % Rx filter
%             fname = strcat(obj.filter_dir,'rxfilter.txt');
%             [g,w] = loadfilter(fname);
             [g,w] = loadfilterstr(filterdata.filterdata);

            obj.filter_g = g;
            obj.filter_w = w;
            obj.filter_b1 = [conj(flipud(w));w(2:end)]./g;
            
            % nfft = length(obj.b1);
            nfft = 1024;
            obj.FFTb1   = abs(fft(obj.filter_b1,nfft))/2;
            obj.fb1  = fs*linspace(0,1,nfft)';
            
            obj.decimation  = (length(obj.filter_b1)+1)/4;
            obj.fsdec       = Constants.wbt.fs/obj.decimation;
            
            % TX amplitude

            ztrd        = obj.trdimpedance;
            obj.txpower = sampledata.transmitpower;
            
            obj.txamp   = sqrt((obj.txpower/4)*(2*ztrd))* (rwbttx+ztrd)/(hv*transf*ztrd);
            
            % Sweep parameters
            obj.sweeprate   = sampledata.sweep;
            obj.slope       = round(configtransducer.slope*10)/10;
            obj.pulselength = sampledata.pulselength;
            obj.fstart      = sampledata.fstart;
            
            obj.fstop   = obj.fstart + obj.sweeprate*obj.pulselength;
            obj.fc      = (obj.fstart + obj.fstop)/2;
            obj.bw      = obj.fstop-obj.fstart;
            
            % Transmit signal

            ytxa    = gentx(obj.fstart,obj.fstop,obj.pulselength,obj.slope)';
            obj.ytx = ytxa*hv*transf*obj.txamp;
            
            nfft        = length(obj.ytx);
            obj.FFTytx  = 2*abs(fft(obj.ytx,nfft)) * sqrt(obj.bw/obj.pulselength) / fs;
            obj.fytx    = fs*linspace(0,1,nfft)';

            % Replicate signal for pulse compression
            
            ytxfil      = conv(obj.ytx,obj.filter_b1);
            ytxcomptmp  = downsample(ytxfil,obj.decimation);
            obj.ytxcomp = ytxcomptmp;

            nfft           = sampledata.count;
            FFTytxcomptmp  = fft(obj.ytxcomp,nfft);
            [fvec, FFTvec ] = freqtransf( FFTytxcomptmp,obj.fsdec,obj.fc );

            obj.FFTytxcomp  = FFTvec;
            obj.fytxcomp    = fvec;
            
            % Gated auto correlation function for replicate signal
            
            obj.ytxcompauto = conv(obj.ytxcomp,flipud(conj(obj.ytxcomp)));
%             figure
%             plot(20*log10(abs(obj.ytxcompauto)))
%             pause
            
            if (obj.slope==0.1)
                obj.ytxcompautogate = obj.ytxcompauto(497:511);
            elseif (obj.slope==0.5)
                obj.ytxcompautogate = obj.ytxcompauto(472:536);
            end
            
        end

    end
    
end

