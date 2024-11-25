%{
Description: Function splits freq. range into individual channels and
             applies the bandpass filtering for each channel using the
             butter function.

Parameters: data                -signal data
            rate                -sampling rate
            N                   -# of channels
            offset_spacing      -spacing for overlap of channels
            order               -order of filter

Output:     filtered_channels   -filtered signal from each channel
            channel_coeffs      -the coeffients returned by butter function
                                 for each channel
%}
function [filtered_channels, channel_coeffs] = bandPassFilter(data, rate, N, offset_spacing, order)
    %frequency range
    min_f = 100;
    max_f = 8000;

    %split frequency range into intervals with logaritmic scale
    freq_log = logspace(log10(min_f), log10(max_f), N+1);
    disp(freq_log);

    filtered_channels = cell(N, 1); %array for filtered signals
    channel_coeffs = cell(N, 2); %array for coeffs. returned by butter()
    bands = zeros(N, 2); %stores the intervals for each channel

    %sets band array with intervals of each channel with overlap to
    %smoothen filtering
    for i = 1:N-1
        bands(i, 1) = freq_log(i) / (1 + offset_spacing);
        bands(i, 2) = freq_log(i+1) * (1 + offset_spacing);
    end

    %for last channel overlap is not applied for upper bound of interval
    bands(N, 1) = freq_log(N);
    bands(N, 2) = freq_log(N+1); 
    disp(bands);

    %normalize frequencies in band array since the butter function
    %only takes the frequency in that form
    norm_freqs = bands / (rate / 2);
        
    for i = 1:N
        cutoff = norm_freqs(i, :);
       
        if cutoff(1) == 0
            disp(cutoff);
            [b, a] = butter(order, cutoff(2), "low"); %low-pass filter for lowest channel
        elseif cutoff(2) >= 1
            disp(cutoff);
            [b, a] = butter(order, cutoff(1), "high"); %high-pass filter for highest channel
        else
            disp(cutoff);
            [b, a] = butter(order, cutoff, "bandpass"); %bandpass filter for all other channels
        end

        filtered_signal = filter(b, a, data); %filters signal using coeffs. from butter()

        %stores coeffs. in array for plotting (output)
        channel_coeffs{i, 1 } = b;
        channel_coeffs{i, 2} = a;
        
        %stores filtered signal for plotting (output)
        filtered_channels{i} = filtered_signal;
    end


end

function [] = signal_filter(input_file, output_file) 
    %read input signal and find sampling data and rate
    [data, rate] = audioread(input_file);
    new_rate = 16000;
    
    %This assumes that all the input signal are above 16KHz.
    %If they are below 16KHz, a new input signal will be chosen
    if(rate > 16000)
        data = resample(data, new_rate, rate);
    end
   
    sz = size(data);

    %if the input signal has 2-channels, convert it into a 1-channel array
    if(sz(2) == 2)
        new_data = mean(data, 2); %takes the average from both channels
    else
        new_data = data;
    end
    
    %play sound of sampling data at new sampling rate
    sound(new_data, new_rate);
    %write signal to output file
    audiowrite(output_file, new_data, new_rate);
    
    %find duration of signal using the data and sampling rate
    time = (0:length(new_data)-1) / new_rate;
    
    %plot sampled data over time
    plot(time, new_data);
    title("Signal Graph at 16KHz");
    xlabel("time(s)");
    ylabel("amplitude");

    %{
    f = 1000; %frequency of cosine signal
    duration = length(new_data) / new_rate; %duration for cosine signal
    t = 0:1/new_rate:duration; %domain and step paramater for time variable
    
    figure;
    y = cos(2 * pi * f .* t(1:end-1)); %cosine signal
    sound(y, new_rate); %sound output from cosine signal
    plot(t(1:2*f/new_rate*new_rate), y(1:2*f/new_rate*new_rate)); %plot for cosine signal at 1KHz for two cycles
    
    title("Cosine Graph of Signal at 1Khz");
    xlabel("time(s)");
    ylabel("Amplitude");
    xlim([0 0.002]);
    grid on;
    %}

    N = 20; %num of channels
    offset = 0.20; %the offset required for overlapping channels
    order = 4; %order of filter

    %filter sound with passband bank 
    [filtered_signals, channel_coeffs] = bandPassFilter(new_data, new_rate, N, offset, order);
    
    %plot of output signal at lowest frequency signal
    figure;
    plot(filtered_signals{1});
    title("Lowest Frequency Channel");
    xlabel("Samples");
    ylabel("Amplitude")

    %plot of output signal at highest frequency signal
    figure;
    plot(filtered_signals{N});
    title("Highest Frequency Channel");
    xlabel("Samples");
    ylabel("Amplitude")

    %plot of frequency response of filter at lowest channel
    [h, w] = freqz(channel_coeffs{1, 1},channel_coeffs{1, 2}); %frequency response
    figure;
    plot(w/pi*(new_rate/2), abs(h));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Frequency Response of Filter of Lowest Channel');

    %plot of frequency response of filter at highest channel
    [h, w] = freqz(channel_coeffs{1, 1},channel_coeffs{1, 2}); %frequency response
    figure;
    plot(w/pi*(new_rate/2), abs(h));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Frequency Response of Filter of highest Channel');

    %plot of amplitude vs freq of filtered signal
    samples = length(filtered_signals{N}); %# of samples in last channel
    frequencies = (0:samples - 1) * (new_rate / samples); %Vector of frequncies
    amp = abs(fft(filtered_signals{N})) / max(abs(filtered_signals{N})); % The FFT normalized using signal from last channel
    figure;
    plot(frequencies(1:samples/2), amp(1:samples/2)); % this only takes positive frequencies
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Amplitude vs Frequency of the Filtered Signal of Highest Channel');

    %envelope extraction step 1: rectify signals
    signals_rectified = cell(N, 1);
    for i = 1:N
        signals_rectified{i} = abs(filtered_signals{i});
    end

    %envelope extraction step 2: apply low-pass filter with 400Hz cutoff
    cutoff = 400 / (new_rate / 2);
    [b, a] = butter(order, cutoff, "low");
    envelopes = cell(N, 1);

    for i = 1:N
        envelopes{i} = filter(b, a, signals_rectified{i});
    end
    
    %plot of envelope at lowest channel
    figure;
    plot(envelopes{1});
    title("Envelope Lowest Frequency Channel");
    xlabel("Samples");
    ylabel("Amplitude")

    %plot of envelope at highest channel
    figure;
    plot(envelopes{N});
    title("Envelope Highest Frequency Channel");
    xlabel("Samples");
    ylabel("Amplitude")
end

signal_filter('Test 1.m4a', 'output.wav')