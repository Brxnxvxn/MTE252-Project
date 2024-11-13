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
end

signal_filter('Test 9.m4a', 'output.wav')