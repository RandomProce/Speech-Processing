clc;
clear all ;
close all;
% Read the signal
%Input Parameters
% sig - Noisy Signal
%babble
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_0dB\0dB\sp03_babble_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_5dB\5dB\sp03_babble_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_10dB\10dB\sp03_babble_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_15dB\15dB\sp03_babble_sn15.wav');
%car
[sig fs1] = audioread('D:\speech project\Noizeus data base\car noise\car_0dB\0dB\sp03_car_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\car noise\car_5dB\5dB\sp03_car_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\car noise\car_10dB\10dB\sp03_car_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\car noise\car_15dB\15dB\sp03_car_sn15.wav');
%station
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Station noise\station_0dB\0dB\sp03_station_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Station noise\station_5dB\5dB\sp03_station_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Station noise\station_10dB\10dB\sp03_station_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Station noise\station_15dB\15dB\sp03_station_sn15.wav');
%Street
% [sig fs1] = audioread('D:\speech project\Noizeus data base\street noise\street_0dB\0dB\sp03_street_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\street noise\street_5dB\5dB\sp03_street_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\street noise\street_10dB\10dB\sp03_street_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\street noise\street_15dB\15dB\sp03_street_sn15.wav');
%Airport
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Airport Noise\airport_0dB\0dB\sp03_airport_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Airport Noise\airport_5dB\5dB\sp03_airport_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Airport Noise\airport_10dB\10dB\sp03_airport_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Airport Noise\airport_15dB\15dB\sp03_airport_sn15.wav');
noisy_signal = sig;
L = length(sig);%Segment Length 128,  50% of overlapping
fr_length = 256;                %Frame Length = 256 Samples
seg_length = 128;               %Segment Length = 128 means that 50% of overlapping
x = L-fr_length;
num_frames = floor(x/seg_length)+1; %Number of Frames calculations
w = 1;
alpha_s = 0.8;
alpha_p = 0.2;
alpha_d = 0.95;
gamma0 = 4.6;
epsalan0 = 1.67;            %Initializations
gamma1 = 3;
% y = num_frames - 96;
b =  hanning(2*w+1);
b = b'./sum(b);
U = 8;
V = 96;
D = 96;
del = 5;
temp = 0;
k = 0;
for i = 1:num_frames
    frame_noisy_signal1(:,i) = (noisy_signal(((0.5*temp*fr_length)+1):(0.5*temp*fr_length)+fr_length));     %Dividing the Noisy signals in to frames with overlapping of 50%
    temp=temp+1;
    frame_noisy_signal(:,i) = hamming(fr_length) .* frame_noisy_signal1(:,i);          %Apply windowing
    power_of_signal(:,i) = ((abs(fft(frame_noisy_signal(:,i),fr_length))).^2);         %FFT of Noisy signal in each frames
%     x1 = (conv(b,power_of_signal(:,i)));                     %Convolving window length of 2*seg_length+1 and FFT of noisy signal
%     Sf(:,i) = x1(2*w+1:end);
    Sf(:,i) = power_of_signal(:,i);
    if i==1
        S(:,i) = alpha_s.*Sf(:,i)+(1-alpha_s).*Sf(:,i);
    else
        S(:,i) = alpha_s.*S(:,i-1)+(1-alpha_s).*Sf(:,i);
    end
        if i==1
            Smin(:,i) = S(:,i);
            Stmp(:,i) = S(:,i);
        else
            if mod(i,D)==0
                Stmp(:,i) = S(:,i);
                Smin(:,i) = min(Stmp(:,i-1),S(:,i));
            else
                Smin(:,i) = min(Smin(:,i-1),S(:,i));
                Stmp(:,i) = min(Stmp(:,i-1),S(:,i));
            end
        end
        for j = 1:fr_length
            Sr(j,i) = S(j,i)./Smin(j,i);
        if Sr(j,i) > del
            I(j,i) = 1;
        else
            I(j,i) = 0;
        end
        
        if i ==1
            p(j,i) = (1-alpha_p)*I(j,i);
        else
            p(j,i) = alpha_p*p(j,i-1)+(1-alpha_p)*I(j,i);
        end
        
        alphad(j,i) = alpha_d+(1-alpha_d)*p(j,i);
        
        if I(j,i) ==0
            if i==1
                lambda_d(j,i) = alphad(j,i)*Sf(j,i)+(1-alphad(j,i))*power_of_signal(j,i);
            else
                lambda_d(j,i) = alphad(j,i)*lambda_d(j,i-1) + (1-alphad(j,i))*power_of_signal(j,i);
            end
        else
            lambda_d(j,i) = lambda_d(j,i-1);
        end
    end
end


freq_bin = 20;
semilogy(power_of_signal(freq_bin,:));
hold on
% semilogy(S(:,freq_bin),'r');
% hold on
semilogy(lambda_d(freq_bin,:),'k','LineWidth',2.5);
% legend({'Power of Signal','Smoothed Signal','Noise estimate'},'Location','best','FontSize',18);
% xlabel('Frame Index','FontSize',20,'FontName','Times');
% ylabel('Power in dB','FontSize',20,'FontName','Times');

seg_SNR = pow2db(mean(mean(power_of_signal))/mean(mean(Smin)))