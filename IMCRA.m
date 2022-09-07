clc;
clear all;
close all;
freq_bin = 20;
%babble
[sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_0dB\0dB\sp03_babble_sn0.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_5dB\5dB\sp03_babble_sn5.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_10dB\10dB\sp03_babble_sn10.wav');
% [sig fs1] = audioread('D:\speech project\Noizeus data base\Babble Noise\babble_15dB\15dB\sp03_babble_sn15.wav');
%car
% [sig fs1] = audioread('D:\speech project\Noizeus data base\car noise\car_0dB\0dB\sp03_car_sn0.wav');
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
% [sig fs1] = audioread('D:\speech project\Noizeus data base\clean\clean\sp03.wav');
% sig = awgn(sig,5,'measured');
fs = 16e3;
y_in_time = resample(sig, fs, fs1);
L = size(y_in_time, 1);
frame_length = 256;
seg_length = 128;
num_frames = floor((L-frame_length)/seg_length)+1;
w = hanning(frame_length);
win = w./sum(w);
f_win_length = 1;
win_freq = hanning(2*f_win_length+1);  % window for frequency smoothing
win_freq = win_freq / sum(win_freq);  % normalize the window function
alpha_eta = 0.92;
alpha_s = 0.6;
alpha_d = 0.85;
beta = 1.47;
eta_min = 0.0158;
gama0 = 4.6;
gamma1 = 3;
zeta0 = 1.67;
Bmin = 1.66;
V = 15;
Nwin = 8;
temp = 0;
y_in_time1 = y_in_time;

for i = 1:num_frames
    frame_noisy_signal1(:,i) = (y_in_time1(((0.5*temp*frame_length)+1):(0.5*temp*frame_length)+frame_length));
    temp = temp+1;
    frame_noisy_signal(:,i) = frame_noisy_signal1(:,i).*win;
    power_of_signal(:,i) = ((abs(fft(frame_noisy_signal(:,i),frame_length))).^2);
    
end

%     Sf1 = conv(win_freq, power_of_signal(i,:));  % frequency smoothing
%     Sf(i,:) = Sf1(2*f_win_length+1:end)
for i = 1:num_frames
    for j = 1:frame_length
        if i ==1 || i==num_frames
            Sf(j,i) = power_of_signal(j,i);
        else
            Sf(j,i) = sum(win_freq'.*[power_of_signal(j,i-1) power_of_signal(j,i) power_of_signal(j,i+1)]);
        end
        if i==1
            S(j,i) = Sf(j,i);  % spec after time smoothing\
            Smin(j,i) = S(j,i);
            gama_min(j,i) = power_of_signal(j,i)/(Bmin*Smin(j,i));
            eta_min(j,i) = S(j,i)/(Bmin*Smin(j,i));
            if gama_min(j,i)<gama0 && eta_min(j,i)<zeta0
                I(j,i) = 1;
            else
                I(j,i) = 0;
            end
            St(j,i) = Sf(j,i);
            Smin(j,i) = Sf(j,i);
            
            lambda_dav(j,i) = power_of_signal(j,i);  % expected noise spectrum
            lambda_dc(j,i) = power_of_signal(j,i);  % modified expected noise spec
            gamma(j,i) = 1;
            Smint(j,i) = Sf(j,i);  % min value get from St
            Smin_sw(j,i) = Sf(j,i);  % auxiliary variable for finding min
            Smint_sw(j,i)= Sf(j,i);
            GH1(j,i) = 1;
            eta(j,i) = 0;
%             eta_2term(i,:)= GH1 .^ 2 .* max(gamma(i,:)-1,0);            
        else
            gamma(j,i) = power_of_signal(j,i)./(lambda_dc(j,i-1));  % update instant SNR
            eta(j,i) = alpha_eta .* GH1(j,i-1) + (1-alpha_eta) .* max(gamma(j,i)-1, 0);  % update smoothed SNR, eq.18, where eta_2term = GH1 .^ 2 .* gamma
            eta(j,i) = (eta(j,i));
            v(j,i) = (gamma(j,i) .* eta(j,i)) ./ (1+eta(j,i));
            GH1(j,i) = (eta(j,i) ./ (1+eta(j,i))).*exp(0.5*expint(v(j,i)));
            S(j,i) = alpha_s .* S(j,i-1) + (1-alpha_s) .* Sf(j,i);
            if mod(i,V)==0
                Smin_sw(j,i) = S(j,i);
                Smin(j,i) = min(Smin_sw(j,i-1),S(j,i));
            else
                Smin(j,i) = min(Smin(j,i-1), S(j,i));
                Smin_sw(j,i) = min(Smin_sw(j,i-1), S(j,i));
            end
        end
        gama_min(j,i) = power_of_signal(j,i) ./ (Bmin.*Smin(j,i));
        zeta(j,i) = S(j,i) ./ (Bmin.*Smin(j,i));
        if (gama_min(j,i)<gama0 && zeta(j,i)<zeta0)
            I_f(j,i) = 1;
        else
            I_f(j,i) = 0;
        end

        if j==1
            temp1(j,i) = I_f(j,i); %sum(I_f(j:j+f_win_length,i));
            temp2(j,i) = I_f(j,i)* Sf(j,i);
        elseif j==frame_length
            temp1(j,i) = I_f(j,i); %sum(I_f(j-f_win_length:j,i));
            temp2(j,i) = I_f(j,i)* Sf(j,i);
        else
            temp1(j,i)=sum(I_f(j-f_win_length:j,i));
            temp2(j,i) = sum(win_freq(1:end-1).*I_f(j-f_win_length:j,i).*power_of_signal(j-f_win_length:j,i));
        end
        
        if temp1(j,i)~=0
            Sf_t(j,i) = temp2(j,i)/temp1(j,i);
        else
            Sf_t(j,i) = St(j,i-1);
        end
        if i==1
            St(j,i) = Sf(j,1);
            Smin_sw(j,i) = Sf(j,i);
            Smint_sw(j,i)= Sf(j,i);
             lambda_dav(j,i) = power_of_signal(j,i);  % expected noise spectrum
            lambda_dc(j,i) = power_of_signal(j,i);  % modified expected noise spec
            gamma(j,i) = 1;
            Smin(j,i) = Sf(j,i);  % noise estimation spec value
            S(j,i) = Sf(j,i);  % spec after time smoothing
            Smin(j,i) = S(j,i);
            St(j,i) = Sf(j,i);  % Sft:smoothing results using speech abscent probability
            GH1(j,i) = 1;
            Smint(j,i) = Sf(j,i);  % min value get from St
            Smin_sw(j,i) = Sf(j,i);  % auxiliary variable for finding min
            Smint_sw(j,i)= Sf(j,i);
            eta(j,i) = 0;
        else
        St(j,i) = alpha_s.*St(j,i-1)+(1-alpha_s).*Sf_t(j,i);
        if mod(i,V)==0
            Smint_sw(j,i) = St(j,i);
            Smint(j,i) = min(Smint_sw(j,i-1),S(j,i));
        else
            Smint(j,i) = min(Smint(j,i-1), St(j,i));
            Smint_sw(j,i) = min(Smint_sw(j,i-1), St(j,i));
        end
        gamma_mint(j,i) = power_of_signal(j,i) ./ (Bmin .*Smint(j,i));
        zetat(j,i) = S(j,i) ./ (Bmin .*Smint(j,i));
%         lambda_d(j,i) = Smin(j,i);
        
        if (gamma_mint(j,i)>=1 & zetat(j,i)<zeta0)
            qhat(j,i) = 1;
            phat(j,i) = 0;
        elseif (1<gamma_mint(j,i)<gamma1 & zetat(j,i)<zeta0)
            qhat(j,i) = (gamma1-gamma_mint(j,i)) / (gamma1-1);
            phat(j,i) = 1 ./ (1+(qhat(j,i)./(1-qhat(j,i))).*(1+eta(j,i)).*exp(-v(j,i)));
        else
            qhat(j,i) = 0;
            phat(j,i) = 1;
        end
        alpha_dt(j,i) = alpha_d + (1-alpha_d) * phat(j,i);
        if I_f(j,i)==1
            lambda_dav(j,i) = alpha_dt(j,i) .* lambda_dav(j,i-1) + (1-alpha_dt(j,i)) .* power_of_signal(j,i);
        else
            lambda_dav(j,i) = lambda_dav(j,i-1);
        end
        lambda_dc(j,i) = lambda_dav(j,i).* beta;
    end
end
end
%         conv_I1 = conv(win_freq, I_f(i,:));  % smoothing
%         conv_I(i,:) = conv_I1(2*f_win_length+1:end);
%         conv_Y1 = conv(win_freq, I_f(i,:).*power_of_signal(i,:));  % eq. 26clc
%         conv_Y(i,:) = conv_Y1(2*f_win_length+1:end);
% %         Ii(i) = sum(I_f(i,:));
%         for j = 1:frame_length
%             if temp1(i,j)~=0
% %                 if Ii(i)~=0
%                 Sf_t(i,j) = conv_Y(i,j)/conv_I(i,j);
%             else
%                 Sf_t(i,j) = St(i-1,j);
%             end
%         St(i,j) = alpha_s.*St(i-1,j)+(1-alpha_s).*Sf_t(i,j);
%         end
%         if mod(i,V)==0
%             Smint_sw(i,:) = St(i,:);
%             Smint(i,:) = min(Smint_sw(i-1,:),S(i,:));
%         else
%             Smint(i,:) = min(Smint(i-1,:), St(i,:));
%             Smint_sw(i,:) = min(Smint_sw(i-1,:), St(i,:));
%         end
%         gamma_mint(i,:) = power_of_signal(i,:) ./ (Bmin .*Smint(i,:));
%         zetat(i,:) = S(i,:) ./ (Bmin .*Smint(i,:));
%         lambda_d(i,:) = Smin(i,:);
%
%         for j = 1:frame_length
%             if (gamma_mint(i,j)>=1 & zetat(i,j)<zeta0)
%                 qhat(i,j) = 1;
%                 phat(i,j) = 0;
%             elseif (1<gamma_mint(i,j)<gamma1 & zetat(i,j)<zeta0)
%                 qhat(i,j) = (gamma1-gamma_mint(i,j)) / (gamma1-1);
%                 phat(i,j) = 1 ./ (1+(qhat(i,j)./(1-qhat(i,j))).*(1+eta(i,j)).*exp(-v(i,j)));
%             else
%                 qhat(i,j) = 0;
%                 phat(i,j) = 1;
%             end
%             alpha_dt(i,j) = alpha_d + (1-alpha_d) * phat(i,j);
%             if I_f(i,j)==1
%                 lambda_dav(i,j) = alpha_dt(i,j) .* lambda_dav(i-1,j) + (1-alpha_dt(i,j)) .* power_of_signal(i,j);
%             else
%                 lambda_dav(i,j) = lambda_dav(i-1,j);
%             end
%             lambda_dc(i,j) = lambda_dav(i,j).* beta;
%         end
%     end
% end
semilogy(power_of_signal(freq_bin,:));
hold on
semilogy(S(freq_bin,:),'r');
hold on
semilogy(lambda_dc(freq_bin,:),'k','LineWidth',2.5);
legend({'Power of Signal','Smoothed Signal','Noise estimate'},'Location','best','FontSize',18);
xlabel('Frame Index','FontSize',20,'FontName','Times');
ylabel('Power in dB','FontSize',20,'FontName','Times');
pow2db(mean(mean(power_of_signal))/mean(mean(abs(lambda_dc))));