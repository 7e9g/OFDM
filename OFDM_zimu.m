clc;
clear;
close all;
%全文原理介绍见：https://zhuanlan.zhihu.com/p/57967971
tic
%% 参数设置

N_sc=64;      %系统子载波数（不包括直流载波）、number of subcarrier
M=4;               %4PSK调制
SNR=-10:1:25;         %仿真信噪比
N_frm=10;            % 每种信噪比下的仿真帧数、frame
Nd=60;               % 每帧包含的OFDM符号数
P_f_inter=6;      %导频间隔
data_station=[];    %导频位置
L=7;                %卷积码约束长度
tblen=6*L;          %Viterbi译码器回溯深度
pilot_interval = 4;       %导频间隔
f_delta = 10e3;     %子载波间隔10KHz
q = 4;
cy_num = 5;%循环次数


is_pilot_k = 1;     %是否插入导频，1为插入（块状导频）
channel_c = 2;%是否通过信道1为通过瑞丽高斯信道2为高斯信道，否则不考虑信道影响
%% 根据采用的子载波数决定IFFT长度 
N_fft=2^ceil(log(N_sc)/log(2));            % FFT 长度
N_cp=N_fft/4;             % 循环前缀长度、Cyclic prefix
N_symbo=N_fft+N_cp;        % 1个完整OFDM符号长度

Band = N_fft * f_delta;                  %基带宽度
fs = q*Band;                                    %数字系统采样率
fc = 10e6;
ts = 1/fs;
fd = 300;                                          %多普勒频偏
pathPower = [-1.0 -1.0 -1.0 0 0 0 -3.0 -5.0 -7.0];
pathDelays = [0 50 120 200 230 500 1600 2300 5000]*1e-9;
% chan = rayleighchan(ts, fd, pathDelays, pathPower);  
rchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
    'MaximumDopplerShift',fd,'FadingTechnique','Sum of sinusoids');

%% 基带数据数据产生
P_data=randi([0 1],1,Nd*N_frm*N_sc*log2(M));

%% 信道编码（卷积码、或交织器）
%卷积码：前向纠错非线性码
%交织：使突发错误最大限度的分散化
trellis = poly2trellis(7,[133 171]);       %(2,1,7)卷积编码
code_data=convenc(P_data,trellis);
data_matintrlv = matintrlv(code_data, log2(M), length(code_data) / log2(M));

%% qpsk调制
data_temp1= reshape(data_matintrlv,[],log2(M));             %以每组2比特进行分组，M=4
data_temp2= bi2de(data_temp1);                             %二进制转化为十进制
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK调制
% figure(1);
scatterplot(modu_data),grid;                  %星座图(也可以取实部用plot函数)

%% 串并转换
data_moded = reshape(modu_data,N_sc,length(modu_data)/N_sc);

%% 补零（使频谱集中）
N_zero = ceil((N_fft-size(data_moded,1))/2);
data_buling = [zeros(N_zero,size(data_moded,2));...
    data_moded;...
    zeros(N_fft-N_zero-size(data_moded,1),size(data_moded,2))];

%% 插入导频
if (is_pilot_k==1)
    pilot_bit_k = randi([0,1],1,log2(M)*N_fft);
    pilot_seq = pskmod(bi2de...
    (reshape(pilot_bit_k,length(pilot_bit_k)/log2(M),log2(M))),M,pi/M);
    pilot_data = insert_pilot_f(data_buling,pilot_seq,pilot_interval,is_pilot_k);
end

%% IFFT
ifft_data = ifft(pilot_data,N_fft).*N_fft;%*sqrt(num_fft);

%% 插入保护间隔、循环前缀
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%把ifft的末尾N_cp个数补充到最前面

%% 并串转换
Tx_data=reshape(Tx_cd,[],1);%由于传输需要

%% 脉冲成型，
sendfir = rcosdesign(0.4,4,fs/Band,'sqrt');
data_upsam = upsample(Tx_data,fs/Band);
data_tx = conv(data_upsam,sendfir,'same');

signal = data_tx;

t1 = ts*(1:length(signal));
% figure;
% plot(t1,real(signal));
%  xlabel('t(s)');
%  ylabel('Magnitude');
%    title('OFDM通信系统发射端时域波形');
% figure;
% plot(t1,imag(signal));
%  xlabel('t(s)');
%  ylabel('Magnitude');
%     title('OFDM通信系统发射端时域与频域波形');
figure;
plot(t1,sqrt(imag(signal).^2+real(signal).^2)); 
 xlabel('t(s)');
 ylabel('Magnitude');
title('OFDM通信系统发射端时域波形');
grid on

% figure;
% fft_y=abs(fft(signal,q*N_fft));
% fft_x=fs*((1:(q*N_fft))/(q*N_fft)-1/2);
% plot(fft_x,20*log10(fftshift(fft_y./max(fft_y))));
%  xlabel('Frequency(Hz)');
%  ylabel('Magnitude');
%  title('OFDM通信系统发射端频域波形');
%  grid on
 
 figure;
%[pxx,f] = pwelch(x,window,noverlap,nfft,fs) pwelch(data_tx,[],[],[],fs,'centered','power') [pxx,f] = 
pwelch(data_tx,[],[],[],fs,'centered','power');
title('OFDM发送信号功率谱');
 
%% 信道（通过多经瑞利信道、或信号经过AWGN信道）
 Ber=zeros(1,length(SNR));
for jj=1:length(SNR)
    for jjj = 1:cy_num
    if(channel_c == 1)%高斯瑞丽
        channel_out = step(rchan,data_tx);
        rx_channel = awgn(channel_out,SNR(jj),'measured');  
        rx_data1 = conv(rx_channel, sendfir, 'same');
        rx_data2 = rx_data1(1:fs/Band:length(rx_data1));
    elseif (channel_c == 2)%高斯
        rx_channel = awgn(data_tx,SNR(jj),'measured');%添加高斯白噪声
        rx_data1 = conv(rx_channel, sendfir, 'same');
        rx_data2 = rx_data1(1:fs/Band:length(rx_data1));
    else %无
        rx_channel = Tx_data;
        rx_data2 = rx_channel;
    end
%% 串并转换

    Rx_data1=reshape(rx_data2,N_fft+N_cp,[]);
    
%% 去掉保护间隔、循环前缀
    Rx_data2=Rx_data1(N_cp+1:end,:);

%% FFT
    fft_data=fft(Rx_data2)./N_fft;

%% 信道校正
% [rx_data_delpilot,H] = get_pilot_f(fft_data,pilot_interval);
[N,NL] = size(fft_data);
H = zeros(N,ceil(NL/(1+pilot_interval)));
output = zeros(N,NL-ceil(NL/(1+pilot_interval)));
jj_g=1;
kk_g=1;
for ii_g = 1:NL
    if mod(ii_g,pilot_interval+1)==1
        H(:,jj_g) = fft_data(:,ii_g);
        jj_g=jj_g+1;
    else
        output(:,kk_g) = fft_data(:,ii_g);
        kk_g=kk_g+1;
    end
end

        rx_data_estimation = chan_estimation_f...
        (output,H,pilot_seq,pilot_interval);

%% 去零
data_aftereq = rx_data_estimation(N_zero+1:N_zero+N_sc,:);
%% 并串转换
    rx_data = reshape(data_aftereq,[],1);
    
%% QPSK解调
    demodulation_data=pskdemod(rx_data,M,pi/M);    
    De_data1 = reshape(demodulation_data,[],1);
    De_data2 = de2bi(De_data1);
    De_Bit = reshape(De_data2,1,[]);

%% （解交织）
	rx_data_jiejiaozi = matdeintrlv(De_Bit, log2(M), length(De_Bit) ./ log2(M));
%% 信道译码（维特比译码）
    rx_c_de = vitdec(rx_data_jiejiaozi,trellis,tblen,'trunc','hard');   %硬判决
%% 计算误码率
    [err, ber] = biterr(rx_c_de(1:length(P_data)),P_data);%译码后的误码率
Ber(jj)=Ber(jj)+ber;
    end
    Ber(jj)=Ber(jj)/cy_num;
end

 

signal1 = rx_channel;
t2 = ts*(1:length(signal1));
% ts*
% figure;
% plot(t2,real(signal1));
%  xlabel('t(s)');
%  ylabel('Magnitude');
%    title('OFDM通信系统接收端时域与频域波形');
% figure;
% plot(t2,imag(signal1));
%  xlabel('t(s)');
%  ylabel('Magnitude');
figure;
plot(t2,sqrt(imag(signal1).^2+real(signal1).^2));
 xlabel('t(s)');
 ylabel('Magnitude');
 title('OFDM通信系统接收端时域波形');
% figure;
% fft_y=abs(fft(signal1,q*N_fft));
% fft_x=fs*((1:(q*N_fft))/(q*N_fft)-1/2);
% plot(fft_x,20*log10(fftshift(fft_y./max(fft_y))));
%  xlabel('Frequency(Hz)');
%  ylabel('Magnitude');
%  title('OFDM通信系统接收端频域波形'); 

figure;
pwelch(rx_channel,[],[],[],fs,'centered','power');
title('OFDM接收信号功率谱');




 figure;
semilogy(SNR,Ber,'r-o');
 %   splot(SNR,Ber,'r-o');
 xlabel('SNR');
 ylabel('BER');
 grid on;
 title('AWGN信道下误比特率曲线');
 
 hold on                                                       %理论高斯信道下BER
a= 4*(1-1/sqrt(M))/log2((M));
k=log2(M);
b= 3*k/(M-1);
ber = a*Q(sqrt(b*10.^(SNR/10)));
semilogy(SNR,ber,'cd-','LineWidth',1);
% axis([1 15 10^-4 1]);
legend('simulation','theory'); 

title('OFDM 16-QAM BER under AWGN channel');
grid on


%   figure;
% %  semilogy(SNR,Ber,'r-o');
%   plot(SNR,Ber,'r-o');
%  xlabel('SNR');
%  ylabel('BER');
%  grid on;
%  title('AWGN信道下误比特率曲线');



 toc
