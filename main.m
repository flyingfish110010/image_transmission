clear;
clear all;
clc;

pic_data = imread("testpic.jpg");% 得到一个三维数组

% 显示图片
% figure(1);
% imshow(pic_data);
% title('Input Image');

% 将图片转换为二进制
pic_data_seq = reshape(pic_data,1,numel(pic_data),1);
trans2bi = int2bit(pic_data_seq,8);
bit_seq = reshape(trans2bi,1,numel(trans2bi));

% 设置参数
fc = 5e6;
fs = 25e6;
Nbit = length(bit_seq);
Nsymb = Nbit / 2;
% Nsymb = 50000;
TB = 1 / 5000000;% 符号周期
SNR = -3 : 2 : 29;
ber = zeros(1,length(SNR));

% data_input = zeros(1,Nbit);
data_input = double(bit_seq);
% data_input = round(rand(1,100000));

% 转为极性码
data_NRZ = 2 .* data_input - 1;
data = reshape(data_NRZ,2,numel(data_NRZ) / 2);
% 分为IQ两路
I = zeros(1,Nsymb);
Q = zeros(1,Nsymb);
I = data(1,:);
Q = data(2,:);

% 插0
zero = fs * TB;
I0 = zeros(1,zero * Nsymb);
Q0 = zeros(1,zero * Nsymb);
index = 1 : zero : Nsymb * zero;
I0(1,index) = I;
Q0(1,index) = Q;

% 平方根升余弦滤波器
NT = 10; % 滤波器阶数
rf = 0.25; % 滚降系数
BW = (1 + rf) / 2 .* (1 / 2 .* TB);% 带宽
psf = rcosfir(rf,NT,zero,1/fs,'sqrt');
Ipulse = conv(I0,psf);
Qpulse = conv(Q0,psf);

% 调制
num = NT * zero * 2;
Ipulse = Ipulse(1,1:end - num);
Qpulse = Qpulse(1,1:end - num);
t = 0 : 1 / fs : Nsymb * zero / fs - 1 / fs;
Imod = Ipulse .* sqrt(2) .* cos(2*pi*fc*t);
Qmod = Qpulse .* sqrt(2) .* sin(2*pi*fc*t);
QPSK = Imod - Qmod;

for k = 1:length(SNR)
    
    % 添加噪声
    QPSK_in = awgn(QPSK,SNR(k),'measured');

    %接收滤波器
    [B,A] = butter(4,[0.08,0.8],'bandpass');
    QPSK_in2 = filter(B,A,QPSK_in);

    % 接收端,解调
    Idem = QPSK_in2 .* sqrt(2) .* cos(2*pi*fc*t);
    Qdem = -QPSK_in2 .* sqrt(2) .* sin(2*pi*fc*t);

    % 低通滤波器
    [B,A] = butter(4,0.95,'low');
    Idem = filter(B,A,Idem);
    Qdem = filter(B,A,Qdem);

    % 匹配滤波
    mtf = rcosfir(rf,NT,zero,fs,'sqrt');
    Idem2 = conv(Idem,mtf);
    Qdem2 = conv(Qdem,mtf);
    Isel = Idem2(1,num+1:end);
    Qsel = Qdem2(1,num+1:end);

    I_sam = Isel(1,index);
    Q_sam = Qsel(1,index);

    threshold = 0.0;
    I_final = sign(I_sam);
    zero_m1 = find(~I_final);
    for i = 1:length(zero_m1)
        I_final(zero_m1(i)) = round(rand*2 - 1);
    end

    Q_final = sign(Q_sam);
    zero_m2 = find(~Q_final);
    for i = 1:length(zero_m2)
        Q_final(zero_m2(i)) = round(rand*2 - 1);
    end

    % 并串变换
    data_dem_parallel(1,:) = I_final;
    data_dem_parallel(2,:) = Q_final;

    data_dem_seq = reshape(data_dem_parallel,1,numel(data_dem_parallel));
    data_output = (data_dem_seq + 1) / 2;
    
    % 误码率
    ber(k) = length(find(data_output ~= data_input)) / Nbit;

    out(:,:,k) = data_output;

end


