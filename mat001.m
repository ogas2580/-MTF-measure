clc;clear;close all;
filename = "Z:\0000083\1.2.156.147522.44.410947.83.5.1.20250105111049.dcm"
dcm = dicomread(filename);
info = dicominfo(filename);
% subplot(1,2,1);imshow(dcm);


%窗宽窗位调整
max_val = max(dcm(:));
min_val = min(dcm(:));
win_width = max_val -min_val;
win_center = (max_val + min_val)/2;
win_width = double(win_width);
win_center = double(win_center);
M=mat2gray(dcm,[win_center-(win_width/2),win_center+(win_width/2)]);
%去除图像边缘黑边
% M = dcm;
i_cut = cutpic(M,0.01);
%i_cut = M;
subplot(1,3,1);imshow(i_cut);

%裁剪出标准的要求区域（求线积分用）
%每个像素对应0.136mm。所以100mm对应735.3像素，算作736pixel；50像素对应367.6像素，算作368pixel
%老师的标准图像：每个像素对应0.188，100mm对应532，50mm对应266
%也就是x,y坐标范围都是从1169-1904，1353-1720
[x_cent,y_cent] = size(i_cut);
xmin = x_cent/2;
ymin = y_cent/2;
xcut = 736;ycut = 368;
img2 = imcrop(i_cut,[(xmin-xcut/2),(ymin-ycut/2),abs(xcut),abs(ycut)]);
subplot(1,3,2);imshow(img2);
%此处得到的img2就是国标要求的ROI区域

%esf函数;img3是对img2作微分得到的"线扩展函数"
[x_size,y_size] = size(img2);
img3 = double(zeros(x_size,y_size));
kernel = [-1, 0, 1];
for x = 1:x_size
    for y = 1:(y_size-1)
        img3(x,:) = conv(img2(x,:), kernel, 'same');
    end
end
% img3_max = max(max(img3));
% img3(find(img3 > img3_max)) = 0;
img3(:,737) = 0;

% for x=1:x_size
%     for y=1:(y_size-1)
%         img3(x,y)=img2(x,y)-img2(x,y+1);
%     end
% end
% 

% 此处是对原图片进行裁切，使得原有的斜边变成竖直的刃边
alpha = 2*deg2rad(1); %2°角
N = 1/tan(alpha);
r = round(N/2);    %r是裁切的半径
for x=1:x_size
    [img3_line_max,img3_max_sit] = max(img3(x,:));
    img4(x,:) = img3(x,(img3_max_sit-r:img3_max_sit+r));
end
%img4是对刃边裁切后的结果，变成竖直边
subplot(1,3,3);imshow(img4);

%求平均ESF信号
[x_size2,y_size2] = size(img4);
ESF_ever = sum(img4,1)/y_size2;


%fft变换（傅里叶变换）https://blog.csdn.net/weixin_42845306/article/details/127062195

% x2 = ESF_ever;
n = y_size2;
LSF = ESF_ever;
% 对LSF进行傅里叶变换得到MTF
MTF = abs(fftshift(fft(LSF)));

% 用MTF零频率的幅度对系数进行归一化处理
MTF = MTF / max(MTF);
% [~,mtfsize] = size(MTF);
% ling = floor(mtfsize/2);
% MTF = MTF /MTF(ling);
% 校正频率轴刻度，单位为lp/mm
% 频率轴的计算基于图像的采样频率和像素间距
N = round(length(MTF)); % MTF的长度
% fs = 1 / (2 * pixel_spacing); % 采样频率，单位为mm^-1

fs = y_size2;
freq_axis_lpmm = fs * (0:(N/2)) / N; % 频率轴，单位为lp/mm

% 绘制MTF曲线
%MTF = flip(MTF);
figure;
plot(freq_axis_lpmm, flip(MTF(1:N/2+1))); % 只绘制正频率部分
xlabel('空间频率 (lp/mm)');
ylabel('MTF值');
title('调制传递函数 (MTF)');
grid on;

% X = abs(fft(x2/(n))); %用fft得出离散傅里叶变换
% X = X(1:floor(n/2+1));
% X(2:end-1) = 2*X(2:end-1);%非0频幅值要乘2
% fs=y_size2;%采样率
% f = fs*(0:(y_size2/2))/y_size2;  %频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
% figure(2);
% plot(f,abs(X));            %画双侧频谱幅度图
%     xlabel("f/Hz")
%     ylabel("幅度")
%     grid on

% function I = cutpic2(I_raw,r)
%     %设置缺省裁减系数r=1
%     if nargin==1
%         r=1;
%     end
%     %裁去纵向黑边
%     % I_raw = int16(I_raw);
%     [m,~] = size(I_raw);bu
%     c = sum(I_raw);
%     z = find(c>m*r);
%     [mm,nn] = size(z);
%     I = I_raw(:,z(mm):z(nn));   
%     %imshow(I);
% end


function I = cutpic(I_raw,r)
    %设置缺省裁减系数r=1
    if nargin==1
        r=1;
    end
    %裁去横向黑边
    [m,n] = size(I_raw);
    b = sum(I_raw,2);
    b = b';
    z = find(b>n*r);
    [mm,nn] = size(z);
    I_raw = I_raw(z(mm):z(nn),:);
    
    %裁去纵向黑边
    [m,n] = size(I_raw);
    c = sum(I_raw);
    z = find(c>m*r);
    [mm,nn] = size(z);
    I = I_raw(:,z(mm):z(nn));   
    %imshow(I);
end