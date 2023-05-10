%%
%һ��RC��ͨ�˲�������
%1��������һ��RC�˲����ķ�ֵ˥�����Ժ���������
%2��������һ��RC�˲�����Ƶ������
%3��ʹ��lsim��ϵͳ���з���
%4��ʹ��FFT��ԭʼ�����źź��˲�������źŽ��з���
%���ݺ�����Gs=100000/(s+100000)
%����C=100nF������fc=1/��2*pi*R*C����fc=15.9kHz��֪��RԼΪ25��
%%
%������źż�����Ƶ�ʵȽ��г�ʼֵ�趨
clc
clear 
close all
A=0.8;                % A ����ֵ
fs=5000000;         % fs ����Ƶ��
F=10000;            % F �ź�Ƶ�ʣ�10kHz
N=10000;            % N ��������
dt=1/fs;            %ʱ����
t=0:dt:(N-1)*dt;    %ʱ������

%�����������źŵ�Ƶ�ʽ����趨
F1=40*10^3;
F2=60*10^3;
F3=200*10^3;

%%
%ѡȡ�����ź�
dataSource=0;       %����ѡ��
switch dataSource
    case 0
        y=A*sin(2*pi*F1*t);
    case 1
        y=A*sin(2*pi*F1*t)+A*sin(2*pi*F2*t)+A*sin(2*pi*F3*t);
    case 2
        y=A*square(2*pi*F1*t,50);
    otherwise
end
%%
%���������źŵ�ʱ��Ƶ����
figure(1);
subplot(2,1,1);
plot(t,y);
title('�����źŵ�ʱ����');   %��ԭʼͼ�����ʱ��ͼ
xlabel('ʱ��/s');
ylabel('��ѹ/v');
h=fft(y,10000);            %���ٸ���Ҷ�任
h_d=abs(fftshift(h));     %ʹƵ��ͼ���м�Ϊ��
f=(-N/2:N/2-1)*fs/N;      %��ȡ��ʱ���ϵĵ�ת��ΪƵ���ϵĵ�
subplot(2,1,2);
plot(f,h_d/10000);         %��ԭʼͼ��Ƶ����ͼ
axis([-200*10^3 200*10^3 0 1]);
title('�����źŵ�Ƶ����');
xlabel('Ƶ��/hz');
ylabel('����')
%%
%�Ե�ͨ�˲������г�ֵ�趨��ͬʱ���Ʒ�Ƶ����
r=25;           %������ֵ������
c=100*10^(-9);       %���ݣ�f��
w=2*pi*f;           %��Ƶ��

Func=tf(1,[r*c,1]);     %ϵͳ�Ĵ��ݺ���
figure(2);
bode(Func);             %ϵͳ�Ĳ���ͼ
title('�˲����Ĳ���ͼ');
%%
%�ź�ͨ���˲�������������źŵ�ʱ��Ƶ����
[yout,tout] = lsim(Func,y,t);           %�˲����ź�ͼ��
figure(3);
subplot(2,1,1);
plot(tout,yout);
title('����źŵ�ʱ����');
xlabel('ʱ��/s');
ylabel('��ѹ/v');
q=fft(yout);
q_d=abs(fftshift(q));
subplot(2,1,2);
plot(f,q_d/10000);
title('����źŵ�Ƶ����');
axis([-200*10^3 200*10^3 0 1]);
xlabel('Ƶ��/hz');
ylabel('����');
%%
%���롢����źŵ�����غ���
figure(4)
subplot(2,1,1);
[Rx,maxlags]=xcorr(y,'unbiased');           %�����źŵ������
if fs>10000                                 %����ʱ���ᵥλ����ǩ,���ڹ۲Ⲩ��
    plot(maxlags/fs*1000,Rx/max(Rx));
else
    plot(maxlags/fs,Rx/max(Rx));
end
title('�����ź�����غ���');
xlabel('ʱ��/s');
ylabel('R(t)');
subplot(2,1,2);
[Rx1,maxlags1]=xcorr(yout,'unbiased');      %����źŵ������
if fs>10000                                 %����ʱ���ᵥλ����ǩ,���ڹ۲Ⲩ��
    plot(maxlags1/fs*1000,Rx1/max(Rx));
else
    plot(maxlags1/fs,Rx1/max(Rx));
end
title('����ź�����غ���');
xlabel('ʱ��/s');
ylabel('R(t)');
%%
%�������롢����źŵĹ�����
figure(5);
subplot(2,1,1);
youtpsd=q_d.*conj(q_d);
plot(f,youtpsd);
title('����źŹ�����');
xlabel('Ƶ��/Hz');
ylabel('W/Hz');
subplot(2,1,2);
ypsd=h_d.*conj(h_d);
plot(f,ypsd);
title('�����źŹ�����');
xlabel('Ƶ��/Hz');
ylabel('W/Hz')