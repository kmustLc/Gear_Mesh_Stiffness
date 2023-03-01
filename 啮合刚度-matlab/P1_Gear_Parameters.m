%%这很神奇
clc
clear
close all

mmm = 1e-3;     % 设置长度单位是m还是mm（1e-3,1→mm）

%%基本参数输入（1是太阳轮、2是行星轮、3是齿圈）

z1 = 41;
z2 = 145;
z3 = 130;

mn = 2;               % 法向模数
An = 20*pi/180;       % 压力角 『弧度制』
B  = 21.56*pi/180;    % 螺旋角 『弧度制』

x1 = 0.345;              % 变位系数
x2 = -0.345;
x3 = 0;


% Ca=0*mmm;           % 齿廓齿顶修形量   (mmm)
% La=2.5*mmm;         % 齿廓齿顶修形量   (mmm)
% Cb=1e-7*mmm;        % 齿向鼓形修形量   (mmm)

b1 = 50*mmm;        % 齿宽            (mmm)
b2 = 45*mmm;        
b3 = 34*mmm;        

bN=1*mmm;           % 切片厚度   (mmm)


r_int1 = 28*mmm;    % 外齿内径，内齿外径   (mmm)
r_int2 = 148*mmm;
r_int3 = 300*mmm;


E=2.07*1e11;    % 弹性模量
v=0.289;        % 泊松比
G=E/(2*(1+v));  % 弹性模量

%%相关参数推导
mt = mn*cos(B);                   % 端面模数
At = atan(tan(An)/cos(B));        % 端面压力角

r1 = z1*mt/2*mmm;               % 分度圆半径 (mmm)
r2 = z2*mt/2*mmm;
r3 = z3*mt/2*mmm; 

rb1 = r1*cos(At);                % 基圆半径 (mmm)
rb2 = r2*cos(At);
rb3 = r3*cos(At);

ra1 = r1+mn*(1+x1)*mmm;         % 齿顶圆半径 (mmm)
ra2 = r2+mn*(1+x2)*mmm;
ra3 = r3-mn*(1+x3)*mmm;

rf1 = r1-mn*(1+0.25-x1)*mmm;    % 齿根圆半径 (mmm)
rf2 = r2-mn*(1+0.25-x2)*mmm;
rf3 = r3+mn*(1+0.25+x3)*mmm;

Bb1 = atan(tan(B)*rb1/r1);       % 基圆螺旋角
Bb2 = atan(tan(B)*rb2/r2);
Bb3 = atan(tan(B)*rb3/r3);

Ba1 = atan(tan(B)*ra1/r1);       % 齿顶圆螺旋角
Ba2 = atan(tan(B)*ra2/r2);
Ba3 = atan(tan(B)*ra3/r3);

Theta = tan(At)-At;        % 分度圆处的展角（是At，不是An）

phif1=pi/z1;               % 齿根圆较y轴的旋转角
phif2=pi/z2;
phif3=pi/z3;

phib1=pi/2/z1+Theta;       % 基圆较y轴的旋转角
phib2=pi/2/z2+Theta;
phib3=pi/2/z3-Theta;

h_fi1=rf1/r_int1;          % 齿轮内径与外径的比值
h_fi2=rf2/r_int2;
h_fi3=r_int3/rf3;

sf1=rf1*2*pi/z1;           % 一个齿所占的弧长
sf2=rf2*2*pi/z2;
sf3=rf3*2*pi/z3;

%%有效啮合长度计算 & 开始结束啮合点压力角
%----------太阳轮和行星轮啮合----------
N1E2 = sqrt(ra1^2-rb1^2);
N2E1 = sqrt(ra2^2-rb2^2);
E1E2 = N1E2+N2E1-(r1+r2)*sin(At); % 有效啮合长度

N1E1 = N1E2-E1E2; 

%------------齿圈和行星轮啮合-----------
N2E3 = sqrt(ra2^2-rb2^2); 
N3E2 = sqrt(ra3^2-rb3^2); 
E2E3 = N2E3-N3E2+(r3-r2)*sin(At); % 有效啮合长度

N2E2 = N2E3-E2E3; 
N3E3 = N3E2+E2E3;

%%基体变形相关参数
Ai=[-5.574*1e-5 60.111*1e-5 -50.952*1e-5 -6.2042*1e-5];
Bi=[-1.9986*1e-3 28.1*1e-3 185.5*1e-3 9.0889*1e-3];
Ci=[-2.3015*1e-4 -83.431*1e-4 0.0538*1e-4 -4.0964*1e-4];
Di=[4.7702*1e-4 -9.9256*1e-3 53.3*1e-3 7.8297*1e-3];
Ei=[0.0271 0.1624 0.2895 -0.1472];
Fi=[6.8045 0.9086 0.9236 0.6904];


%% 一个切片的宽度
L1  = bN/cos(Bb1);
L2  = bN/cos(Bb2);
L3  = bN/cos(Bb3);