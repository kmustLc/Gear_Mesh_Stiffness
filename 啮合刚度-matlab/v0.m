%%
clc;
clear all;
close all;
%% 定义参数
syms r1 r2 rb1 rb2 rf1 rf2 ra1 ra2 rint1 rint2; %分度圆、基圆、齿根圆、齿顶圆半径、齿轮内径
syms zn1 zn2 z1 z2 x1 x2 b1 b2 mn12 mt12; %主从动齿轮当量齿数、齿数、变位系数、齿宽、法向模数、端面模数
syms pn12 pbn12 pbt12; %法面齿距、法向基圆齿距、端面基圆齿距
syms B12 Bb12 An12 At12 lamda12 lamdaa12 lamdab12 hfp12; %分度圆螺旋角、基圆螺旋角、法向压力角、端面压力角、总重合度、端面重合度、轴向重合度、齿根高
syms E G v; %材料弹性模量、剪切模量、泊松比
syms Kb1 Ks1 Kr1 Kf1 Kb2 Ks2 Kr2 Kf2 Kh Kt Kti angle;
syms w; %定义一个齿自啮合点开始，齿轮转过的角度
syms y1 y2; %计算一个轴向切片上的刚度时微分的横向截面离啮合点的距离
syms lk1 lk2 ly1 ly2; %一个轴向切片上任一啮合点K距齿根的距离、任一微分点Y距齿根的距离
syms Ak1 Ak2 Ay1 Ay2; %一个轴向切片上任一啮合点K的压力角、任一微分点Y的压力角
syms theta thetak1 thetak2 thetaf1 thetaf2 thetay1 thetay2; %分度圆处的展角、一个轴向切片上任一啮合点K的展角、齿根处的展角、任一微分点Y的展角
syms phi1 phi2 phik1 phik2 phiy1 phiy2 phif1 phif2; %主被动齿轮基圆处较y轴的旋转角、一个轴向切片上任一啮合点K较y轴的旋转角、任一微分点Y较y轴的旋转角、齿根较y轴的旋转角
syms rk1 rk2 ry1 ry2; %一个轴向切片上任一啮合点K的向径、任一微分点Y的向径
syms sk1 sk2 sy1 sy2; %一个轴向切片上任一啮合点K距y轴的距离、任一微分点Y距y轴的距离
syms N1 N2 bN; %齿轮沿轴向分成多少个切片、每片的厚度
syms Iy1 Iy2 AAy1 AAy2; %微分截面上惯性矩和截面积
syms uf1 uf2 sf1 sf2; %啮合力作用方向与轮齿中线交点到齿根圆最高点距离、齿根圆上一个轮齿所占的弧长
syms Ai Bi Ci Di Ei Fi LL1 LL2 MM1 MM2 PP1 PP2 QQ1 QQ2 Kf1cell Kf2cell; %算轮体变形刚度时用的系数
syms N1N2 N1B2 N2B1 B1B2; %啮合线上的各种长度
% syms A1start A1stop; %主动轮端面上开始啮合点压力角与结束啮合点压力角
syms Ca La Cb; %修形参数

%% 设置修形参数
Ca=0; %齿廓齿顶修形量 m
La=2.5e-3; %齿廓齿顶修形长度 m
Cb=1e-10; %齿向鼓形修形量 m

%% 齿轮副参数
B12=29/180*pi; %分度圆螺旋角
mn12=1.51; %法向模数 mm
mt12=mn12/cos(B12); %端面模数 mm
An12=18/180*pi; %法向压力角
At12=atan(tan(An12)/cos(B12)); %端面压力角
Bb12=atan(tan(B12)*cos(At12)); %基圆螺旋角
pn12=pi*mn12*1e-3; %法面齿距 m
pbn12=pn12*cos(An12); %法向基圆齿距 m
pbt12=pbn12/cos(B12); %端面基圆齿距 m
theta=tan(At12)-At12; %分度圆处的展角
% lamda12=4.9281; %总重合度
hfp12=1.25*mn12*1e-3; %齿根高
E=2.07*1e11; %弹性模量 pa
G=8.17*1e10; %剪切模量 pa
v=0.289; %泊松比
bN=1e-3; %设置齿轮切片的厚度
Ai=[-5.574*1e-5 60.111*1e-5 -50.952*1e-5 -6.2042*1e-5];
Bi=[-1.9986*1e-3 28.1*1e-3 185.5*1e-3 9.0889*1e-3];
Ci=[-2.3015*1e-4 -83.431*1e-4 0.0538*1e-4 -4.0964*1e-4];
Di=[4.7702*1e-4 -9.9256*1e-3 53.3*1e-3 7.8297*1e-3];
Ei=[0.0271 0.1624 0.2895 -0.1472];
Fi=[6.8045 0.9086 0.9236 0.6904];

%% 主动轮的参数
r1=0.022444; %主动轮分度圆半径 m
rint1=0.01; %主动轮内径 m
z1=26; %主动轮齿数
b1=0.034; %主动轮宽度 m
x1=0.1527; %主动轮变位系数
ra1=r1+mn12*(1+x1)*1e-3; %主动轮齿顶圆半径 m
rf1=r1-mn12*(1+0.25-x1)*1e-3; %主动轮齿根圆半径 m
rb1=r1*cos(At12); %主动轮基圆半径 m
% zn1=z1/(cos(B12))^3; %主动齿轮当量齿数
sf1=2*rf1*pi/z1; %主动轮齿根圆上一个轮齿所占弧长
phif1=pi/z1; %主动轮齿根处较y轴的旋转角
phi1=pi/2/z1+theta; %主动齿轮基圆处较y轴的旋转角
N1=b1/bN;
L1=ra1-rf1*cos(phif1); %%% 齿全高
s1=rb1/cos(At12)*sin(pi/2/z1);

%% 从动轮的参数
r2=0.0578365; %从动轮分度圆半径 m
rint2=0.02; %从动轮内径 m
z2=67; %从动轮齿数
b2=0.031; %从动轮宽度 m
x2=-0.3361; %从动轮变位系数
ra2=r2+mn12*(1+x2)*1e-3; %从动轮齿顶圆半径 m
rf2=r2-mn12*(1+0.25-x2)*1e-3; %从动轮齿根圆半径 m
rb2=r2*cos(At12); %从动轮基圆半径 m
zn2=z2/(cos(B12))^3; %从动齿轮当量齿数
sf2=2*rf2*pi/z2; %从动轮齿根圆上一个轮齿所占弧长
phif2=pi/z2; %从动轮齿根处较y轴的旋转角
phi2=pi/2/z2+theta; %从动齿轮基圆处较y轴的旋转角
N2=b2/bN;
L2=ra2-rf2*cos(phif2);
s2=rb2/cos(At12)*sin(pi/2/z2);

%% 补充齿轮副参数
% lamdab12=min(b1,b2)*tan(Bb12)/pbt12; %轴向重合度
% lamdaa12=lamda12-lamdab12; %端面重合度
% A1start=2*At12-acos(rb2/ra2);
% A1stop=acos(rb1/ra1);

N1N2=((r1+r2)^2-(rb1+rb2)^2)^0.5; %理论啮合线总长度
N1B2=(ra1^2-rb1^2)^0.5; %理论啮合线左端到开始啮合点的长度
N2B1=(ra2^2-rb2^2)^0.5; %理论啮合线右端到啮合结束点的长度
B1B2=N1B2+N2B1-N1N2; %实际啮合线长度

%% 循环计算一个轴向切片直齿轮转一个周期的啮合刚度
Kt=[];
angle=[];
outsidestep=bN*tan(Bb12)/rb1; %计算的啮合点角度间隔
insidestep=0.002; %计算的微分角度间隔
for w=0:outsidestep:B1B2/rb1
    
    % 计算主动齿轮一个轴向切片直齿轮上某一啮合位置处的单齿啮合刚度
    
    Ak1=atan((N1B2-B1B2+w*rb1)/rb1);
    thetak1=tan(Ak1)-Ak1;
    thetay1=tan(Ay1)-Ay1;
    phik1=phi1-thetak1;
    phiy1=phi1-thetay1;
    rk1=rb1/cos(Ak1);
    ry1=rb1/cos(Ay1);
    sk1=rk1*sin(phik1);
    lk1=rk1*cos(phik1)-rf1*cos(phif1);
%     sy1=ry1*sin(phiy1);
    ly1=ry1*cos(phiy1)-rf1*cos(phif1);
    Ly1=L1-ly1;
    y1=lk1-ly1;
    sy1=ry1*sin(phiy1)-Ca*((La-Ly1)/La)^2;
%     if La-Ly1>0
%         sy1=sy1-Ca*((La-Ly1)/La)^2;
%     else
%         sy1=sy1-0;
%     end
    Iy1=2*sy1^3*(bN/cos(Bb12));
    AAy1=2*sy1*(bN/cos(Bb12));
    fKb1=(y1*cos(Ak1)-sy1*sin(Ak1))^2/(E*Iy1)*diff(y1,Ay1);
    fKs1=1.2*cos(Ak1)^2/(G*AAy1)*diff(y1,Ay1);
    fKr1=sin(Ak1)^2/(E*AAy1)*diff(y1,Ay1);
    Kh=pi*E*(bN/cos(Bb12))/(4*(1-v^2));
    fKb1m=0;
    for Ay1=-phif1:insidestep:Ak1
        fKb1m=fKb1m+eval(fKb1)*(-0.001);
    end
    Kb1=1/fKb1m;
    fKs1m=0;
    for Ay1=-phif1:insidestep:Ak1
        fKs1m=fKs1m+eval(fKs1)*(-0.001);
    end
    Ks1=1/fKs1m;
    fKr1m=0;
    for Ay1=-phif1:insidestep:Ak1
        fKr1m=fKr1m+eval(fKr1)*(-0.001);
    end
    Kr1=1/fKr1m;

    Kf1cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(pi/z1)^2;(rf1/rint1)^2;(rf1/rint1)/(pi/z1);1/(pi/z1);(rf1/rint1);1]);
    [LL1,MM1,PP1,QQ1]=deal(Kf1cell{:});
    uf1=rb1/cos(Ak1)-rf1;
    Kf1=(bN/cos(Bb12))*E/(LL1*(uf1/sf1)^2+MM1*(uf1/sf1)+PP1*(1+QQ1*(tan(Ak1))^2));

    % 计算从动齿轮一个轴向切片直齿轮上某一啮合位置处的单齿啮合刚度
    Ak2=atan((N2B1-w*rb1)/rb2);
    thetak2=tan(Ak2)-Ak2;
    thetay2=tan(Ay2)-Ay2;
    phik2=phi2-thetak2;
    phiy2=phi2-thetay2;
    rk2=rb2/cos(Ak2);
    ry2=rb2/cos(Ay2);
    sk2=rk2*sin(phik2);
    lk2=rk2*cos(phik2)-rf2*cos(phif2);
%     sy2=ry2*sin(phiy2);
    ly2=ry2*cos(phiy2)-rf2*cos(phif2);
    Ly2=L2-ly2;
    y2=lk2-ly2;
    sy2=ry2*sin(phiy2)-Ca*((La-Ly2)/La)^2;
%     if La-Ly2>0
%         sy2=sy2-Ca((La-Ly2)/La)^2;
%     else
%         sy2=sy2-0;
%     end
    Iy2=2*sy2^3*(bN/cos(Bb12));
    AAy2=2*sy2*(bN/cos(Bb12));
    fKb2=(y2*cos(Ak2)-sy2*sin(Ak2))^2/(E*Iy2)*diff(y2,Ay2);
    fKs2=1.2*cos(Ak2)^2/(G*AAy2)*diff(y2,Ay2);
    fKr2=sin(Ak2)^2/(E*AAy2)*diff(y2,Ay2);
    fKb2m=0;
    for Ay2=-phif2:insidestep:Ak2
        fKb2m=fKb2m+eval(fKb2)*(-0.001);
    end
    Kb2=1/fKb2m;
    fKs2m=0;
    for Ay2=-phif2:insidestep:Ak2
        fKs2m=fKs2m+eval(fKs2)*(-0.001);
    end
    Ks2=1/fKs2m;
    fKr2m=0;
    for Ay2=-phif2:insidestep:Ak2
        fKr2m=fKr2m+eval(fKr2)*(-0.001);
    end
    Kr2=1/fKr2m;
    Kf2cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(pi/z2)^2;(rf2/rint2)^2;(rf2/rint2)/(pi/z2);1/(pi/z2);(rf2/rint2);1]);
    [LL2,MM2,PP2,QQ2]=deal(Kf2cell{:});
    uf2=rb2/cos(Ak2)-rf2;
    Kf2=(bN/cos(Bb12))*E/(LL2*(uf2/sf2)^2+MM2*(uf2/sf2)+PP2*(1+QQ2*(tan(Ak2))^2));

    % 计算一个轴向切片直齿轮上某一啮合位置时的单齿啮合刚度
    Kti=1/(1/Kb1+1/Ks1+1/Kr1+1/Kf1+1/Kb2+1/Ks2+1/Kr2+1/Kf2+1/Kh);
    Kt=[Kt Kti];
    angle=[angle w];
    syms Ay1 Ay2;
end
Kt=Kt*cos(Bb12)^2;
 
%% 画端面一个轴向切片直齿轮的啮合刚度变化曲线
figure(1)
set(gcf,'position',[0 0 1500 1000]);
plot(angle/pi*180,Kt,'*');
 set(gca,'fontsize',30);
       xlabel('主动齿轮从啮合起始点开始转过的角度（°）','fontsize',30);
       ylabel('啮合刚度（N/m）','fontsize',30);
       title('轴向薄片齿轮上一个轮齿的时变啮合刚度变化曲线','fontsize',30);
       
%%
angle_t=0:outsidestep:2*pi/z1*(z1+3)-outsidestep;
Ktt=zeros(1,length(angle_t));
for i=0:1:z1
    Ktt(round(length(angle_t)/(z1+3))*i+1:round(length(angle_t)/(z1+3))*i+length(angle)) ...
    =Ktt(round(length(angle_t)/(z1+3))*i+1:round(length(angle_t)/(z1+3))*i+length(angle))+Kt;
end
Ktt=Ktt(round(length(angle_t)/(z1+3)):round(length(angle_t)/(z1+3))*(z1+1)-1);
angle_t=0:outsidestep:(length(Ktt)-1)*outsidestep;

%% 画端面直齿轮的啮合刚度变化曲线
figure(2)
set(gcf,'position',[0 0 1500 1000]);
plot(angle_t(1:30)/pi*180,Ktt(1:30),'linewidth',5);
 set(gca,'fontsize',30);
       xlabel('主动齿轮从啮合起始点开始转过的角度（°）','fontsize',30);
       ylabel('啮合刚度（N/m）','fontsize',30);
       title('单个轴向薄片齿轮的综合时变啮合刚度变化曲线','fontsize',30);

%% 计算出整个斜齿轮的总啮合刚度
delta=bN*tan(Bb12)/rb1;
Kt_total=Ktt;
Kt_slide=Ktt;
for i=1:1:min(N1,N2)-1
    x=(i-1)*bN;
    R=Cb/2+(min(b1,b2))^2/8/Cb;
    Ci=R-sqrt(R^2-(min(b1,b2)/2-x)^2);
    Kt_slide(1:i)=Ktt(length(Ktt)-i+1:length(Ktt));
    Kt_slide(i+1:length(Ktt))=Ktt(1:length(Ktt)-i);
    Kt_total=Kt_total+Kt_slide*(s1-Ci)/s1;
end

%% 画整个斜齿轮的总啮合刚度
figure(3)
set(gcf,'position',[0 0 1500 1000]);
plot(angle_t(9:19)/pi*180,Kt_total(9:19),'linewidth',5);
 set(gca,'fontsize',30);
       xlabel('主动齿轮从啮合起始点开始转过的角度（°）','fontsize',30);
       ylabel('啮合刚度（N/m）','fontsize',30);
       title('斜齿轮总体的啮合刚度变化曲线','fontsize',30);    

%%
% figure(4)
% set(gcf,'position',[0 0 1500 1000]);
% plot(angle_t(1:294)/pi*180,Kt_total(1:294),'b','linewidth',5);hold on;
% plot(kisssoft.ANGLE+24.4179,kisssoft.STIFFNESS*31*1e6,'r','linewidth',5);
% legend('本文计算方法','KISSsoft软件计算方法');
%  set(gca,'fontsize',30);
%  ylim([5.2e8 5.7e8]);
%        xlabel('主动齿轮从啮合起始点开始转过的角度（°）','fontsize',30);
%        ylabel('啮合刚度（N/m）','fontsize',30);
%        title('啮合刚度变化曲线对比','fontsize',30);  

















