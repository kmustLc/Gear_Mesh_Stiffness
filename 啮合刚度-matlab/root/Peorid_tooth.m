function [ww,kk]=Peorid_tooth(w,k,z)
% 将一个齿的啮合刚度进行平移叠加，得到整个齿轮的时变啮合刚度
% 返回叠加后的刚度随转角的变化曲线
% w为转角矩阵（弧度）
% k为刚度
% z为齿数
    
    w=w*180/pi;
    step=360/z; % 一个齿所占的弧长
    
    % 将原始数据进行插值
    num = 1000;
    x0=linspace(w(1),w(end),num);
    y0=spline(w,k,x0);
    
    y=y0;
    
    f=round(step/w(end)*num);

    for k=f:num
        k0=k-f+1;
        y(k)=y0(k0)+y0(k);
    end
    y1 = y0;
    for k=1:num-f
        k0=k+f;
        y(k)=y1(k)+y1(k0);
    end
    ww = x0;
    kk = y;
end