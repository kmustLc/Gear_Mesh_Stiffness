%% 齿圈和行星轮啮合刚度计算

step = 0.01;
w23=0:step:E2E3/rb2; % w的范围
k23=zeros(size(w23));

kks2=k23;
kkb2=k23;
kkr2=k23;
kkf2=k23;

kks3=k23;
kkb3=k23;
kkb=k23;
kkr3=k23;
kkf3=k23;

for k = 1:length(w23)
        w=w23(k);

    %%  赫兹接触刚度
        kh23 = pi*E*L2/4/(1-v^2);

    %% 『行星轮』任意点K参数计算
        Ak2  = atan((N2E2+w*rb2)/rb2);
        tk2  = tan(Ak2)-phib2;

        yk2 = rb2*(tk2+phib2)*sin(tk2)+rb2*cos(tk2);
        xk2 = rb2*(tk2+phib2)*cos(tk2)-rb1*sin(tk2);

    % 『行星轮』任意点Y参数计算
        syms ty2; 

        yy2 = rb2*(ty2+phib2)*sin(ty2)+rb1*cos(ty2);
        xy2 = rb1*(ty2+phib2)*cos(ty2)-rb1*sin(ty2);        

    % 其他相关参数
        y2  = yk2-yy2;
        Iy2 = 2/3*xy2^3*L2;
        Sy2 = 2*xy2*L2;

        t_start  = -phib2;                    % 基圆处t的值
        t_end    = t_start+tan(Ak2);          % 啮合点k处t的值
        dy_dt_2  = rb2*(ty2+phib2)*cos(ty2);  % yy2对ty2的求导

    
    % 『行星轮』弯曲变形刚度
        fkb2 = (y2*cos(tk2)-xk2*sin(tk2))^2/(E*Iy2)*dy_dt_2;
        kb2  = 1/integral(matlabFunction(fkb2),t_start,t_end);
        
        kkb2(k)=kb2;

    % 『行星轮』剪切变形刚度
        fks2 = 1.2*cos(tk2)^2/G/Sy2*dy_dt_2;
        ks2  = 1/integral(matlabFunction(fks2),t_start,t_end);
        
        kks2(k)=ks2;

    % 『行星轮』径向压缩刚度
        fkr2 = sin(tk2)^2*dy_dt_2;
        kr2  = 1/integral(matlabFunction(fkr2),t_start,t_end);
        
        kkr2(k)=kr2;

    % 『行星轮』基体变形引起的刚度计算
        uf2 = rb2/cos(tk2)-rf2;
        kf2cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif2)^2; h_fi2^2; h_fi2/phif2; 1/phif2; h_fi2; 1]);
        [LL2,MM2,PP2,QQ2]=deal(kf2cell{:});
        kf2=1/(cos(At)^2/(E*L2)*(LL2*(uf2/sf2)^2+MM2*(uf2/sf2)+PP2*(1+QQ2*(tan(Ak2))^2)));

        kkf2(k)= kf2;

    %% 『齿圈』任意点K参数计算
        Ak3 = atan((N3E2+w*rb2)/rb3);
        tk3 = tan(Ak3)+phib3;

        yk3 = rb3*(tk3-phib3)*sin(tk3)+rb3*cos(tk3);
        xk3 = rb3*sin(tk3)-rb3*(tk3-phib3)*cos(tk3);

    % 『齿圈』任意点Y参数计算
        syms ty3;

        yy3 = rb3*(ty3-phib3)*sin(ty3)+rb3*cos(ty3);
        xy3 = rb3*sin(ty3)-rb3*(ty3-phib3)*cos(ty3);
  
        % 其他相关参数
        y3  = yy3-yk3;
        Iy3 = 2/3*xy3^3*L3;
        Sy3 = 2*xy3*L3;

        t_start  = (N3E2+w*rb2)/rb3+phib3;
        Af = acos(rb3/rf3);

        t_end    = tan(Af)+phib3;

        dy_dt_3  = rb3*(ty3-phib3)*cos(ty3);

    % 『齿圈』弯曲变形刚度
        fkb3 = (y3*cos(tk3)-xk3*sin(tk3))^2/(E*Iy3)*dy_dt_3;
        kb3  = 1/integral(matlabFunction(fkb3),t_start,t_end);

        kkb3(k)=kb3;
        
    % 『齿圈』剪切变形刚度
        fks3 = 1.2*cos(tk3)^2/G/Sy3*dy_dt_3;
        ks3  = 1/integral(matlabFunction(fks3),t_start,t_end);
        
        kks3(k)=ks3;

    % 『齿圈』径向压缩刚度
        fkr3 = sin(tk3)^2*dy_dt_3;
        kr3  = 1/integral(matlabFunction(fkr3),t_start,t_end);

        kkr3(k)=kr3;

    % 『齿圈』基体变形引起的刚度计算
       uf3 = rf3-rb3/cos(tk3);
       kf3cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif3)^2; h_fi3^2; h_fi3/phif3; 1/phif3; h_fi3; 1]);
       [LL3,MM3,PP3,QQ3]=deal(kf3cell{:});
       kf3=L3*E/(LL3*(uf3/sf3)^2+MM3*(uf3/sf3)+PP3*(1+QQ3*(tan(Ak3))^2));

       kkf3(k)=kf3;

    %% 刚度汇总
        kk=1/(1/kb3+1/ks3+1/kr3+1/kf3+1/kb2+1/ks2+1/kr2+1/kf2+1/kh23);
        k23(k)=kk;
end