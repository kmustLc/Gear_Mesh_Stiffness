%% 太阳轮和行星轮啮合刚度计算

step = 0.01;
w12=0:step:E1E2/rb1; % w的范围
k12=zeros(size(w12));

kks1=k12;
kkb1=k12;
kkr1=k12;
kkf1=k12;

kks2=k12;
kkb2=k12;
kkr2=k12;
kkf2=k12;

for k = 1:length(w12)

        w=w12(k);

    %%  赫兹接触刚度
        kh12 = pi*E*L1/4/(1-v^2);

    %% 『太阳轮』任意点K参数计算
        Ak1 = atan((N1E1+w*rb1)/rb1);
        tk1 = tan(Ak1)-phib1; 

        yk1 = rb1*(tk1+phib1)*sin(tk1)+rb1*cos(tk1);
        xk1 = rb1*(tk1+phib1)*cos(tk1)-rb1*sin(tk1);

    % 『太阳轮』任意点Y参数计算
        syms ty1;

        yy1 = rb1*(ty1+phib1)*sin(ty1)+rb1*cos(ty1);
        xy1 = rb1*(ty1+phib1)*cos(ty1)-rb1*sin(ty1);

    % 其他相关参数
        y1  = yk1-yy1;
        Iy1 = 2/3*xy1^3*L1;
        Sy1 = 2*xy1*L1;

        t_start  = -phib1;                    % 基圆处t的值
        t_end    = t_start+tan(Ak1);          % 啮合点k处t的值
        dy_dt_1  = rb1*(ty1+phib1)*cos(ty1);   % yy1对ty1的求导

    % 『太阳轮』弯曲变形刚度
        fkb1 = (y1*cos(tk1)-xk1*sin(tk1))^2/(E*Iy1)*dy_dt_1;
        kb1  = 1/integral(matlabFunction(fkb1),t_start,t_end);
        
        kkb1(k)=kb1;

    % 『太阳轮』剪切变形刚度
        fks1 = 1.2*cos(tk1)^2/G/Sy1*dy_dt_1;
        ks1 = 1/integral(matlabFunction(fks1),t_start,t_end);

        kks1(k)=ks1;

    % 『太阳轮』径向压缩刚度
        fkr1 = sin(tk1)^2*dy_dt_1;
        kr1 = 1/integral(matlabFunction(fkr1),t_start,t_end);

        kkr1(k)=kr1;

    % 『太阳轮』基体变形引起的刚度计算
       uf1 = rb1/cos(tk1)-rf1;
       kf1cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif1)^2; h_fi1^2; h_fi1/phif1; 1/phif1; h_fi1; 1]);
       [LL1,MM1,PP1,QQ1]=deal(kf1cell{:});
       kf1=1/(cos(At)^2/(E*L2)*(LL1*(uf1/sf1)^2+MM1*(uf1/sf1)+PP1*(1+QQ1*(tan(Ak1))^2)));
       
       kkf1(k)=kf1;

    %% 『行星轮』任意点K参数计算
        Ak2 = atan((N2E1-w*rb1)/rb2);
        tk2 = tan(Ak2)-phib2;

        yk2 = rb2*(tk2+phib2)*sin(tk2)+rb2*cos(tk2);
        xk2 = rb2*(tk2+phib2)*cos(tk2)-rb2*sin(tk2);

    % 『行星轮』任意点Y参数计算
        syms ty2;
        
        yy2 = rb2*(ty2+phib2)*sin(ty2)+rb2*cos(ty2);
        xy2 = rb2*(ty2+phib2)*cos(ty2)-rb2*sin(ty2);

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

        kkf2(k)=kf2;
    %% 刚度汇总
        kk=1/(1/kb1+1/ks1+1/kr1+1/kf1+1/kb2+1/ks2+1/kr2+1/kf2+1/kh12);
        k12(k)=kk;
end