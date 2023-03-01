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
        Ak1 = atan((N1E1_12+w*rb1)/rb1);
        Thetak1 = tan(Ak1)-Ak1;
        tk1 = tan(Ak1)-phib1; % 力F在点K处要分解的力的角
        phik1 = phib1-Thetak1;
        rk1 = rb1/cos(Ak1);
        yk1 = rk1*cos(phik1);
        xk1 = rk1*sin(phik1);

    % 『太阳轮』任意点Y参数计算
        syms Ay1;

        Thetay1 = tan(Ay1)-Ay1;
        phiy1 = phib1-Thetay1;
        ry1 = rb1/cos(Ay1);
        yy1 = ry1*cos(phiy1);
        xy1 = ry1*sin(phiy1);

        y1  = yk1-yy1;

        Iy1 = 2/3*xy1^3*L1;
        Sy1 = 2*xy1*L1;

    % 『太阳轮』弯曲变形刚度

        fkb01 = (y1*cos(tk1)-xk1*sin(tk1))^2/(E*Iy1)*diff(yy1,Ay1);
        fkb01 = matlabFunction(fkb01);
        fkb1 = -integral(fkb01,0,Ak1);
        kb1=1/fkb1;
        
        kkb1(k)=kb1;

    % 『太阳轮』剪切变形刚度
        fks01 = 1.2*cos(Ak1)^2/G/Sy1*diff(y1,Ay1);
        fks01 = matlabFunction(fks01);
        fks1 = -integral(fks01,0,Ak1);
        ks1=1/fks1;

        kks1(k)=ks1;

    % 『太阳轮』径向压缩刚度
        fkr01 = sin(Ak1)^2*diff(y1,Ay1);
        fkr01 = matlabFunction(fkr01);
        fkr1 = -integral(fkr01,0,Ak1);
        kr1=1/fkr1;

        kkr1(k)=kr1;

    % 『太阳轮』基体变形引起的刚度计算
       uf1 = rb1/cos(Ak1-phik1)-rf1;
       kf1cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif1)^2; h_fi1^2; h_fi1/phif1; 1/phif1; h_fi1; 1]);
       [LL1,MM1,PP1,QQ1]=deal(kf1cell{:});
       kf1=L1*E/(LL1*(uf1/sf1)^2+MM1*(uf1/sf1)+PP1*(1+QQ1*(tan(Ak1))^2));
        
       kkf1(k)=kf1;

    %% 『行星轮』任意点K参数计算
        Ak2 = atan((N2E1_12-w*rb1)/rb2);
        Thetak2 = tan(Ak2)-Ak2;
        A2 = Ak2-Thetak2;
        phik2 = phib2-Thetak2;
        rk2 = rb2/cos(Ak2);
        yk2 = rk2*cos(phik2);
        xk2 = rk2*sin(phik2);
        lk2 = yk2-rf2*cos(phif2);
    % 『行星轮』任意点Y参数计算
        syms Ay2;
        Thetay2 = tan(Ay2)-Ay2;
        phiy2 = phib2-Thetay2;
        ry2 = rb2/cos(Ay2);
        yy2 = ry2*cos(phiy2);
        xy2 = ry2*sin(phiy2);
        ly2 = yy2-rf2*cos(phif2);

        y2  = yk2-yy2;

        Iy2 = 2/3*xy2^3*L2;
        Sy2 = 2*xy2*L2;
    % 『行星轮』弯曲变形刚度
        fkb02 = (y2*cos(Ak2)-xk2*sin(Ak2))^2/(E*Iy2)*diff(y2,Ay2);
        fkb02 = matlabFunction(fkb02);
        fkb2 = -integral(fkb02,0,Ak2);
        kb2=1/fkb2;

        kkb2(k)=kb2;

    % 『行星轮』剪切变形刚度
        fks02 = 1.2*cos(Ak2)^2/G/Sy2*diff(y2,Ay2);
        fks02 = matlabFunction(fks02);
        fks2 = -integral(fks02,0,Ak2);
        ks2=1/fks2;

        kks2(k)=ks2;

    % 『行星轮』径向压缩刚度
        fkr02 = sin(Ak2)^2*diff(y2,Ay2);
        fkr02 = matlabFunction(fkr02);
        fkr2 = -integral(fkr02,0,Ak2);
        kr2=1/fkr2;
        
        kkr2(k)=kr2;

    % 『行星轮』基体变形引起的刚度计算
        uf2 = rb2/cos(Ak2-phik2)-rf2;
        kf2cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif2)^2; h_fi2^2; h_fi2/phif2; 1/phif2; h_fi2; 1]);
        [LL2,MM2,PP2,QQ2]=deal(kf2cell{:});
        kf2=L2*E/(LL2*(uf2/sf2)^2+MM2*(uf2/sf2)+PP2*(1+QQ2*(tan(Ak2))^2));

    %% 刚度汇总
        kk=1/(1/kb1+1/ks1+1/kr1+1/kf1+1/kb2+1/ks2+1/kr2+1/kf2+1/kh12);
        k12(k)=kk;
end