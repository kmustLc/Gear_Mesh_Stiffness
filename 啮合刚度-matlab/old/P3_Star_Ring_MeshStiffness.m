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
kkr3=k23;
kkf3=k23;

for k = 1:length(w23)
        w=w23(k);
    %%  赫兹接触刚度
        kh23 = pi*E*L2/4/(1-v^2);

    %% 『行星轮』任意点K参数计算
        Ak22     = atan((N2E2_23+w*rb2)/rb2);
        Thetak22 = tan(Ak22)-Ak22;
        phik22   = phib2-Thetak22;
        rk22     = rb2/cos(Ak22);
        yk22     = rk22*cos(phik22);
        xk22     = rk22*sin(phik22);

    % 『行星轮』任意点Y参数计算
        syms Ay22;
        Thetay22 = tan(Ay22)-Ay22;
        phiy22   = phib2-Thetay22;
        ry22     = rb2/cos(Ay22);
        yy22     = ry22*cos(phiy22);
        xy22     = ry22*sin(phiy22);

        y22  = yk22-yy22;

        Iy22 = 2/3*xy22^3*L2;
        Sy22 = 2*xy22*L2;
    % 『行星轮』弯曲变形刚度
        fkb02 = (y22*cos(Ak22)-xk22*sin(Ak22))^2/(E*Iy22)*diff(y22,Ay22);
        fkb02 = matlabFunction(fkb02);
        fkb2 = -integral(fkb02,0,Ak22);
        kb2=1/fkb2;
        
        kkb2(k)=kb2;

    % 『行星轮』剪切变形刚度
        fks02 = 1.2*cos(Ak22)^2/G/Sy22*diff(y22,Ay22);
        fks02 = matlabFunction(fks02);
        fks2 = -integral(fks02,0,Ak22);
        ks2=1/fks2;
        
        kks2(k)=ks2;

    % 『行星轮』径向压缩刚度
        fkr02 = sin(Ak22)^2*diff(y22,Ay22);
        fkr02 = matlabFunction(fkr02);
        fkr2 = -integral(fkr02,0,Ak22);
        kr2=1/fkr2;
        
        kkr2(k)=kr2;

    % 『行星轮』基体变形引起的刚度计算
        uf22 = rb2/cos(Ak22-phik22)-rf2;
        kf2cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif2)^2; h_fi2^2; h_fi2/phif2; 1/phif2; h_fi2; 1]);
        [LL2,MM2,PP2,QQ2]=deal(kf2cell{:});
        kf2=L2*E/(LL2*(uf22/sf2)^2+MM2*(uf22/sf2)+PP2*(1+QQ2*(tan(Ak22))^2));
        
        kkf2(k)=kf2;

    %% 『齿圈』任意点K参数计算
        Ak3 = atan((N3E2_23+w*rb2)/rb3);
        Thetak3 = tan(Ak3)-Ak3;
        phik3 = phib3+Thetak3;
        rk3 = rb3/cos(Ak3);
        yk3 = rk3*cos(phik3);
        xk3 = rk3*sin(phik3);

    % 『齿圈』任意点Y参数计算
        syms Ay3;
        Thetay3 = tan(Ay3)-Ay3;
        phiy3 = phib3+Thetay3;
        ry3 = rb3/cos(Ay3);
        yy3 = ry3*cos(phiy3);
        xy3 = ry3*sin(phiy3);
  
        y3  = yy3-yk3;
    
        Iy3 = 2/3*xy3^3*L3;
        Sy3 = 2*xy3*L3;
        
    % 『齿圈』弯曲变形刚度
        fkb03 = (y3*cos(Ak3)-xk3*sin(Ak3))^2/(E*Iy3)*diff(y3,Ay3);
        
        fkb03 = matlabFunction(fkb03);
        fkb3 = integral(fkb03,Ak3,A3_m);
        kb3=1/fkb3;

        kkb3(k)=kb3;

    % 『齿圈』剪切变形刚度
        fks03 = 1.2*cos(Ak3)^2/G/Sy3*diff(y3,Ay3);
        fks03 = matlabFunction(fks03);
        fks3 = integral(fks03,Ak3,A3_m);
        ks3=1/(fks3);
        
        kks3(k)=ks3;

    % 『齿圈』径向压缩刚度
        fkr03 = sin(Ak3)^2*diff(y3,Ay3);
        fkr03 = matlabFunction(fkr03);
        fkr3 = integral(fkr03,Ak3,A3_m);
        kr3=1/fkr3;

        kkr3(k)=kr3;

    % 『齿圈』基体变形引起的刚度计算
       uf3 = rf3-rb3/cos(Ak3+phik3);
       kf3cell=num2cell([Ai;Bi;Ci;Di;Ei;Fi]'*[1/(phif3)^2; h_fi3^2; h_fi3/phif3; 1/phif3; h_fi3; 1]);
       [LL3,MM3,PP3,QQ3]=deal(kf3cell{:});
       kf3=L3*E/(LL3*(uf3/sf3)^2+MM3*(uf3/sf3)+PP3*(1+QQ3*(tan(Ak3))^2));

       kkf3(k)=kf3;

    %% 刚度汇总
        kk=1/(1/kf3+1/kf2+1/kh23);
%         kk=1/(1/kb3+1/ks3+1/kr3+1/kf3+1/kb2+1/ks2+1/kr2+1/kf2+1/kh23);
        k23(k)=kk;
end