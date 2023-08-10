clc;     %% 清除不了全局变量，只能清除普通变量
clear all; %% 清除所有的变量，包括全局变量global
close all;   %% 关闭所有窗口
%%
step = 0.0001; %% 步长
n =200000;    %% 取点个数
t = 0:step:step*n; %% 总仿真时长
%% 初始条件
x1=0.5;
x2=-1;
hat_theta_1=0;
hat_theta_2=0;
%% design parameters
gama1=0.5;gama2 =0.5;sigma1=1;sigma2 =1;
k1=5;k2=10;T=2;
eta0=1;
eta =0.08;
%% neural  network node
b =1;hidden1 = 89;

for i=1:1:n
    tt= i*step;
    %% 扰动
        d1 =-0.5*exp(-x1^2);
        d2 = 0.02*cos(x2*tt);
    %% 跟踪轨迹
    xd =0.5*sin(tt);
    xds(i) =xd;
    dot_xd =0.5*cos(tt);
    %% Scalar time-varying function 
	if tt<T
        beta =(eta0-eta)*((T-tt)/T)^3+eta;
        dot_beta=-3*(1-eta)*(T-tt)^2/(T^3); 
    else 
        beta=eta;
        dot_bata=0;
    end;
    beta1s(i)=beta;
        beta2s(i)=-beta;
%% error transfer
    e = x1-xd;
    e1s(i)=e;
%%    dot_e = dot_x1 - dot_xd;
%%    dot_zeta_1 =P *(alfa_1+zeta_2+d1-(dot_beta/beta)*e);  
   %% neural networks 1    
      Z1 = [x1;dot_beta;beta;xd];
      dd1 = [-11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 10.25 10.5 10.75 11;
           -11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 10.25 10.5 10.75 11;
           -11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 10.25 10.5 10.75 11;
           -11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 10.25 10.5 10.75 11];
for j = 1:hidden1
    varphi_1(j) = exp(-((norm(Z1-dd1(:,j)))^2)/(2*b^2));
end
Phi_1=norm(varphi_1)^2+1;
 %% define P
 P=((1+tan((pi*e)/(2*beta)))^2)*(pi/(2*beta));  
 zeta_1=tan((pi*e)/(2*beta));
%% alfa1
alfa_1=-(1/P)*k1*zeta_1-zeta_1*hat_theta_1*Phi_1;
dot_hat_theta_1 = gama1*P*zeta_1^2*Phi_1-sigma1*hat_theta_1;
alfa_1s(i) = alfa_1;
%% step2
 zeta_2 = x2-alfa_1;

 %% neural networks 2   
      Z2 = [x1;x2;dot_hat_theta_1];
      dd2 = [-11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 1.25 10.5 10.75 11;
                    -11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 1.25 10.5 10.75 11;
          -11 -10.75 -10.5 -10.25 -10 -9.75 -9.5 -9.25 -9 -8.75 -8.5 -8.25 -8 -7.75 -7.5 -7.25 -7 -6.75 -6.5 -6.25 -6 -5.75 -5.5 -5.25 -5 -4.75 -4.5 -4.25 -4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5  9.75 10 1.25 10.5 10.75 11];
for j = 1:hidden1
    varphi_2(j) = exp(-((norm(Z2-dd2(:,j)))^2)/(2*b^2));
end
Phi_2 =norm(varphi_2)^2+1;
%% u
u =-k2*zeta_2-zeta_2*hat_theta_2*Phi_2-P*zeta_1;
dot_hat_theta_2 = gama2*zeta_2^2*Phi_2-sigma2*hat_theta_2;
us(i)= u;
%% update law
hat_theta1s(i)=hat_theta_1;
hat_theta2s(i)=hat_theta_2;
hat_theta_1=hat_theta_1+dot_hat_theta_1*step;
hat_theta_2=hat_theta_2+dot_hat_theta_2*step;
%% plant
x1 = x1 +(cos(x1)+x2+d1)*step;
x2 = x2 +(cos(x1*x2)+u+d2)*step;
x1s(i)=x1;
x2s(i)=x2;
end
figure(1);
plot(t(1:n),e1s,'r',t(1:n),beta1s,'g--',t(1:n),beta2s,'b--','linewidth',2);
legend('e_{1}','boundary curve','boundary curve');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$e_1$}','Interpreter','latex','FontSize',14);


figure(2);
plot(t(1:n), us,'b','linewidth',1.5);
legend('u_{1}');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$u_1$}','Interpreter','latex','FontSize',14);


figure(3);
plot(t(1:n), hat_theta1s,'r', 'linewidth',2);
legend({'$\hat\theta_1$'},'Interpreter','latex','FontSize',14)
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$||\hat \theta_1||$}','Interpreter','latex','FontSize',14);


figure(4);
plot(t(1:n), x1s,'r', t(1:n), x2s,'b--', 'linewidth',2);
legend('x_{1}','x_{2}');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$x_1,x_2$}','Interpreter','latex','FontSize',14);


figure(5);
plot(t(1:n),x1s,'r',t(1:n),xds,'b--','linewidth',2);
legend('x_{1}','x_{d}');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$x_1,x_d$}','Interpreter','latex','FontSize',14);


figure(6);
plot(t(1:n),  hat_theta2s,'b--', 'linewidth',2);
legend({'$\hat\theta_2$'},'Interpreter','latex','FontSize',14)
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
ylabel('{$||\hat \theta_2||$}','Interpreter','latex','FontSize',14);

