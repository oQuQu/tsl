clear all;
clc;
step=0.0001;
n=400000;
t=0:step:step*n;
% system parameters
p1=1;
p2=1;
% initial conditions
x1=0.2;
x2=-0.3;
hat_b1=0;
hat_b2=0;
% design parameters
c1=15;
c2=20;
gamma1=0.01;
gamma2=0.01;
sigma1=0.1;
sigma2=0.1;

T0=10;

for i=1:1:n
    tt=i*step;
    % time varying constraint functions
    F11=0.6-0.8*sin(tt);
    F12=0.9+0.1*cos(tt);
    F21=1-0.3*sin(tt);
    F22=0.4+0.1*cos(tt);
    % the derivative of time-varying constraint functions
    dF11=-0.8*cos(tt);
    dF12=-0.1*sin(tt);
    ddF11=0.8*sin(tt);
    ddF12=-0.1*cos(tt);
    dF21=-0.2*cos(tt);
    dF22=-0.1*sin(tt);
    
    F11s(i)=-F11;
    F12s(i)=F12;
    F21s(i)=-F21;
    F22s(i)=F22;
        %% �Ŷ�
        d1 =-0.5*exp(-x1^2);
        d2 = 0.02*cos(x2*tt);
   
    % the boundary of time-varying constraints
    uF11=-0.3;
    uF12=0.5;
    uF21=0.4;
    uF22=0.2;
    
    % desired signal
    yd=0.5*sin(tt);
    yds(i)=yd;
    dyd=0.5*cos(tt);
    ddyd=-0.5*sin(tt);
    
    % nonlinear transformation
    zeta1=(uF11+x1)/(F11+x1)+(x1-uF12)/(F12-x1);
    zeta2=(uF21+x2)/(F21+x2)+(x2-uF22)/(F22-x2);   
    alpha0=(uF11+yd)/(F11+yd)+(yd-uF12)/(F12-yd);
    
    % virtual error
    z1=zeta1-alpha0;
    z1s(i)=z1;
    % tracking error
    e=x1-yd;
    es(i)=e;
    
    % define some variables
    eta01=(F11-uF11)/((F11+yd)^2)+(F12-uF12)/((F12-yd)^2);
    eta02=-dF11*(uF11+yd)/((F11+yd)^2)-dF12*(yd-uF12)/((F12-yd)^2);
    eta11=(F11-uF11)/((F11+x1)^2)+(F12-uF12)/((F12-x1)^2);
    eta12=-dF11*(uF11+x1)/((F11+x1)^2)-dF12*(x1-uF12)/((F12-x1)^2);
    eta21=(F21-uF21)/((F21+x2)^2)+(F22-uF22)/((F22-x2)^2);
    eta22=-dF21*(uF21+x2)/((F21+x2)^2)-dF22*(x2-uF22)/((F22-x2)^2);
    zeta23=(uF21*F22-F21*uF22)/(F21-uF21+F22-uF22);
    ell1=eta12-eta01*dyd-eta02;
    
     % core functions
    varphi1=2*x1^2+0.25;
    varphi2=x1^2*x2^2+x2^2+0.5;
    Phi11=eta11^2*varphi1^2;
    Phi12=eta11^2*zeta23^2;
    Phi13=ell1^2;    
    Phi1=Phi11+Phi12+Phi13;
   
    % partial derivatives of Phi1
    pPhi11eta11=2*eta11*varphi1^2;
    pPhi11varphi1=2*varphi1*eta11^2;
    pvarphi1x1=4*x1;
    pPhi12eta11=2*eta11*zeta23^2;
    pPhi12zeta23=2*zeta23*eta11^2;
    pPhi13ell1=2*ell1;
    
    % virtual controller 
    Gamma1=c1*z1+hat_b1*z1*Phi1;
    alpha1=-Gamma1/eta11;
    alpha1s(i)=alpha1;
    % adpative law for hat_b1
    hat_b1s(i)=hat_b1;
    dot_hatb1=gamma1*z1^2*Phi1-sigma1*hat_b1;
    hat_b1=hat_b1+dot_hatb1*step;
    
    
    % partial derivatives
    pa1eta11=Gamma1/(eta11^2);
    pa1Gam1=-1/eta11;
    pGamma1Phi1=hat_b1*z1;
    pGamma1z1=c1+hat_b1*Phi1;
    pz1zeta1=1;
    pz1alpha0=-1;
    %partial{zeta1}
    pzeta1x1=(F11-uF11)/((F11+x1)^2)+(F12-uF12)/((F12-x1)^2);
    pzeta1F11=-(uF11+x1)/((F11+x1)^2);
    pzeta1F12=-(x1-uF12)/((F12-x1)^2);
   % partial{alpha0}
    palpha0yd=(F11-uF11)/((F11+yd)^2)+(F12-uF12)/((F12-yd)^2);
    palpha0F11=-(uF11+yd)/((F11+yd)^2);
    palpha0F12=-(yd-uF12)/((F12-yd)^2);
   % partial{eta01}
    peta01yd=2*(F12-uF12)/((F12-yd)^3)-2*(F11-uF11)/((F11+yd)^3);
    peta01F11=(F11+yd-2*(F11-uF11))/((F11+yd)^3);
    peta01F12=(F12-yd-2*(F12-uF12))/((F12-yd)^3);
  %  partial{eta02}
    peta02yd=-(dF11*(F11+yd)-2*dF11*(uF11+yd))/((F11+yd)^3)-(dF12*(F12-yd)+2*dF12*(yd-uF12))/((F12-yd)^3);
    peta02F11=2*dF11*(uF11+yd)/((F11+yd)^3);
    peta02F12=2*dF12*(yd-uF12)/((F12-yd)^3);
    peta02dF11=-(uF11+yd)/((F11+yd)^2);
    peta02dF12=-(yd-uF12)/((F12-yd)^2);       
 %   partial{eta11}
    peta11x1=2*(F12-uF12)/((F12-x1)^3)-2*(F11-uF11)/((F11+x1)^3);
    peta11F11=(F11+x1-2*(F11-uF11))/((F11+x1)^3);
    peta11F12=(F12-x1-2*(F12-uF12))/((F12-x1)^3);
 %   partial{eta12}    
    peta12x1=-(dF11*(F11+x1)-2*dF11*(uF11+x1))/((F11+x1)^3)-(dF12*(F12-x1)+2*dF12*(x1-uF12))/((F12-x1)^3);
    peta12F11=2*dF11*(uF11+x1)/((F11+x1)^3);
    peta12F12=2*dF12*(x1-uF12)/((F12-x1)^3);
    peta12dF11=-(uF11+x1)/((F11+x1)^2);
    peta12dF12=-(x1-uF12)/((F12-x1)^2);   
    
    % partial{Phi1}/partial{x1}
    pPhi11x1=pPhi11eta11*peta11x1+pPhi11varphi1*pvarphi1x1;
    pPhi12x1=pPhi12eta11*peta11x1;
    pPhi13x1=pPhi13ell1*peta12x1;
    pPhi1x1=pPhi11x1+pPhi12x1+pPhi13x1;
    pz1x1=pz1zeta1*pzeta1x1;
    pGamma1x1=pGamma1z1*pz1x1+pGamma1Phi1*pPhi1x1;
     % partial{alpha1}/partial{x1}   
    palpha1x1=pa1eta11*peta11x1+pa1Gam1*pGamma1x1;
    
    % partial{alpha1}/partial{yd}
    pz1yd=-palpha0yd;
    pell1yd=-dyd*peta01yd-peta02yd;
    pPhi13yd=pPhi13ell1*pell1yd;
    pPhi1yd=pPhi13yd;
    pGamma1yd=pGamma1z1*pz1yd+pGamma1Phi1*pPhi1yd;
    
    palpha1yd=pa1Gam1*pGamma1yd;

    % partial{alpha1}/partial{F11}
    pPhi11F11=pPhi11eta11*peta11F11;
    pPhi12F11=pPhi12eta11*peta11F11;
    pell1F11=peta12F11-dyd*peta01F11-peta02F11;
    pPhi13F11=pPhi13ell1*pell1F11;
    pPhi1F11=pPhi11F11+pPhi12F11+pPhi13F11;
    pz1F11=pz1zeta1*pzeta1F11+pz1alpha0*palpha0F11;
    pGamma1F11=pGamma1z1*pz1F11+pGamma1Phi1*pPhi1F11;
    
    palpha1F11=pa1eta11*peta11F11+pa1Gam1*pGamma1F11;
    
      % partial{alpha1}/partial{F12}
    pPhi11F12=pPhi11eta11*peta11F12;   
    pPhi12F12=pPhi12eta11*peta11F12;   
    pell1F12=peta12F12-dyd*peta01F12-peta02F12;  
    pPhi13F12=pPhi13ell1*pell1F12;
    pPhi1F12=pPhi11F12+pPhi12F12+pPhi13F12;  
    pz1F12=pz1zeta1*pzeta1F12+pz1alpha0*palpha0F12;  
    pGamma1F12=pGamma1z1*pz1F12+pGamma1Phi1*pPhi1F12;  
    
    palpha1F12=pa1eta11*peta11F12+pa1Gam1*pGamma1F12;
    
    % partial{alpha1}/partial{F21}
     pzeta23F21=(-F22*(F21-uF21+F22-uF22)-(uF21*F22-F21*uF22))/((F21-uF21+F22-uF22)^2);
     pzeta23F22=(uF21*(F21-uF21+F22-uF22)-(uF21*F22-F21*uF22))/((F21-uF21+F22-uF22)^2);
     pPhi12zeta23=2*zeta23*eta11^2;
     pPhi1F21=pPhi12zeta23*pzeta23F21;
     pGamma1F21=pGamma1Phi1*pPhi1F21;
     
     palpha1F21=pa1Gam1*pGamma1F21;
     
     % partial{alpha1}/partial{F22}
     pPhi1F22=pPhi12zeta23*pzeta23F22;
     pGamma1F22=pGamma1Phi1*pPhi1F22;
     
     palpha1F22=pa1Gam1*pGamma1F22;
     
     % partial{alpha1}/partial{dyd}
     pell1dyd=-eta01;
     pPhi1dyd=pPhi13ell1*pell1dyd;
     pGamma1dyd=pGamma1Phi1*pPhi1dyd;
     
     palpha1dyd=pa1Gam1*pGamma1dyd;
     
     % partial{alpha1}/partial{dF11}
     pell1dF11=peta12dF11-peta02dF11;
     pPhi1dF11=pPhi13ell1*pell1dF11;
     pGamma1dF11=pGamma1Phi1*pPhi1dF11;
     
     palpha1dF11=pa1Gam1*pGamma1dF11;

     % partial{alpha1}/partial{dF12}
     pell1dF12=peta12dF12-peta02dF12;
     pPhi1dF12=pPhi13ell1*pell1dF12;
     pGamma1dF12=pGamma1Phi1*pPhi1dF12;
     
     palpha1dF12=pa1Gam1*pGamma1dF12;
     
     % partial{alpha1}/partial{hat_b1}
     palpha1hatb1=-z1*Phi1/eta11;
     
     % beta1
     beta1=palpha1yd*dyd+palpha1dyd*ddyd+palpha1hatb1*dot_hatb1+palpha1F11*dF11+palpha1dF11*ddF11+palpha1F12*dF12+palpha1dF12*ddF12+palpha1F21*dF21+palpha1F22*dF22;
     % ell2
     ell2=eta22-beta1;
     % Phi2
     Phi2=eta21^2*varphi2^2+ell2^2+(palpha1x1*varphi1)^2+(palpha1x1*x2)^2;
     % z2
     z2=zeta2-alpha1;
     z2s(i)=z2;
     % controller u
     u=-(c2*z2+hat_b2*z2*Phi2)/eta21;
     us(i)=u;
     %  adaptive law for hat_b2
     dot_hatb2=gamma2*z2^2*Phi2-sigma2*hat_b2;
     hat_b2s(i)=hat_b2;
     hat_b2=hat_b2+dot_hatb2*step;
     
     f1=cos(x1)+d1;
     f2=cos(x1*x2)+d2;
     
     dx1=x2+f1;
     dx2=u+f2;
     
     x1s(i)=x1;
     x2s(i)=x2;
     
     x1=x1+dx1*step;
     x2=x2+dx2*step;
     
end
% 
figure(1);
plot(t(1:n), es,'g','linewidth',2);
legend('z_{1}');
xlabel('Time(s)','FontName','Times New Roman','FontSize',14);
ylabel('{$z_1$}','Interpreter','latex','FontSize',14);
% 
figure(2);
plot(t(1:n), us,'g','linewidth',2);
legend('u_{1}');
xlabel('Time(s)','FontName','Times New Roman','FontSize',14);
ylabel('{$u_1$}','Interpreter','latex','FontSize',14);


figure(3);
plot(t(1:n),x1s,'r--','linewidth',2);
legend('x_{1}');
xlabel('Time(s)','FontName','Times New Roman','FontSize',14);
ylabel('{$x_1$}}','Interpreter','latex','FontSize',14);

% figure(1);
% clear all;
% load('v1');
% E=plot(t(1:n),us,'r','linewidth',2);
% E=plot(t(1:n),x1s,'r',t(1:n),xds,'b--','linewidth',2);
% hold on
% clear all;
% load('zhao');
% E=plot(t(1:n),x1s,'g','linewidth',2);
% legend('$x_{1}$','$x_{d}$','[26]','Interpreter','latex');
% hold on
% xlabel('Time(sec)','Interpreter','latex','FontName','Times New Roman','FontSize',16);
% ylabel('Tracking performance','Interpreter','latex','FontName','Times New Roman','FontSize',16);
% set(get(gca,'YLabel'),'FontName','Times New Roman','FontSize',16);
% set(get(gca,'TITLE'),'FontName','Times New Roman','FontSize',16);
% set(gca,'fontsize',16);
% 
% axes('Position',[0.2;0.25;0.2;0.25]); % 生成子图,1200�?600分别对应6秒和8秒时�?
% clear all;
% load('v1');
% plot(t(250:2000),us(250:2000),'g','LineWidth',2)
% hold on
% clear all;
% load('zhaokai');
% plot(t(250:2000),us(250:2000),'b--','LineWidth',2)
% set(gca,'YLim',[-50 200]);
% set(gca,'XLim',[0.25 2]);
% grid on

figure(1);
plot(t(1:n),F11s,'r',t(1:n),F12s,'g',t(1:n),x1s,'b',t(1:n),yds,'m--','linewidth',2);
legend({'$-{F}_{11}(t)$','$F_{12}(t)$','$x_1(t)$','$y_d(t)$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);

figure(2);
plot(t(1:n),F21s,'r',t(1:n),F22s,'g',t(1:n),x2s,'b','linewidth',2);
legend({'$-{F}_{21}(t)$','$F_{22}(t)$','$x_2(t)$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);

figure(3);
subplot(211);
plot(t(1:n),alpha1s,'r','linewidth',2);
legend({'$\alpha_1$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
subplot(212);
plot(t(1:n),us,'r','linewidth',2);
legend({'$u$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);

figure(4);
plot(t(1:n),es,'r','linewidth',2);
legend({'$e=x_1-y_d$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);

figure(5);
plot(t(1:n),hat_b1s,'r',t(1:n),hat_b2s,'b--','linewidth',2);
legend({'$\hat b_1$','$\hat b_2$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);

figure(6);
subplot(211);
plot(t(1:n),z1s,'r','linewidth',2);
legend({'$z_1$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
subplot(212);
plot(t(1:n),z2s,'r','linewidth',2);
legend({'$z_2$'},'Interpreter','latex');
xlabel('Time(sec)','FontName','Times New Roman','FontSize',14);
