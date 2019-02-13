%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% last update 12Feb2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Here, you have to choose your material among the following %%%%%%%%%% 
%%%%%% Take care, not all the models are available for all the materials! %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Material='GaN';
%Material='AlN';
%Material='InN';
%Material='ZnO';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=300;                  % Temperature [Kelvin]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library
ExtractParameters_WZ

Nk=100;                   %% number of k-points for the dispersion
[k_list,k]=kWZ_f(Nk,a);   %% function to compute the k-vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Activate the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kp_6x6_Luttinger   = 1;  % EC, HH, LH, SO (The most used model)
kp_8x8_Luttinger   = 1;  % EC, HH, LH, SO (The most used model)

plot_Luttinger_parabole_small_k=0;
plot_Luttinger_parabole_large_k=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add the strain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exx=0.000; %% exx = (a0-a)/a0
ezz=-2*c13/c33*exx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;FS=20;
c=[
0 1 0
0 0 1
1 0 0
1 0 1
0 1 1
1 1 0
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kp_6x6_Luttinger  == 1
  i=i+1;
  E{i}=kp_6bands_DKK_WZ_strain_f(k_list, Dcr, Dso, AA, DD, exx, ezz);
  %E{i}=kp_6bands_WZ_strain_f(k_list, Dcr, Dso, AA, DD, exx, ezz);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 6x6: Luttinger');
end
if kp_8x8_Luttinger  == 1
  i=i+1;
  E{i}=kp_8bands_DKK_WZ_strain_f(k_list, Eg, EP, Dcr, Dso, me_z, me_xy, AA, ac1, ac2, DD, exx, ezz);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 8x8: Luttinger');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Luttinger_parabole_small_k == 1 || plot_Luttinger_parabole_large_k == 1;
  i=i+1;
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{0 0 0}Luttinger parabole');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[50 200 1000 900]);

FS=20;
LW=2;

%xscale=[-2 2];
%yscale=[-2 4];
xscale=[-1.201 1.201];
yscale=[-0.101 0.0601];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on;

xlabel('Wavevector (nm-1)')
ylabel('Energy (eV)')
xlim(xscale)
ylim(yscale)
set(gca,'ytick',-1:0.02:1)

text(0.75*(xscale(2)-xscale(1))+xscale(1),0.95*(yscale(2)-yscale(1))+yscale(1),s);
text(-0.03*(xscale(2)-xscale(1))+xscale(1),-0.06*(yscale(2)-yscale(1))+yscale(1),{'kx:[100]'},'fontsize',FS);
text(0.92*(xscale(2)-xscale(1))+xscale(1),-0.06*(yscale(2)-yscale(1))+yscale(1),{'kz:[001]'},'fontsize',FS);

title(strcat(M{1},' Wurtzite bandstructure @T=',num2str(T),'K'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(E)
  
  plot(k*1e-9,E{j},'linewidth',LW,'color',c(j,:))
  
end

if exx==0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}unstrained'));
elseif exx>0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}exx=',num2str(exx*100),'% tensiled'));
elseif exx<0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}exx=',num2str(exx*100),'% compressed'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Luttinger_parabole_small_k == 1 || plot_Luttinger_parabole_large_k == 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  h=6.62606896E-34;               %% Planck constant [J.s]
  hbar=h/(2*pi);
  e=1.602176487E-19;              %% electron charge [Coulomb]
  m0=9.10938188E-31;              %% electron mass [kg]
  H0=hbar^2/(2*m0) ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Ehh=Dcr+Dso/3
  Elh=(Dcr-Dso/3)/2 + sqrt( ((Dcr-Dso/3)/2)^2 + 2*(Dso/3)^2 )
  Ech=(Dcr-Dso/3)/2 - sqrt( ((Dcr-Dso/3)/2)^2 + 2*(Dso/3)^2 )

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  A1=AA(1);A2=AA(2);A3=AA(3);A4=AA(4);A5=AA(5);A6=AA(6);A7=AA(7)*e*1e-10;
  
  if plot_Luttinger_parabole_small_k == 1
    mhh001 = 1 / -(A1+A3) ;
    mhh100 = 1 / -(A2+A4) ;
    
    mlh001 = 1 / -(A1+(Elh/(Elh-Ech))*A3) ;
    mlh100 = 1 / -(A2+(Elh/(Elh-Ech))*A4) ;
    
    mch001 = 1 / -(A1+(Ech/(Ech-Elh))*A3) ;
    mch100 = 1 / -(A2+(Ech/(Ech-Elh))*A4) ;
  end
  if plot_Luttinger_parabole_large_k == 1
    mhh001 = 1 / -(A1+A3) ;
    mhh100 = 1 / -(A2+A4-A5) ;
    
    mlh001 = 1 / -(A1+A3) ;
    mlh100 = 1 / -(A2+A4+A5) ;
    
    mch001 = 1 / -A1 ;
    mch100 = 1 / -A2 ;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  f = 1/e*H0*k(1:Nk).^2;
  plot(k(1:Nk)*1e-9 ,  f/me_xy+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mhh100+Ehh ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mlh100+Elh ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mch100+Ech ,'linewidth',LW,'color','k','linestyle','--')
  
  f = 1/e*H0*k(Nk:end).^2;
  plot(k(Nk:end)*1e-9 ,  f/me_z+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mhh001+Ehh ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mlh001+Elh ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mch001+Ech ,'linewidth',LW,'color','k','linestyle','--')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%