function[E]=kp_8bands_DKK_WZ_f(k_list, Eg, EP, Dcr, Dso, me_z, me_xy, AA)

% Stefan Birner (Nextnano)
% PhD thesis: "Modeling of semiconductor nanostructures and semiconductor-electrolyte interfaces" (2011)
% Chapter3, page 36: "Multi-band k.p envelope function approximation"
% Download:
% https://mediatum.ub.tum.de/doc/1084806/1084806.pdf
% https://www.nextnano.com/downloads/publications/PhD_thesis_Stefan_Birner_TUM_2011_WSIBook.pdf

% S. L. Chuang et al. PRB, 54, 2491 (1996)
% "k.p method for strained wurtzite semiconductors"

% M. Kumagai, S. L. Chuang et al. PRB, 57, 15304 (1998)
% "Analytical solutions of the block-diagonalized Hamiltonian for strained wurtzite semiconductors"

% Seoung-Hwan Park and S. L. Chuang, PRB, 59, 4726 (1999)
% "Crystal-orientation effects on the piezoelectric field and electronic properties of strained wurtzite semiconductors"

% D. J. Dugdale, S. Brand, and R. A. Abram, PRB, 61, 12934 (2000)
% "Direct calculation of k"p parameters for wurtzite AlN, GaN, and InN"

% W. J. Fan et al. J. Crystal Growth, 287, 28 (2006)
% "Electronic structures of wurtzite ZnO and ZnO/MgZnO quantum well"

% P. Rinke, PRB, 77, 075202 (2008)
% "Consistent set of band parameters for the group-III nitrides AlN, GaN, and InN"

% http://cmt.dur.ac.uk/ssterratum/paper/node4.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dcr = Dcr*e;
Dso = Dso*e;
Eg  = Eg*e;
EP1z  = EP(1)*e;
EP2xy = EP(2)*e;
P1    = sqrt(EP1z *hbar^2/(2*m0)) ;
P2    = sqrt(EP2xy*hbar^2/(2*m0)) ;

me_z  = me_z*m0;
me_xy = me_xy*m0;

A1=AA(1);A2=AA(2);A3=AA(3);A4=AA(4);A5=AA(5);A6=AA(6);A7=AA(7)*e*1e-10;

% m{//}  = mz  parallel to the c-axis (same as z)
% m{_|_} = mxy perpendicular to the c-axis

A1z  = hbar^2/2*(1/me_z-1/m0)-P1^2/Eg;    % mass of electron in CB
A2xy = hbar^2/2*(1/me_xy-1/m0)-P2^2/Eg;  % mass of electron in CB

L1 = H0*(A2+A4+A5-1)+P1^2/Eg;
L2 = H0*(A1-1)+P2^2/Eg;

M1 = H0*(A2+A4-A5-1);
M2 = H0*(A1+A3-1);
M3 = H0*(A2-1);

N1 = H0*2*A5       + P1^2/Eg;
N2 = H0*sqrt(2)*A6 + P1*P2/Eg;

N3 = 1i*sqrt(2)*A7;

B1=0;B2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
Hdiag = H0*k^2*[1 1 1 1]  +  [Eg 0 0  0 ];

H4=[

  0      1i*P2*kx    1i*P2*ky   1i*P1*kz
  0         0           0          0
  0         0           0          0
  0         0           0          0

];

HH4 = H4' + H4 + diag(Hdiag); 

H_DKK=[

A2xy*(kx^2+ky^2)+A1z*kz^2       B2*ky*kz                B2*kx*kz             B1*kx*ky
        B2*ky*kz        L1*kx^2+M1*ky^2+M2*kz^2         N1*kx*ky             N2*kx*kz
        B2*kx*kz                N1*kx*ky         M1*kx^2+L1*ky^2+M2*kz^2     N2*ky*kz
        B1*kx*ky                N2*kx*kz                N2*ky*kz        M3*(kx^2+ky^2)+L2*kz^2

];


Hso=[

 0   0   0   0   0   0   0   0
 0   0   1   0   0   0   0   1i
 0  -1   0   0   0   0   0   1
 0   0   0   0   0  -1i -1   0
 0   0   0   0   0   0   0   0
 0   0   0  -1i  0   0  -1   0
 0   0   0   1   0   1   0   0
 0   1i -1   0   0   0   0   0

];

Hcr=[0   Dcr   Dcr   0   0   Dcr   Dcr   0 ];
Hcr=diag(Hcr);

H=[HH4  zeros(4,4) ; zeros(4,4)  HH4] + [H_DKK  zeros(4,4) ; zeros(4,4)  H_DKK] + Hso*Dso/(3i) + Hcr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end