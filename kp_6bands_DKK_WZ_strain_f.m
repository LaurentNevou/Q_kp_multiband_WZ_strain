function[E]=kp_6bands_DKK_WZ_strain_f(k_list, Dcr, Dso, AA, DD, exx, ezz)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dcr = Dcr*e;
Dso = Dso*e;

A1=AA(1);A2=AA(2);A3=AA(3);A4=AA(4);A5=AA(5);A6=AA(6);A7=AA(7)*e*1e-10;
D1=DD(1)*e;D2=DD(2)*e;D3=DD(3)*e;D4=DD(4)*e;D5=DD(5)*e;D6=DD(6)*e;

% here is the modification of the parameters due to the 6 bands instead of 8 bands model
%L1 = H0*(A2+A4+A5-1)+P1^2/Eg;
%L2 = H0*(A1-1)+P2^2/Eg;
L1 = H0*(A2+A4+A5-1);   
L2 = H0*(A1-1);

M1 = H0*(A2+A4-A5-1);
M2 = H0*(A1+A3-1);
M3 = H0*(A2-1);

% here is the modification of the parameters due to the 6 bands instead of 8 bands model
%N1 = H0*2*A5       + P1^2/Eg;
%N2 = H0*sqrt(2)*A6 + P1*P2/Eg;
N1 = H0*2*A5       ;
N2 = H0*sqrt(2)*A6 ;

N3 = 1i*sqrt(2)*A7;

B1=0;
B2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
exy = 0; eyx=0;
ezx = 0; exz=0;
eyz = 0; ezy=0;
ee  = exx+eyy+ezz;

l1 = D2+D4+D5;
l2  = D1;
m1 = D2+D4-D5;
m2 = D1+D3;
m3 = D2;
n1 = 2*D5;
n2 = sqrt(2)*D6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Strained Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% since it does not depend on k, it can be outside the loop %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hst=[

l1*exx+m1*eyy+m2*ezz          n1*exy                  n2*exz
       n1*exy          m1*exx+l1*eyy+m2*ezz           n2*eyz
       n2*exz                 n2*eyz             m3*(exx+eyy)+l2*ezz 

];

Hst = [Hst  zeros(3,3) ; zeros(3,3)  Hst];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

Hdiag = H0*k^2  *ones(1,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DKK Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_DKK=[

    L1*kx^2+M1*ky^2+M2*kz^2         N1*kx*ky             N2*kx*kz
           N1*kx*ky         M1*kx^2+L1*ky^2+M2*kz^2      N2*ky*kz
           N2*kx*kz                 N2*ky*kz        M3*(kx^2+ky^2)+L2*kz^2
];

H_DKK = [H_DKK  zeros(3,3) ; zeros(3,3)  H_DKK];

%%%%%%%%%%%%%%%%%%%%%%%%% Spin-orbit Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hso=[

 0   1   0   0   0   1i
-1   0   0   0   0   1
 0   0   0  -1i -1   0
 0   0  -1i  0  -1   0
 0   0   1   1   0   0
 1i -1   0   0   0   0

];

%%%%%%%%%%%%%%%%%%%%%%% Crystal field Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hcr=[   Dcr   Dcr   0      Dcr   Dcr   0 ];
Hcr=diag(Hcr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = diag(Hdiag) + H_DKK + Hso*Dso/(3i) + Hcr + Hst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end