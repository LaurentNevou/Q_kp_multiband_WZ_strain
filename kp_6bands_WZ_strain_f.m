function[E]=kp_6bands_WZ_strain_f(k_list, Dcr, Dso, AA, DD, exx, ezz)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
exy = 0; eyx=0;
ezx = 0; exz=0;
eyz = 0; ezy=0;
ee  = exx+eyy+ezz;

%l1 = D2+D4+D5;
%l2  = D1;
%m1 = D2+D4-D5;
%m2 = D1+D3;
%m3 = D2;
%n1 = 2*D5;
%n2 = sqrt(2)*D6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

% keep in mind
% Delta1=Dcr
% Delta2=Delta3=Dso/3

lambda = H0*( A1*kz^2 + A2*(kx^2+ky^2) ) + D1*ezz + D2*(exx+eyy);
theta  = H0*( A3*kz^2 + A4*(kx^2+ky^2) ) + D3*ezz + D4*(exx+eyy);
F      = Dcr + Dso/3 + lambda + theta ;
G      = Dcr - Dso/3 + lambda + theta ;
K      = H0*A5*( kx + 1i*ky )^2    + D5*( exx-eyy + 2i*exy);
H      = H0*A6*( kx + 1i*ky ) * kz + D6*( ezx + 1i*eyz );
Delta  = sqrt(2)*Dso/3 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hamiltonian WZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HH=[

    F      -K'      -H'       0       0      0
   -K       G        H        0       0    Delta
   -H       H'     lambda     0     Delta    0
    0       0        0        F      -K      H
    0       0      Delta     -K'      G     -H'
    0     Delta      0        H'     -H    lambda
    
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(HH)/e ;

end

end