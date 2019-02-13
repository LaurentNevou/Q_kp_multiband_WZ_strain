function[k_list,k]=kWZ_f(Nk,a)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Brillouin zone vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% x-valley
kL=linspace(-1,0,Nk);
k1 = [ kL'  0*kL'  0*kL' ];

%%% z-valley
kx=linspace(0,1,Nk);
k2 = [ 0*kx'  0*kx'  1*kx' ];


k_list=[ k1 ; k2  ]*2*pi/a * 0.1;

k=sqrt( k_list(:,1).^2 + k_list(:,2).^2 + k_list(:,3).^2 );
k=[ -k(1:Nk) ; k(Nk+1:end)  ];

end