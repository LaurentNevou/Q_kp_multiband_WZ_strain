%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Extract general parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lattice parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'a');
  if idx==1
    a=M{2}(i-1)*1e-10;
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP1z and EP2xy from the Kane model

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'EP1z');
  if idx==1
    EP = [M{2}(i-1) M{2}(i)];
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6c (the band gap in a direct band gap semiconductor)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg6cG');
  if idx==1
    Eg6c = M{2}(i-1);
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "alpha" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'alphaG');
  if idx==1
    alphaG = M{2}(i-1);
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "beta" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'betaG');
  if idx==1
    betaG = M{2}(i-1);
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eg   = Eg6c - (alphaG*T^2) ./ (T+betaG);   %Eg = Eg0 - (a*T.^2)./(T + b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crystal field splitting DeltaCR

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Dcr');
  if idx==1
    Dcr(1) = M{2}(i-1);
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split-off band energy Delta-so

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Dso');
  if idx==1
    Dso(1) = M{2}(i-1);
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the longitudinal electrons (z)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'me_z');
  if idx==1
    me_z = [M{2}(i-1) ];
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the transverse electrons (xy)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'me_xy');
  if idx==1
    me_xy = [M{2}(i-1) ];
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the holes

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'A1');
  if idx==1
    AA = [M{2}(i-1:i+5) ];
    break  
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%