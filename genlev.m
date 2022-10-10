function lev = genlev(which,nz,dz,dzmin,a)

switch which

  case 'uniform'
    disp('Uniform level distribution');
    if (nargin < 3)
      error('need to supply nz,dz as arguments');
    end
    if (nargin > 3)
      disp('NOTE: additional arguments ignored');
    end
    z=(0:nz-1)'*dz;
    dz=z(2:end)-z(1:end-1);

    lev.description='Uniform level distribution';
    lev.formula='dz(k)=dz';
    lev.dz=dz;

  case 'cosmo2'
    disp('COSMO-2 level distribution');
    if (nargin > 1)
      disp('NOTE: additional arguments ignored');
    end
    z=importdata('levels_60.2.txt');
    dz=z(2:end)-z(1:end-1);
    nz=numel(z);

    lev.description='COSMO-2 level distribution';
    lev.formula='z(z)=vcoord(k)';
    
  case 'cosmo-de'
    disp('COSMO-DE level distribution');
    if (nargin > 1)
      disp('NOTE: additional arguments ignored');
    end
    z=importdata('levels_cosmode.txt');
    dz=z(2:end)-z(1:end-1);
    nz=numel(z);

    lev.description='COSMO-DE level distribution';
    lev.formula='z(z)=vcoord(k)';

  case 'cosmo-de-65'
    disp('COSMO-DE level distribution');
    if (nargin > 1)
      disp('NOTE: additional arguments ignored');
    end
    z=importdata('levels_cosmode_new1.txt');
    dz=z(2:end)-z(1:end-1);
    nz=numel(z);

    lev.description='COSMO-DE level distribution';
    lev.formula='z(z)=vcoord(k)';

  case 'cosmo-de-80'
    disp('COSMO-DE level distribution');
    if (nargin > 1)
      disp('NOTE: additional arguments ignored');
    end
    z=importdata('levels_cosmode_new2.txt');
    dz=z(2:end)-z(1:end-1);
    nz=numel(z);

    lev.description='COSMO-DE level distribution';
    lev.formula='z(z)=vcoord(k)';

  case 'L60.1'
    disp('L60.1 level distribution');
    if (nargin > 1)
      disp('NOTE: additional arguments ignored');
    end
    z=importdata('levels_60.1.txt');
    dz=z(2:end)-z(1:end-1);
    nz=numel(z);

    lev.description='L60.1 level distribution';
    lev.formula='z(z)=vcoord(k)';

  case 'power'
    disp('Power law level distribution');
    if (nargin < 5)
      error('need to supply nz,dz,dzmin,a as arguments');
    end
    if (nargin > 5)
      disp('NOTE: additional arguments ignored');
    end
    dz=dz+(dz-dzmin)*(2*(0:nz-2)'/(nz-2)-1).^a;
    z=cumsum([0; dz]);

    lev.deescription='Power law level distribution';
    lev.formula='dz(z)=dz+(dz-dzmin)*(2*(k-1)/(nk-1)-1).^a';
    lev.dz=dz;
    lev.dzmin=dzmin;

  case 'tanh'
    disp('Tanh level distribution');
    if (nargin < 5)
      error('need to supply nz,dz,dzmin,a as arguments');
    end
    if (nargin > 5)
      disp('NOTE: additional arguments ignored');
    end
    dz=dz+(dzmin-dz)/tanh(2*a)*tanh(2*a*((1:nz-1)'-nz/2)/(1-nz/2));
    z=cumsum([0; dz]);

    lev.deescription='Tanh level distribution';
    lev.formula='dz(k)=dz+(dzmin-dz)/tanh(2*a)*tanh(2*a*((k-(nk+1)/2)/(1-(nk+1)/2))';
    lev.dz=dz;
    lev.dzmin=dzmin;
    lev.a=a;

  case 'dwd'
    disp('DWD level distribution');
    if (nargin < 3)
      error('need to supply nz,a as arguments');
    end
    if (nargin > 3)
      disp('NOTE: additional arguments ignored');
    end
    a = dz;
    dz=(0:(nz-2))'/(nz-2);
    dz = dz.^a;
    dz = dz./sum(dz)*(22000.0-(nz-1)*20.0)+20.0;
    z=cumsum([0; dz]);

    lev.deescription='DWD level distribution';
    lev.formula='';
    lev.dz=dz;
    lev.a=a;

  otherwise
    error('undefined level distribution chosen');

end

% fill standard values
lev.name=which;
lev.nz=nz;
lev.dz=dz;
lev.z=z;
lev.k=1:lev.nz;

end
