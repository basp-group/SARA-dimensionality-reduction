function [u,v,w,R]=util_uv_read(stringname, nbr)
% from Arwa
% uvread reads realistic uv coverages
% Input:
% stringname   - u,v,w, cooordinates file
% nber            - number of uv points to keep
% Output:
% u,v,w            - u,v,w coordinates


c=3E+08;
freq=1.01458E+09;% -- MeerKAT simulation
lambda=c/freq;
h=0; % angular hour -- MeerKAT simulation
dec=-40*pi/180; %declination -- MeerKAT simulation
COV = importdata(stringname);
xyz=COV.data; % coordinates in units of meter
x = xyz(2:end,1);
y = xyz(2:end,2);
z = xyz(2:end,3);

% get rid of 0 components
xdummy=x(abs(x)+abs(y)>0);
ydummy=y(abs(x)+abs(y)>0);
zdummy=z(abs(x)+abs(y)>0); 
xyz=[xdummy ydummy zdummy];
sz = size(zdummy, 1);
if nargin==2
   index=randi([1 sz],1,nbr);
   x=xyz(index,1);
   y=xyz(index,2);
   z=xyz(index,3);    
end

xyz=[x'; y';z'];
rot=(1/lambda).*[sin(h) cos(h) 0;...
    -sin(dec)*cos(h) sin(dec)*sin(h) cos(dec);...
    cos(dec)*cos(h) -cos(dec)*sin(h) sin(dec)];
% convert baselines from meters to units of wavelengths
uvw=rot*xyz;
bmax=max(sqrt(uvw(1,:).^2+uvw(2,:).^2)); %maximum baseline
uvw=pi*uvw./bmax;
u=uvw(1,:)';
v=uvw(2,:)';
w=uvw(3,:)';
u=[u;-u];
v=[v; -v];
w=[w;-w];
  
R=length(u);