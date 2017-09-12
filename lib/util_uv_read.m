function [u,v,w,R]=util_uv_read(stringname, nbr)
% from Arwa
% uvread reads realistic uv coverages
% Input:
% stringname   - u,v,w, cooordinates file
% nber            - number of uv points to keep
% Output:
% u,v,w            - u,v,w coordinates


c=299792458;
freq=720*1e6;% -- MeerKAT simulation
lambda=c/freq;
h=0; % angular hour -- MeerKAT simulation
dec=-30*pi/180; %declination -- MeerKAT simulation
uvw = load(stringname);
xyz = uvw.uvw;
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

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
% uvw=rot*xyz;

uvw = xyz;
bmax=max(sqrt(uvw(1,:).^2+uvw(2,:).^2)); %maximum baseline
% bmax = max(max(abs(uvw(1,:))),max(abs(uvw(2,:))));
uvw=pi*uvw./(bmax);
u=uvw(1,:)';
v=uvw(2,:)';
w=uvw(3,:)';
% u=[u;-u];
% v=[v;-v];
% w=[w;-w];
  
R=length(u);
end