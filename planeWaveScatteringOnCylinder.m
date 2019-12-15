%function planeWaveScatteringOnCylinder(varargin)
%

tic % measure execution time

% parameters
% ( d/2 * resolution/zoom ) should be integer!!!
d = 0.5; % [0.2 1 4 30]
zoom = 25; % [25 250]
resolution = 100;

R = d/2; % radius of the cylinder
frequency = 1e+09; % 1GHz -> lambda ~= 0.3m

disp(sprintf(' d = %.2f \n zoom = %d ', d, zoom));

% physical constants
C = 2.997924580003452e+08; % physconst('LightSpeed')
MU_O = 4*pi*1e-7;
EPS_O = 8.8541878176e-12;

% wave vector: k0 = 2πf√με
k0 = 2*pi * frequency * sqrt(MU_O*EPS_O);

%  coordinates
x = -zoom:zoom/resolution:zoom;
y = -zoom:zoom/resolution:zoom;
[xg,yg] = meshgrid(x,y);
zg = xg + i*yg;
r = abs(zg);
fi = angle(zg);

% calculate the field
result = exp(-i*k0*xg); % init with the incident wave
% calculate the zero component
result = result + ... 
      -((-i)^0)*besselj(0,k0*R)/besselh(0,2,k0*R) * ...
      besselh(0,2,k0*r).*exp(i*0*(fi));
n = 1;
while true
  result = result + ... 
      -((-i)^n)*besselj(n,k0*R)/besselh(n,2,k0*R) * ...
      besselh(n,2,k0*r).*exp(i*n*(fi));
  n = -n;
  result = result + ... 
      -((-i)^n)*besselj(n,k0*R)/besselh(n,2,k0*R) * ...
      besselh(n,2,k0*r).*exp(i*n*(fi));
  n = -n+1;
  if (abs(result(resolution,R*resolution/zoom+resolution)) < 0.1)
      break
  end
end

% mask out the cylinder
result = result.*(r>R);

% rendering the result
figure();
imagesc(x,y,abs(result));
hold on;
title(sprintf('f = %g Hz \n d = {%0.2f}m ({%.2f}\\lambda)', frequency, d, d/(C/frequency)));
xlabel('x [m]');
ylabel('y [m]');
axis('square');
%caxis([0 2]);
colorbar;
% draw the cylinder
fill(R*cos(theta), R*sin(theta),'w','LineStyle','none');
% export as image
%saveas(gcf, sprintf('plot_d%.2f_w%d.png', d, zoom), 'png');

toc % measure execution time

%end % function