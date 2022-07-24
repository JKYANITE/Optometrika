
function TransComp()
%
%  analysis of light transmittance of glass fibre of cylindrical section
%
%  Zheng Ji, 11.07.2022

clear all
close all
clc

startup;

% initialisation

origin_glass_fibre = 50; % origin of cylinder
diameter_glass_fibre = 10; % diameter of glass fibre in (mu)metre for macroscopic demonstration

length_glass_fibre = 20; % cylinder height
entry_polymer_matrix = 'pmma'; % refractive index of polymer matrix
entry_glass_fibre = 'bk7';   % refractive index of glass fibre

diameter_glass_fibre_iba = diameter_glass_fibre*1e-6; % diameter of glass fibre in micrometre
radius_glass_fibre_iba = diameter_glass_fibre_iba/2;
width_model = 5e-6;

number_glass_fibre = 9
%volume_fraction = number_glass_fibre*pi*(refrindx( wavelength,entry_glass_fibre ))^2/(w*h)

% definition of spectrum range
wavelength_start = 380;
wavelength_end = 800;

% create a container for optical elements (Bench class)
bench = Bench;
ref_bench = Bench;

%subfunction of fibre/matrix alignment
polymer_fibre_setup(origin_glass_fibre, diameter_glass_fibre, entry_polymer_matrix, bench, length_glass_fibre, entry_glass_fibre); 

% screen
screen = ScreenGeneric( [ origin_glass_fibre+4*diameter_glass_fibre 0 0 ], 10 * diameter_glass_fibre, 10 * diameter_glass_fibre, 256, 256, 'wf' );
bench.append( screen );
ref_bench.append( screen );

% create divergent rays
nrays = 100;
d = -5; % distance between rays of different wavelengths
rays_width = 30; % width of collimated rayNs
rays_dir = [1 0 0]; %rays incidence angle



democase_rays = ['spectrum'];
switch(democase_rays)

    % demonstrate spectrum of 380nm-700nm in the unit of 1nm
    case 'spectrum'
        for wavelength = wavelength_start:wavelength_end
            
            n_f = 1.4794+3920.*wavelength^(-2)-808.*wavelength^(-4)+0.1.*wavelength^(-6);
           

            rays_bundle(1,wavelength-(wavelength_start-1))=...
                Rays( nrays, 'collimated', [ 0 0 -10+0.05*(wavelength-379) ], rays_dir, rays_width, 'linear', 'air', ...
                      wavelength*1e-9, 1/255.*wavelength22color(wavelength, 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb') );
            
            wavelength_local = wavelength.*1e-9;
            %phase_lag = 4.*pi.*radius_glass_fibre_iba.*abs(refrindx(wavelength_local,entry_glass_fibre)-refrindx(wavelength_local,entry_polymer_matrix))/(wavelength_local);
            phase_lag = 4.*pi.*radius_glass_fibre_iba.*abs(n_f-refrindx(wavelength_local,entry_polymer_matrix))/(wavelength_local);

            fun = @(gamma) 2.*phase_lag.*sin(phase_lag.*cos(gamma)).*(sin(gamma)).^2;
            efficiency_factor = integral(fun,0,pi/2);

            geometrical_shadow_ratio = 2.*radius_glass_fibre_iba/width_model; % Geomatrical shadow ratio
            transmittance = 1-geometrical_shadow_ratio.*efficiency_factor; % Transmittance of single-fibre composite

            W (1,wavelength-(wavelength_start-1)) = wavelength;
            W (2,wavelength-(wavelength_start-1)) = refrindx(wavelength_local,entry_polymer_matrix);
            W (3,wavelength-(wavelength_start-1)) = n_f;
          % W (3,wavelength-(wavelength_start-1)) = refrindx(wavelength_local,entry_glass_fibre);
            W (4,wavelength-(wavelength_start-1)) = transmittance;
            W (5,wavelength-(wavelength_start-1)) = efficiency_factor;
            W (6,wavelength-(wavelength_start-1)) = phase_lag;
        end

    % demonstrate spectrum of 380nm-700nm in the unit of 10nm
    case 'spectrum/10'
        for wavelength = wavelength_start/10:wavelength_end/10
            rays_bundle(1,wavelength-37)=...
                Rays( nrays, 'collimated', [ 0 -10 -10+0.8*(wavelength-37) ], rays_dir, 30, 'linear', 'air', ...
                     wavelength*1e-8, 1/255.*wavelength22color(wavelength*10, 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb') );


        end

    % manual ray configuration
    case 'manual config'
        rays_bundle = [
              Rays( nrays, 'collimated', [ 0 0  0   ], rays_dir, rays_width, 'linear', 'air', 700e-9, 1/255.*[ 255 0 0 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -d   ], rays_dir, rays_width, 'linear', 'air', 607e-9, 1/255.*[ 255 165 0 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -2*d ], rays_dir, rays_width, 'linear', 'air', 580e-9, 1/255.*[ 255 255 0 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -3*d ], rays_dir, rays_width, 'linear', 'air', 510e-9, 1/255.*[ 0 255 0 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -4*d ], rays_dir, rays_width, 'linear', 'air', 490e-9, 1/255.*[ 0 255 255 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -5*d ], rays_dir, rays_width, 'linear', 'air', 440e-9, 1/255.*[ 0 0 255 ] ) ...
              Rays( nrays, 'collimated', [ 0 0 -6*d ], rays_dir, rays_width, 'linear', 'air', 380e-9, 1/255.*[ 97 0 97 ] )
            ];
end

sz_rays_bundle = size(rays_bundle);

for i = 1:sz_rays_bundle(1,2)

    cache = [bench.trace(rays_bundle(i), 0)];
    cachesize = size(cache);
    stack(1, 1+cachesize(1,2)*(i-1):cachesize(1,2)*(i-1)+cachesize(1,2)) = bench.trace(rays_bundle(i), 0);
    % automatically expand Rays matrix for newly added rays for rays_through
    cache2 = [ref_bench.trace(rays_bundle(i), 0)];
    cache2size = size(cache2);
    stack2(1, 1+cache2size(1,2)*(i-1):cache2size(1,2)*(i-1)+cache2size(1,2)) = ref_bench.trace(rays_bundle(i), 0);
    
    wavelength_rays_bundle(i,:) = stack(end).w;
    offset_rays_bundle(i,:) = stack(end).r(:,2);
    offset_ref(i,:) = stack2(end).r(:,2);
end

rays_through = stack;

% trace the rays, enable tracing rays that miss some bench elements by setting the second input parameter to 0

%ref
rays2_through = stack2;
%rays_all_through = [rays_through rays2_through];


% draw bench elements and draw rays as arrows
%bench.draw( rays_through, 'lines' , 0.1);  % display everything, the other draw option is 'lines'
%bench.draw( rays2_through, 'lines' , 0.1);







%
deviated_angle_matrix = [abs(rad2deg(atan((offset_rays_bundle-offset_ref)/(origin_glass_fibre+4*diameter_glass_fibre-(origin_glass_fibre-2*diameter_glass_fibre)))))]; %output of deviation in degree in abs

deviated_rays = zeros(1,sz_rays_bundle(1,2));

for counter = 1:sz_rays_bundle(1,2)
for z = 1:nrays

    if deviated_angle_matrix(counter,z) > 2.5
        deviated_rays(1,counter) = deviated_rays(1,counter)+1;
    end
    
end
end

wavelength_rays_bundle(:,1) = wavelength_rays_bundle(:,1).*1e9;
Haze = deviated_rays./nrays;
wa(:,1) = wavelength_rays_bundle(:,1);

%
for counter = 1:sz_rays_bundle(1,2)
color(counter,:) = 1/255.*wavelength22color( wa(counter,1), 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb');
end
colormap(flipud(color))

rainbowplot(wa',Haze);
xlabel('wavelength(μm)');
ylabel('Haze(%)');


figure();
subplot(2,2,1)
x = W(1,:);
y1 = W(2,:);
y2 = W (3,:);
plot(x,y1,'-',x,y2,'-');
legend('RI of polymer matrix','RI of glass fibre')
xlabel('wavelength(μm)')
ylabel('Refractive Index n')

subplot(2,2,2)
x2 = W (6,:);
y3 = W (5,:);
plot(x2,y3);
xlabel('Phase lag \rho')
ylabel('Efficiency factor Q_{ext}(\rho)')


subplot(2,2,3)
y3 = W (4,:);
plot(x,y3);
xlabel('wavelength(μm)')
ylabel('Transmittance')

subplot(2,2,4)
y4 = W (6,:);
plot(x,y4);
xlabel('wavelength(μm)')


%phi_z = asin(RAY(1,:).n(:,3)./RAY(1,:).n(:,1))
%
% 
% plot(phi_y,off-off_ref,'.');
% xlabel('Beam Angle Y (rad)');
% ylabel('Offset (mm)');
% %zlabel('Offset (mm)');
% title('Beam Offset on screen (referenced)');






end
%% 

function colorCode = wavelength22color(wavelength, varargin)

  % default arguments
  maxIntensity = 1;
  gammaVal = 0.8;
  colorSpace = 'rgb';

  for iargin=1:2:(nargin-1)
    switch varargin{iargin}
      case 'maxIntensity' 
        maxIntensity = varargin{iargin + 1};
      case 'gammaVal'
        gammaVal = varargin{iargin + 1};
      case 'colorSpace'
        switch varargin{iargin + 1}
          case 'rgb'
            colorSpace = 'rgb';
          case 'hsv'
            colorSpace = 'hsv';
            otherwise
            error('Invalid colorspace defined');
        end
        otherwise
        error('Invalid argument passed');
    end
  end

	function outputVal = adjust(inputVal, factor)

		if (inputVal == 0)
	  	outputVal = 0;
  	  else
			outputVal = (inputVal * factor)^gammaVal;
  	  end

    end

	if (wavelength >= 380) && (wavelength < 440)
		r = -(wavelength - 440) / (440 - 380);
    g = 0;
    b = 1;
	elseif (wavelength >= 440) && (wavelength < 490)
		r = 0;
 		g = (wavelength - 440) / (490 - 440);
    b = 1;
  elseif (wavelength >= 490) && (wavelength < 510)
  	r = 0;
    g = 1;
    b = -(wavelength - 510) / (510 - 490);
  elseif (wavelength >= 510) && (wavelength < 580)
  	r = (wavelength - 510) / (580 - 510);
    g = 1;
    b = 0;
  elseif (wavelength >= 580) && (wavelength < 645)
    r = 1;
    g = -(wavelength - 645) / (645 - 580);
    b = 0;
  elseif (wavelength >= 645) && (wavelength < 780)
    r = 1;
    g = 0;
    b = 0;
    else
  	r = 0;
    g = 0;
    b = 0;
    end
    
  if (wavelength >= 380) && (wavelength < 420)
  	factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
  elseif (wavelength >=  420) && (wavelength < 700)
    factor = 1;
  elseif (wavelength >= 700) && (wavelength < 780)
  	factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
  else
    factor = 0;
  end

  r = adjust(r, factor);
  g = adjust(g, factor);
  b = adjust(b, factor);

  rgbCode = [r, g, b];

  switch colorSpace
    case 'rgb'
      colorCode = rgbCode;
    case 'hsv'
      colorCode = rgb2hsv(rgbCode);
  end

  colorCode = colorCode * maxIntensity;

end
% A simple tool to convert a wavelength in nm to an RGB, hexadecimal or HSL colour.
% https://academo.org/demos/wavelength-to-colour-relationship/

function polymer_fibre_setup(O, D, entry_polymer_matrix, bench, h, entry_glass_fibre)

%%
% add optical elements in the order they are encountered by light rays


upperbound = Plane( [ O-2*D 0 0 ], 200, 200, { 'air', entry_polymer_matrix } );
bench.append( upperbound );

% NOTE THAT THE CYLIDRICAL LENS IS DOUBLED HERE, BECAUSE RAYS ENCOUNTER IT TWICE!!!

% cylindrical lens front surface
lens1a = CylinderLens( [ O-D 0 2*h ], D, 4*h, { entry_polymer_matrix ,entry_glass_fibre });
lens1a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens1a );
% cylindrical lens back surface
lens1b = CylinderLens( [ O-D 0 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens1b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens1b );



% cylindrical lens front surface
lens2a = CylinderLens( [ O-D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens2a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens2a );
% cylindrical lens back surface
lens2b = CylinderLens( [ O-D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens2b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens2b );
% cylindrical lens front surface
lens3a = CylinderLens( [ O-D/2 D/2*sqrt(3) 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens3a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens3a );
% cylindrical lens back surface
lens3b = CylinderLens( [ O-D/2 D/2*sqrt(3) 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens3b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens3b );



% cylindrical lens front surface
lens4a = CylinderLens( [ O 0 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens4a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens4a );
% cylindrical lens back surface
lens4b = CylinderLens( [ O 0 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens4b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens4b );


% cylindrical lens front surface
lens5a = CylinderLens( [ O+D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens5a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens5a );
% cylindrical lens back surface
lens5b = CylinderLens( [ O+D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens5b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens5b );
% cylindrical lens front surface
lens6a = CylinderLens( [ O+D/2 D/2*sqrt(3) 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens6a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens6a );
% cylindrical lens back surface
lens6b = CylinderLens( [ O+D/2 D/2*sqrt(3) 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens6b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens6b );


% cylindrical lens front surface
lens7a = CylinderLens( [ O+D 0 2*h ], D, 4*h, { entry_polymer_matrix,entry_glass_fibre } );
lens7a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens7a );
% cylindrical lens back surface
lens7b = CylinderLens( [ O+D 0 2*h ], D, 4*h, { entry_glass_fibre,entry_polymer_matrix} );
lens7b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens7b );

% polymer matrix lower bound
lowerbound = Plane( [ O+2*D 0 0 ], 200, 200, { entry_polymer_matrix, 'air' } );
bench.append( lowerbound );
end
% model construct

function rainbowplot(x, y)
%------------------------------------------------------------
%  RAINBOWPLOT Colorful linear 2-D plot
%  This function plots a line colored like a rainbow. 
%  The line is defined by x versus y pairs, which is the same
%  as a regular 2-D plot.
%
%   Usage Examples,
%
%   x = 1:100; y = randn(1,100);  
%   rainbowplot(x, y);
%
%   Kun Liu 
%   Version 1.00
%   June, 2006
%------------------------------------------------------------

if size(x, 1)~=1 || size(y, 1)~=1 
    error('x and y must be one dimensional vector...');    
end
if size(x, 2) ~= size(y, 2)
    error('x and y must have the same number of elements...');
end

length = size(x, 2);
d = length:-1:1;
wavelength = patch([x nan],[y nan], [d nan], 'EdgeColor', 'interp');
end 
%rainbowplot