function TransComp()
%
%  analysis of light transmittance of glass fibre of cylindrical section
%
%  Zheng Ji, 11.07.2022



startup;

O = 50; % origin of cylinder
D = 10; % cylinder diameter
h = 20; % cylinder height
RI_matrix = 'pmma'; % refractive index of polymer matrix
RI_fibre = 'bk7';   % refractive index of glass fibre

% create a container for optical elements (Bench class)
bench = Bench;
ref_bench = Bench;

polymer_fibre_setup(O, D, RI_matrix, bench, h, RI_fibre);

% screen
screen = ScreenGeneric( [ O+4*D 0 0 ], 10 * D, 10 * D, 256, 256, 'wf' );
bench.append( screen );
ref_bench.append( screen );

% create divergent rays
nrays = 31;
d = -5; % distance between rays of different wavelengths

rays_dir = [1 0 0];


democase = 'spectrum/10';
switch(democase)

    % demonstrate spectrum of 380nm-700nm in the unit of 1nm
    case 'spectrum'
        for p = 380:700
            
            stack3(1,p-379)=...
                Rays( nrays, 'collimated', [ 0 0 -10+0.5*(p-379) ], rays_dir, 30., 'linear', 'air', ...
                      p*1e-9, 1/255.*wavelength22color(p, 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb') );
        end
        RAY = stack3;

    % demonstrate spectrum of 380nm-700nm in the unit of 10nm
    case 'spectrum/10'
        for p = 38:70
            stack3(1,p-37)=...
                Rays( nrays, 'collimated', [ 0 0 -10+0.8*(p-37) ], rays_dir, 30., 'linear', 'air', ...
                     p*1e-8, 1/255.*wavelength22color(p*10, 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb') );


        end
        RAY = stack3;

    % manual ray configuration
    case 'manual config'
        RAY = [
            Rays( nrays, 'collimated', [ 0 0  0   ], rays_dir, 30., 'linear', 'air', 700e-9, 1/255.*[ 255 0 0 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -d   ], rays_dir, 30., 'linear', 'air', 607e-9, 1/255.*[ 255 165 0 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -2*d ], rays_dir, 30., 'linear', 'air', 580e-9, 1/255.*[ 255 255 0 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -3*d ], rays_dir, 30., 'linear', 'air', 510e-9, 1/255.*[ 0 255 0 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -4*d ], rays_dir, 30., 'linear', 'air', 490e-9, 1/255.*[ 0 255 255 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -5*d ], rays_dir, 30., 'linear', 'air', 440e-9, 1/255.*[ 0 0 255 ] ) ...
             Rays( nrays, 'collimated', [ 0 0 -6*d ], rays_dir, 30., 'linear', 'air', 380e-9, 1/255.*[ 97 0 97 ] )
            ];
end


N = size(RAY);


for i = 1:N(1,2)

    cache = [bench.trace(RAY(i), 0)];
    cachesize = size(cache);
    stack(1, 1+cachesize(1,2)*(i-1):cachesize(1,2)*(i-1)+cachesize(1,2)) = bench.trace(RAY(i), 0);
    % automatically expand Rays matrix for newly added rays for rays_through
    cache2 = [ref_bench.trace(RAY(i), 0)];
    cache2size = size(cache2);
    stack2(1, 1+cache2size(1,2)*(i-1):cache2size(1,2)*(i-1)+cache2size(1,2)) = ref_bench.trace(RAY(i), 0);
end

rays_through = stack;

% trace the rays, enable tracing rays that miss some bench elements by setting the second input parameter to 0

%ref
rays2_through = stack2;
%rays_all_through = [rays_through rays2_through];


% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'lines' , 0.1);  % display everything, the other draw option is 'lines'
%bench.draw( rays2_through, 'lines' , 0.1);
 
% figure();
% %opl = rays_through(end).opl;
% %opl_ref = rays2_through(end).opl;
% 
% off = vecnorm(rays_through(end).r(:,2:3)')
% off_ref = vecnorm(rays2_through(end).r(:,2:3)')
% 
% off-off_ref;
% off./off_ref;
% 
% phi_y = asin(RAY(1,:).n(:,2)./RAY(1,:).n(:,1));
% phi_z = asin(RAY(1,:).n(:,3)./RAY(1,:).n(:,1));
% 
% 
% plot(phi_y,off-off_ref,'.');
% xlabel('Beam Angle Y (rad)');
% ylabel('Offset (mm)');
% %zlabel('Offset (mm)');
% title('Beam Offset on screen (referenced)');








end

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

function polymer_fibre_setup(O, D, RI_matrix, bench, h, RI_fibre)

%%
% add optical elements in the order they are encountered by light rays


upperbound = Plane( [ O-2*D 0 0 ], 200, 200, { 'air', RI_matrix } );
bench.append( upperbound );

% NOTE THAT THE CYLIDRICAL LENS IS DOUBLED HERE, BECAUSE RAYS ENCOUNTER IT TWICE!!!

% cylindrical lens front surface
lens1a = CylinderLens( [ O-D 0 2*h ], D, 4*h, { RI_matrix ,RI_fibre });
lens1a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens1a );
% cylindrical lens back surface
lens1b = CylinderLens( [ O-D 0 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens1b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens1b );



% cylindrical lens front surface
lens2a = CylinderLens( [ O-D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens2a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens2a );
% cylindrical lens back surface
lens2b = CylinderLens( [ O-D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens2b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens2b );
% cylindrical lens front surface
lens3a = CylinderLens( [ O-D/2 D/2*sqrt(3) 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens3a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens3a );
% cylindrical lens back surface
lens3b = CylinderLens( [ O-D/2 D/2*sqrt(3) 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens3b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens3b );



% cylindrical lens front surface
lens4a = CylinderLens( [ O 0 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens4a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens4a );
% cylindrical lens back surface
lens4b = CylinderLens( [ O 0 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens4b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens4b );


% cylindrical lens front surface
lens5a = CylinderLens( [ O+D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens5a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens5a );
% cylindrical lens back surface
lens5b = CylinderLens( [ O+D/2 -D/2*sqrt(3) 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens5b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens5b );
% cylindrical lens front surface
lens6a = CylinderLens( [ O+D/2 D/2*sqrt(3) 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens6a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens6a );
% cylindrical lens back surface
lens6b = CylinderLens( [ O+D/2 D/2*sqrt(3) 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens6b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens6b );


% cylindrical lens front surface
lens7a = CylinderLens( [ O+D 0 2*h ], D, 4*h, { RI_matrix,RI_fibre } );
lens7a.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens7a );
% cylindrical lens back surface
lens7b = CylinderLens( [ O+D 0 2*h ], D, 4*h, { RI_fibre,RI_matrix} );
lens7b.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens7b );

% polymer matrix lower bound
lowerbound = Plane( [ O+2*D 0 0 ], 200, 200, { RI_matrix, 'air' } );
bench.append( lowerbound );
end
% model construct