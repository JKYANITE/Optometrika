clear all
close all
clc

% Initialisation
N = 9; % Number of fibres
R_f = 7500e-9; % Radius of cylindrical fibre
w_s = 10000e-9
w = 60000e-9;
h = 60000e-9;
f = N*pi*(R_f)^2/(w*h)  % Volume fraction
 
% dS = 1; % Unit surface area
% I = 1; % Intensity of transmitted light
% I_0 = 1; % Uniform intensity of incident light flux


% T = I/I_0; % Transmittance of the incident light in monochromatic light


theta = 0; % Scattering angle - angle between the incident light and scattered light



wavelength_start = 400;
wavelength_end = 800;

democase = ['single'];

switch(democase)
    case 'single'
for lambda = wavelength_start:wavelength_end
% rho = 2*k*R_f*abs(n_f-n_m); % Phase lag sustained by the centre ray
%lambda = 589e-9;

n_f = 1.49+2000.*lambda^(-2)-808.*lambda^(-4)+0.1.*lambda^(-6);
n_m = 1.4794+4557.*lambda^(-2)+903.*lambda^(-4)+0.1.*lambda^(-6);
% %  a_1 = 0.751;
% %  a_2 = 1.16*10^4;
% %  n_m = sqrt(1/(a_1-a_2/lambda^2)+1);


lambda_local = lambda.*1e-9;

% n_a = 1;

rho = 4.*pi.*R_f.*abs(n_f-n_m)/(lambda_local);
% rho = 4.*pi.*R_f.*0.1/(lambda_local);

% c_ext = 4/k*Re*A(0); % Effective extinction width - energy incident on a strip with width
% Q_ext = c_ext/(2*R_f); % Efficiency factor
fun = @(gamma) 2.*rho.*sin(rho.*cos(gamma)).*(sin(gamma)).^2;
Q_ext = integral(fun,0,pi/2);
% Q_ext = (2*rho^2)/3; % Approximation

% gamma = 0; % Angle between incident light direction and surface of the fibre
% delta = 2*R_f*cos(gamma); % Light path travelled in the fibre


g = 2.*R_f/w_s; % Geomatrical shadow ratio
T_s = 1-g.*Q_ext; % Transmittance of single-fibre composite

%T_m=(1-((n_m-n_a)/(n_m+n_a))^2)^2*exp(0);

W (1,lambda-(wavelength_start-1)) = lambda;
W (2,lambda-(wavelength_start-1)) = n_m;
W (3,lambda-(wavelength_start-1)) = n_f;
W (4,lambda-(wavelength_start-1)) = T_s;
W (5,lambda-(wavelength_start-1)) = Q_ext;
W (6,lambda-(wavelength_start-1)) = rho;

end

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

    case 'multi'
        d = 5*1e-3;
        G_f = 2*(f/pi)^(1/2); % geometric shadow fraction
        p = 2*f*d/(pi*R_f); % expected number of times through multiplication of the probability

        for lambda = wavelength_start:wavelength_end
            
            n_f = 1.4792+4500.*lambda^(-2)-808.*lambda^(-4)+0.1.*lambda^(-6);
            n_m = 1.4794+4557.*lambda^(-2)+903.*lambda^(-4)+0.1.*lambda^(-6);
            
            lambda_local = lambda.*1e-9;
            
            rho = 4.*pi.*R_f.*abs(n_f-n_m)/(lambda_local);
            fun = @(gamma) 2.*rho.*sin(rho.*cos(gamma)).*(sin(gamma)).^2;
            Q_ext = integral(fun,0,pi/2);
            
            T = (1-2.*Q_ext.*(f./pi).^0.5).^((d./R_f)*(f./pi).^0.5);

            W (1,lambda-(wavelength_start-1)) = lambda;
            W (2,lambda-(wavelength_start-1)) = n_m;
            W (3,lambda-(wavelength_start-1)) = n_f;
            W (4,lambda-(wavelength_start-1)) = T;
            W (5,lambda-(wavelength_start-1)) = Q_ext;
            W (6,lambda-(wavelength_start-1)) = rho;
        end
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


    case'Beer-Lambert'
        d = 5*1e-3;
        G_f = 2*(f/pi)^(1/2); % geometric shadow fraction
        p = 2*f*d/(pi*R_f); % expected number of times through multiplication of the probability

        for lambda = wavelength_start:wavelength_end
            
            n_f = 1.486+2500.*lambda^(-2)-808.*lambda^(-4)+0.1.*lambda^(-6);
            n_m = 1.4794+4557.*lambda^(-2)+903.*lambda^(-4)+0.1.*lambda^(-6);

%            n_f = 1.4794+3920.*lambda^(-2)-808.*lambda^(-4)+0.1.*lambda^(-6);
%            n_m = 1.4794+4557.*lambda^(-2)+903.*lambda^(-4)+0.1.*lambda^(-6);
            n_a = 1;
            

         
            
            fun_fresnel = @(A) 1./(2*(sin(A)).^2.*...
                        (cos(n_f.*sin(A)./n_m)).^2./(sin(A+n_f.*sin(A)./n_m)).^2.*...
                        (1+1./(cos(A-n_f.*sin(A)./n_m)).^2));
            
            T_fresnel = integral(fun_fresnel,-pi/2,pi/2)
            
            
            T_m1 = (1-((n_m-n_a)/(n_m+n_a))^2)^2*exp(0);

            T_m2 = (1-((n_f-n_m)/(n_f+n_m))^2)^2*exp(0)*T_m1;

            T_m3 = (1-((n_m-n_f)/(n_m+n_f))^2)^2*exp(0)*T_m2;

            T_m4 = (1-((n_a-n_m)/(n_a+n_m))^2)^2*exp(0)*T_m3;

            lambda_local = lambda.*1e-9;
            
            rho = 4.*pi.*R_f.*abs(n_f-n_m)/(lambda_local);
            fun = @(gamma) 2.*rho.*sin(rho.*cos(gamma)).*(sin(gamma)).^2;
            Q_ext = integral(fun,0,pi/2);
            
            T = (1-2.*Q_ext.*(f/pi).^0.5).^((d/R_f)*(f/pi)^0.5)*T_m4;

            W (1,lambda-(wavelength_start-1)) = lambda;
            W (2,lambda-(wavelength_start-1)) = n_m;
            W (3,lambda-(wavelength_start-1)) = n_f;
            W (4,lambda-(wavelength_start-1)) = T_m4;
            W (5,lambda-(wavelength_start-1)) = Q_ext;
            W (6,lambda-(wavelength_start-1)) = T_fresnel;
        end
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

       
end

% T_fresnel = 2*(sin(A))^2*...
%                 (cos(n_f*sin(A)/n_m))^2/sin(A+n_f*sin(A)/n_m)^2*...
%                 (1+1/(cos(A-n_f*sin(A)/n_m))^2)
