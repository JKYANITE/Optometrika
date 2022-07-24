clear all
close all
clc


A = 0:0.01:pi/2
B = 1.6*sin(A)/1.61
fun_fresnel = @(B) 2*(sin(A)).^2.*...
                        (((cos(B)))./sin(A+(asin(B)))).^2.*...
                        (1+1./(cos(A-asin(B))).^2);
fun_fresnel(2.26841)
rad2deg(2.26841)

colormap summer

rainbowplot(A/pi*2*90,fun_fresnel(B)), grid on

   %         T_fresnel = integral(fu1.5resnel,0,pi/2)


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
p = patch([x nan],[y nan], [d nan], 'EdgeColor', 'interp');
end