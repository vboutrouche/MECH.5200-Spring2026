function [ x_up, y_up, x_lo, y_lo, x_cl, y_cl ] = get_NACA(NACA_airfoil, normalized_x_vals, params)
%=======================================================================
% get_NACA_4_digits.m
%=======================================================================
%
% Inputs and Outputs
% ==================
%  i* NACA_airfoil: This input currently is a 4-digit number related to 
%       the NACA-4-digit-airfoil classification. eg: 2412 which has the
%       following properties:
%           o Percent Camber: 2/100
%           o Position of max camber: 4/10*chord
%           o Percent thickness:12
%           
%  i* normalized_x_vals: This is a vector of x-values at which the airfoil
%       coordinates are to be computed. This should always be in the 
%       interval x = [0,1]. the leading edge will be at x=0, and the t.e. 
%       at x = 1.0.
%           eg: normalized_x_vals = 1-sin(0:0.05:pi/2); 
%
%  i* params: The params input allows us to hand params to this function
%       which may apply to the overall solution process. The parameter of
%       interest in this case is the plot functionality. 
%           eg: params.plot_on = 1   ==> plots the airfoil 
%           eg: params.plot_on = 0   ==> does not plot the airfoil 
%
%
%  o* x_up
%
%  o* y_up
%
%  o* x_lo
%
%  o* x_lo
%
%  o* x_cl
%
%  o* x_cl
%
%=======================================================================
% Example useage:
% 
%    params.plot_on = 1;
%    NACA = 2412;
%    xval = 1-sin(0:pi/100:pi/2);
%    [ x_up, y_up, x_lo, y_lo, x_cl, y_cl ] = get_NACA_4(NACA, xval, params);
%=======================================================================
%
% Author(s)
% =========
% This code was written by David J. Willis, 2007.
% 
%=======================================================================
% Copyright (C) MIT, 2007, all rights reserved.
%=======================================================================
%
% References
% ==========
%
% this code was written based on information in the following sources:
%
% Article: 232 of sci.physics.computational.fluid-dynamics
% From: hulburt@leland.Stanford.EDU (Greg Payne)
% Subject: How to calculate NACA 4- and 5-digit sections<-Here it is!
% Organization: Stanford University
% Date: Sat, 7 May 1994 04:37:33 GMT
% 
% Abbot, Ira H., and von Doenhoff, Albert E.  Theory of 
% Wing Sections.  copyright 1959, Dover Publications, Inc.  
% (ISBN 0-486-60586-8)
% 
% Reigels, Dr. Friedrich W.  Aerofoil Sections.  copyright 
% 1961, Butterworth & Co. (translated from German by D. G.
% Randall )
%=======================================================================

airfoil_4_digits = NACA_airfoil;

%===
% Initialization of the airfoil
%===
xval = normalized_x_vals; 

% percent max camber
camb = floor(airfoil_4_digits/1000);
maxc = camb/100;

% Position of max camber
cpos = floor((airfoil_4_digits-camb*1000)/100);
xmxc = cpos/10;   

% thickness of airfoil
thck = ((airfoil_4_digits-camb*1000 - cpos*100))/100;

%===
% Compute the airfoil Thickness
%===
yt_c = 5.*thck.*(   0.29690.*xval.^0.5 - ...
                    0.12600.*xval.^1.0 - ...
                    0.35160.*xval.^2.0 + ...
                    0.28430.*xval.^3.0 - ...
                    0.10150.*xval.^4.0 );
%===
% Compute the camber lines
%===
for(i=1:length(xval))
    if(xval(i) == 0)
        yc_c(i) = 0;
        thet(i) = pi/2;   

    elseif(xval(i) <= xmxc)
        yc_c(i) =  (maxc).*(1.0./(xmxc.^2.0)).*...
                (2.0.*xmxc.*(xval(i)) - (xval(i)).^2.0); 
        
        thet(i) = atan((maxc).*(1.0./(xmxc.^2.0)).*...
                (2.0.*xmxc - 2.0.*(xval(i)).^1.0));    
                    
    else
        yc_c(i) =  (maxc).*(1.0./(1.0-xmxc).^2).*...
                ((1.0-2.0.*xmxc)+2.0.*xmxc.*(xval(i))-(xval(i)).^2.0);
        
        thet(i) = atan((maxc).*(1.0./(1.0-xmxc).^2).*...
                ((1.0-2.0.*xmxc)+2.0.*xmxc.*(1)-2.0.*(xval(i)).^1.0));    
    end
end

%===
% Combining the airfoil camber an thickness
%===
for(i=1:length(xval))
    x_up(i) = xval(i) - yt_c(i)*sin(thet(i));  
    y_up(i) = yc_c(i) + yt_c(i)*cos(thet(i));
    x_lo(i) = xval(i) + yt_c(i)*sin(thet(i)); 
    y_lo(i) = yc_c(i) - yt_c(i)*cos(thet(i));
end
x_cl = xval;
y_cl = yc_c;


%===
% Plotting Functions
%===
if(params.plot_on)
    
    figure
    plot(x_up, y_up,'r','LineWidth',2)
    hold on
    plot(x_lo, y_lo,'b','LineWidth',2)
    plot(x_cl, y_cl,'k--')

    plot(x_up, y_up,'ro','MarkerSize',2,'MarkerFaceColor','y')
    plot(x_lo, y_lo,'bo','MarkerSize',2,'MarkerFaceColor','c')
    axis equal
    grid on
    xlabel('x/c')
    ylabel('y/c')
    title(['A plot of the NACA ',num2str(NACA_airfoil),' airfoil section'])
end

