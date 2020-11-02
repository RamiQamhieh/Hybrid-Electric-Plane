function [t,p,rho] = T_P_rho(alt,units)
%function [t,p,rho]=statm2(alt,units)
%Gives the temperature, density and pressure
%for a given altitude. Works in SI or US units.
%Accepts a vector for altitude.

%Your function should look like
%[t, p, rho]=statm(alt,units)

%Determine the constants
if (units =='SI')
    R=287;
    g=9.81;
    
    h1 = 11000 ; h2=25000;
    h3 = 47000 ; h4 = 53000;
    h5 = 79000 ; h6 = 90000;
    
    t0=288.16;t1=216.66;
    t2 = 282.66 ; t3 = 165.66;
    
    a1=-6.5E-3;a2=3E-3;
    a3 = -4.5E-3;
    
   p0=1.013e5;
   p1=2.2346e4;
   p2=p1*exp(-g/R/t1*(h2-h1));
   p3=p2*(t2/t1)^(-g/a2/R);
   p4=p3*exp(-g/R/t2*(h4-h3));
   p5=p4*(t3/t2)^(-g/a3/R);
   p6=p5*exp(-g/R/t3*(h6-h5));
   

   rho0=1.226;
   rho1=3.5932e-1;
   rho2=rho1*exp(-g/R/t1*(h2-h1));
   rho3=rho2*(t2/t1)^(-g/a2/R-1);
   rho4=rho3*exp(-g/R/t2*(h4-h3));
   rho5=rho4*(t3/t2)^(-g/a3/R-1);
   rho6=rho5*exp(-g/R/t3*(h6-h5));
    
    %include additional values for h,t,a here
    %you will also need p and rho at cutoff points
    
elseif (units=='US')
    R=1716;
    
    h1=36100;h2=82000;h3=154200;
    h4=173900;h5=259000;h6=295000;
    
    a1=-3.6E-3;a2=1.65E-3; a3=-2.5E-3;
    
    t0=518.67;t1=388.70999;
    t2=508.79;t3=298.19;
    
    g=32.2;
    
    p0=2116.2;
    p1=470.5716;
    p2=p1*exp(-g/R/t1*(h2-h1));
    p3=p2*(t2/t1)^(-g/a2/R);
    p4=p3*exp(-g/R/t2*(h4-h3));
    p5=p4*(t3/t2)^(-g/a3/R);
    
    rho0=2.3769e-3;
    rho1=7.05253667e-4;
    rho2=rho1*exp(-g/R/t1*(h2-h1));
    rho3=rho2*(t2/t1)^(-g/a2/R-1);
    rho4=rho3*exp(-g/R/t2*(h4-h3));
    rho5=rho4*(t3/t2)^(-g/a3/R-1);

    %include additional values for h,t,a here
    %you will also need p and rho at cutoff points

end

i=1;%initialize index for t
for h=alt
    if h<=h1
        t(i)=t0+a1*h;
        p(i)=p0*(t(i)/t0)^(-g/a1/R);
        rho(i)=rho0*(t(i)/t0)^(-g/a1/R-1);
        
    elseif h<h2
        t(i)=t1;
        p(i)=p1*exp(-g/R/t1*(h-h1));
        rho(i)=rho1*exp(-g/R/t1*(h-h1));
        
    elseif h<=h3
        t(i)=t1+a2*(h-h2);
        p(i)=p2*(t(i)/t1)^(-g/a2/R);
        rho(i)=rho2*(t(i)/t1)^(-g/a2/R-1);
        
    elseif h<h4
        t(i)=t2;
        p(i)=p3*exp(-g/R/t2*(h-h3));
        rho(i)=rho3*exp(-g/R/t2*(h-h3));
        
    elseif h<=h5
        t(i)=t2+a3*(h-h4);
        p(i)=p4*(t(i)/t2)^(-g/a3/R);
        rho(i)=rho4*(t(i)/t2)^(-g/a3/R-1);
        
    elseif h<=h6
        t(i)=t3;
        p(i)=p5*exp(-g/R/t3*(h-h5));
        rho(i)=rho5*exp(-g/R/t3*(h-h5));
                
    %continue with elseif's for higher altitudes 
    end
    i=i+1;  %increment the index for t
end
t=t';
p=p';
rho=rho';
