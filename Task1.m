%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clc;
wsectionz;
%~~~~~~~~~~~~~~~~~~~~~~~Data Input~~~~~~~~~~~~~~~~~~~%
%grade of steel in KSI
fy=50;   
%length of column in inchs
%l=50*12;
for i=2:1:10
    l=5*i*12;
    5*i
%young's modulus of elasticity for steel in KSI 
E=29000;
%End condition parameter
K=2;

phi=0.85;
jj=get_section_property(W,'W24x94');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
jh = W{jj,2};
d = jh(1);
b_f = jh(2);
t_w = jh(3);
t_f =jh(4);
A=jh(5);
a_total=2*t_f*b_f + t_w*(d-2*t_f);
I_xx= jh(6);
r_xx=(I_xx/A)^0.5;
I_x= (1/12*(d-2*t_f)^3*t_w  +  2*(b_f*t_f*((d-t_f)/2)^2+1/12*b_f*t_f^3));
r_x=(I_x/a_total)^0.5;
I_yy=jh(7);J=jh(8); weight=jh(9);
I_y=2*1/12*t_f*b_f^3 + 1/12*t_w^3*(d-2*t_f);
r_yy=(I_yy/A)^0.5;
r_y=(I_y/(2*b_f*t_f + t_w*(d-2*t_f)))^0.5;
% %..............Properties about X axis....................................
% 
 fprintf('Calculated Moment of Inertia about X axis is : %0.4f\n',I_x);
 fprintf('Moment of Inertia about X axis as given in table(LRFD) : %0.4f\n',I_xx);
%......................Properties about Y axis.............................

 fprintf('Calculated Moment of Inertia about Y axis is : %0.4f\n',I_y);
 fprintf('Moment of Inertia about Y axis as given in (LRFD Manual): %0.4f\n',I_yy);
bloadx = (pi*pi * E * I_xx)/(K*l*K*l);
bloady = (pi*pi * E * I_yy)/(K*l*K*l);

fprintf('Buckling load of the column in x direction is : %0.3f Kips\n',bloadx);
fprintf('Buckling load of the column in y direction is : %0.3f Kips\n',bloady);
bload=min(bloadx,bloady);

if(K*l/r_xx>=K*l/r_yy)
   lambda=K*l*sqrt(fy/E)/(r_xx *pi);
   if (lambda<=1.5) 
       fcr=fy*(0.658)^(lambda^2);
   else
    fcr=fy*(0.877/(lambda^2));
   end
else
   lambda=K*l*sqrt(fy/E)/(r_yy*pi);
   if (lambda<=1.5) 
       fcr=fy*(0.658)^(lambda^2);
   else
       fcr=fy*(0.877/(lambda^2));
   end
end
fprintf('--------------------------------------------------\n')
fprintf('Design strength of the column is : %0.3f Kips\n',phi*A*fcr);
fprintf('--------------------------------------------------\n')
end



