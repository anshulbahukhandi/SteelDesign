clc;
wsectionz;
jj=get_section_property(W,'W14x132');
jh = W{jj,2};
d = jh(1); b_f = jh(2); t_w = jh(3); t_f =jh(4); A=jh(5); II_xx= jh(6);II_yy=jh(7);J=jh(8); weight=jh(9);Zx=jh(10);Zy=jh(11);
Fr=10;
Fy=36;
P=800;
Mx=200;
My=0;
M2=Mx;
M1=-Mx;
Lb=15;
L=12*Lb;
lamdaf=b_f/(2*t_f);
lamdaw=d/t_w;
lamdapf=65/sqrt(Fy);
lamdapw=640/sqrt(Fy);
lamdarf=141/sqrt(Fy-Fr);
lamdarw=970/sqrt(Fy);
Ryy = (II_yy/A)^(1/2);
Rxx = (II_xx/A)^(1/2);
Sxx = II_xx/(d/2);
de=d-t_w;
Cw = (de^2)*((b_f)^3)*t_f/24;
Lp=300*Ryy/(sqrt(Fy)*12);
X1=40014.98*sqrt(J*A)/Sxx;
X2=(4*Cw/II_yy)*((Sxx/(11200*J))^2);
Lr=((Ryy*X1)*sqrt(1+sqrt(1+X2*((Fy-Fr)^2))))/((Fy-Fr)*12);
% Cb=12.5*Mmax/(2.5*Mmax+3*Ma+4*Mb+3*Mc);
pi=3.1415;
phi=0.85;
Cb=1.75+1.05*(M1/M2)+0.3*((M1/M2)^2);
Cm=0.6-0.4*(M1/M2);

Mp=Fy*Zx/12;
%young's modulus of elasticity for steel in KSI 
E=29000;
%End condition parameter
K=1;

bloadx = (pi*pi * E * II_xx)/(K*L*K*L);
B1=Cm/(1-(P/bloadx));
if(B1>1)
    Mx=B1*Mx;
end
if(K*L/Rxx>=K*L/Ryy)
   lambda=K*L*sqrt(Fy/E)/(Rxx *pi);
   if (lambda<=1.5) 
       Fcr=Fy*(0.658)^(lambda^2);
   else
       Fcr=Fy*(0.877/(lambda^2));
   end
else
   lambda=K*L*sqrt(Fy/E)/(Ryy*pi);
   if (lambda<=1.5) 
       Fcr=Fy*(0.658)^(lambda^2);
   else
       Fcr=Fy*(0.877/(lambda^2));
   end
end
Pn=A*Fcr;
fprintf('--------------------------------------------------\n')
fprintf('Design strength of the Column is : %0.3f Kips\n',Pn);
fprintf('--------------------------------------------------\n')
if(lamdaf<lamdapf)
    if(Lb<=Lp)
        Mn=Mp;
    elseif((Lp<Lb)&&(Lb<Lr))
        Mr=(Fy-Fr)*Sxx/12;
        Mn=Cb*(Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp)));
        if(Mn>Mp)
            Mn=Mp;
        end
    else
        Me=(Sxx*X1*sqrt(2)*Ryy/(12*Lb))*(sqrt(1+(X1*X1*X2*Ryy*Ryy*0.5/(12*12*Lb*Lb))));
        Mn=Me/12;
    end
elseif((lamdapf<lamdaf)&&(lamdaf<lamdarf))||((lamdapw<lamdaw)&&(lamdaw<lamdarw))
    if(lamdaf<=lamdapf)
        Mn=Mp;
    elseif((lamdapf<lamdaf)&&(lamdaf<lamdarf))
        Mr=(Fy-Fr)*Sxx/12;
        Mn=Cb*(Mp-(Mp-Mr)*((lamdaf-lamdapf)/(lamdarf-lamdapf)));
        if(Mn>Mp)
            Mn=Mp;
        end
    end
    if(lamdaw<=lamdapw)
        Mn=Mp;
    elseif((lamdapw<lamdaw)&&(lamdaw<lamdarw))
        Mr=(Fy(i)-Fr)*Sxx/12;
        Mn=Cb*(Mp-(Mp-Mr)*((lamdaw-lamdapw)/(lamdarw-lamdapw)));
        if(Mn>Mp)
            Mn=Mp;
        end
    end
    if(Lb<Lp)
        Mn=Mp;
    elseif((Lp<Lb)&&(Lb<Lr))
       Mr=(Fy-Fr)*Sxx/12;
       Mn=Cb*(Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp)));
       if(Mn>Mp)
          Mn=Mp;
       end 
    else
        Me=(Sxx*X1*sqrt(2)*Ryy/(12*Lb))*(sqrt(1+(X1*X1*X2*Ryy*Ryy*0.5/(12*12*Lb*Lb))));
        Mn=Me/12;
    end
end
fprintf('--------------------------------------------------\n')
fprintf('Design strength of the Beam is : %0.3f Kips-ft\n',Mn);
fprintf('--------------------------------------------------\n')
Pe=P/(phi*Pn);
if(Pe>=0.2)
   if(Pe+(8/9)*((Mx/(0.9*Mn))+(My/(0.9*Mn)))<=1)
      fprintf('%s Section is acceptable to LRFD\n',W{jj,1});
   else
      fprintf('%s Section is not acceptable to LRFD\n',W{jj,1});
   end
end
if(Pe<0.2)
    if(Pe/2+(Mx/(0.9*Mn)+(My/(0.9*Mn)))<=1)
      fprintf('%s Section is acceptable to LRFD\n',W{jj,1});
   else
      fprintf('%s Section is not acceptable to LRFD\n',W{jj,1});
   end
end
        