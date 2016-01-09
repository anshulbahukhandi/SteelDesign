clc;
wsectionz;
jj=get_section_property(W,'W24x62');
jh = W{jj,2};
d(1) = jh(1); b_f(1) = jh(2); t_w(1) = jh(3); t_f(1) =jh(4); A(1)=jh(5); II_xx(1)= jh(6);II_yy(1)=jh(7);J(1)=jh(8); weight(1)=jh(9);Zx(1)=jh(10)
Zy(1)=jh(11);
jj1=get_section_property(W,'W24x94');
jh = W{jj1,2};
d(2) = jh(1); b_f(2) = jh(2); t_w(2) = jh(3); t_f(2) =jh(4); A(2)=jh(5); II_xx(2)= jh(6);II_yy(2)=jh(7);J(2)=jh(8); weight(2)=jh(9);Zx(2)=jh(10)
Zy(2)=jh(11);
Fr=10;
n=2;
Fy=50;
M=zeros(40,1);
for i=1:n
wdl=500;
wll=1000;
w=(1.2*wdl+1.6*wll)/1000;
Mb=w*Lb*Lb/8;
Mmax=Mb;
Ma=3*w*Lb*Lb/32;
Mc=Ma;
lamdaf=b_f(i)/(2*t_f(i));
lamdaw=d(i)/t_w(i);
lamdapf=65/sqrt(Fy);
lamdapw=640/sqrt(Fy);
lamdarf=141/sqrt(Fy-Fr);
lamdarw=970/sqrt(Fy);
Ryy = (II_yy(i)/A(i))^(1/2);
Sxx = II_xx(i)/(d(i)/2);
de=d(i)-t_w(i);
Cw = (de^2)*((b_f(i))^3)*t_f(i)/24;
Lp=300*Ryy/(sqrt(Fy)*12)
X1=40014.98*sqrt(J(i)*A(i))/Sxx;
X2=(4*Cw/II_yy(i))*((Sxx/11200*J(i))^2);
Lr=((Ryy*X1)*sqrt(1+sqrt(1+X2*((Fy-Fr)^2))))/((Fy-Fr)*12)
Lc=0;Lcp=0;Lnc=0;Lncp=0;Lnci=0;Lci=0;
% Cb=12.5*Mmax/(2.5*Mmax+3*Ma+4*Mb+3*Mc);
Cb=1;
Mp=Fy*Zx(i)/12;
for Lb=1:1:40
if(lamdaf<lamdapf)
    if(Lb<=Lp)
        Mn=Mp;
    elseif((Lp<Lb)&&(Lb<Lr))
        Mr=(Fy-Fr)*Sxx/12;
        Mn=Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp));
        if(Mn>Mp)
            Mn=Cb*(Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp)));
        end
    else
        Me=(Sxx*X1*sqrt(2)*Ryy/(12*Lb))*(sqrt(1+(X1*X1*X2*Ryy*Ryy*0.5/(12*12*Lb*Lb))));
        Mn=Me/12;
    end
elseif((lamdapf<lamdaf)&&(lamdaf<lamdarf))||((lamdapw<lamdaw)&&(lamdaw<lamdarw))
    Lnc=Lb;
    if(lamdaf<=lamdapf)
        Mn=Mp;
    elseif((lamdapf<lamdaf)&&(lamdaf<lamdarf))
        Mr=(Fy-Fr)*Sxx/12;
        Mn=Mp-(Mp-Mr)*((lamdaf-lamdapf)/(lamdarf-lamdapf));
        if(Mn>Mp)
            Mn=Cb*Mn;
        end
    end
    if(lamdaw<=lamdapw)
        Mn=Mp;
    elseif((lamdapw<lamdaw)&&(lamdaw<lamdarw))
        Mr=(Fy-Fr)*Sxx/12;
        Mn=Mp-(Mp-Mr)*((lamdaw-lamdapw)/(lamdarw-lamdapw));
        if(Mn>Mp)
            Mn=Cb*Mn;
        end
    end
    if(Lb<Lp)
        Mn=Mp;
    elseif((Lp<Lb)&&(Lb<Lr))
       Mr=(Fy-Fr)*Sxx/12;
       Mn=Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp));
       if(Mn>Mp)
          Mn=Cb*Mn;
       end 
    else
        Me=(Sxx*X1*sqrt(2)*Ryy/(12*Lb))*(sqrt(1+(X1*X1*X2*Ryy*Ryy*0.5/(12*12*Lb*Lb))));
        Mn=Me/12;
    end
end
M(Lb,i)=Mn;
end
end
L=zeros(40,1);
for Lb=1:1:40
L(Lb,1)=Lb;
end
fig=figure(1); clf; grid on; axis square; hold on;
xlabel('Length (ft)'); ylabel('Design Strength (k-ft)'); title('Effect of Non Compact');
hold on;
plot(L,M);
hold off;