clc;
wsectionz;
jj=get_section_property(W,'W18x192');
jh = W{jj,2};
d = jh(1); b_f = jh(2); t_w = jh(3); t_f =jh(4); A=jh(5); II_xx= jh(6);II_yy=jh(7);J=jh(8); weight=jh(9);Zx=jh(10);Zy=jh(11);
Fr=10;
n=3;
Fy(1)=36;Fy(2)=50;Fy(3)=60;
M=zeros(40,1);
for i=1:n
P=800;
Mx=200;
My=0;
M2=Mx;
M1=-Mx;
lamdaf=b_f/(2*t_f);
lamdaw=d/t_w;
lamdapf=65/sqrt(Fy(i));
lamdapw=640/sqrt(Fy(i));
lamdarf=141/sqrt(Fy(i)-Fr);
lamdarw=970/sqrt(Fy(i));
Ryy = (II_yy/A)^(1/2);
Sxx = II_xx/(d/2);
de=d-t_w;
Cw = (de^2)*((b_f)^3)*t_f/24;
Lp=300*Ryy/(sqrt(Fy(i))*12)
X1=40014.98*sqrt(J*A)/Sxx;
X2=(4*Cw/II_yy)*((Sxx/(11200*J))^2);
Lr=((Ryy*X1)*sqrt(1+sqrt(1+X2*((Fy(i)-Fr)^2))))/((Fy(i)-Fr)*12)
Cb=1.75+1.05*(M1/M2)+0.3*((M1/M2)^2);
% Cb=1;
Mp=Fy(i)*Zx/12
for Lb=1:1:40

if(lamdaf<lamdapf)
    if(Lb<=Lp)
        Mn=Mp;
    elseif((Lp<Lb)&&(Lb<Lr))
        Mr=(Fy(i)-Fr)*Sxx/12;
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
        Mr=(Fy(i)-Fr)*Sxx/12;
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
       Mr=(Fy(i)-Fr)*Sxx/12;
       Mn=Cb*(Mp-(Mp-Mr)*((Lb-Lp)/(Lr-Lp)));
       if(Mn>Mp)
          Mn=Mp;
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
fig=figure(2); clf; grid on; axis square; hold on;
xlabel('Length (ft)'); ylabel('Design Strength (k-ft)'); title('Effect of Non Compact');
hold on;
plot(L,M);
hold off;