clear all;
clc;
Alpha_c=0.9; Tao_g=0.95; Alpha_f=0.6;Alpha_b=0.8;Mf=20;
Am=0.605; Ab=2; Acond=1.5*0.25;
Arm=0.25; Arc=0.75; Aam=0.5; Aac=1.5;
Ac=1.395;
hi_dash=5.8;
Beta_c=0.25;
X=0.33;
etao=0.15;
Beta=0.89;
Alpha_p=0.80;
As=1; 
g=9.8;
Vf=0.01;
Neta_m=0.13;
Apha_f=0.6;
Tao_c=0.95;
Li=0.1;
Ki=0.166;
Di=0.0125;
Neta_o=0.15;
Beta_o=0.0045;
beta=0.89;
M_dot_f=0.04;
LT=2;
N=4;
V=7.11;
Delta_t=3600;
ho = 9.5;
%ho=5.7+3.8*V;
h_i = 5.7;
%h_i=2.8+3*V;
Lg=0.0035;
Kg=0.816;
%Kg=0.780;
econd=0.95;
ef=0.95;
UA_SL=0;
hi_dash=5.8;
F_dash=0.968;
M_dot_cond=0.001;
sigma = (5.67).*(10^(-8));

KP=3000; dP=(1.6)*10^(3);
CP=796;phiP = 0.01; 

Kp=6;  Kb=0.035; Lb=0.005; Kcond =210;Lcond=0.05;

% %  May. Month

Ib=[288.64	471.72	622.98	731.31	782.58	773.48	710.86	595.71	426.51	239.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%incident solar intensity on solar still (W/m2)
Ie=[572.876	757.4147	861.6151	894.4769	853.1213	748.7437	599.0055	411.9275	214.2461	129.1556 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];%incident solar intensity on solar still (W/m2)
Iw=[157.5443	331.8115	524.2421	707.2609	853.1213	946.7553	983.9697	949.3991	827.9625	632.8997 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%incident solar intensity on solar still (W/m2)
Ta=[30.8	30.8	30.8	30.1	30.6	31.8	33.8	35.3	36.6	37.6	38.5	40.3	40	38.6	38	38	37.3	36	34.1	31	28.1	27	27.6	28]; %ambient temp. in centrigrde   
  

for i=1:24
    Tf0di(1)=26.8;
    Tf(1) = 27;   
    Tcond(1)=26;
    Tf0d(1)=Tf0di(1);
    hba=((Li/Ki)+(1/h_i))^(-1);
    It(i)=(Ie(i)+Iw(i))/2;
    Pf(i)=exp((25.317)-(5144./(273+Tf(i))));                % partial vapour pressure at condesing cover
    Pcond(i)=exp((25.317)-(5144./(273+Tcond(i))));          % partial vapour pressure at fluid temp.
    eff=(((1/ef)+(1/econd)-1)^(-1));                        % from 'internal HTC from fluid surface to internal surface of the glass cover'

    % Total Internal Energy inside Solar Still
    hc(i) = (0.884).*((Tf(i)-Tcond(i))+(((Pf(i)-Pcond(i)).*(Tf(i)+273))./((268.9).*(10^(3))-Pf(i))))^(1/3);         % Convective HTC
    he(i)=(hc(i)).*(0.016273).*(Pf(i)-Pcond(i))./(Tf(i)-Tcond(i));                                                  % Evaporate HTC
    hr(i)=sigma.*(eff).*((((Tf(i)+273)^(2))+(Tcond(i)+273)^(2))).*(Tf(i)+Tcond(i)+546);                             % Radioactive HTC
    h_1(i) = hc(i) + he(i) + hr(i);                                                                                 % Total Internal Energy


    ho=5.7+3.8*7;
    h_i=2.8+3*4;
    hi_dash=2.8+3*7;
    %ho = 9.5;
    %h_i = 5.7;
    %hi_dash = 5.8;
    U_tcond_a=(((1./ho)+(Kcond./Lcond))^(-1));
    Tv(i)=(Tf(i)+Tcond(i))./2;

    % Latent heat of vaporization of fluid
    %if Tv(i)>70
    if Tv(i)>70
    %       Doubt
        Lv(i) = 3.1625.*(10^3)+(1-(7.616).*(10^(-4)).*Tv(i)); 
    %else Tv(i)<70
    else Tv(i)<70
        Lv(i) = (2.4935.*(10^(3))).*((1-(9.4779.*(10^(-4)).*Tv(i)))+(1.3132.*(10^(-7)).*((Tv(i))^2)) - (4.7974.*(10^(-3)).*(Tv(i)^3)));
    end
    % Specific Heat of base fluid
    Cf(i) = 4.217 - 0.00561*(Tf0d(i)) + 0.00129*((Tf0d(i))^(1.5)) - (0.000115*(Tf0d(i))^2) + 4.149*10^(-6)*(Tf0d(i))^(2.5);
    Cf(i) = Cf(i)*1000;
    nanonanoCf(i) = phiP*CP + (1-phiP)*Cf(i);
    % Thermal Conductivity of basefluid
    Kf(i) = 0.565 +0.00263*(Tf0d(i)) - 0.000125*((Tf0d(i))^(1.5)) - 1.515*10^(-6)*((Tf0d(i))^2) - 0.000941*((Tf0d(i))^(0.5));
    nanoKf(i) = Kf(i)*...
            (1 - phiP + (2*phiP*(KP/(KP - Kf(i))))*log((KP + Kf(i))/2*Kf(i)))/...
            (1 - phiP + (2*phiP*(Kf(i)/(KP - Kf(i))))*log((KP + Kf(i))/2*Kf(i)));


    Tcond(i) = (M_dot_f*Tf(i)-M_dot_cond*Lv(i)+U_tcond_a*Ta(i)*Acond)./(M_dot_f*nanonanoCf(i)+U_tcond_a*Acond);

    Tcond(i+1) = Tcond(i);
    U_ba = (((1./h_i)+(Kb./Lb))^(-1));
    % Viscosity of nanafluid
    Mu_f(i)=1./(557.82-19.408.*(Tf0d(i))+((0.136).*(Tf0d(i)^2))-((3.116).*(Tf0d(i)^3).*(10^(-4))));  

    nanouf(i) = Mu_f(i)*(1-phiP)^(-2.5); 
    
    % Density of BaseFluid
    Rho_f(i)=999.79+0.0683.*(Tf0d(i))-0.0107.*(Tf0d(i)^2)+0.00082.*(Tf0d(i)^(2.5))-2.303*(10^(-5)).*(Tf0d(i)^3);
    nanodf(i) = phiP*dP + (1-phiP )*Rho_f(i);
    Beta_f(i)=1./(Tf0d(i)+273.15);

    %           Doubt

    %Grf(i)=((1).*(g.*(Beta_f(i)).*(Rho_f(i))^(2)).*(X^3).*(dT(i)))./((Mu_f(i))^2);

    % Reynold Number
    Ref(i)=(nanodf(i).*Vf.*X)./nanouf(i);
    % Prandlit Number
    Prf(i)=((nanouf(i)).*(nanonanoCf(i)))./(nanoKf(i));
    n=1/4;
    c=0.54;
    % Prandlit Number
    Nuf(i)=c.*((Ref(i).*(Prf(i)))^(n));

    h_bf(i)=(Nuf(i)*nanoKf(i))./X;


    Tb(i) = ((1-Alpha_f)*(Alpha_b)*((Tao_g)^2)*(1-Beta_c)*It(i)*Am +(h_bf(i)*Tf(i)+U_ba*Ta(i))*Ab)./(h_bf(i)*Ab+U_ba*Ab);

    U_tc_a=(((1/ho) + (Lg/Kg))^(-1));
    %U_bc_f=(((1./h_i)+(Kg./Lg))^(-1));
    
    %           Doubt
    U_bc_f=(((1/h_i)+(Lg/Kg))^(-1));



    Neta_c(1)=Neta_m./(Tao_g.*Beta_c);


    Tc(i)=((Alpha_c*Tao_g*Beta_c*It(i)*Am)-(Neta_c(i)*Tao_g*Beta_c*It(i)*Am)+(U_tc_a*Ta(i)+ U_bc_f*Tf(i))*Am)./((U_tc_a+U_bc_f)*Am);

    Neta_c(i+1)=Neta_o*(1-Beta_o*(Tc(i)-Ta(i)));
    % Grf(i)=((1).*(g.*(Beta_f(i)).*(Rho_f(i)^(2)).*(X^3).*(dT(i))))./((Mu_f(i))^2);
    
    % Prandlt number
    Prf(i)=((nanouf(i)).*(nanonanoCf(i)))./(nanoKf(i));

    % Nusslet Number 
    Nuf(i)=c.*(((Ref(i)).*(Prf(i)))^(n));
    hcf(i)=((Nuf(i)).*(X))./nanoKf(i);
    
    % HTC (convective CPC   ) 
    Kfc(i)=0.565+0.00263*(Tf0di(i))-0.000125*(Tf0di(i)^(1.5))-1.515*(10^(-6))*(Tf0di(i)^2)-0.000941*(Tf0di(i)^(0.5));
    nanoKfc(i) = Kfc(i)*...
            (1 - phiP + (2*phiP*(KP/(KP - Kfc(i))))*log((KP + Kfc(i))/2*Kfc(i)))/...
            (1 - phiP + (2*phiP*(Kfc(i)/(KP - Kfc(i))))*log((KP + Kfc(i))/2*Kfc(i)));

    Mu_fc(i)=1./(557.82-19.408.*(Tf0d(i))+((0.136).*(Tf0d(i)^2))-((3.116).*(Tf0d(i)^3).*(10^(-4))));  
    nanoufc(i) =  Mu_fc(i)*(1-phiP)^(-2.5);
    
    Rho_fc(i)=999.79+0.0683.*(Tf0d(i))-0.0107.*(Tf0d(i)^2)+0.00082.*(Tf0d(i)^(2.5))-2.303*(10^(-5)).*(Tf0d(i)^3);
    nanodf(i) = phiP*dP + (1-phiP )*Rho_fc(i);
    
    C_fc(i)=(10^(3)).*(4.217- 0.00561.*(Tf0di(i))+(0.00129).*(Tf0di(i)^(1.5))-(0.000115).*(Tf0di(i)^(2))+(4.149).*((10)^(-6)).*(Tf0di(i)^(2.5)));   %SpecifIb Heat
    nanonanoCfc(i) = phiP*CP + (1-phiP)* C_fc(i);
    Beta_fc(i)=1./(Tf0di(i)+273.15);
    
    Gzm(i)=(M_dot_f*nanonanoCfc(i))./(nanoKfc(i).*LT);  
    Grm(i)=(6.*(Di^2).*g)./(nanoufc(i)^2);                  % Delta(T) is taken to be 6
    Prnf(i)=(nanoufc(i).*nanonanoCfc(i))./nanoKfc(i);
    Renf(i)=4*M_dot_f/(3.14*Di*nanoufc(i));
    f=(0.79*(log(Renf(i))-1.64))^(-2);
    %Nuc(i)=(1.75).*(Gzm(i)+0.0083.*((Grm(i).*Prm(i))^(0.75)))^(1/3);
    Nuc(i)=(f/8)*(Renf(i)-1000)*Prnf(i)./(1+12.7*((f/8)^.5)*(((Prnf(i))^(2/3))-1));
    h_p_f(i)=(Nuc(i).*nanoKfc(i))./(Di);
    h_p_f(i)=h_bf(i);

    U_ba=(((1./h_i)+(Kb./Lb))^(-1));
    U_bc_f=(((1./h_i)+(Kg./Lg))^(-1));
    U_tc_p=((Lg/Kg)+(1/h_i))^(-1); 
    U_tc_a=((Lg/Kg)+(1/ho))^(-1);
    U_tp_a(i) = ((1./h_p_f(i))+(1./(hi_dash))+(Li./Ki))^(-1) + ((1./(U_tc_a))+(1./(U_tc_p)))^(-1);
    %U_tp_a(i)=((((1/ho)+(1/h_i)+(Lg/Kg))^(-1))+((1./h_p_f(i))+(Li./Ki)+(1./hi_dash))^(-1));
    PF_c(i)=h_p_f(i)./(F_dash*h_p_f(i)+U_tp_a(i));
    PF_1=U_tc_p./(U_tc_a+U_tc_p);
    U_l_1=((U_tc_a.*U_tc_p)./(U_tc_a+U_tc_p));
    U_l_2(i)=U_l_1+U_tp_a(i);
    PF_2(i)=h_p_f(i)./(F_dash.*h_p_f(i)+U_l_2(i));
    U_l_m(i)=((U_l_2(i)*h_p_f(i))./(U_l_2(i)+F_dash*h_p_f(i)));
    U_l_c(i)=((U_tp_a(i).*h_p_f(i))./(U_tp_a(i)+F_dash.*h_p_f(i)));

    AlphaTao_2eff=(Tao_g^2).*Alpha_p.*(1-Beta).*(Aam./Arm);
    AlphaTao_ceff(i)=(PF_c(i).*Alpha_p.*Tao_g.*(Arc./Aac));
    AcFrc(i)=((M_dot_f*nanonanoCf(i))./(U_l_c(i)))*(1-exp(-(F_dash*(U_l_c(i))*Ac)/(M_dot_f*nanonanoCf(i))));
    AmFrm(i)=(((M_dot_f*nanonanoCf(i))./U_l_m(i))*(1-exp(-(F_dash*U_l_m(i)*Am./(M_dot_f*nanonanoCf(i))))));
    ArFrUl_1(i)=(AcFrc(i)*U_l_c(i))+(AmFrm(i)*U_l_m(i)*(1-(AcFrc(i)*U_l_c(i))./(M_dot_f*nanonanoCf(i))));
    AFrUlm(i)=AmFrm(i)*U_l_m(i); 


    K_k(i)=1-((ArFrUl_1(i))./(M_dot_f*nanonanoCf(i)));
    Km(i)=(1-((AFrUlm(i))/(M_dot_f*nanonanoCf(i))));


    A1(i)=((1/(U_tc_a+U_tc_p)).*It(i));         
    A2(i)=((Alpha_c*Tao_g*beta)+(U_tc_p/(U_l_2(i)+h_p_f(i)))*(AlphaTao_2eff+(PF_1*Alpha_c*Tao_g*beta)));
    A3=((U_tc_p*h_p_f(i))/(2*M_dot_f*nanonanoCf(i)*(U_l_2(i)+h_p_f(i))));
    A4(i)=(((PF_2(i)*AlphaTao_2eff)+(PF_2(i)*PF_1*Alpha_c*Tao_g*beta))*AmFrm(i));
    A5(i)=((1+Km(i))*((AcFrc(i)*AlphaTao_ceff(i))+(((PF_2(i)*AlphaTao_2eff)+(PF_2(i)*PF_1*Alpha_c*Tao_g*beta))*AmFrm(i)*(1-((AcFrc(i)*U_l_c(i))/(M_dot_f*nanonanoCf(i)))).*((1-(K_k(i))^(N-1))/(1-K_k(i))))));
    A6(i)=A3.*(A4(i)+A5(i));
    A(i)=A1(i)*(A2(i)+A6(i));
    B(i)=((1/(U_tc_a+U_tc_p))*((U_tc_a+((U_tc_p*U_l_2(i))/(U_l_2(i)+h_p_f(i)))+((U_tc_p*h_p_f(i))/(2*M_dot_f*nanonanoCf(i)*(U_l_2(i)+h_p_f(i))))*(AFrUlm(i)+((1+Km(i))*ArFrUl_1(i).*((1-(K_k(i))^(N-1))/(1-K_k(i)))))).*Ta(i)));
    C(i)=((1/(U_tc_a+U_tc_p))*((((U_tc_p*h_p_f(i))/(2*(U_l_2(i)+h_p_f(i))))*(1+Km(i)).*(K_k(i)^(N-1))).*Tf0di(i)));
    D(i)=(Beta_o.*((A(i)+B(i)+C(i))-25));
    Y(i)=(etao.*(1-D(i)));
    Z1(i)=((etao*Beta_o.*It(i))/(U_tc_a+U_tc_p));
    Z2(i)=((Tao_g*beta)+((U_tc_p*PF_1*Tao_g*beta)/(U_l_2(i)+h_p_f(i))));
    Z3(i)=(((U_tc_p*h_p_f(i))/(2*M_dot_f*nanonanoCf(i)*(U_l_2(i)+h_p_f(i))))*((PF_2(i)*PF_1*Tao_g*beta*AmFrm(i))+((1+Km(i))*PF_2(i)*PF_1*Tao_g*beta*AmFrm(i)*(1-((AcFrc(i)*U_l_c(i))/(M_dot_f*nanonanoCf(i)))).*((1-K_k(i)^(N-1))/(1-K_k(i))))));
    Z(i)=(1-(Z1(i).*(Z2(i)+Z3(i))));
    Netac(i)=(Y(i)./Z(i));


    AlphaTao_1eff(i)=Tao_g*Beta_c*(Alpha_c-Netac(i))*(Aam./Arm);
    AlphaTao_ceff(i)=(PF_c(i)*Alpha_p*Tao_g*(Arc./Aac));
    AlphaTao_meff(i)=AlphaTao_2eff+PF_1*AlphaTao_1eff(i);

    AFrAplhaTao_1(i)=(AcFrc(i).*AlphaTao_ceff(i)+ (PF_2(i).*AlphaTao_meff(i)+AmFrm(i))*(1-(AcFrc(i)*U_l_c(i))./(M_dot_f*nanonanoCf(i))));
    Alpha_T_Eff(i)=AFrAplhaTao_1(i)*(1-((K_k(i))^N))./(1-K_k(i));
    Frm=M_dot_f*nanonanoCf(i).*(1-exp((-1)*F_dash*U_l_m(i)*Arm./(M_dot_f*nanonanoCf(i))))./(U_l_m(i)*Arm);


    UA_eff(i)=ArFrUl_1(i)*(1-((K_k(i))^N))./(1-K_k(i));

    A1(i) = ((M_dot_f*nanonanoCf(i)*U_tcond_a*Acond)./((M_dot_f*nanonanoCf(i))+U_tcond_a*Acond));
    A2(i) = (h_bf(i)*U_ba*Ab)./(h_bf(i)+U_ba);
    A3(i) = (UA_SL)*As + UA_eff(i);
    A4(i) = (h_1(i)*Am*U_tc_a)./(U_tc_a+U_bc_f);
    a(i) = (A1(i)+A2(i)+A3(i)+A4(i))./(M_dot_f*nanonanoCf(i));

    F1(i) = h_bf(i).*((((1-Alpha_f).*Alpha_b.*(Tao_g^2).*(1-Beta_c).*It(i).*Am)+(Ab*Ta(i).*U_ba))./(h_bf(i)+U_ba));
    F2(i) = Alpha_f.*(Tao_g^2).*(1-Beta_c).*It(i)*Am+UA_SL.*Ta(i)*As+Alpha_T_Eff(i)*Ib(i)+UA_eff(i)*Ta(i);
    F3(i) = (M_dot_f*nanonanoCf(i)).*((U_tcond_a*Ta(i)*Acond - M_dot_cond*Lv(i))./(M_dot_f*nanonanoCf(i)+U_tcond_a*Acond));
    F4(i) = (h_1(i)*Am).*(((Alpha_c-Neta_c(i)).*Tao_g.*Beta_c.*It(i)+U_tc_a*Ta(i))./(U_tc_a+U_bc_f));
    f_t(i) = (F1(i)+F2(i)+F3(i)+F4(i))./(M_dot_f*nanonanoCf(i));

    Tf(i)=(((f_t(i)/a(i)).*(1-(exp(-a(i).*(3600)))))+(Tf0d(i).*(exp((-a(i)).*3600))));          % Here Delta(T) is 3600
    Tf(i+1) = Tf(i);
    TT1(i)=(AFrAplhaTao_1(i)).*(1-K_k(i)^N).*(Ib(i))./((1-K_k(i))*(M_dot_f*nanonanoCf(i)));
    TT2(i)=(ArFrUl_1(i))*(Ta(i))*(1-K_k(i)^N)./((M_dot_f*nanonanoCf(i)).*(1-K_k(i)));
    TFNN(i)=TT1(i)+TT2(i)+(K_k(i)^N).*Tf0di(i);
    Tf0d(i+1)=TFNN(i);
    Tf0di(i+1)=Tf(i);



    % Thermal Gain
    c1(i) = (It(i).*Beta_c).*(Alpha_c.*Tao_g.*Am-Neta_c(i).*Tao_g)./(U_tc_a+U_bc_f);
    Neta_gth(i) = ((h_1(i))./(It(i).*Am)).*(((Am.*U_tc_a.*(Tf(i)-Ta(i)))./(U_tc_a+U_bc_f))-c1(i));
    

    % Thermal Loss
    Neta_cth(i) = ((M_dot_f*nanonanoCf(i))./(It(i))).*((f_t(i)./a(i))-Tf0di(i)).*(1-exp(-a(i).*(3600)));           % Delta(T) is taken to be 3600



end



Tf
TFNN
x = 1:23;

y = Tf(x);
plot1 = plot(x,y)


hold on

y = Tc(x);
plot2 = plot(x,y)

hold on

y = Tb(x);
plot3 = plot(x,y)

ylabel("Temperture,(°C)")
xlabel('Time(Hours)')
title('Plot of Variation of Different Temp with Time')
legend([plot1,plot2,plot3],'Base Fluid','Solar Cell','Basin');




y = Rho_f(x);
figure;
plot(x,y,'-o','MarkerIndices',1:1:length(y))
title('Variation of Density of Basefluid')
xlabel('Time(Hours)')
ylabel("Density,kg/m3")

y = Kf(x);
figure;
plot(x,y,'-o','MarkerIndices',1:1:length(y))
title('Thermal Conductivity of Basefluid')
xlabel('Time(Hours)')
ylabel("Thermal Conductivity,(W/mK)")

y = Cf(x);
figure;
plot(x,y,'-o','MarkerIndices',1:1:length(y))
title('Specific Heat of Basefluid')
xlabel('Time(Hours)')
ylabel("Specific Heat,J/(kg °C)")



Neta_gth
Neta_cth

i = 1:23;
y = Neta_gth(i);
x = (Tf0d(i)-Ta(i))./It(i);
plot1 = plot(x,y)

hold on

y = Neta_cth(i);
plot2 = plot(x,y)
title("Thermal Gain & Loss of Base Fluid")
xlabel('(Tf0-Ta)/It')
ylabel("Thermal Gain & Loss")
legend([plot1,plot2],'Thermal Gain','Thermal Loss')






