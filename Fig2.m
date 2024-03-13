clear,clc
close all;
syms al be r q Nreal
cc=0.005;
N=10;
p=0.2;
%% pi/4
t=pi/4;
al=cos(t/2);
be=sin(t/2);

%% UPS_1/4
z=0;
for r=0:cc:0.99
    z=z+1;
    jj=0;    
    
    for q=0:0.02:1
        jj=jj+1;
        s=1-r;
        up=al^4*s+be^4*s^(N-1)*(1-q)+2*al^2*be^2*s^(N/2)*(1-q)^(1/2);
        %Success probability of UPS
        P=al^2*(1-q*r)+be^2*(1-q);  %Success probability of UPS defined in Eq.5
        %Fidelity of UPs
        Fid_UPS(jj)=up/P;  %Fidelity of UPS defined in Eq.6
        up2=4*(abs(al*be))^2*(N^2)*(1-q)*(s^N);
        kh=al^2*(1-r)+be^2*(1-q)*s^(N-1);
        QFI_UPS(jj)=up2/(P*kh);
    end
    
    [Max_f,I_f]=max(Fid_UPS);   
    FidUPS4(z)=Max_f;  
    qq_UPSf(z)=(I_f-1)*0.02;
    P_fidp4(z)=al^2*(1- qq_UPSf(z)*r)+be^2*(1- qq_UPSf(z));
end
%% WMRPS_p1_p2_1/4
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.2;
        p_2=0.2;
        
    %Success probability of WMRPS defined in Eq.8
    P_MRp4(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N;  

    %Fidelity of WMRPS defined in Eq.8
    FidMRp4(z)=(1/P_MRp4(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);
end
%% DONOTHING_1/4

z=0;
for r=0:cc:0.99
    z=z+1;

    %Fidelity of Do Nothing defined in Eq.7
    FidDNp4(z)=((abs(al))^4)+((abs(al*be))^2)*(r^N)+2*(abs(al*be))^2*(1-r)^(N/2)+((abs(be))^4)*(1-r)^(N);
end
%% pi/2
t=pi/2;
al=cos(t/2);
be=sin(t/2);
%% UPS_1/2
z=0;
for r=0:cc:0.99
    z=z+1;
    jj=0;    

    for q=0:0.02:1
        jj=jj+1;
        s=1-r;
        up=al^4*s+be^4*s^(N-1)*(1-q)+2*al^2*be^2*s^(N/2)*(1-q)^(1/2);
        %Success probability of UPS
        P=al^2*(1-q*r)+be^2*(1-q);  %Success probability of UPS defined in Eq.5
        %Fidelity of UPs
        Fid_UPS(jj)=up/P;  %Fidelity of UPS defined in Eq.6
        up2=4*(abs(al*be))^2*(N^2)*(1-q)*(s^N);
        kh=al^2*(1-r)+be^2*(1-q)*s^(N-1);
        QFI_NM(jj)=up2/(P*kh);
    end
    
    [Max_f,I_f]=max(Fid_UPS);   
    FidUPS2(z)=Max_f;  
    qq_UPSf(z)=(I_f-1)*0.02;
    P_fidp2(z)=al^2*(1- qq_UPSf(z)*r)+be^2*(1- qq_UPSf(z));
end
%% WMRPS_p1_p2_1/2
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.2;
        p_2=0.2;

    %Success probability of WMRPS defined in Eq.8
    P_MRp2(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N;  

    %Fidelity of WMRPS defined in Eq.8
    FidMRp2(z)=(1/P_MRp2(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);
end
%% DONOTHING_1/2
z=0;
for r=0:cc:0.99
    z=z+1;

    %Fidelity of Do Nothing defined in Eq.7
    FidDNp2(z)=((abs(al))^4)+((abs(al*be))^2)*(r^N)+2*(abs(al*be))^2*(1-r)^(N/2)+((abs(be))^4)*(1-r)^(N);
end
%% 3pi/4
t=3*pi/4;
al=cos(t/2);
be=sin(t/2);
%% UPS_3/4
z=0;
for r=0:cc:0.99
    z=z+1;
    jj=0;    

    for q=0:0.02:1
        jj=jj+1;
        s=1-r;
        up=al^4*s+be^4*s^(N-1)*(1-q)+2*al^2*be^2*s^(N/2)*(1-q)^(1/2);
        %Success probability of UPS
        P=al^2*(1-q*r)+be^2*(1-q);  %Success probability of UPS defined in Eq.5
        %Fidelity of UPs
        Fid_UPS(jj)=up/P;  %Fidelity of UPS defined in Eq.6
        up2=4*(abs(al*be))^2*(N^2)*(1-q)*(s^N);
        kh=al^2*(1-r)+be^2*(1-q)*s^(N-1);
        QFI_NM(jj)=up2/(P*kh);
    end

    [Max_f,I_f]=max(Fid_UPS);   %find the optimal strength
    FidUPS3p4(z)=Max_f;  
    qq_UPSf(z)=(I_f-1)*0.02;
    P_fid3p4(z)=al^2*(1- qq_UPSf(z)*r)+be^2*(1- qq_UPSf(z));
end
%% WMRPS_p1_p2_3/4
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
    p_1=0.2;
    p_2=0.2;

    %Success probability of WMRPS defined in Eq.8
    P_MR3p4(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N; 

    %Fidelity of WMRPS defined in Eq.8
    FidMR3p4(z)=(1/P_MR3p4(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);
end
%% DONOTHING_3/4
z=0;
for r=0:cc:0.99
    z=z+1;

    %Fidelity of Do Nothing defined in Eq.7
    FidDN3p4(z)=((abs(al))^4)+((abs(al*be))^2)*(r^N)+2*(abs(al*be))^2*(1-r)^(N/2)+((abs(be))^4)*(1-r)^(N);
end
%% PLOT
figure(1)
tf=0:cc:0.99;
L1=plot(tf,P_fidp4,'m','LineStyle','-.','LineWidth',1.5); hold on
L2=plot(tf, P_MRp4,'b','LineStyle','-.','LineWidth',1.5); hold on
L3=plot(tf,P_fidp2,'m','LineWidth',1.5); hold on
L4=plot(tf, P_MRp2,'b','LineWidth',1.5); hold on
L5=plot(tf,P_fid3p4,'m--','LineWidth',1.5); hold on
L6=plot(tf, P_MR3p4,'b--','LineWidth',1.5); hold on
legend([L1,L2,L3,L4,L5,L6],'(\eta=\pi/4)P_{UPS}','(\eta=\pi/4)P_{WMRPS}','(\eta=\pi/2)P_{UPS}','(\eta=\pi/2)P_{WMRPS}','(\eta=3\pi/4)P_{UPS}','(\eta=3\pi/4)P_{WMRPS}');
xlabel('r')
ylabel('Success probability')
ylim([0 1]);
grid on

figure(2)
tf=0:cc:0.99;
L7=plot(tf, FidUPS4,'m','LineStyle','-.','LineWidth',1.5); hold on
L8=plot(tf, FidMRp4,'b','LineStyle','-.','LineWidth',1.5); hold on
L9=plot(tf, FidDNp4,'k','LineStyle','-.','LineWidth',1.5); hold on
L10=plot(tf, FidUPS2,'m','LineWidth',1.5); hold on
L11=plot(tf, FidMRp2,'b','LineWidth',1.5); hold on
L12=plot(tf, FidDNp2,'k','LineWidth',1.5); hold on
L13=plot(tf, FidUPS3p4,'m--','LineWidth',1.5); hold on
L14=plot(tf, FidMR3p4,'b--','LineWidth',1.5); hold on
L15=plot(tf, FidDN3p4,'k--','LineWidth',1.5); hold on

legend([L13,L14,L15],'Fid_{UPS}','Fid_{WMRPS}','Fid_{DN}');

% Draw a circle to indicate the grouping for t=\pi/4
ellipsePosition1 = [0.444642857142857,0.678571428571429,0.1,0.16];
annotation('ellipse', ellipsePosition1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
annotation('textbox', 'String', '\eta=\pi/4', 'Position', [ellipsePosition1(1)+ellipsePosition1(3)+0.02, ellipsePosition1(2), 0.1, 0.1], 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

% Draw a circle to indicate the grouping for t=\pi/2
ellipsePosition2 = [0.5,0.2678571428571429,0.15,0.3];
annotation('ellipse', ellipsePosition2, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
annotation('textbox', 'String', '\eta=\pi/2', 'Position', [ellipsePosition2(1)-ellipsePosition1(3)+0.02, ellipsePosition2(2)+ellipsePosition2(4)-0.1, 0.1, 0.1], 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

% Draw a circle to indicate the grouping for t=\3pi/4
ellipsePosition3 = [0.28,0.15,0.1,0.18];
annotation('ellipse', ellipsePosition3, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
annotation('textbox', 'String', '\eta=3\pi/4', 'Position', [ellipsePosition3(1)-ellipsePosition3(3)+0.02, ellipsePosition3(2)-0.01, 0.1, 0.1], 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');


xlabel('r')
ylabel('Fidelty')
ylim([0 1]);
grid on
