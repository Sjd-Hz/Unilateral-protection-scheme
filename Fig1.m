clear,clc
close all;

cc=0.005;
N=10;
t=pi/2;
al=cos(t/2);
be=sin(t/2);
%% Do Nothing
z=0;
for r=0:cc:0.99
    z=z+1;
    
    %Fidelity of Do Nothing 
    FidDN(z)=((abs(al))^4)+((abs(al*be))^2)*(r^N)+2*(abs(al*be))^2*(1-r)^(N/2)+((abs(be))^4)*(1-r)^(N); %Fidelity of Do Nothing defined in Eq.7 
end
%% WMRPS_p1_p2=0.1
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.1;
        p_2=0.1;
    %Success probability of WMRPS when p1=p2=0.1
    P_MR1(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N; %Success probability of WMRPS defined in Eq.8
    %Fidelity of WMRPS when p1=p2=0.1
    FidMR1(z)=(1/P_MR1(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);    %Fidelity of WMRPS defined in Eq.8
end
%% WMRPS_p1_p2=0.2
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.2;
        p_2=0.2;
    %Success probability of WMRPS when p1=p2=0.2
    P_MR2(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N;  %Success probability of WMRPS defined in Eq.8
    %Fidelity of WMRPS when p1=p2=0.2
    FidMR2(z)=(1/P_MR2(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);    %Fidelity of WMRPS defined in Eq.8
end
%% WMRPS_p1_p2=0.6
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.6;
        p_2=0.6;
    %success probability of WMRPS when p1=p2=0.6
    P_MR6(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N;  %Success probability of WMRPS defined in Eq.8
    %fidelity of WMRPS when p1=p2=0.6
    FidMR6(z)=(1/P_MR6(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);  %Fidelity of WMRPS defined in Eq.8
end

%% WMRPS_p1_p2=0.9
z=0;
for r=0:cc:0.99
    z=z+1;
    s=1-r;
        p_1=0.9;
        p_2=0.9;
    %Success probability of WMRPS when p1=p2=0.9
    P_MR9(z)=abs(al)^2 * (1 - p_2)^N + abs(be)^2 * (1 - p_1)^N * (r * (1 - p_2) + s)^N;  %Success probability of WMRPS defined in Eq.8
    %Fidelity of WMRPS when p1=p2=0.9
    FidMR9(z)=(1/P_MR9(z))*(abs(al)^4 * (1 - p_2)^N + ...
                    abs(al*be)^2 * (1 - p_1)^N * r^N * (1 - p_2)^N + ...
                    2 * abs(al*be)^2 * (1 - p_1)^(N/2) * s^(N/2) * (1 - p_2)^(N/2) + ...
                    abs(be)^4 * (1 - p_1)^N * s^N);  %Fidelity of WMRPS defined in Eq.8
end
%% UPS
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
        Fid_NM(jj)=up/P;  %Fidelity of UPS defined in Eq.6
        up2=4*(abs(al*be))^2*(N^2)*(1-q)*(s^N);
        kh=al^2*(1-r)+be^2*(1-q)*s^(N-1);
        QFI_NM(jj)=up2/(P*kh);
    end
    
    [Max_f,I_f]=max(Fid_NM); 
    [Max_Q,I_Q]=max(QFI_NM); 
   
    FidNM(z)=Max_f;
    FishNM(z)=Max_Q;
   
    qq_NMf(z)=(I_f-1)*0.02;
    qq_NMR(z)=(I_Q-1)*0.02;
    P_fid(z)=al^2*(1- qq_NMf(z)*r)+be^2*(1- qq_NMf(z));
end

%% PLOT
figure
tf=0:cc:0.99;

% Plot FidMR1\FidMR2\FidMR9 curves
L1 = plot(tf, FidMR1, 'r--', 'LineWidth', 2); hold on  %p1=p2=0.1 fidelity of WMPRS
%L11=plot(tf,FidNM,'m','LineWidth',1.5); hold on
L2 = plot(tf, FidMR2, 'b:', 'LineWidth', 2); hold on  %p1=p2=0.2 fidelity of WMPRS
L3 = plot(tf, FidMR9, 'm-.', 'LineWidth', 2); hold on  %p1=p2=0.9 fidelity of WMPRS
%L4 = plot(tf, FidMR6, 'b:', 'LineWidth', 2); hold on
LD = plot(tf, FidDN, 'k', 'LineWidth',1.5); hold on  %fidelity of Do Nothing

grid on
% Fill the area between FidMR1 and FidMR9
fill([tf, fliplr(tf)], [FidMR1, fliplr(FidMR9)], 'g', 'FaceAlpha', 0.05);

xlabel('r')
ylabel('Fidelity')
axis tight;
ylim([0 1]);

legend([L1, L2, L3,LD], 'Fid_{WMRPS_{0.1}}', 'Fid_{WMRPS_{0.2}}','Fid_{WMRPS_{0.9}}', 'Fid_{DN}');

figure(2)
% Plot P_MR1 P_MR2 P_MR9 curves (Success proobability)
L1 = plot(tf, P_MR1, 'r--', 'LineWidth', 2); hold on  %p1=p2=0.1 success probability of WMPRS
L2 = plot(tf, P_MR2, 'b:', 'LineWidth', 2); hold on   %p1=p2=0.1 success probability of WMPRS
L3 = plot(tf, P_MR9, 'm-.', 'LineWidth', 2); hold on   %p1=p2=0.1 success probability of WMPRS
%L4 = plot(tf, P_MR6, 'b:', 'LineWidth', 2); hold on
% Fill the area between P_MR1 and P_MR9
fill([tf, fliplr(tf)], [P_MR1, fliplr(P_MR9)], 'g', 'FaceAlpha', 0.05);

ylim([0 1]);
xlabel('r')
ylabel('Success probability')

grid on
legend([L1, L2,L3], 'P_{WMRPS_{0.1}}','P_{WMRPS_{0.2}}','P_{WMRPS_{0.9}}');
