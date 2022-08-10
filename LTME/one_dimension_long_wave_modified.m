%REMARK LONG NONLINEAR 1-D WAVE MODEL PROPAGATION OVER 
%CONSTANT DEPTH CENTERED DIFFERENCES STAGGERED GRID

% DT,DX = TIME AND SPACE DISCRETIZATION STEPS 
% IM,NM = MAX VALUES OF SPACE(I) AND TIME(N) STEPS 
% Z0,PER = AMPLITUDE AND PERIOD OF INCOMING WAVES 
% Z = FREE SURFACE ELEVATION WITH RESPECT TO THE SWL 
% H0 = INITIAL WATER DEPTH 
% H,HN = PAST AND PRESENT WATER DEPTHS 
% U,UN = PAST AND PRESENT VELOCITY VALUES

clear all


%fid_read=fopen('input.dat', 'r');

%celldata=textscan(fid_read, '%f', 6)
%data=celldata{1}
%DT=data(1);
%DX=data(2);
%IM=data(3);
%NM=data(4);
%Z0=data(5);
%PER=data(6);

DT=1
DX=300
IM=5
NM=120
Z0=100
PER=12

IX=IM;
%IX=3;
%celldata2=textscan(fid_read, '%f');
%H0(1:51)=celldata2{1}
H0(1:51)=240
H0(2:51)=470
for i=1:IX
    H(i)=H0(i)
end
T=0.0
G=9.81
C=sqrt(G*H(1))
L=C*PER
N=1  
U(1:IX)=0
Z(1:IX)=0
Z(2:IX)=150
HN(1:IX)=0

%%if(mod(n,20)/=0)
    

while (mod(N,20)||(N<=NM))
    T=T+DT
    
    if((N-1)*DT > DX/C)
        Z1=Z(1)-Z0*sin(2*pi*(N)*DT/PER)
        Z2=Z(2)-Z0*sin(2*pi*(N)*DT/PER-DX/L)
        Z1=Z1+DT/DX*C*(Z2-Z1)
        Z(1)=Z0*sin(2*pi*(N+1)*DT/PER)+Z1-3
        HN(1)=H0(1)+Z(1)
    else
        Z1=0
        Z(1)=Z0*sin(2*pi*(N+1)*DT/PER)+Z1-3
        HN(1)=H0(1)+Z(1)
    end
    for i=2:(IX-1)
        HN(i)=H(i)-DT/DX/2*((H(i)+H(i+1))*U(i+1)-(H(i)+H(i-1))*U(i)) % CONTINUITY EQUATION
        Z(i)=HN(i)-H0(i)
        UN(i)=U(i)-DT/DX/8*((U(i+1)+U(i)).^2-(U(i)+U(i-1)).^2)-G*DT/DX*(Z(i)-Z(i-1))  %EQUILIBRIUM EQUATION     
    end
    
    UN(1)=UN(2)
	UN(IX)=Z(IX-1)*sqrt(G/HN(IX-1))
    
    for i=1:IX
        H(i)=HN(i)
        U(i)=UN(i)
    end
    for i=1:IX
        finalu(N,i)=UN(i)
        finalh(N,i)=HN(i)
    end 
    N=N+1
end
figure;
plot(1:24, finalh(1:24,1))
title('model (1)')
xlabel('time')
ylabel('tide(m)')
figure;
%finalh(1:24,2)=finalh(1:24,2)*2
plot(1:24, finalh(1:24,2))
title('model(2)')
xlabel('time')
ylabel('tide(m)')

%fclose(fid_read);
