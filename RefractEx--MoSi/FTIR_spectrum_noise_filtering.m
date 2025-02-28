clc;        clear;      close all;      
root = cd;  format long   
% path1='FTIR_MoSi_Spectra\Transmission';
% path1='FTIR_MoSi_Spectra\Ref12';
% path1='FTIR_MoSi_Spectra\Veemax\R30';

load('FTIR_MoSi_Spectra.mat');  
wv1=linspace(MoSi15_caf2_IR.wv(1),MoSi15_caf2_IR.wv(end),300); 
nCaF2=interp1(nk_CaF2(1,:),real(nk_CaF2(3,:)),wv1)+1i*interp1(nk_CaF2(1,:),imag(nk_CaF2(3,:)),wv1); 
T0=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.tran,wv1);          
R12=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref12,wv1); 
R30=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref30,wv1);     
R35=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref35,wv1);       
R40=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref40,wv1);     
R45=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref45,wv1);  
R50=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref50,wv1);     
R55=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref55,wv1);  
R60=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref60,wv1);  

% ind1=find(MoSi15_caf2_IR.wv>=16.5);       
% ext4=[zeros(1,length(MoSi15_caf2_IR.wv)-length(ind1)) -5*ones(1,length(ind1))];       
% MoSi15_caf2_IR.ref60=smoothdata(MoSi15_caf2_IR.ref60-ext4,'gaussian',80);
% MoSi15_caf2_IR.ref12=MoSi15_caf2_IR.ref12-22;

%% Load FTIR data
path1='FTIR_MoSi_Spectra\Ref12';
[wvlen1,Spec1,~]=LoadSpectra (path1,'CaF2.spa',1,root);
[wvlen,Spec]=noise_MoSi(wvlen1,Spec1);
figure(1);  
plot(wvlen,Spec,'--k','linewidth',1.1);    
hold on;        plot(wv1,R12,'--b','linewidth',1.1);

function [wv,spec]=noise_MoSi(wvlen,Spec)
    ind1=find(wvlen>=7.8 & wvlen<=9);       wvlen(ind1)=[];     Spec(ind1)=[];
    ind1=find(wvlen>=10.3 & wvlen<=12.3);   wvlen(ind1)=[];     Spec(ind1)=[];
    ind1=find(wvlen>=16.2 & wvlen<=19);     wvlen(ind1)=[];     Spec(ind1)=[];
    wv=linspace(wvlen(1),wvlen(end),200);
    spec=smoothdata(interp1(wvlen,Spec,wv),'gaussian',10);
end