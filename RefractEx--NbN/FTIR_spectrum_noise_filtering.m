clc;        clear;      close all;      
root = cd;  format long   
% path1='FTIR_NbN_Spectra\Transmission';
% path1='FTIR_NbN_Spectra\Ref12';
% path1='FTIR_NbN_Spectra\Veemax\R30';

%% Load FTIR data
path1='FTIR_NbN_Spectra\Ref12';
[wvlen1,Spec1,~]=LoadSpectra (path1,'NbN.spa',1,root);
[wvlen,Spec]=noise_NbN(wvlen1,Spec1);
figure(1);  
plot(wvlen1,Spec1,'-b',wvlen,Spec,'--k','linewidth',1.1);    

function [wv,spec]=noise_NbN(wvlen,Spec)
    ind1=find(wvlen>=7.8 & wvlen<=9);       wvlen(ind1)=[];     Spec(ind1)=[];
    ind1=find(wvlen>=10.3 & wvlen<=12.3);   wvlen(ind1)=[];     Spec(ind1)=[];
    ind1=find(wvlen>=16.2 & wvlen<=19);     wvlen(ind1)=[];     Spec(ind1)=[];
    wv=linspace(wvlen(1),wvlen(end),200);
    spec=smoothdata(interp1(wvlen,Spec,wv),'gaussian',10);
end