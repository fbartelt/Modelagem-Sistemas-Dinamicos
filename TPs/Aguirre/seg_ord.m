function [numd,dend]=seg_ord(kdc,qsi,wn);
% function [numd,dend]=seg_ord(kdc,qsi,wn)
% Monta a funcao de transferencia de um sistema
% de segunda ordem dados o kdc, qsi e o wn

% Luis A. Aguirre - PPGEE
% B.H. 13/07/90

[a1,b1,c1,d1]=ord2(wn,qsi);
[numd,dend]=ss2tf(a1,b1,c1,d1,1);
ji=sqrt(-1);
numd=numd*(kdc*dend(length(dend))/(numd(length(numd))));
dend=dend;


