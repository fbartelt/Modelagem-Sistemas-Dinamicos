function [num,den]=delay(tau,order);

% [num,den]=delay(tau,order) determines a proper T.F.
% of specified order which approximates a given pure delay
% of tau seconds

% Luis A. Aguirre - Sheffield 19/06/91
% Taken from:
% Glader et al. IJC, Vol. 53, no. 2, 369-390, 1991.


for i=0:order

% denominator determination
  den(1+i)=fat(order+i)*tau^(order-i)/(fat(i)*fat(order-i));
% numerator determination
  num(1+i)=fat(order+i)*(-tau)^(order-i)/(fat(i)*fat(order-i));

end;
