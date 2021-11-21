function x=fat(i)

% fat(i) computes de factorial of i

if i == 0
  x=1;
elseif i == 1
  x=1;
else
  x=i;
  for j=i-1:-1:1
    x=x*j;
  end;
end;
