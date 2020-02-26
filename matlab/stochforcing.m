function ret=stochforcing(nrows, massarray) 

l=0.125;
C=1 ;
nx=sqrt(nrows/2);

for col=0:1
  etax = C*col;
  etay = C*(1-col);
  for i = 1:2:nrows
    j = i / 2;
    x = mod(j,nx) / (nx - 1);
    y = floor(j / nx) /(nx - 1);
    ret(i,1+col)= exp(-2. * ((x-0.5) * (x-0.5) + (y-0.5) * (y-0.5)) / (4. * l * l)) * ...
       (etay - 2. * etay * x + 2. * etax * y - etax) / (2. * l * l);
    if(massarray(i) == 0.0) ret(i,1+col)=0.0D0; end;
    ret(i+1,1+col)=0.0;
  end
end 

end 
