function symtest(V,pm)
  [n,m]=size(V);
  m
  nx=sqrt(n/2);
  for offset=[0:1]
    offset
    for i=1:m
      ind=(1:2:n) +offset;
      norm(V(ind,i)+pm*reshape(fliplr(reshape(V(ind,i),nx,nx)),nx^2,1))
    end
  end
end
