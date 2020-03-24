n  = 40;
c  = 4;
b  = 8;

P2 = speye(n);

ord1 = [];
ord2 = [];
for i = 1:n/b
    ord1 = [ord1, (i-1)*b+1:(i-1)*b+c];
    ord2 = [ord2, (i-1)*b+c+1:i*b];
end

P2 = P2([ord1,ord2],:);

P2 * (1:n)'