function [a,ad,q,p,nn] = operators(qc)

I = ([1 0; 0 1]);

H = ((1/sqrt(2))*[1 1;1 -1]);

Z = ([1 0; 0 -1]);
X = ([0 1; 1 0]);
iY = (1i*[0 -1i; 1i 0]);  

minus = (1/2)*(X-iY);
plus = (1/2)*(X+iY);

n = (1/2)*(I-Z);

dim = 2^qc;

temp = (zeros(dim));

a = temp;
ad = temp;

for j = 0:1:qc-1
    if j ~= 0
        eyeFill = (eye(2^(j)));
    else
        eyeFill = 1;
    end
    
    dop = kron(eyeFill,minus);
    op = kron(eyeFill,plus);
    
    for k = j+1:1:qc-1
        dop = kron(dop,plus);
        op = kron(op,minus);
    end
    
    nopSum = temp;
    for k = 0:1:j-1
        if k ~= 0
            eyeFill = (eye(2^(k)));
        else
            eyeFill = (1);
        end
        nop = kron(eyeFill,n);
        if k ~= qc-1
            eyeFill = eye(2^(qc-k-1));
        else
            eyeFill = (1);
        end
        nop = kron(nop,eyeFill);
        nop = nop/2^k;
        
        nopSum = nopSum+nop;
    end
    
    nopSum = nopSum+(eye(dim)/2^j);
    nopSum = sqrtm(nopSum);
    
    a = a + nopSum*op;
    ad = ad + nopSum*dop;
end

weight = (2^(0.5*(qc-1)));
a = a*weight;
ad = ad*weight;

%% Position-momentum
q = (1/sqrt(2))*(ad+a);
p = (1i/sqrt(2))*(ad-a);

nn = zeros(2^qc,2^qc);
nn(2^qc,2^qc) = 1;

end
