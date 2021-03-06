function u = invariant_probability_distribution_silagadze(Y,m,evenly_extend)
% u = invariant_probability_distribution_silagadze(Y,m)
% Inputs
%   Y: Vector to Smooth (Y must be positif and without NaN)
%   m: Span
% 
% Output
%   u: Probability distribution of the peaks
% By Y.Poltera 2015
% From ZK Silagadze - ?1995 http://arxiv.org/pdf/hep-ex/9506013.pdf

if nargin==2
    evenly_extend = false;
end

if evenly_extend
    Y =  [Y(m+1:-1:2);Y;Y(end-1:-1:end-m)];% evenly extend Y
end

P = zeros(length(Y),length(Y));
A = zeros(length(Y),1);
for i=2:length(Y)-1
    m_sup = min(length(Y)-i,m);
    m_sub = min(i-1,m);
    e_sup = sum(exp((Y(i+(1:m_sup))-Y(i))./sqrt(Y(i+(1:m_sup))+Y(i))));
    e_sub = sum(exp((Y(i-(1:m_sub))-Y(i))./sqrt(Y(i-(1:m_sub))+Y(i))));
    
    A(i) = 1/(e_sup+e_sub);
    
    P(i,i+1) = A(i) * e_sup;
    P(i,i-1) = A(i) * e_sub;
end

P(1,2) = 1;
P(end,end-1) = 1;

p_upper_diag = diag(P,1);
p_lower_diag = diag(P,-1);
factorU = zeros(length(Y),1);
for i=2:length(Y)
    factorU(i) = prod(p_upper_diag(1:i-1)./p_lower_diag(1:i-1));
end
u = zeros(length(Y),1);
u(1) = 1/(1+sum(factorU));
for i=2:length(Y)
    u(i) = factorU(i)*u(1);
end

if evenly_extend
    u = u(m+1:end-m);
else
    u(1:m) = NaN;
    u(end-m+1:end) = NaN;
end

end