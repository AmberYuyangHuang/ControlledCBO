function y = polycomb(varargin)

% Simplified Version of Combvec (Mark Beale, 12-15-93)

y = varargin{1};
for i=2:length(varargin)
    z = varargin{i};
    y = [copy_blocked(y,size(z,2)); copy_interleaved(z,size(y,2))];
end

%=========================================================
function b = copy_blocked(m,n)

[mr,mc] = size(m);
b = zeros(mr,mc*n);
ind = 1:mc;
for i=[0:(n-1)]*mc
  b(:,ind+i) = m;
end
%=========================================================

function b = copy_interleaved(m,n)

[mr,mc] = size(m);
b = zeros(mr*n,mc);
ind = 1:mr;
for i=[0:(n-1)]*mr
  b(ind+i,:) = m;
end
b = reshape(b,mr,n*mc);