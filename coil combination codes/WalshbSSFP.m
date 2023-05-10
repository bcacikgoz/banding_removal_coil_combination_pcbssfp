function [combined, m] = WalshbSSFP(profiles_selected_slice)

% input dimensions are Nx x Ny x Npc x Ncoil

ncoils = size(profiles_selected_slice, 4);


mag = sqrt(sum(profiles_selected_slice .* conj(profiles_selected_slice),4));
s_raw=profiles_selected_slice ./ repmat(mag + eps,[1 1 1 ncoils]); clear mag; 
s_raw = profiles_selected_slice;
R = ismrm_correlation_matrix(s_raw);
m = ismrm_eig_power(R);

[Nx, Ny, Npc, Ncoil] = size(profiles_selected_slice);
combined = size(Nx, Ny, Npc);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Npc
        p = squeeze(profiles_selected_slice(i,j,k,:));
        mp = squeeze(m(i, j, :));
%         mp = (abs(mp)/sum(abs(mp))).*angle(mp);
        mp = mp/sum(abs(mp));
        combined(i,j,k) = mp.'*p;
        end
    end
end


end

function [Rs]=ismrm_correlation_matrix(s)
% function [Rs]=correlation_matrix(s);
%
% function correlation_matrix calculates the sample correlation matrix (Rs) for
% each pixel of a multi-coil image s(y,x,coil)
%
% input:
%    s   complex multi-coil image s(y,x,coil)
% output:
%    Rs  complex sample correlation matrices, Rs(y,x,coil,coil)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

[rows,cols,~,ncoils]=size(s);
disp(ncoils)
Rs=zeros(rows,cols,ncoils); % initialize sample correlation matrix to zero
for x = 1:rows
    for y = 1:cols
        for i=1:ncoils
            for j=1:i-1
		        Rs(x,y,i,j)=squeeze(s(x,y,:,i))'*squeeze(s(x,y,:,j));
		        Rs(x,y,j,i)=conj(Rs(x,y,i,j)); % using conjugate symmetry of Rs
            end
            Rs(x,y,i,i)=squeeze(s(x,y,:,i))'*squeeze(s(x,y,:,i));
        end
    end
    if mod(x, round(rows/100))==0
        disp(strcat("Estimating Sense Maps: %", num2str(round(x*100/rows))));
    end
end

end



function [v,d]=ismrm_eig_power(R)
% function [v,d]=eig_power(R);
%
% vectorized method for calculating the dominant eigenvector based on
% power method. Input, R, is an image of sample correlation matrices
% where: R(y,x,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%
% v is the dominant eigenvector
% d is the dominant (maximum) eigenvalue

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

rows=size(R,1);cols=size(R,2);ncoils=size(R,3);
N_iterations=2;
v=ones(rows,cols,ncoils); % initialize e.v.

d=zeros(rows,cols);
for i=1:N_iterations
    v=squeeze(sum(R.*repmat(v,[1 1 1 ncoils]),3));
	d=ismrm_rss(v);
    d( d <= eps) = eps;
	v=v./repmat(d,[1 1 ncoils]);
end

p1=angle(conj(v(:,:,1)));
% (optionally) normalize output to coil 1 phase
v=v.*repmat(exp(sqrt(-1)*p1),[1 1 ncoils]);
v=conj(v);

end



function y = ismrm_rss(x,dim)
%
%   [mag] = ismrm_rss(samples, dim)
%
%   Computes root-sum-of-squares along a single dimension.
%
%
%   INPUT:
%     - x   : multi-dimensional array of samples
%     - dim : dimension of reduction; defaults to last dimension
%
%   OUTPUT:
%
%     - y       : root sum of squares result
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%


if nargin==1
    dim=ndims(x);
else
    if isempty(dim); dim=ndims(x); end
end

y = squeeze(sqrt(sum(real(x).^2 + imag(x).^2,dim)));
end