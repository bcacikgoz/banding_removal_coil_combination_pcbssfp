function combined_image = BandingRemoval(profiles, method, p)

if nargin == 2
    p = 2;
end

if method == "GeometricSolutionStable"
    combined_image = LinearizedStableGeometricSolution_inscript(profiles);
elseif method == "DensityWeightedSum"
    combined_image = DensityWeightedSum_inscript(profiles, p);
else

for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        prof = squeeze(profiles(x,y,:));
        if method == "ComplexSum"
            combined_image(x,y) = mean(squeeze(profiles(x,y,:)));

        elseif method == "GeometricSolution"
            combined_image(x,y) = GeometricSolution_inscript(squeeze(profiles(x,y,:)));

        elseif method == "MaxIntensity"
            combined_image(x,y) = max(abs(profiles(x,y,:)));

        elseif method == "p-Norm"
            combined_image(x,y) = norm(squeeze(profiles(x,y,:)), p);

        end

    end
end
end

end

function GSw = LinearizedStableGeometricSolution_inscript(profiles, patch_size)

if nargin == 1
    patch_size = 5;
end
for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        p = squeeze(profiles(x,y,:));
        GS = GeometricSolution(p);
        if abs(GS)>max(abs(p))
            GS = mean(p);
        end
        GSmap(x,y) = GS;
    end
end


stride = (patch_size-1)/2;
[Nx,Ny,Npc] = size(profiles);
GSw = zeros(Nx,Ny);
for x = (1+stride):(Nx-stride)
    for y = (1+stride):(Ny-stride)
        sum_num = 0;
        sum_denom = 0; 
        for xx = (x-stride):(x+stride)
            for yy = (y-stride):(y+stride)
                pp = squeeze(profiles(xx,yy,:));
                pp_gs = GSmap(xx,yy);
                sum_num = sum_num + 2*real(conj(pp(3)-pp_gs)*(pp(3)-pp(1)));
                sum_denom = sum_denom + 2*((abs(pp(1)-pp(3)))^2);
            end
        end
        w = sum_num/sum_denom;
        p = squeeze(profiles(x,y,:));
        GSw(x,y) = w*p(1) + (1-w)*p(3);
    end
end

for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        p = squeeze(profiles(x,y,:));
        GS = GSw(x,y);
        if abs(GS)>max(abs(p))
            GS = mean(p);
        end
        GSw(x,y) = GS;
    end
end
end
                    

function M = GeometricSolution_inscript(profile)
if length(profile)==4
I1 = profile(1); I2 = profile(2); 
I3 = profile(3); I4 = profile(4);
x1 = real(I1); y1 = imag(I1);
x2 = real(I2); y2 = imag(I2);
x3 = real(I3); y3 = imag(I3);
x4 = real(I4); y4 = imag(I4);
M = ((x1*y3-x3*y1)*(I2-I4) - (x2*y4-x4*y2)*(I1-I3))/...
    ((x1-x3)*(y2-y4) + (x2-x4)*(y3-y1));
if var(real(profile))<max(real(profile))*2e-4
    M = mean(profile);
end
else
    Npc = length(profile);
    angle_step = 360/Npc;
    sample_step = 90/angle_step;
    for i = 1:(Npc/4)
        I1 = profile(i); I2 = profile(i+sample_step); 
        I3 = profile(i+2*sample_step); I4 = profile(i+3*sample_step);   
        x1 = real(I1); y1 = imag(I1);
        x2 = real(I2); y2 = imag(I2);
        x3 = real(I3); y3 = imag(I3);
        x4 = real(I4); y4 = imag(I4);  
        M(i) = ((x1*y3-x3*y1)*(I2-I4) - (x2*y4-x4*y2)*(I1-I3))/...
        ((x1-x3)*(y2-y4) + (x2-x4)*(y3-y1));
    end
    M = mean(M);
end
end

function dw = DensityWeightedSum_inscript(profiles, c)

if nargin==1
    c = 15;
end

for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        p = squeeze(profiles(x,y,:));
        D = abs(repmat(p, [1, length(p)]) - repmat(p.', [length(p),1]));
        meandist = mean(mean(D));
%         for pcs = 1:length(p)
%             s = normpdf(D(pcs,:), 0, meandist);
%             s = s/max(s);
%             density_arr(pcs) = sum(s);
%         end
        w = 1./sum(D,1);
        w = w.^c;
        w = w/sum(w);
        dw(x,y) = w*p;
    end
end

end
