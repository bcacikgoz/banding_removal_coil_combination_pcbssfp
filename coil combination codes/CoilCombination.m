function combined_profiles = CoilCombination(uncombined_images, method,p)

[Nx,Ny,Npc,Nc] = size(uncombined_images);
%%%%% Pre-rotation of ellipses
uncombined_images = uncombined_images.*exp(-1i*angle(repmat(mean(uncombined_images,3), [1, 1, Npc, 1])));

switch method
    case "ModWalsh"
        combined_profiles = WalshbSSFP(uncombined_images);
    case "Walsh"
        [~, weights] = openadapt(permute(squeeze(squeeze(sum(uncombined_images, 3))), [3 1 2]), false, eye(size(uncombined_images,4)));
        weights = permute(weights, [2 3 1]);
        reweighted_uncombined_images = uncombined_images.*repmat(permute(weights, [1 2 4 3]), [1 1 Npc 1]);
        combined_profiles = sum(reweighted_uncombined_images, 4);
    case "p-Norm"
        weights = vecnorm(repmat(abs(sum(uncombined_images,3)),[1 1 Npc 1]), p, 4);
        reweighted_uncombined_images = uncombined_images.*weights;
        combined_profiles = sum(reweighted_uncombined_images, 4);
    case "PancakeStack"
        combined_profiles = sum(uncombined_images,4);
end

end




