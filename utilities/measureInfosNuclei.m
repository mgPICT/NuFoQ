function [nuclei] = measureInfosNuclei(dec, dec_mod, seg_dec, spheres, spots, saveImgs, seg_name)

    spheres = double(spheres);
    seg_dec = double(seg_dec);
    
	nbNuclei = max(seg_dec(:));
	nuclei = zeros(nbNuclei, 4);
    
    dims = size(dec);
    
    seg_nospots = seg_dec .* (spheres==0);

	fprintf('Measure information from nuclei\n');
	reps=15;
    fprintf(['\t[',repmat(' ',1,reps),']']) %make the "background"
    did = 0;
	for i_nuc = 1:nbNuclei
		% Get nuclear intensity
		nuclei(i_nuc, 1) = sum(dec(seg_dec == i_nuc));
        
		% Get nuclear background
		nuclei(i_nuc, 2) = sum(dec(seg_nospots == i_nuc));
        
		% Get number of spots
		nuclei(i_nuc, 3) = sum(spots(:,6) == i_nuc);
        
        % Get volume
        nuclei(i_nuc, 4) = sum(seg_dec(:) == i_nuc);
		
		if sum(i_nuc == round(1:(nbNuclei/reps):nbNuclei)) == 1
            fprintf(repmat('\b',1,reps+1-did)) %send the cursor back to the start
            did = did +1;
            fprintf('-');
            fprintf(repmat(' ',1,reps-did));
            fprintf(']');
        end
	end

end
