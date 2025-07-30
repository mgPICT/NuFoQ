function [Iout,whatScale,Voutx,Vouty,Voutz]=HoG_blobs3D(I,options)
%
% This function HoG_blobs3D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to blobs, according
% to the method (and code) described by Frangi
%
% [J,Scale,Vx,Vy,Vz] = HoG_blobs3D(I, Options)
%
% inputs,
%   I : The input image volume (vessel volume)
%   Options : Struct with input options,
%       .HoGScaleRange : The range of sigmas used, default [1 8]
%       .HoGScaleRatio : Step size between sigmas, default 2
%       .HoGAlpha : HoG vesselness constant, treshold on Lambda2/Lambda3
%					   determines if its a line(vessel) or a plane like structure
%					   default .5;
%       .HoGBeta  : HoG vesselness constant, which determines the deviation
%					   from a blob like structure, default .5;
%       .HoGC     : HoG vesselness constant which gives
%					   the threshold between eigenvalues of noise and 
%					   vessel structure. A thumb rule is dividing the 
%					   the greyvalues of the vessels by 4 till 6, default 500;
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
%       .Output : select the type of output image ('sum', 'max', 'all')
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Vx,Vy,Vz: Matrices with the direction of the smallest eigenvector, pointing
%				in the direction of the line/vessel.
%
% Literature, 
%	Manniesing et al. "Multiscale Vessel Enhancing Diffusion in 
%		CT Angiography Noise Filtering"
%
% Example,
%   % compile needed mex file
%   mex eig3volume.c
%
%   load('ExampleVolumeStent');
%   
%   % HoG Filter the stent volume
%   options.BlackWhite=false;
%   options.HoGScaleRange=[1 1];
%   Vfiltered=HoGFilter3D(V,options);
%
%   % Show maximum intensity plots of input and result
%   figure, 
%   subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
%   subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
%   subplot(2,2,3), imshow(V(:,:,100),[])
%   subplot(2,2,4), imshow(Vfiltered(:,:,100),[])
%
% Written by D.Kroon University of Twente (May 2009)

% Constants vesselness function

defaultoptions = struct('HoGScaleRange', [1 10], 'HoGScaleRatio', 2, 'HoGAlpha', 0.5, 'HoGBeta', 0.5, 'HoGC', 500, 'verbose',true,'BlackWhite',true,'Output','max');

% Process inputs
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('HoGFilter3D:unknownoption','unknown options found');
    end
end

% Use single or double for calculations
if(~isa(I,'double')), I=single(I); end

sigmas=options.HoGScaleRange(1):options.HoGScaleRatio:options.HoGScaleRange(2);
sigmas = sort(sigmas, 'ascend');

% HoG filter for all sigmas
for i = 1:length(sigmas),
    % Show progress
    if(options.verbose)
        disp(['Current HoG Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigmas(i));

    if(sigmas(i)>0)
        % Correct for scaling
        c=(sigmas(i)^2);
        Dxx = c*Dxx; Dxy = c*Dxy;
        Dxz = c*Dxz; Dyy = c*Dyy;
        Dyz = c*Dyz; Dzz = c*Dzz;
    end
    
    % Calculate eigen values
    if(nargout>2)
        [Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    else
        [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    end
    
    % Free memory
    clear Dxx Dyy  Dzz Dxy  Dxz Dyz;

    % Calculate absolute values of eigen values
%     LambdaAbs1=abs(Lambda1);
%     LambdaAbs2=abs(Lambda2);
%     LambdaAbs3=abs(Lambda3);

    % The blobness Features
    if options.BlackWhite
        tra = Lambda1+Lambda2+Lambda3;
        det = Lambda1.*Lambda2.*Lambda3;
    else
        tra = imcomplement(Lambda1+Lambda2+Lambda3);
        det = imcomplement(Lambda1.*Lambda2.*Lambda3);
    end
    tra(tra<0) = 0;
    det(det<0) = 0;

    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D
%     S = sqrt(LambdaAbs1.^2+LambdaAbs2.^2+LambdaAbs3.^2);
	
    % Free memory
    clear LambdaAbs1 LambdaAbs2 LambdaAbs3

    %Compute blobness function
    tra = tra./3;
    det = det.^(1/3);
    
    A = range(tra(:)).*options.HoGAlpha;
    B = range(det(:)).*options.HoGBeta;
    
    expTra = (1-exp(-(tra./A)));
    expDet = (1-exp(-(det./B)));

%     keyboard
    % Free memory
%     clear S A B C Ra Rb
    clear tra det

    %Compute blobness function
    Voxel_data = expTra.* expDet;
    
    % Free memory
    clear expTra expDet;
    
%     if(options.BlackWhite)
%         Voxel_data(Lambda2 < 0)=0; Voxel_data(Lambda3 < 0)=0;
%     else
%         Voxel_data(Lambda2 > 0)=0; Voxel_data(Lambda3 > 0)=0;
%     end
        
    % Remove NaN values
    Voxel_data(~isfinite(Voxel_data))=0;

    % Add result of this scale to output
    if(i==1)
        Iout=Voxel_data;
        if(nargout>1)
            whatScale = ones(size(I),class(Iout));
        end
        if(nargout>2)
            Voutx=Vx; Vouty=Vy; Voutz=Vz;
        end
    else
        if(nargout>1)
            whatScale(Voxel_data>Iout)=i;
        end
        if(nargout>2)
            Voutx(Voxel_data>Iout)=Vx(Voxel_data>Iout);
            Vouty(Voxel_data>Iout)=Vy(Voxel_data>Iout);
            Voutz(Voxel_data>Iout)=Vz(Voxel_data>Iout);
        end
        % Keep maximum filter response
        switch options.Output
            case 'max'
                Iout=max(Iout,Voxel_data);
            case 'sum'
                Iout = Iout + Voxel_data;
            case 'all'
                Iout = cat(4, Iout, Voxel_data);
        end
    end
end

