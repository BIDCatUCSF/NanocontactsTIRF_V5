function watershedMaskStack = watershedMask( maskStack )
%watershedMask Creates a stack of the mask after watershed transform
%  
% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% May 2017

    %%
    % Find the size of the mask
    dimensions = size(maskStack);
    
    % preallocate for speed
    dStack = zeros(dimensions);
    LdStack = zeros(dimensions);
    bw2Stack = zeros(dimensions);
    mStack = zeros(dimensions);
    d2Stack = zeros(dimensions);
    Ld2Stack = zeros(dimensions);
    watershedMaskStack = zeros(dimensions);
    % Perform the watershed as per the mathworks formula
    for ii = 1:dimensions(3)
        dStack(:,:,ii) = -bwdist(~maskStack(:,:,ii));
        LdStack(:,:,ii) = watershed(dStack(:,:,ii));
        bw2Stack(:,:,ii) = maskStack(:,:,ii);
        for jj = 1:dimensions(1)
            for kk = 1:dimensions(2)
                if LdStack(jj,kk,ii) == 0
                    bw2Stack(jj,kk,ii) = 0;
                end % if
            end % for
        end % for
        mStack(:,:,ii) = imextendedmin(dStack(:,:,ii),1);
        d2Stack(:,:,ii) = imimposemin(dStack(:,:,ii),mStack(:,:,ii));
        Ld2Stack(:,:,ii) = watershed(d2Stack(:,:,ii));
        watershedMaskStack(:,:,ii) = maskStack(:,:,ii);
        for jj = 1:dimensions(1)
            for kk = 1:dimensions(2)
                if Ld2Stack(jj,kk,ii) == 0
                    watershedMaskStack(jj,kk,ii) = 0;
                end % if
            end % for
        end % for
    end % for

% http://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
end % watershedMask

