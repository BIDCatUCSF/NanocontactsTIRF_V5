function outStack = makeBinary(inStack,threshold)
%makeBinary Creates a stack of ones and zeros from the input stack using
%the threshold.

% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% May 2017

    %% 
    % Find the size of the stack to be stabilized
    dimensions = size(inStack);
    c = 'uint8';

    % If greater or equal to threshold, 1, else zero.
    outStack = zeros(dimensions,c);
    for ii = 1:dimensions(3)
        for jj = 1:dimensions(1)
            for kk = 1:dimensions(2)
                if inStack(jj,kk,ii) >= threshold
                    outStack(jj,kk,ii) = 1;
                end
            end
        end
    end

end % makeBinary

