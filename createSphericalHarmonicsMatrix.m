%% Create spherical harmonics matrix
%
% Note that we use real-valued N3D-normalized spherical harmonics.
%
%   Y =
%                        n = 0                      n = 1                                    n = N
%                    /------------\  /---------------------------------------------\  ... --------------\
%
%       gamma_1   /  Y_0^0(gamma_1)  Y_1^-1(gamma_1)  Y_1^0(gamma_1)  Y_1^1(gamma_1)  ...  Y_N^N(gamma_1)  \
%       gamma_2   |  Y_0^0(gamma_2)  Y_1^-1(gamma_2)  Y_1^0(gamma_1)  Y_1^1(gamma_1)  ...  Y_N^N(gamma_2)  |
%       gamma_3   |  Y_0^0(gamma_3)  Y_1^-1(gamma_3)  Y_1^0(gamma_1)  Y_1^1(gamma_1)  ...  Y_N^N(gamma_3)  |
%         ...     |     ...             ...              ...             ...          ...     ...          |
%       gamma_Q   \  Y_0^0(gamma_Q)  Y_1^-1(gamma_Q)  Y_1^0(gamma_1)  Y_1^1(gamma_1)  ...  Y_N^N(gamma_Q)  /
%
%
%
% Directions: vector theta and vector phi of length Q
% Spherical harmonics order: N
function Y = createSphericalHarmonicsMatrix(thetaVec, phiVec, N)
    
    Q = numel(thetaVec);
    Y = zeros(Q,(N+1)^2);
    
    for n=0:N
        %% normalization constant
        F = sqrt( (2*n+1)/(4*pi) );
        
        %% Legendre part
        P_pos = bsxfun(@times, legendre(n,cos(thetaVec(:))), ...
            sqrt(factorial(n-(0:n)')./factorial(n+(0:n)')) );
        
        % MATLAB uses the Condon-Shortley-Phase in legendre!
        % Cancel the Condon-Shortley-Phase which is a term (-1)^m
        % More details: http://mathworld.wolfram.com/Condon-ShortleyPhase.html
        P_pos = P_pos .* ( (-1).^(0:n)'*ones(1,Q) );
        % negative m
        P_neg = P_pos(end:-1:2,:) ;
        % all together
        P = [P_neg; P_pos];
        
        %% trigonometric part
        mx = (-n:n)';
        R =   (mx>0)  * ones(1,Q) .* cos(abs(mx)*phiVec(:)') * sqrt(2) ...
            + (mx<0)  * ones(1,Q) .* sin(abs(mx)*phiVec(:)') * sqrt(2) ...
            + (mx==0) * ones(1,Q) ;
        
        %% total spherical harmonics
        acnx = n.^2+n+(-n:n);
        Y(:,acnx+1) = F .* P' .* R';
    end
end
