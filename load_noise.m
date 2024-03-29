if ~isfile('noise.mat')
    dim.nnoise = 5; % no of realizations
    dim.tunp = 2; % no of parameters subjected to noise 
    Nobs = 50;      % No. of more observations to perform
    %% zero noise

    zero_noise_variance = 0;
    zero_noise_a = zeros(dim.tunp,Nobs,dim.nnoise);
    zero_noise_b = zeros(dim.tunp,Nobs,dim.nnoise);
    
    for i = 1:dim.nnoise
        zero_noise_a(:, :, i) = zero_noise_variance*randn(dim.tunp, Nobs);
    end

    zero_noise_c = cat(1, zero_noise_a, zero_noise_b);
    zero_noise = 1 + zero_noise_c;
    %% Low noise
    
    low_noise_variance = 0.01;
    low_noise_a = zeros(dim.tunp,Nobs,dim.nnoise);
    low_noise_b = zeros(dim.tunp,Nobs,dim.nnoise);
    
    for i = 1:dim.nnoise
        low_noise_a(:, :, i) = low_noise_variance*randn(dim.tunp, Nobs);
    end

    low_noise_c = cat(1, low_noise_a, low_noise_b);
    low_noise = 1 + low_noise_c;
    
    %% medium noise
    
    medium_noise_variance =  0.05;
    medium_noise_a = zeros(dim.tunp,Nobs,dim.nnoise);
    medium_noise_b = zeros(dim.tunp,Nobs,dim.nnoise);
    
    for i = 1:dim.nnoise
        medium_noise_a(:, :, i) = medium_noise_variance*randn(dim.tunp, Nobs);
    end

    medium_noise_c = cat(1, medium_noise_a, medium_noise_b);
    medium_noise = 1 + medium_noise_c;
    
    %% high noise
    
    high_noise_variance = 0.15;
    high_noise_a = zeros(dim.tunp,Nobs,dim.nnoise);
    high_noise_b = zeros(dim.tunp,Nobs,dim.nnoise);
    
    for i = 1:dim.nnoise
        high_noise_a(:, :, i) = high_noise_variance*randn(dim.tunp, Nobs);
    end

    high_noise_c = cat(1, high_noise_a, high_noise_b);
    high_noise =  1 + high_noise_c;

    save noise.mat zero_noise_variance zero_noise low_noise_variance low_noise medium_noise_variance medium_noise high_noise_variance high_noise
else
    load noise.mat
end