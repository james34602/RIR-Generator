c = 340;                    % Sound velocity (m/s)
fs = 48000;                 % Sample frequency (samples/s)
r = [2 1.5 2 ; 1 1.5 2 ; 2.5 3.5 1.2 ; 2.2 4.4 1.2];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [12 12 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)  //Error test  [36 46 66]; 
beta = 0.4;                 % Reverberation time (s)
n = 8192;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);
h=h';
csvwrite('impExamp5.csv',h);