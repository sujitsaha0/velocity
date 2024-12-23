clc
clear
%% parameters
EI=0.01;
alpha_E =-0.5;
gta =-1;
bta = 0.0;
alpha_e =0;
alpha_s =0;
fi=0.03;
Ha=0;
Du =0;
S =0;
Da =50;
k =20;
t=5;
efs_ratio = 80 / (6.95 * 10^(-10));
sigma_ratio = 25000 / 0.1;
row_ratio = 5.18 / 1.025;
P = 1;
%% Initialize arrays
z = linspace(0,1,100);
u = zeros(100, 1);
%% Main loop to compute u
for i = 1:100
    y = z(i);
    sum = 0;
    for n = 1:500
        equation = @(lamda) (cos(lamda) - bta * lamda * sin(lamda));
        lamda_guess = (2 * n - 1) * pi / 2;
        lamda = Newton_Raphson(equation, lamda_guess);
         gta_a = gta * (1 + bta * k * (sinh(gta) / gta));
         alpha_1 = alpha_e / ((1 + alpha_e * alpha_s)^2 + alpha_e^2);
         alpha_2 = (1 + alpha_e * alpha_s) / ((1 + alpha_e * alpha_s)^2 + alpha_e^2);
         efs_r = 1 + ((2 * (efs_ratio - 1) * fi) / (((efs_ratio + 2) - (efs_ratio - 1) * fi)));
         row_r = (1 - fi) + (row_ratio) * fi;
        
         mu_r = 1 / ((1 - fi)^2.5);
         sigma_r = 1 + ((3 * (sigma_ratio - 1) * fi) / ((sigma_ratio + 2) - (sigma_ratio - 1) * fi));
        
        
         c2 = sqrt((sigma_r / mu_r) * alpha_2 * Ha^2 + (1 / Da));
         c3 = row_r / mu_r;
        c4 = ((gta_a * k * sqrt(efs_r)) / (k^2 - c2^2 * efs_r)) * (cosh(c2) * tanh(k / sqrt(efs_r)));
        c41 = -(P * c4) / (mu_r * c2^2 * cosh(c2) + mu_r * bta * c2^3 * sinh(c2));
        c42 = -(c4 * alpha_1 * sigma_r * Ha^2 * S) / (mu_r * c2^2 * cosh(c2) + mu_r * bta * c2^3 * sinh(c2));
        c43 = -(c4 * k^2 * efs_r * gta_a) / ((k^2 - c2^2 * efs_r) * (cosh(c2) + bta * c2 * sinh(c2)));
        c44 = -(c4 * bta * k^3 * sqrt(efs_r) * gta_a * tanh(k / sqrt(efs_r))) / ((k^2 - c2^2 * efs_r) * (cosh(c2) + bta * c2 * sinh(c2)));
        c5 = (c2 * efs_r * gta_a * sinh(c2)) / ((k^2 - c2^2 * efs_r));
        c51 = (P * c5) / (mu_r * c2^2 * cosh(c2) + mu_r * bta * c2^3 * sinh(c2));
        c52 = (c5 * alpha_1 * sigma_r * Ha^2 * S) / (mu_r * c2^2 * cosh(c2) + mu_r * bta * c2^3 * sinh(c2));
        c53 = (c5 * k^2 * efs_r * gta_a) / ((k^2 - c2^2 * efs_r) * (cosh(c2) + bta * c2 * sinh(c2)));
        c54 = (c5 * bta * k^2 * sqrt(efs_r) * gta_a * tanh(k / sqrt(efs_r))) / ((k^2 - c2^2 * efs_r) * (cosh(c2) + bta * c2 * sinh(c2)));
        c6 = (gta_a * sqrt(efs_r) * tanh(k / sqrt(efs_r))) / (c2^2 * k);
        c61 = (P * c6) / mu_r;
        c62 = (c6 * alpha_1 * sigma_r * Ha^2 * S) / mu_r;
        c7 = ((k^2 * efs_r * gta_a^2) / (2 * mu_r * (k^2 - c2^2 * efs_r))) * (1 / (cosh(k / sqrt(efs_r))^2));
        c8 = (k * efs_r^1.5 * gta_a^2 * tanh(k / sqrt(efs_r))) / (2 * mu_r * (k^2 - c2^2 * efs_r));
        E_s = (c41 + c51 + c61) / (alpha_E + alpha_E * Du - c42 - c43 - c44 - c52 - c53 - c54 - c62 - c7 - c8);
        c1 = (P + alpha_1 * sigma_r * E_s * Ha^2 * S) / mu_r;
        S1 = (EI * c2^2 + c3) / (2 * c3 * EI);
        S2 = (c2^2 + lamda^2) / (2 * c3 * EI);
        mu_n = sqrt(2 * S2 - S1^2);
     
        D = (1 / (cosh(c2) + bta * c2 * sinh(c2))) * ((-c1 / c2^2) - ((bta * k^3 * sqrt(efs_r) * E_s * gta_a) / (mu_r * (k^2 - c2^2 * efs_r))) * tanh(k / sqrt(efs_r)) - (k^2 * efs_r * E_s * gta_a) / (mu_r * (k^2 - c2^2 * efs_r)));
        A1 = (D * lamda * cosh(c2) * sin(lamda)) / (lamda^2 + c2^2) + D * (c2 / (lamda^2 + c2^2)) * (sinh(c2) * cos(lamda));
        A2 = ((k^2 * efs_r * E_s * gta_a) / (mu_r * cosh(k / sqrt(efs_r)) * (k^2 - c2^2 * efs_r))) * (((lamda * efs_r * cosh(k / sqrt(efs_r)) * sin(lamda)) / (efs_r * lamda^2 + k^2)) + ((sqrt(efs_r) * k) / (efs_r * lamda^2 + k^2)) * (sinh(k / sqrt(efs_r)) * cos(lamda)));
        E_n = (-1) * ((4 * lamda) / (2 * lamda + sin(2 * lamda))) * (A1 + A2 + (c1 / (c2^2 * lamda)) * sin(lamda));
        S_N = E_n * cos(lamda * y) * exp(-S1 * t) * (cos(mu_n * t) + (S1 / mu_n) * sin(mu_n * t));
        sum = sum + S_N;
    end
    v1 = D * cosh(c2 * y) + (c1 / c2^2) + ((k^2 * efs_r * E_s * gta_a) / (mu_r * (k^2 - c2^2 * efs_r))) * (cosh((k / sqrt(efs_r)) * y) / cosh(k / sqrt(efs_r)));
    R1 = v1 + sum;
    u(i, 1) = R1;
end
hold on, plot(z,u)