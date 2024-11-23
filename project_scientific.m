clear;
close all;

%Initialization of parameters
dx = 1/40; %change for different spacial steps
dt = 10*(dx)^2;
mu = 10;
x0 = -1;
xJ = 1;
t0 = 0;
tN = 0.5;
N = round((tN - t0) / dt) + 1;
J = round((xJ - x0) / dx) + 1;

%Stopping criterias
maxIterNum = 10^5;
tol = min([dt^2, dx^2]);

%Initialization of IC 
x = linspace(x0, xJ, J); % Use linspace to generate evenly spaced points
tt = linspace(t0, tN, N); % Similarly for the time grid
u = zeros(J, N);

for i = 1 : J
    u(i,1) = initCon(x(i));
end

% LU-decomposition
subD = -mu * ones(J-3,1);
leadD = (2 + 2*mu) * ones(J-2,1);
supD = -mu * ones(J-3,1);

C_mat = [mu, 2 - 2 * mu, mu];
D_Vec = zeros(J-2,1);

%Boundary conditions:
u(1,:) = exactSol(x0, tt, tol, maxIterNum);
u(J,:) = exactSol(xJ, tt, tol, maxIterNum);


for n = 2 : N
    u_Old = u(2:J-1, n-1);
    t = (n-1)*dt;
    bc_left = mu*(u(1,n-1) + u(1,n));
    bc_right = mu*(u(J,n-1) + u(J,n));
    
    D_Vec(1) = bc_left + C_mat(2) * u_Old(1) + C_mat(3) * u_Old(2);
    for j = 2 : J-3
        D_Vec(j) = C_mat(1) * u_Old(j-1) + C_mat(2) * u_Old(j) + ...
            C_mat(3) * u_Old(j+1);
    end
    D_Vec(J-2) =  C_mat(1) * u_Old(J-3) + C_mat(2) * u_Old(J-2) + bc_right;

    u(2:J-1, n) = tridiagonal_vector(subD, leadD, supD, D_Vec);

end


fig1 = figure(1);
colormap("winter")
surf(tt,x,u,'EdgeColor','none');
xlabel('t')
ylabel('x')
zlabel('u approx')
colorbar
fontsize(fig1, 15, "points")

uExact = zeros(J,N);
for n = 1 : N
    t = (n-1)*dt;
    uExact(:,n) = exactSol(x, t, tol, maxIterNum);
end

fig2 = figure(2);
colormap("spring")
surf(tt,x,uExact,'EdgeColor','none');
xlabel('t')
ylabel('x')
zlabel('u exact')
colorbar
fontsize(fig2, 15, "points")


L2_error = sqrt(sum( (uExact - u).^2, "all")/(N*J));
L_inf_error = max(abs(uExact - u),[],"all");

% Display errors
fprintf('L2 Error: %.6f\n', L2_error);
fprintf('L_inf Error: %.6f\n', L_inf_error);

function u = exactSol(x, t, tol, maxIterNum)
    uval = 1/2;
    temp = 10;
    n = 0;
    while (any(abs(temp) > tol) && n < maxIterNum)
        temp = 2*(-1)^n * cos(pi*(2*n+1)*x)/(pi*(2*n+1))*exp(-pi^2*(2*n+1)^2*t);
        uval = uval + temp;
        n = n + 1;
    end
    u=uval;
end

function u = initCon(x)
    eps = 1.e-15;
    if (abs(x) < 1/2)
        u = 1;
    elseif (abs(x-1/2)<eps)
        u = 1/2;
    else 
        u = 0;
    end
end