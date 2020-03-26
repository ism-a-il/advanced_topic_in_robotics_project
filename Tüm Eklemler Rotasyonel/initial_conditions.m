clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial_conditions deðiþkenleri
global n x y z a R
global Phi H Phit J l l_s hh h_s
global h1 h2 h3 h4 h5 h6 h7
global l12 l23 l34 l45 l56 l67 l7t lsum
global l_s_7t
global H1 H2 H3 H4 H5 H6 H7
global preH
global sizelnn sizelnn_s sizehn_s sizePhinn
global templ temph
global sampling_time
global theta_dots
global counter_forward
global counter_inverse
global output_matrix
global output_matrix_forward i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deðiþkenler tanýmlanýyor.
n = 7;
a = 2;
R = zeros(3,3*n);
counter_forward = 0;
counter_inverse = 0;
sampling_time = 0.01;
H = zeros(6*n,n);
Phi = zeros(6*n,6*n);
l = zeros(3*n,n);
sizelnn = [3,1];
l_s = zeros(3*n,3*n);
sizelnn_s = [3 3];
hh = zeros(3,n);
h_s = zeros(3,3*n);
sizehn_s = [3 3];
sizePhinn = [6 6];
Phit = zeros(6,6*n);
theta_dots = [0;0;0;0;0;0;0];
output_matrix = zeros(31,7);
output_matrix_forward = zeros(6,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Koordinat eksen matrisleri oluþturuluyor.
x = zeros(3,n);
x(1,:) = 1;
y = zeros(3,n);
y(2,:) = 1;
z = zeros(3,n);
z(3,:) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h matrisi oluþturuluyor.
h1 = y(:,1);
h2 = z(:,2);
h3 = x(:,3);
h4 = z(:,4);
h5 = x(:,5);
h6 = x(:,6);
h7 = y(:,7);
hh = [h1 h2 h3 h4 h5 h6 h7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l matrisi oluþturuluyor.
l12 = a*y(:,1);   l(1:sizelnn(1),2) = l12;
l23 = a*y(:,2);   l(1+sizelnn(1):2*sizelnn(1),3) = l23;
l34 = a*y(:,3);   l(1+2*sizelnn(1):3*sizelnn(1),4) = l34;
l45 = a*y(:,4);   l(1+3*sizelnn(1):4*sizelnn(1),5) = l45;
l56 = a*y(:,5);   l(1+4*sizelnn(1):5*sizelnn(1),6) = l56;
l67 = a*y(:,6);   l(1+5*sizelnn(1):6*sizelnn(1),7) = l67;
l7t = a*y(:,7);
lsum = l12 + l23 + l34 + l45 + l56 + l67 + l7t;

for i = 1:n
    for j = 1:n
        
       if(i >= j || j == i+1)
           continue
       end
       
       l(1+(i-1)*sizelnn(1):i*sizelnn(1),j) = l(1+(i-1)*sizelnn(1):i*sizelnn(1),j-1) + ...
            l(1+(j-2)*sizelnn(1):(j-1)*sizelnn(1),j);
       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l_s matrisi oluþturuluyor.
for i = 1:n
    for j = 1:n
        
       if(i >= j)
           continue
       end
       
       templ = l(1+(i-1)*sizelnn(1):i*sizelnn(1),j);
       l_s(1+(i-1)*sizelnn_s(1):i*sizelnn_s(1),1+(j-1)*sizelnn_s(2):j*sizelnn_s(2)) = [0 -templ(3) templ(2);
                                                                                       templ(3) 0 -templ(1);
                                                                                      -templ(2) templ(1) 0];       
    end
end

l_s_7t = [0 -l7t(3) l7t(2); l7t(3) 0 -l7t(1); -l7t(2) l7t(1) 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phi matrisi oluþturuluyor.
for i = 1:n
    for j = 1:n
        
        if(i < j)
            continue
        elseif(i == j)
            Phi(1+(i-1)*sizePhinn(1):i*sizePhinn(1),1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = eye(6);
        else
            Phi(1+(i-1)*sizePhinn(1):i*sizePhinn(1),1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = [eye(3) zeros(3);
                -l_s(1+(j-1)*sizelnn_s(1):j*sizelnn_s(1),1+(i-1)*sizelnn_s(2):i*sizelnn_s(2)) eye(3)];
        end   
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_s matrisi oluþturuluyor.
for j = 1:n
    
    temph = hh(:,j);
    h_s(:,1+(j-1)*sizehn_s(2):j*sizehn_s(2)) = [0 -temph(3) temph(2);
                                                temph(3) 0 -temph(1);
                                               -temph(2) temph(1) 0];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preH matrisi oluþturuluyor.
H1 = [h1; zeros(3,1)];
H2 = [h2; zeros(3,1)];
H3 = [h3; zeros(3,1)];
H4 = [h4; zeros(3,1)];
H5 = [h5; zeros(3,1)];
H6 = [h6; zeros(3,1)];
H7 = [h7; zeros(3,1)];
preH = [H1 H2 H3 H4 H5 H6 H7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H matrisi oluþturuluyor.
for i = 0:(n-1)
    
    H((n-1)*i+1:(n-1)*(i+1), i+1) = preH(:,i+1);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phit matrisi oluþturuluyor.
for j = 1:n
        
        if(j <= n-1)
            Phit(:,1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = zeros(6);
        else
            Phit(:,1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = [eye(3) zeros(3); -l_s_7t eye(3)];
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J matrisi oluþturuluyor.
J = Phit*Phi*H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
