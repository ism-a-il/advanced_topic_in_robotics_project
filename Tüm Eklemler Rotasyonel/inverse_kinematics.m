function output_matrix = inverse_kinematics( Vt )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial_conditions de�i�kenleri
global n x y z a R
global Phi H Phit J l_s hh h_s
global h1 h2 h3 h4 h5 h6 h7
global l12 l23 l34 l45 l56 l67 l7t
global l_s_7t
global H1 H2 H3 H4 H5 H6 H7
global preH
global sizelnn sizelnn_s sizehn_s sizePhinn
global templ temph
global sampling_time
global theta_dots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward_kinematics de�i�kenleri
global sizeRn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse_kinematics de�i�kenleri
global invJlineer Jlineer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% De�i�kenler tan�mlan�yor.
sizeRn = 3;
l = zeros(3*n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R matrisi olu�turuluyor.
for j = 1:n
    
    R(:,1+(j-1)*sizeRn:j*sizeRn) = eye(3) + ...
        h_s(:,1+(j-1)*sizehn_s(2):j*sizehn_s(2))*sin(theta_dots(j)*sampling_time) + ...
        (h_s(:,1+(j-1)*sizehn_s(2):j*sizehn_s(2))^2)*(1-cos(theta_dots(j)*sampling_time));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x(:,7) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*R(:,19:21)*x(:,7);
y(:,7) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*R(:,19:21)*y(:,7);
z(:,7) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*R(:,19:21)*z(:,7);

x(:,6) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*x(:,6);
y(:,6) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*y(:,6);
z(:,6) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*R(:,16:18)*z(:,6);

x(:,5) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*x(:,5);
y(:,5) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*y(:,5);
z(:,5) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*R(:,13:15)*z(:,5);

x(:,4) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*x(:,4);
y(:,4) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*y(:,4);
z(:,4) = R(:,1:3)*R(:,4:6)*R(:,7:9)*R(:,10:12)*z(:,4);

x(:,3) = R(:,1:3)*R(:,4:6)*R(:,7:9)*x(:,3);
y(:,3) = R(:,1:3)*R(:,4:6)*R(:,7:9)*y(:,3);
z(:,3) = R(:,1:3)*R(:,4:6)*R(:,7:9)*z(:,3);

x(:,2) = R(:,1:3)*R(:,4:6)*x(:,2);
y(:,2) = R(:,1:3)*R(:,4:6)*y(:,2);
z(:,2) = R(:,1:3)*R(:,4:6)*z(:,2);

x(:,1) = R(:,1:3)*x(:,1);
y(:,1) = R(:,1:3)*y(:,1);
z(:,1) = R(:,1:3)*z(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h matrisi g�ncelleniyor.
h1 = y(:,1);
h2 = z(:,2);
h3 = x(:,3);
h4 = z(:,4);
h5 = x(:,5);
h6 = x(:,6);
h7 = y(:,7);
hh = [h1 h2 h3 h4 h5 h6 h7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l matrisi g�ncelleniyor.
l12 = a*y(:,1);   l(1:sizelnn(1),2) = l12;
l23 = a*y(:,2);   l(1+sizelnn(1):2*sizelnn(1),3) = l23;
l34 = a*y(:,3);   l(1+2*sizelnn(1):3*sizelnn(1),4) = l34;
l45 = a*y(:,4);   l(1+3*sizelnn(1):4*sizelnn(1),5) = l45;
l56 = a*y(:,5);   l(1+4*sizelnn(1):5*sizelnn(1),6) = l56;
l67 = a*y(:,6);   l(1+5*sizelnn(1):6*sizelnn(1),7) = l67;
l7t = a*y(:,7);

for i = 1:n
    for j = 1:n
        
       if(i >= j || j == i+1)
           continue
       end
       
       l(1+(i-1)*sizelnn(1):i*sizelnn(1),j) = l(1+(i-1)*sizelnn(1):i*sizelnn(1),j-1) + ...
            l(1+(j-2)*sizelnn(1):(j-1)*sizelnn(1),j);
       
    end
end
output_matrix(1:21,1:7) = l;
output_matrix(22:24,1) = l7t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l_s matrisi g�ncelleniyor.
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
% Phi matrisi g�ncelleniyor.
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
% h_s matrisi g�ncelleniyor.
for j = 1:n
    
    temph = hh(:,j);
    h_s(:,1+(j-1)*sizehn_s(2):j*sizehn_s(2)) = [0 -temph(3) temph(2);
                                                temph(3) 0 -temph(1);
                                               -temph(2) temph(1) 0];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preH matrisi g�ncelleniyor.
H1 = [h1; zeros(3,1)];
H2 = [h2; zeros(3,1)];
H3 = [h3; zeros(3,1)];
H4 = [h4; zeros(3,1)];
H5 = [h5; zeros(3,1)];
H6 = [h6; zeros(3,1)];
H7 = [h7; zeros(3,1)];
preH = [H1 H2 H3 H4 H5 H6 H7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H matrisi g�ncelleniyor.
for i = 0:(n-1)
    
    H((n-1)*i+1:(n-1)*(i+1), i+1) = preH(:,i+1);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phit matrisi g�ncelleniyor.
for j = 1:n
        
        if(j <= n-1)
            Phit(:,1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = zeros(6);
        else
            Phit(:,1+(j-1)*sizePhinn(2):j*sizePhinn(2)) = [eye(3) zeros(3); -l_s_7t eye(3)];
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J matrisi g�ncelleniyor.
J = Phit*Phi*H;
Jlineer = J(4:6,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% invJ matrisi olu�turuluyor.
invJlineer = pinv(Jlineer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta_dots hesaplan�yor.
theta_dots = invJlineer*Vt;
output_matrix(25:31,1) = theta_dots;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

