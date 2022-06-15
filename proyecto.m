%rutina para elementos finitos 2D
%Condicion inicial
Ini=randn(4,5);
[a b]=size(Ini);
%matriz que tiene |# nodo|coordenadas del nodo (x,y) cada nodo habla de una
%funcion base y viceversa
Nodos=zeros(a*b,3)
Nodos(:,1)=[1:a*b]' 
for i=0:b-1
Nodos(a*i+1 : a*(i+1), 2)=(i+1)*ones(a,1)
Nodos(a*i+1 : a*(i+1), 3)=[1:a]' 
end
%matriz que tiene |# elemento|nodos del elemento| (sentido antihorario) 
%empezando el de más abajo- luego derecha. Si hay a*b nodos hay (a-1)(b-1)
%elementos
Elementos=zeros((a-1)*(b-1),5)
Elementos(:,1)=[1:(a-1)*(b-1)]' 
for i=0:(b-2) 
    Elementos((a-1)*i+1 : (a-1)*(i+1),2)=[a*i+2 : a*(i+1)]'
    Elementos((a-1)*i+1 : (a-1)*(i+1),3)=[a*i+2 : a*(i+1)]'+a
    Elementos((a-1)*i+1 : (a-1)*(i+1),4)=[a*i+1 : a*(i+1)-1]'+a
    Elementos((a-1)*i+1 : (a-1)*(i+1),5)=[a*i+1 : a*(i+1)-1]'
end
%matriz que tiene |# de esquina |nodo 1 |nodo 2| flujo|
Newmann=zeros(2*(a+b)-4,4)
Newmann(:,1)=[1:2*(a+b)-4]'
Newmann(1:a-1,2)=[1:a-1]'
Newmann(1:a-1,3)=[2:a]'
Newmann(a:a+b-1,2)=[a:a:a*b]
Newmann(a:a+b-1,3)=Newmann(a+1:a+b,2)
Newmann(a+b:2*a+b-2,2)=[a*b-1:-1:a*(b-1)+1]
Newmann(a+b-1:2*a+b-3,3)=Newmann(a+b:2*a+b-2,2)
Newmann(2*a+b-1:2*(a+b)-4,2)=[a*(b-2)+1:-a:a+1]
Newmann(2*a+b-2:2*(a+b)-4,3)=[Newmann(2*a+b-1:2*(a+b)-4,2); Newmann(1,2)]
%matriz que tiene |# de esquina |nodo 1 |nodo 2| 
%Dirichlet
%Tiempo
T=[0:0.1:10];

Solv=zeros([Nodos(:,1) length(T)])
Sol=zeros([a b length(T)])
Solv(:,1)=Ini(:);

%Ensamble de la matriz de rigidez
%Cuadratura de 9 puntos en [0,1]
puntos=([-sqrt(3/5) 0 sqrt(3/5)]+1)*0.5
[x,y]=meshgrid(puntos,puntos)
w=[5/9 8/9 5/9]
W=w'*w
%Rigidez local (valor func (en malla: dos espacios)| posicion en el grad
%(1-2)| subindice 
d_phi(:,:,1,1)= -1+y
d_phi(:,:,2,1)= -1+x
d_phi(:,:,2,2)= -x
d_phi(:,:,1,2)= 1
d_phi(:,:,1,3)= y
d_phi(:,:,2,3)=x
d_phi(:,:,1,4)= -y
d_phi(:,:,2,4)= 1


%Matriz global
A=sparse(Nodos(end:1));
B=sparse(Nodos(end:1));
%Ensamble de A y B
for k=1:length(Nodos(end:1))
    %coeficiente difusion K(x,y)
    K=ones(size(x))
    %Matriz de rigidez local armado
    M=zeros(4)
    for i=1:4
        for j=1:4
            M(i,j)=0.25*sum(sum((W.*K).*(d_phi(:,:,1,i).*d_phi(:,:,1,j)+d_phi(:,:,2,i).*d_phi(:,:,2,j)))); M(i,j)
        end
    end
    %verdadero ensamble
    A(Elementos(k,:),Elementos(k,:))=M;
    B(Elementos(k,:),Elementos(k,:))=[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;
end
%Euler implícito
for t=2:length(T)
    Solv(:,t)=(A+(T(t)-T(t-1))*B)\(B*Solv(:,t-1));
end
%Poniendo como imágenes
for t=1:length(T)
    Sol(:,:,t)=reshape(Solv(:,t) , size(Ini));
end
%Graficando
figure
for t=1:length(T)
    imagesc(Sol(:,:,t))
    colorbar
    pause(0.1)
end




