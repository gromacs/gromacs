% Matlab code to compute the inverse matrix of poly coefficient aux matrix
BLOCK_SIZE_MIN=3;
BLOCK_SIZE_MAX_1D=4096;
BLOCK_SIZE_MAX_2D=64;
BLOCK_SIZE_MAX_3D=16;

x=generate_coef_inverse_1(6);
x=generate_coef_inverse_2(6,6);
disp(x);
x=generate_coef_inverse_3(6,6,6);
disp(x);

fileID = fopen('PolyRegressionCoefAux1D.f32','w');
for i=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_1D
    Ai=generate_coef_inverse_1(i);
    fwrite(fileID,[i],"single");
    fwrite(fileID,Ai,'single');
end
fclose(fileID);

fileID = fopen('PolyRegressionCoefAux2D.f32','w');
for i=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_2D
    for j=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_2D
        Ai=generate_coef_inverse_2(i,j);
        fwrite(fileID,[i,j],"single");
        fwrite(fileID,Ai,'single');
    end
end
fclose(fileID);

fileID = fopen('PolyRegressionCoefAux3D.f32','w');
for i=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_3D
    for j=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_3D
        for k=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX_3D
            Ai=generate_coef_inverse_3(i,j,k);
            fwrite(fileID,[i,j,k],"single");
            fwrite(fileID,Ai,'single');
        end
    end
end
fclose(fileID);

function Ai = generate_coef_inverse_1(n1)
    M=3;
    A=zeros(M);
    for i = 0:n1-1
        C=[1, i, i * i];
        A=A+C.'*C;
    end
    Ai=inv(A);
end

function Ai = generate_coef_inverse_2(n1,n2)
    M=6;
    A=zeros(M);
    for i = 0:n1-1
        for j=0:n2-1
            C=[1, i, j, i * i, i * j, j * j];
            A=A+C.'*C;
        end
    end
    Ai=inv(A);
end

function Ai = generate_coef_inverse_3(n1,n2,n3)
    M=10;
    A=zeros(M);
    for i = 0:n1-1
        for j=0:n2-1
            for k=0:n3-1
                C=[1, i, j, k, i * i, i * j, i * k, j * j, j * k, k * k];
                A=A+C.'*C;
            end
        end
    end
    Ai=inv(A);
end

