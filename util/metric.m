function diff = metric(A,B,C,X,Y,Z)
diff = 0;
N = 0;
for i = 1:size(A,2)
    for j = 1:size(A{i},3)
        diff = diff + norm(A{i}(:,:,j)*X*B{i}(:,:,j)-Y*C{i}(:,:,j)*Z, 'fro');
        N = N+1;
    end
end
diff = diff/N;