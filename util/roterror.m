function err=roterror(X1,X2)
err = norm( so3_vec( skewlog( ( X1(1:3,1:3)'*X2(1:3,1:3) ) ) ) );
end