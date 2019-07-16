function F = GetFreq(X,Y,FrameRate,XC,YC)
    X = X - XC;
    Y = Y - YC;
    theta = unwrap(atan2(Y,X));
    F = -FrameRate*diff(theta)/(2*pi);
end