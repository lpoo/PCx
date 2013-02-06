c = [-1; -2];
b = [2; 7 ;3];
A = [-2 1;-1 2;1 2];

lbPrimal = zeros(size(c));
[xprimal, fprimal, flag, output, Y] = PCx(c,sparse(A),b,[],[],lbPrimal, ...
					  [],[]); %,optimset('MaxIter',200,'Display','iter'));
yprimal = Y.ineqlin;

ubDual = zeros(size(b));
[ydual, fdual, flag, output, X] = PCx(-b,sparse(A'),c,[],[],[],ubDual,[]);%, ...
				      %optimset('MaxIter',200,'Display', 'iter'));
xdual = -X.ineqlin;

%**************************************************************************
% primal
%     max [24  14]'*X
%       s.t.   |3   2|      |1200|
%              |4   1| X <= |1000|
%              |2   1|      | 700|
%              X>=0

% c = -[24;14];
% b = [1200;1000;700];
% A = [3 2;4 1;2 1];
% 
% lbPrimal = zeros(size(c));
% [xprimal, fprimal, flag, output, Y] = ...
%     linprog(c,sparse(A),b,[],[],lbPrimal,[],[],optimset('MaxIter',200,'Display','iter'));
% [xprimal, fprimal, flag, output, Y] = ...
%     PCx(c,sparse(A),b,[],[],lbPrimal,[],[],optimset('MaxIter',200,'Display','iter'));
% yprimal = Y.ineqlin;
% yprimal(abs(yprimal)<1e-8)=0;
% 
% lbDual = zeros(size(b));
% [ydual, fdual, flag, output, X] = ...
%     linprog(b,sparse(-A'),c,[],[],lbDual,[],[],optimset('MaxIter',200,'Display','iter'));
% [ydual, fdual, flag, output, X] = ...
%     PCx(b,sparse(-A'),c,[],[],lbDual,[],[],optimset('MaxIter',200,'Display','iter'));
% xdual = X.ineqlin;
% ydual(abs(ydual)<1e-8)=0;

fprintf('          Primal              Dual\n');
fprintf('  OBJ     %f                  %f\n\n', fprimal, fdual);
fprintf('  x1      %f                  %f\n', xprimal(1), xdual(1));
fprintf('  x2      %f                  %f\n\n', xprimal(2), xdual(2));
fprintf('  y1      %f                  %f\n', yprimal(1), ydual(1));
fprintf('  y2      %f                  %f\n', yprimal(2), ydual(2));
fprintf('  y3      %f                  %f\n\n', yprimal(3), ydual(3));

z1 = b - A*xprimal;
fprintf(' compl. slackness for primal: %.6e\n', norm(z1.*yprimal));
fprintf('                              %.6e\n', norm(xprimal.*(c-A'*yprimal)));


z2 = c - (A')*ydual;
fprintf(' compl. slackness for dual: %.6e\n',norm(z2.*xdual));
fprintf('                              %.6e\n', norm(ydual.*(b-A*xdual)));
