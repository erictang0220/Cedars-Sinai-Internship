% solve current for shim 
function [current] = solveShimCurrent(bodyBzEachCoil, bodyB0, mask, DClimit)
    
% DClimit=7;
[~, ~, ~, coilNum] = size(bodyBzEachCoil);
lb0 = -ones(coilNum,1)*DClimit;
ub0 = ones(coilNum,1)*DClimit;
X0 = zeros(coilNum,1); % initial val

% field offset
lb0(length(lb0)+1) = -10000;
ub0(length(ub0)+1) = 10000;
X0 = zeros(coilNum+1,1);

B0f = bodyB0(mask); % get B0 of the prostate region 

Bzf = zeros(length(B0f),coilNum); % dummy initialization
for i = 1:coilNum        
    Bz_tmp = bodyBzEachCoil(:,:,:,i);
    Bzf(:,i) = Bz_tmp(mask); % get 1A Bz of the prostate region
end

% field offset
Bzf(:,coilNum+1)=0;

% convert from Hz to ppm (ppm is used in Chris's code)
B0f = B0f./123242249;
Bzf = Bzf./123242249;


options7 = optimset('Algorithm','trust-region-reflective','MaxFunEvals',1e14, 'MaxIter',...
    1e12,'TolFun',1e-30,'TolX',1e-16,'Disp','Iter');
[current, ~, ~] = lsqlin(Bzf,B0f,[],[],[],[],lb0,ub0,X0,options7);

end
