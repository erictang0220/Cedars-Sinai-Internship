function [sumBz] = BzCalc(bodyB0, bodyBzEachCoil, current, prostateMask)
nCoils = size(bodyBzEachCoil, 4);
% orginal 
sumBz = zeros(192,192,96);

% physical phantom
% sumBz = zeros(112,84,36);

% digital phantom
% sumBz = zeros(512,512,40);

% digital phantom 3
% sumBz = zeros(512,512,50);

for k = 1:nCoils
    tmp = bodyBzEachCoil(:,:,:,k)*current(k);
    sumBz = sumBz + tmp;
end
% figure;hist(sumBz(prostateMask),100)
% figure;hist(bodyB0(prostateMask)-sumBz(prostateMask),100)
end
