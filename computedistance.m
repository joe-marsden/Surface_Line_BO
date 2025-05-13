function dist = computedistance(known, unknown)
% computedistance - calculate distance between two lines
%   summation of the route square difference between corresponding points
%   along the two lines

currentdist = 0;

for i = 1:length(known)
    currentdist = currentdist + sqrt((unknown(i,1)-known(i,1))^2);
end

dist = currentdist;

end