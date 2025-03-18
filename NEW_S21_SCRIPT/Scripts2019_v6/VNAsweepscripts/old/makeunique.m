function good=makeunique(bad)
%removes double values in first col
[bla,b]=unique(bad(:,1));
good=bad(b,:);