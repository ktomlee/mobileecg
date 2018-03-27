acceptable = load('RECORDS-acceptable.txt');
unacceptable = load('RECORDS-unacceptable.txt');

TotalECG = {};

i = 1;
while i < 774
    TotalECG{1,i} = load([num2str(acceptable(i,1)), '.txt']);
    if i < 226
        TotalECG{2,i} = load([num2str(unacceptable(i,1)), '.txt']);
    end
    i = i + 1;
end