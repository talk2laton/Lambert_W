tic
for i = 1:1
    Lambert_W(1, rand(1, 2000000));
end
time1 = toc;

tic
for i = 1:1
    lambertw(1, rand(1, 2000000));
end
time2 = toc;