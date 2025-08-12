clear
T = readtable('sdp_Mittelmann_mosek1121.xlsx');
for i = 1:size(T,1)
out  = struct();
out.name = table2array(T(i,1));
out.name = out.name{1};
out.m = table2array(T(i,2));
out.n = table2array(T(i,3));
out.Xrank =  table2array(T(i,4));
out.Srank =  table2array(T(i,5));
out.pinf = table2array(T(i,6));
out.dinf = table2array(T(i,7));
out.gap = table2array(T(i,8));
out.iter = table2array(T(i,9));
out.time = table2array(T(i,10));
filename = strcat(out.name,'_res.mat')
save(filename,'out')
1;
end
1;