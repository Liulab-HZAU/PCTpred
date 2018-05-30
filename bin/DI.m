addpath(genpath('../softwares/calculate_evolutionary_constraints_v2.0'))
[ex] = textread('new_protein.txt','%s');
for i = 1:1
    a = ex{i};
    disp(i);
    disp(a);
    path_in = strcat ('../result_files/co-evolution/',a,'.msa');
    path_out = strcat ('../result_files/co-evolution/',a,'.txt');
    calculate_evolutionary_constraints(path_in,'Query_1',path_out)
end
