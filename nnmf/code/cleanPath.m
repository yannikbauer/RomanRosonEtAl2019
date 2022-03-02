% clean matlab path
p = path;
cellpath = strsplit(p,':')';

cc = 10;
for i=1:numel(cellpath)
    if mod(i,floor(numel(cellpath)/10)) == 0 % prints the progress
        fprintf('... %d%% completed\n',cc)
        cc = cc+10;
    end
    
    if contains(cellpath{i},'.git')
        rmpath(cellpath{i})
    end    
end




















