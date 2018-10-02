clear

processedDataFile = [cd, '/SiouxFalls6180_pp.mat'];
load(processedDataFile, 'pathSinks', 'pathSources');


%%
OD_set=[pathSources(1),pathSinks(1)];
ODpath_set=cell(1,1);
ODpath_set{1,1}=1;

for i=2:length(pathSources)
    match=0;
    dummy=OD_set==[pathSources(i),pathSinks(i)];
    if sum(sum(dummy,2)==2)==1
        match=1;
        location=find(sum(dummy,2)==2);
        ODpath_set{location,1}=[ODpath_set{location,1},i];
    end
    
    
    if match==0
        OD_set=[OD_set;[pathSources(i),pathSinks(i)]];
        ODpath_set=[ODpath_set;cell(1,1)];
        ODpath_set{end,1}=i;
    end
end

NumPaths_OD=zeros(size(ODpath_set,1),1);
for i=1:size(ODpath_set,1)
    NumPaths_OD(i)=length(ODpath_set{i,1});
end

OD_set=[OD_set,NumPaths_OD];
save('OD_info.mat', 'ODpath_set', 'OD_set');
            
           