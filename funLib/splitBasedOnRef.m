function splittedVec=splitBasedOnRef(splitThisVec,referenceVec)
refPosIdx=ismember(splitThisVec,referenceVec);
splitIdx=find(refPosIdx);
splittedVec={};
for i=1:2:numel(splitIdx)-1
    startIdx=splitIdx(i);
    endIdx=splitIdx(i+1);
    splittedVec{end+1}=splitThisVec(startIdx:endIdx);
end