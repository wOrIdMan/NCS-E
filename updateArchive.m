function [archiveX, archiveObj, archiveG, archiveH, archiveLen] = updateArchive(archiveX, archiveObj, archiveG, archiveH, offspringX, offspringObj, offspringG, offspringH, archiveLen, archiveMaxLen)
m = size(offspringX,1);
for i = 1:m
    if archiveLen<archiveMaxLen
        archiveLen = archiveLen+1;
        replaceedI = archiveLen;
    else
        replaceedI = randi(archiveMaxLen);
    end
    archiveX(replaceedI,:) = offspringX(i,:);
    archiveObj(replaceedI) = offspringObj(i);
    archiveG(replaceedI,:) = offspringG(i,:);
    archiveH(replaceedI,:) = offspringH(i,:);
end

end