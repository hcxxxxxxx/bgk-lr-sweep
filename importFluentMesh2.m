function [blks]=importFluentMesh2(str,dms_factor)



fid=fopen(str);
txt = fscanf(fid,'%c');
fclose(fid);
blks = readBlocks(txt);
blks.ncomment = 0;
blks.nbcfaces = [];
while ~isempty(txt)
    [curBlk,txt] = nxtBlk(txt);
    [type,dim] = determineType(curBlk);
    switch type
        case 0
%             blks.ncomment = blks.ncomment+1;
%             blks.comments = readComment(curBlk,blks.comments);
        case 1
            blks.meshGenComments = readMeshGenComment(curBlk, blks.meshGenComments);
        case 2
            blks.dimensions = dim;
            if dim~=2
                error('only 2d meshes are acceptale');
            end
        case 10
            [hdr,data] = readNodes(curBlk);
            if isempty(data)
                blks.meshInfo(1) = hdr(3);
            else
                blks.nodeZones = [blks.nodeZones;hdr];
                blks.points = [blks.points;round(data,10)];
            end
        case 12
            [hdr,data] = readCells(curBlk);
            if isempty(data)
                blks.meshInfo(3) = hdr(3);
            else
                blks.cellZones = [blks.cellZones;hdr];
                blks.celltype = [blks.celltype; data];
            end
        case 13
            [hdr,data] = readFaces(curBlk);
            if isempty(data)
                blks.meshInfo(2) = hdr(3);
            else
                blks.nbcfaces = [blks.nbcfaces,hdr(2),hdr(3)];
                blks.faceZones = [blks.faceZones;hdr];
                blks.faces = [blks.faces;data];
            end
        case 45
            [hdr,data] = readZonesInfo(curBlk);
    end
end
fprintf(' ------------------------------------------------------------\n');
fprintf(' Reading %s...\n',str);
fprintf(' Buffering for file scan...\n\n');

blks.points = blks.points * dms_factor;


% fprintf('\t“%s \n',blks.comments(3));
% fprintf('\t“%s \n',blks.comments(8));
% % fprintf('\t“%s \n',blks.comments(9));
% for i=12:15
%     fprintf('\t“%s \n',blks.comments(i));
% end

end


function [curBlk,txt] = nxtBlk(txt)

o = strfind(txt,'(');
c = strfind(txt,')');
val = ones(size(c));
val = [val,-val];
[l,Index] = sort([o,c]);
val = cumsum(val(Index));
i1 = l(find(val == 1, 1 ));
i2 = l(find(val == 0, 1 ));
curBlk = strtrim(txt(i1+1:i2-1));
txt = strtrim(txt(i2+1:end));
end

function [type,dim] = determineType(curBlk)

type = sscanf(curBlk,'%d');
dim = 0;
if type(1) == 2
    dim = type(2);
    type = type(1);
end
end

function comments = readMeshGenComment(curBlk,comments)

loc = strfind(curBlk,'"');
comments = [comments;string(strtrim(curBlk(loc(1)+1:loc(end)-1)))];
end

function blks = readBlocks(txt)

blks.comments = [];
blks.meshGenComments = [];
blks.dimensions = [];
blks.meshInfo = [];
blks.points = [];
blks.nodeZones = [];
blks.cellZones = [];
blks.cells = [];
blks.faceZones = [];
blks.faces = [];
blks.celltype = [];
end

function comments = readComment(curBlk,comments)

loc = strfind(curBlk,'"');
comments = [comments;string(strtrim( curBlk(loc(1)+1:loc(end)-1)))]; 
end

function [hdr,data] = readHeader(str)
data = [];
oLocs = strfind(str,'(');
cLocs = strfind(str,')');
hdr = sscanf(str(oLocs(1)+1:cLocs(1)-1),'%x');
hdr = hdr';
if numel(cLocs) > 1
    data = str(oLocs(2):cLocs(end));
end
end

function [hdr,data] = readNodes(curBlk)

[hdr,data] = readHeader(curBlk);
if hdr(1) ~= 0
    data = str2num(data(2:end-1));
end
end

function [hdr,data] = readCells(curBlk)

[hdr,data] = readHeader(curBlk);
% if ~ismember(hdr(1),[0,6])
if ~ismember(hdr(1),0)
    if ~isempty(data)
    data = sscanf(data(2:end-1),'%d');
    else
        data = ones(hdr(3) - hdr(2) + 1, 1) * hdr(5);
    end
end
end

function [hdr,data] = readFaces(curBlk)

[hdr,data] = readHeader(curBlk);
if hdr(1) ~= 0
    data = sscanf(data(2:end-1),'%x');
end
if numel(hdr) == 4
    return;
end
if hdr(5) == 0
    data = reshape(data,5,hdr(3)-hdr(2)+1);
    data(1,:) = [];
else
    data = reshape(data,4,hdr(3)-hdr(2)+1);
end
data = data';
end

function [hdr,data] = readZonesInfo(curBlk)
hdr = [];
data = [];
end
