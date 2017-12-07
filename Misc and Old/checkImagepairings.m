%check to make sure all images have their pairing e.g. novel and respective
%replaced image with same naming

cd('C:\Users\seth.koenig\Desktop\All moved\')

list = ls;

for l = 1:length(list)
    bmp = strfind(list(l,:),'.bmp');
    if ~isempty(bmp)
        if ~strcmp(list(l,bmp-1),'m') %then novel image
            if ~exist([list(l,1:bmp-1) 'm.bmp'],'file')
               disp(list(l,:))
            end
        end
    end
end

%%
%for REPEATS ONLY make sure none of the pixels have been manipulated

cd('C:\Users\seth.koenig\Desktop\All repeats\')

list = ls;

for l = 1:length(list)
    bmp = strfind(list(l,:),'.bmp');
    if ~isempty(bmp)
        if ~strcmp(list(l,bmp-1),'p') %then novel image
            img1 = imread(list(l,:));
            img2 = imread([list(l,1:bmp-1) 'p.bmp']);
            if any(any(any(img1 ~= img2))) %a pixel has changed
               disp(list(l,:))
            end
        end
    end
end