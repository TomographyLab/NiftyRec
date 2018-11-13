function h = tt_imageIn3D(im,origin,xend,yend)

%Jamie McClelland
%Copyright 2010-2011 CMIC-UCL.
%Gower Street, London, UK

%function to draw a 2D image in 3D space
%im is the 2D image
%origin is the 3D location of 'outside corner' of the 1st pixel in im
%xend is the 3D location of the 'outside corner' of the last pixel in the
%first row of im (i.e. im(end,1)
%yend is the 3D lcoation of the 'outside corner' of the last pixel in the
%first column of im (i.e. im(1,end)

imsize = size(im);

x_vect = (xend-origin)./imsize(1);
y_vect = (yend-origin)./imsize(2);

verts = zeros(prod(imsize+[1 1]),3);
verts(1,:) = origin;
v = 2;
for x = 1:imsize(1)+1
    for y = 1:imsize(2)
        verts(v,:) = verts(v-1,:) + y_vect;
        v = v + 1;
    end
    verts(v,:) = verts(v-1-imsize(2),:) + x_vect;
    v = v + 1;
end

faces = zeros(prod(imsize),4);
cdata = zeros(prod(imsize),1);
f = 1;
for x = 1:imsize(1)
    for y = 1:imsize(2)
        faces(f,:) = ((x-1)*(imsize(2)+1)+y)+[0 1 imsize(2)+2 imsize(2)+1];
        cdata(f) = im(x,y);
        f = f + 1;
    end
end

h = patch('vertices',verts,'faces',faces,'facevertexcdata',cdata,'edgecolor','none','facecolor','flat');
