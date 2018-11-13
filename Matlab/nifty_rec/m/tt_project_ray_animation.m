function tt_project_ray_animation(v,det_size,source,vsize_mm,scale,trans,rot,step,proj)

%Jamie McClelland
%Copyright 2010-2011 CMIC-UCL.
%Gower Street, London, UK

if nargin < 9
    proj = tt_project_ray_mex(v,det_size,source,vsize_mm,scale,trans,rot,step);
end

%scale proj so that maximum proj value is same as maximum v value
proj = proj*double(max(v(:)))/max(proj(:));

%make new plot for each projection
for p = 1:size(proj,3)
    
    %first find out where detector is and plot image there
    origin = [0 0 0 1]';
    xend = [scale(p,1) 0 0 1]';
    yend = [0 scale(p,2) 0 1]';
    
    x_rot = formAffineMatrix([0 0 0 -180*rot(p,1)/pi 0 0]);
    y_rot = formAffineMatrix([0 0 0 0 -180*rot(p,2)/pi 0]);
    z_rot = formAffineMatrix([0 0 0 0 0 -180*rot(p,3)/pi]);
    aff_mat = formAffineMatrix([trans(p,:) 0 0 0])*(z_rot*(y_rot*x_rot));
    origin = aff_mat*origin;
    xend = aff_mat*xend;
    yend = aff_mat*yend;
    
    cla;
    h_proj = tt_imageIn3D(proj(:,:,p),origin(1:3)',xend(1:3)',yend(1:3)');
    %set(h_proj,'facealpha',0.5);
    hold on;
    
    %draw box around volume and plot 3 central slices
    verts = [0 0 0;
        vsize_mm(1) 0 0;
        vsize_mm(1) vsize_mm(2) 0;
        0 vsize_mm(2) 0;
        0 0 vsize_mm(3);
        vsize_mm(1) 0 vsize_mm(3);
        vsize_mm(1) vsize_mm(2) vsize_mm(3);
        0 vsize_mm(2) vsize_mm(3)];
    faces = [1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        4 1 5 8;
        1 4 3 2;
        5 8 7 6];
    patch('vertices',verts,'faces',faces,'edgecolor','k','facecolor','none');
    
    vsize_vox = size(v);
    vdim = vsize_mm./vsize_vox;
    xslice_vox = round(vsize_vox(1)./2);
    xslice_mm = (xslice_vox-0.5)*vdim(1);
    yslice_vox = round(vsize_vox(2)./2);
    yslice_mm = (yslice_vox-0.5)*vdim(2);
    zslice_vox = round(vsize_vox(3)./2);
    zslice_mm = (zslice_vox-0.5)*vdim(3);
%     h_xslice = tt_imageIn3D(squeeze(v(xslice_vox,:,:)),[xslice_mm 0 0],[xslice_mm vsize_mm(2) 0],[xslice_mm 0 vsize_mm(3)]);
%     h_yslice = tt_imageIn3D(squeeze(v(:,yslice_vox,:)),[0 yslice_mm 0],[vsize_mm(1) yslice_mm 0],[0 yslice_mm vsize_mm(3)]);
%     h_zslice = tt_imageIn3D(v(:,:,zslice_vox),[0 0 zslice_mm],[vsize_mm(1) 0 zslice_mm],[0 vsize_mm(2) zslice_mm]);
%     set(h_xslice,'facealpha',0.5);
%     set(h_yslice,'facealpha',0.5);
%     set(h_zslice,'facealpha',0.5);


    %plot source and projection lines from source to corners of detector
    plot3(source(p,1),source(p,2),source(p,3),'o','markersize',10,'markerfacecolor','r','markeredgecolor','r');
    xyend = yend + xend - origin;
%     s_to_o_step = step*(origin(1:3)' - source(p,:))/sqrt(sum((origin(1:3)' - source(p,:)).^2));
%     s_to_x_step = step*(xend(1:3)' - source(p,:))/sqrt(sum((xend(1:3)' - source(p,:)).^2));
%     s_to_y_step = step*(yend(1:3)' - source(p,:))/sqrt(sum((yend(1:3)' - source(p,:)).^2));
%     s_to_xy_step = step*(xyend(1:3)' - source(p,:))/sqrt(sum((xyend(1:3)' - source(p,:)).^2));
%     plot3(source(p,1):s_to_o_step(1):origin(1),source(p,2):s_to_o_step(2):origin(2),source(p,3):s_to_o_step(3):origin(3),'.g');
%     plot3(source(p,1):s_to_x_step(1):xend(1),source(p,2):s_to_x_step(2):xend(2),source(p,3):s_to_x_step(3):xend(3),'.g');
%     plot3(source(p,1):s_to_y_step(1):yend(1),source(p,2):s_to_y_step(2):yend(2),source(p,3):s_to_y_step(3):yend(3),'.g');
%     plot3(source(p,1):s_to_xy_step(1):xyend(1),source(p,2):s_to_xy_step(2):xyend(2),source(p,3):s_to_xy_step(3):xyend(3),'.g');
    plot3([source(p,1) origin(1)],[source(p,2) origin(2)],[source(p,3) origin(3)],'g');
    plot3([source(p,1) xend(1)],[source(p,2) xend(2)],[source(p,3) xend(3)],'g');
    plot3([source(p,1) yend(1)],[source(p,2) yend(2)],[source(p,3) yend(3)],'g');
    plot3([source(p,1) xyend(1)],[source(p,2) xyend(2)],[source(p,3) xyend(3)],'g');
    
    
    
    daspect([1 1 1]);
    view(3)
    lims = axis;
    axis([min(lims(1),lims(3)) max(lims(2),lims(4)) min(lims(1),lims(3)) max(lims(2),lims(4)) lims(5) lims(6)]);
    axis vis3d;
    colormap(gray);
    
    drawnow;
    
end
