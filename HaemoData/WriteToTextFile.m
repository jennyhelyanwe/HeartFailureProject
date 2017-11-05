function WriteToTextFile(x, filename, pwp_toggle)
fw = fopen(filename, 'w');
% ED frame. 
frame_num = 1;
fprintf(fw, '1\t');
for i = 1:length(x.ds_ed(end,:))
    fprintf(fw, '%f\t', x.ds_ed(end, i));
end
fprintf(fw, '\n');
% ED to end IVC frames. 
for i = 1:length(x.ed_eivc(:,1));
    fprintf(fw, '%d\t', frame_num + i);
    for j = 1:length(x.ed_eivc(1,:));
        fprintf(fw, '%f\t', x.ed_eivc(i,j));
    end
    fprintf(fw, '\n');
end
frame_num = frame_num + i;
% End IVC to ES frames. 
for i = 1:length(x.eivc_es(:,1));
    fprintf(fw, '%d\t', frame_num + i);
    for j = 1:length(x.eivc_es(1,:));
        fprintf(fw, '%f\t', x.eivc_es(i,j));
    end
    fprintf(fw, '\n');
end
frame_num = frame_num + i;
if pwp_toggle == 1
    % ES to end IVR frames. 
    for i = 1:length(x.es_eivr(:,1));
        fprintf(fw, '%d\t', frame_num + i);
        for j = 1:length(x.es_eivr(1,:));
            fprintf(fw, '%f\t', x.es_eivr(i,j));
        end
        fprintf(fw, '\n');
    end
    frame_num = frame_num + i;
    % End IVR to DS frames.
    for i = 1:length(x.eivr_ds(:,1));
        fprintf(fw, '%d\t', frame_num + i);
        for j = 1:length(x.eivr_ds(1,:));
            fprintf(fw, '%f\t', x.eivr_ds(i,j));
        end
        fprintf(fw, '\n');
    end
    frame_num = frame_num + i;
else
    % ES to DS frames. 
     for i = 1:length(x.es_ds(:,1));
        fprintf(fw, '%d\t', frame_num + i);
        for j = 1:length(x.es_ds(1,:));
            fprintf(fw, '%f\t', x.es_ds(i,j));
        end
        fprintf(fw, '\n');
    end
    frame_num = frame_num + i;
end
% DS to ED frames. 
for i = 1:length(x.ds_ed(1:(end-1),1));
    fprintf(fw, '%d\t', frame_num + i);
    for j = 1:length(x.ds_ed(1,:));
        fprintf(fw, '%f\t', x.ds_ed(i,j));
    end
    fprintf(fw, '\n');
end
fclose(fw);
end