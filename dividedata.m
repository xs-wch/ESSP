function [target_data,face_data,face_invert_data,object_data] = dividedata(mycnt,ana_begin,ana_end)
% input: 
% mycnt: the data structure loaded from .cnt file
% ana_begin: the begining of the ERP
% ana_end: the end of the ERP

% output: 
% format:channel*sample*trail
% target_data: the response of stimulation 'target'
% face_data: the response of stimulation 'face'
% face_invert_data: the response of stimulation 'inverted face'
% object_data: the response of stimulation 'object'

n=length(mycnt.event);
target_num = 0;face_num = 0;face_invert_num = 0;object_num = 0;

for i=1:n
    switch mycnt.event(i).stimtype
        case 1 % target
           target_num = target_num+1; 
        case 2 % face
            face_num = face_num+1;
        case 3 % inverted face
            face_invert_num = face_invert_num+1;
        case 4 % object
            object_num = object_num+1;
    end
end

[channum,len] = size(mycnt.data);
clear len
target_data=zeros(channum,ana_end-ana_begin,target_num);
face_data=zeros(channum,ana_end-ana_begin,face_num);
face_invert_data=zeros(channum,ana_end-ana_begin,face_invert_num);
object_data = zeros(channum,ana_end-ana_begin,object_num);

i1=1;i2=1;i3=1;i4=1;

for i=1:n
    switch mycnt.event(i).stimtype
        case 1
            target_data(:,:,i1) = mycnt.data(:,mycnt.event(i).offset+ana_begin+1:mycnt.event(i).offset+ana_end);
            i1=i1+1;
        case 2
            face_data(:,:,i2) = mycnt.data(:,mycnt.event(i).offset+ana_begin+1:mycnt.event(i).offset+ana_end);
            i2=i2+1;
        case 3
            face_invert_data(:,:,i3) = mycnt.data(:,mycnt.event(i).offset+ana_begin+1:mycnt.event(i).offset+ana_end);
            i3=i3+1;
        case 4
            object_data(:,:,i4) = mycnt.data(:,mycnt.event(i).offset+ana_begin+1:mycnt.event(i).offset+ana_end);
            i4=i4+1;
    end
end
