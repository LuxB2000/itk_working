%
% Create the test date using Matlab
%
function [] = create_data()

    function v = create_sq(v,bor,c,val,verbose)
        
        for x=c(1)-bor(1):c(1)+bor(1)
            for y=c(2)-bor(2):c(2)+bor(2)
                for z=c(3)-bor(3):c(3)+bor(3)
                    v(x,y,z) = val;
                end
            end
        end
        
        if(verbose)
            for z=1:size(v,3)
                figure, imshow(v(:,:,z),[]), title(num2str(z));
            end
        end
        
    end
    
    function [fid] = write_mha_header(info,fid)
        fprintf(fid,strcat('ObjectType = ', 32, info.objectType,'\n'));
        fprintf(fid,strcat('NDims = ', 32, num2str(info.numberOfDimensions),'\n'));
        fprintf(fid,strcat('BinaryData = ', 32, info.binaryData,'\n'));
        fprintf(fid,strcat('BinaryDataByteOrderMSB = ', 32, info.binaryDataByteOrderMSB,'\n'));
        fprintf(fid,strcat('CompressedData = ', 32, info.compressedData,'\n'));
        fprintf(fid,strcat('TransformMatrix = ', 32 ...
            , num2str(info.transformMatrix(1)), 32 ,num2str(info.transformMatrix(2)), 32 ,num2str(info.transformMatrix(3)), 32  ...
            , num2str(info.transformMatrix(4)), 32 ,num2str(info.transformMatrix(5)), 32 ,num2str(info.transformMatrix(6)), 32  ...
            , num2str(info.transformMatrix(7)), 32 ,num2str(info.transformMatrix(8)), 32 ,num2str(info.transformMatrix(9)) ...
            , '\n'));
        fprintf(fid,strcat('Offset = ', 32 ...
            , num2str(info.offset(1)), 32, num2str(info.offset(2)), 32, num2str(info.offset(3)) ...
            , '\n'));
        fprintf(fid,strcat('CenterOfRotation = ', 32 ...
            , num2str(info.centerOfRotation(1)), 32 ...
            , num2str(info.centerOfRotation(2)), 32 ...
            , num2str(info.centerOfRotation(3)) ...
            , '\n'));
        fprintf(fid,strcat('AnatomicalOrientation = ', 32, info.anatomicalOrientation, '\n' ));
        
        fprintf(fid,strcat('ElementSpacing = ', 32 ...
            , num2str(info.elementSpacing(1)), 32 ...
            , num2str(info.elementSpacing(2)), 32 ...
            , num2str(info.elementSpacing (3)) ...
            , '\n'));
        fprintf(fid,strcat('DimSize = ', 32 ...
            , num2str(info.elementSize(1)), 32 ...
            , num2str(info.elementSize(2)), 32 ...
            , num2str(info.elementSize (3)) ...
            , '\n'));
        fprintf(fid,'ElementNumberOfChannels = %i\n',info.ElementNumberOfChannels);
        fprintf(fid,strcat('ElementType = ', 32, info.elementType, '\n'));
        fprintf(fid,strcat('ElementDataFile = ', 32, info.elementDataFile,'\n'));
    end

    function [] = write_mha_volume(v,spacing, fileName)
        info.fileName = strcat(fileName, '.mha' );
        info.numberOfDimensions=3;
        info.objectType='Image';
        info.binaryData='True';
        info.binaryDataByteOrderMSB='False';
        info.compressedData='False';
        info.transformMatrix=[1 0 0 0 1 0 0 0 1];
        info.offset=[0 0 0];
        info.centerOfRotation=[0,0,0];
        info.anatomicalOrientation='RAI';
        info.elementSpacing=spacing;
        info.elementSize=size(v); %- ones(1,3);
        info.ElementNumberOfChannels = 1;
        info.elementType='MET_USHORT';
        info.elementDataFile='LOCAL';

        % open the file, if exist replace it.
        fid=fopen(info.fileName,'w');
        if(fid<0)
            fprintf('could not open file %s\n',info.fileName);
            return
        end
        % write the header
        write_mha_header(info,fid);
        fclose(fid);
        
        % write the image
        fid=fopen(info.fileName,'a','native');
        fwrite(fid,uint16(v),'uint16');
        fclose(fid);
        
%         sz = size(v);
%         for x=1:sz(1)
%             for y=1:sz(2)
%                 for z=1:sz(3)
%                     fprintf(fid,'(v(x,y,z)));
%                 end
%             end
%         end
    end

V = zeros(20,30,10);
spacing = [0.5,0.5,0.35];

V_fixed = create_sq( V, [3,5,3], [12,10,5], 125, 0 );
V_fixed = create_sq( V_fixed, [3,3,3], [12,14,5], 255, 0 );

write_mha_volume( V_fixed, spacing, 'fixed_volume' );

V_moving = create_sq( V, [3,5,3], [8,10,5], 125, 0 );
V_moving = create_sq( V_moving, [3,3,3], [8,14,5], 255, 0 );

write_mha_volume( V_moving, spacing, 'moving_volume' );

V_atlas = create_sq( V, [3,5,3], [15,20,4], 1, 0 );
V_atlas = create_sq( V_atlas, [3,3,3], [15,24,4], 2, 0 );

write_mha_volume( V_atlas, spacing, 'atlas_volume' );


end
