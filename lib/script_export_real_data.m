
visibility_file_name = 'data/at166B.3C129.c0';
filename_stokesI = [visibility_file_name,'I.vis']; % stokes I
filename_stokesV = [visibility_file_name,'V.vis']; % stokes V


format='%f %f %f %f %f'; % we have 5 columns, [u v y_r y_i sigma]

% read filename_stokesV

FID = -1;
t = 0;
while FID<1 && t<100
     t=t+1;
    FID = fopen(filename_stokesV)
end
if FID <1
     disp('Error opening the file');
end

outV = textscan(FID, format);
fclose(FID);
dataV=cell2mat(outV);

% noise std estimate
y_V = (dataV(:,3)./dataV(:,5)) + 1i*(dataV(:,4)./dataV(:,5)); %whitened data V

% read filename_stokesI
FID = -1;
t = 0;
while FID<1 && t<100
     t=t+1;
    FID = fopen(filename_stokesI);
end
if FID <1
     disp('Error opening the file');
end

outI = textscan(FID, format);
fclose(FID);
dataI=cell2mat(outI);

y_I = (dataI(:,3)+1i*dataI(:,4))./dataI(:,5); %whitened data I

uvw = zeros(length(y_I),3);
uvw(:,1) = dataI(:,1);
uvw(:,2) = dataI(:,2);

save('data.mat','uvw','y_I','y_V','-v7.3');
