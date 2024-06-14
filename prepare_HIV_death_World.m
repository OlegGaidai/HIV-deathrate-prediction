clear all
clc
close all
cd('F:\HIV')
fid=fopen('HIV_death_data.txt');
fidC = fopen('Countries.txt'); Z_raw=[]; Z_Country=[];
tline = fgetl(fid);  tlineC = fgetl(fidC);
while ischar(tline)
    numb = str2num(tline);
    if (numb>0)&(numb<1e5)&(length(tlineC)==3)
        Z_raw=[Z_raw; numb];  Z_Country=[Z_Country; tlineC];
    end
    tline = fgetl(fid); tlineC = fgetl(fidC);
end
fclose(fid); fclose(fidC);
Z_raw = Z_raw/1e3;

figure; plot(Z_raw,'s-'); set(gcf,'color','white'); axis tight;
grid on; xlabel('j'); ylabel('R_j','Rotation',90); title('World');
N_tot=length(Z_raw);

%%
IND_C =[];  country_name_running='XXX';
for j=1:N_tot
    if ~strcmp(Z_Country(j,:),country_name_running)
        country_name_running=Z_Country(j,:); IND_C = [IND_C, j];
    end
end

Z_ind=[]; 
N_countries = length(IND_C);

for jj=1:N_countries
    Z_ind = [Z_ind, IND_C(jj)]; 
end
%%

length_R = 30; Z=zeros(1,N_countries*length_R); SS = zeros(length_R, N_countries); 
for j=1:length_R                                                           % day of the year
    shift_j = N_countries*(j-1);
    for jj=1:N_countries
        Z(jj+shift_j) = Z_raw(Z_ind(jj)+j-1); SS(j,jj) = Z_raw(Z_ind(jj)+j-1);
    end
end
figure; plot(Z,'s-'); set(gcf,'color','white'); axis tight;
grid on; xlabel('j'); ylabel('R_j','Rotation',90); title('World');

% fid=fopen('World_R.txt','w'); fprintf(fid,'%e\n',Z); fclose(fid);

[X,Y] = meshgrid(1:length_R, 1:N_countries);
%%
figure; set(gcf,'color','white'); 
s=surf(X,Y,SS','FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
% daspect([5 5 1])
axis tight
view(-50,30)
camlight left
xlabel('years'); ylabel('countries')
s.EdgeColor = 'none';
colorbar
