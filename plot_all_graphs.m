function [fn_plot_graphs]=plot_all_graphs()
%% Use Thomas' colormaps for comparison
grassmap = uint8([240 200 60; 230 200 60; 220 200 60; 210 200 60; 200 200 60; 180 200 60; 160 200 60; 140 200 60; 120 200 60]);
watermap = uint8([60 200 240; 50 180 220; 40 160 200; 30 140 180; 20 120 180; 10 100 160; 0 80 140]);
nitrogenmap = uint8([50 30 20; 80 40 18; 110 50 16; 140 60 14; 170 70 12; 200 80 10; 230 90 10]);

% FUNCTION - PLOT GRAPHS - a messy business in matlab - best to keep instructions neatly out of everyones way
%% --------------------------------------------------------------------------
%define parameters 
global time fieldsize species Bmax raindata QV
% calculated variables and arrays
global field_species mid_resource deep_resource 
global loop counter     
global time_series_plant time_series_resource plant_transect resource_transect max_connected_cells ave_connected_cells  
% functions
global fn_plot_graphs
%This if loop is set up to plot graphs as data vs year when 312 yrs of tree-ring data are used
if (time>=312)
    year=loop+1658; %use this line when running with 312 years data
else
    year=loop; %i.e. when anything other than the tree-ring rainfall is used
end

%keyboard;
%% --------------------plot field maps-------------------  
nu_fieldsize=fieldsize; %Plot data that is removed from any suggestion of an edge effect (cheating - but ask me why)
plot_field=zeros(nu_fieldsize,nu_fieldsize,species+1);
for h=1:species
    for i=1:nu_fieldsize
        for j=1:nu_fieldsize
            %plot_field(i,j,h)=field_species(i+9,j+9,h);
            plot_field(i,j,h)=field_species(i,j,h);
        end
    end
end
plot_field(:,:,species+1)=plot_field(:,:,1)+plot_field(:,:,2);%Calculate combined biomass of both species
figure;
%colormap(bone);
colormap(grassmap);
imagesc(plot_field(:,:,3));
clim([0 300]);
hold off, colorbar; %hold command may or maynot be required depending on computer
title(['Combined biomass density map for year ',num2str(loop)])
xlabel('Distance across field (m)'), ylabel('Distance downslope (m)')

%figure;
%%colormap(bone);
%colormap(grassmap);

% imagesc(plot_field(:,:,1));
% hold off,colorbar;
% title(['Grass density map for year ',num2str(year)])
% xlabel('Distance across field (m)'), ylabel('Distance downslope (m)')
% 
% figure;
% %colormap(bone);
% colormap(grassmap);
% imagesc(plot_field(:,:,2));
% hold off,colorbar;
% title(['Shrub density map for year ',num2str(year)])
% xlabel('Distance across field (m)'), ylabel('Distance downslope (m)')
% counter=0;

%% --------------------plot transects---------------------------------
% 
% %initialise arrays containing transect data for graphs - transect is taken
% %down center line 25 of the grid (arbitrary choice of line - any can be used
% plant_transect=zeros(nu_fieldsize,species+1);resource_transect=zeros(nu_fieldsize,5); %store for biomass per cell of each species
% for i=1:nu_fieldsize
%     plant_transect(i,1)=i; resource_transect(i,1)=i;       
%     for j=1:species        
%         %plant_transect(i,j+1)=field_species(i,25,j); %Plot actual biomass
%         plant_transect(i,j+1)=field_species(i+9,25,j)/Bmax(j); %Plot normalised biomass
%     end
%     %Use the following 5 lines to plot ACTUAL resource levels
%     %for j=1:randp, %plot actual resource levels
%     %    resource_transect(i,j+1)=top_resource(i,25,j);
%     %    resource_transect(i,j+3)=mid_resource(i,25,j);              
%     %    resource_transect(i,j+5)=deep_resource(i,25,j);
%     %end 
%     %Use the following 9 lines to plot resource data that is normalised against input values data
%     resource_transect(i,3)=mid_resource(i+9,25,2)/QV(2); %Normalise nitrogen - actual nitrogen divided by Nitrogen added to cell
%     resource_transect(i,5)=deep_resource(i+9,25,2)/QV(2); %Normalise nitrogen
%     if (loop<=312) %use to normalise long term data
%         resource_transect(i,2)=mid_resource(i+9,25,1)/raindata(loop,2); %divide water in cell by water added to cell       
%         resource_transect(i,4)=deep_resource(i+9,25,1)/raindata(loop,2);
%     else 
%         resource_transect(i,2)=mid_resource(i+9,25,1)/243; %use to normalise constant rainfall data  (see begining of loop in main)    
%         resource_transect(i,4)=deep_resource(i+9,25,1)/243;
%     end        
% end
% %plot graphs of biomass and resource on transect
% figure;plot(plant_transect(:,1),plant_transect(:,2),'-+',plant_transect(:,1),plant_transect(:,3),'-*');
% legend ('grass','shrub');title(['Normalised species density transect of column 25 for year ',num2str(year)])
% xlabel('Distance across field (m)'), ylabel('Normalised Biomass')
% figure;
% plot(resource_transect(:,1),resource_transect(:,2),'-+',resource_transect(:,1),resource_transect(:,3),'-*',resource_transect(:,1),resource_transect(:,4),'-o',resource_transect(:,1),resource_transect(:,5),'-x');
% legend ('mid water', 'mid nitrogen','deep water','deep nitrogen');
% title(['Normalised resource density transect of column 25 for year ',num2str(year)])   
% xlabel('Distance across field (m)'), ylabel('Normalised Resources')

%% -----GRAPHS AND FINAL OUTPUT FILES----------------------------------------

if (loop==time)
    figure;     %Plot graph of change in biomass with time
    plot(time_series_plant(:,1),time_series_plant(:,2),'-+',time_series_plant(:,1),time_series_plant(:,3),'-*');
    legend ('grass', 'shrub');title ('species density');
    xlabel('Time (years)'), ylabel('Average species density (g/m^2)');
end

%% Plot data in csv file, if required (I used this for debugging)
%    csvwrite('biomass.dat',time_series_plant)   %write biomass/time data to csv file
%    csvwrite('resources.dat',time_series_resource) %write resource/time data fo csv file