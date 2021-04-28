function [file_names_final,data] = tecmagImport(folder,nuc)

%This function will take in tecmag exported spectra, contained in a file folder,
%and return a Nx3 array, comprising the frequency/ppm range, and the real and imaginary parts of
%the spectra, respectively. Meant to be used with XTALFIT2
%Data must have Nucleus, and no unnecessary/superfluous files


        file_array = dir(folder);
        
        %initialize file name array
        file_names = cell(1,length(file_array));
        
        for k = 1:1:length(file_array)
           str = file_array(k).name;
           
           if contains(str,'txt') && contains(str,nuc)
               file_names(k) = {str}; 
           end
        end
        
        %clean up and sort array of empty entries
        file_names_final= file_names(~cellfun('isempty',file_names));
        file_names_final = natsortfiles(file_names_final);
        %begin importing data using array of file names
        
        data = cell(1,26);
        y = cell(1,length(file_names_final));
        for t = 1:1:length(file_names_final)
            y{1,t} = table2array(scImport(file_names_final{t})); %Import function
        end
        
        angle = 0;
        
        for k = 1:2:26
            data{1,k} = angle; data{1,k+1} = angle;
            angle = angle + 15;
        end
        data = cell2mat(data);
        z = [];
        for m = 1:length(y)
            temp = y{1,m}; temp(:,2) = []; temp = fliplr(temp); 
            z = [z, temp];
        end
           
        %put it all together
        data = [data;z];                    
end
        

