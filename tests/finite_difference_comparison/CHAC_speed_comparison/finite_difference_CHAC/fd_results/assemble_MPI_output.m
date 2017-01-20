function array_global = assemble_MPI_output(prda,path,run_num)

domain_info = importdata(['domain_info' run_num '.frt']);

xg = domain_info(1);
yg = domain_info(2);
zg = domain_info(3);

num_x_subdomain = domain_info(4);
num_y_subdomain = domain_info(5);
num_z_subdomain = domain_info(6);

plot_start_x = 1;
plot_end_x = xg;
plot_start_y = 1;
plot_end_y = yg;
plot_start_z = 1;
plot_end_z = zg;

% Cycle through each subdomain
for x_subdomain = 0:num_x_subdomain-1
    for y_subdomain = 0:num_y_subdomain-1
        for z_subdomain = 0:num_z_subdomain-1
            
            if x_subdomain < num_x_subdomain-1
                px = floor(xg/num_x_subdomain);
            else
                px = xg - (num_x_subdomain-1)*floor(xg/num_x_subdomain);
            end 
            
            if y_subdomain < num_y_subdomain-1
                py = floor(yg/num_y_subdomain);
            else
                py = yg - (num_y_subdomain-1)*floor(yg/num_y_subdomain);
            end 
            
            if z_subdomain < num_z_subdomain-1
                pz = floor(zg/num_z_subdomain);
            else
                pz = zg - (num_z_subdomain-1)*floor(zg/num_z_subdomain);
            end 
            
            Is1 = x_subdomain*floor(xg/num_x_subdomain) + 1;
            Ie1 = Is1 + px - 1;
            
            Is2 = y_subdomain*floor(yg/num_y_subdomain) + 1;
            Ie2 = Is2 + py - 1;
            
            Is3 = z_subdomain*floor(zg/num_z_subdomain) + 1;
            Ie3 = Is3 + pz - 1;
            
            subdomain_index = strcat(sprintf('%02d', x_subdomain), sprintf('%02d', y_subdomain), sprintf('%02d', z_subdomain));

            ext = '.frt';  

            %% domain parameter
            fname = [path prda subdomain_index ext];
            fid = fopen(fname);
            skip1 = fread(fid,1,'int32');
            BD = fread(fid,inf,'double');
            fclose(fid);
            [r1 c1] = size(BD);
            [r1 px*py*pz];
            array_local = reshape(BD(1:r1),px,py,pz);
            
            array_global(Is1:Ie1,Is2:Ie2,Is3:Ie3) = array_local;
            
        end
    end
end
