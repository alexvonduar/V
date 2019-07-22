function run
if ~isdeployed
    addpath util
end

%% todo options

todo.save_results_image = 0;
todo.benchmark = 0;
todo.calibrate = 0;
todo.ortho_rectify = 0;
todo.save_ortho_images = 0;
todo.debugfile = 1;

%% plot options (0 if not plotted, figure number otherwise)
plots.hvps = 1;%1;
plots.z = 1;%1;
plots.hl = 1;%1;
plots.gthl = 0;%1;
plots.benchmark = 0;%2; % display precision curve and AUC for each dataset (the figure number will be plots.benchmark*100+#dataset)
plots.manhattan = 0;%1;
plots.orthorectify = 0; %3; rectified images (the figure number will be plots.orthorectify*100+#plane)

%% path(s) to the image (set(s))
imgDir = cell(1,1);

% some test images provided with this code
imgDir{1,1} = './test_images/'; 

% the following datasets, as well as the ground truth horizon lines in a
% unified format can be obtained from our webpage:
% https://members.loria.fr/GSimon/software/v/

% imgDir{1,end+1} = '/home/gsimon/ownCloud/Matlab/fastvp2/YorkUrbanDB/';
% imgDir{1,end+1} = '/home/gsimon/ownCloud/Matlab/fastvp2/EurasianCitiesBase/';
% imgDir{1,end+1} = '/home/gsimon/Documents/MATLAB/HLW dataset/images/Tests/';

%% for each image set...
for is = 1:size(imgDir,2)
    imageList = glob([imgDir{1,is}, '*.jpg'])
    nImages = numel(imageList);
    params = default_params();
    params.include_infinite_hvps = 1;  % includes the detection of infinite 
    % horizontal VPs. TO BE REMOVED if one want to get the same results as
    % in our ECCV'2018 paper
    hl_error = [];
    %     scene = strsplit(imgDir{1,i}, '/');
    %     scene = scene{end-1};
    for i = 1:nImages
        close all

        % read the image
        fprintf('%d / %d\n', i, nImages);
        im = imread(imageList{i});
        width = size(im, 2);
        height = size(im, 1);
        focal = max(width, height) / 2; % fake focal length

        if todo.debugfile
            [pathstr, name, ext] = fileparts(imageList{i});
            debug_fname = strcat(pathstr,'/',name,'_mat_debug.txt');
            params.debug_fileid = fopen(debug_fname, 'w');
        end

        % call V
        [hl, hvps, hvp_groups, z, z_group, ls] = V(im, width, height, focal, params);
        
        if params.debug_fileid > 2
            fclose(params.debug_fileid);
        end
       
        % plot the results
        if plots.hvps
            figure(plots.hvps); clf(); imshow(im); hold on;
            ax = gca;
            for j = 1:numel(hvp_groups)
                hg = hvp_groups{j};
                for k = 1:length(hg)
                    %                     ax.ColorOrderIndex = j;
                    %                     plot([hvps(j,1);(ls(hg(k),1)+ls(hg(k),3))/2],[hvps(j,2);(ls(hg(k),2)+ls(hg(k),4))/2], '--', 'linewidth', 1);
                    ax.ColorOrderIndex = j;
                    plot([ls(hg(k),1);ls(hg(k),3)],[ls(hg(k),2);ls(hg(k),4)], '-', 'linewidth', 2);
                end
            end
            drawnow();
        end
        if plots.z
            figure(plots.z);
            if plots.z ~= plots.hvps
                clf(); imshow(im); hold on;
                ax = gca;
            end
            zg = z_group;
            for k = 1:length(zg)
                %                 ax.ColorOrderIndex = numel(hvp_groups)+1;
                %                 plot([z(1);(ls(zg(k),1)+ls(zg(k),3))/2],[z(2);(ls(zg(k),2)+ls(zg(k),4))/2], '--', 'linewidth', 1);
                ax.ColorOrderIndex = numel(hvp_groups)+1;
                plot([ls(zg(k),1);ls(zg(k),3)],[ls(zg(k),2);ls(zg(k),4)], '-', 'linewidth', 2);
            end
            drawnow();
        end
        if plots.hl
            figure(plots.hl);
            if plots.hl ~= plots.hvps && plots.hl ~= plots.z
                clf(); imshow(im); hold on;
            end
            plot(hl(1:2,1), hl(1:2,2), '-c', 'linewidth', 3);
            drawnow();
        end
        if plots.gthl || todo.benchmark
            [pathstr, name, ext] = fileparts(imageList{i});
            gthl_fname = strcat(pathstr,'/',name,'hor.mat');
            if exist(gthl_fname, 'file')
                m = load(gthl_fname);
                gthl_homo = m.horizon;
                gthl_a = gthl_homo(1);
                gthl_b = gthl_homo(2);
                gthl_c = gthl_homo(3);
                gthl(1,1) = 0;
                gthl(1,2) = -gthl_c/gthl_b;
                gthl(2,1) = width;
                gthl(2,2) = -(gthl_a*width+gthl_c)/gthl_b;
            end
        end
        if plots.gthl && exist(gthl_fname, 'file')
            figure(plots.gthl);
            if plots.gthl ~= plots.hl
                clf(); imshow(im); hold on;
            end
            plot(gthl(1:2,1), gthl(1:2,2), '--y', 'linewidth', 3);
            drawnow();
        end
        
       % benchmark
        if todo.benchmark && exist(gthl_fname, 'file')
            dist1 = abs(gthl(1,2) - hl(1,2));
            dist2 = abs(gthl(2,2) - hl(2,2));
            hl_error(i) = max(dist1,dist2)/height;
            if plots.benchmark
                figure(plots.benchmark*100+is); clf(); hold on;
                auc = calc_auc(hl_error',1,'AUC Simon et al.',1);
                drawnow();
            else
                auc = calc_auc(hl_error',0,'',0);
            end
        end
        
        % find Manhattan directions and calibrate the camera
        
        if todo.calibrate
            [focal, manh_vps, confident] = calibrate(z, hvps, width, height);
            if plots.manhattan && confident >= 0
                figure(plots.manhattan);
                u0 = width/2;
                v0 = height/2;
                if z(2) > v0
                    posy = 0;
                else
                    posy = height;
                end
                switch confident
                    case 3
                        plot([u0; manh_vps(1,1)], [posy; manh_vps(1,2)], '-r', 'LineWidth', 3);
                        plot([u0; manh_vps(2,1)], [posy; manh_vps(2,2)], '-g', 'LineWidth', 3);
                        plot([u0; manh_vps(3,1)], [posy; manh_vps(3,2)], '-b', 'LineWidth', 3);
                    case 2
                        plot([u0; manh_vps(1,1)], [posy; manh_vps(1,2)], '--r', 'LineWidth', 3);
                        plot([u0; manh_vps(2,1)], [posy; manh_vps(2,2)], '-g', 'LineWidth', 3);
                        plot([u0; manh_vps(3,1)], [posy; manh_vps(3,2)], '-b', 'LineWidth', 3);
                    case 1
                        plot([u0; manh_vps(1,1)], [posy; manh_vps(1,2)], '--r', 'LineWidth', 3);
                        plot([u0; manh_vps(2,1)], [posy; manh_vps(2,2)], '--g', 'LineWidth', 3);
                        plot([u0; manh_vps(3,1)], [posy; manh_vps(3,2)], '-b', 'LineWidth', 3);
                    otherwise
                end
                drawnow();
            end
        end
        
        % save the results image
        if plots.hvps && todo.save_results_image
            [pathstr, name, ext] = fileparts(imageList{i});
            fname = strcat(pathstr,'/',name,'res.png');
            saveas(figure(plots.hvps),fname,'png');
        end
                
        % ortho-rectify all vertical planes
        
        if todo.ortho_rectify 
            K = [];
            if focal > 0
                K = [focal 0 width/2; 0 focal height/2; 0 0 1];
            end
            hl_homo = line_hmg_from_two_points([hl(1,1) hl(1,2)],[hl(2,1) hl(2,2)]);
            [imR,maskR,transform] = orthorectify_from_vps_and_lines(im, hvps, hvp_groups, z, z_group, ls, 4, K, hl_homo, 0);
            if todo.save_ortho_images
                root = strsplit(imageList{i},'.jpg');
                root = root{1};
                for j = 1:numel(imR)
                    if ~isempty(imR{j})
                        imwrite(imR{j},strcat(root,'_R_',num2str(j),'.jpg'));
                        %                         imwrite(maskR{j},strcat(root,'_Rmask_',num2str(j),'.jpg'));
                        %                         tform = transform{j};
                        %                         save(strcat(root,'_tform_',num2str(j),'.mat'),'tform');
                    end
                end
            end
            if plots.orthorectify
                for j = 1:numel(imR)
                    if ~isempty(transform{j})
                        figure(plots.orthorectify*100+j); clf();
                        imshow(imR{j}); drawnow();
                    end
                end
            end
        end
    end
end
