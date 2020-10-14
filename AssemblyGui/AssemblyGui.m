function varargout = AssemblyGui(varargin)
% ASSEMBLYGUI MATLAB code for AssemblyGui.fig
%      ASSEMBLYGUI, by itself, creates a new ASSEMBLYGUI or raises the existing
%      singleton*.
%
%      H = ASSEMBLYGUI returns the handle to a new ASSEMBLYGUI or the handle to
%      the existing singleton*.
%
%      ASSEMBLYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSEMBLYGUI.M with the given input arguments.
%
%      ASSEMBLYGUI('Property','Value',...) creates a new ASSEMBLYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AssemblyGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AssemblyGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AssemblyGui

% Last Modified by GUIDE v2.5 02-Mar-2018 03:50:39

% Begin initialization code - DO NOT EDIT

% Reference paper: Herzog et al. 2020 "Scalable and accurate automated method 
% for neuronal ensemble detection in spiking neural networks"
% https://www.biorxiv.org/content/10.1101/2020.10.12.335901v1
% Rubén Herzog October 2020
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AssemblyGui_OpeningFcn, ...
    'gui_OutputFcn',  @AssemblyGui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AssemblyGui is made visible.
function AssemblyGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AssemblyGui (see VARARGIN)

% checking if parallel comp is present
pkgs=ver;
if sum(strcmp('Parallel Computing Toolbox',{pkgs(:).Name}))==1
    % if present, start parallel pool
    gcp;
end
handles.output = hObject;

% default parameters
handles.pars.filename = '';
handles.pars.sr = 20000;
handles.pars.bin = 0.02;

handles.pars.npcs = 3;
handles.pars.dc = 0.01;
handles.pars.minspk = 3;
handles.pars.minsize = 3;
handles.pars.cent_thr = 99.9;
handles.pars.nsur = 1000;
handles.pars.prct = 99.9;
handles.pars.inner_corr = 5;

% variables to store the last value
handles.pars.old_npcs = nan;
handles.pars.old_dc = nan;
handles.pars.old_minspk = nan;
handles.pars.old_minsize = nan;
handles.pars.old_cent_thr = nan;
handles.pars.old_nsur = nan;
handles.pars.old_prct = nan;
handles.pars.old_inner_corr = nan;

% Choose default command line output for AssemblyGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AssemblyGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AssemblyGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------- LOADING DATA --------------------------------------------
% --- Executes on button press in load_butt.
function load_butt_Callback(hObject, eventdata, handles)
% hObject    handle to load_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open file dialog
cla(handles.raster_axes)
cla(handles.delta_rho_axes)
cla(handles.pca_axes)
cla(handles.ens_ras_axes)
cla(handles.corr_axes)
cla(handles.core_axes)
cla(handles.inner_corr_axes)
cla(handles.npcs_axes)
try
    if isfield(handles,'filepath')
        [handles.pars.filename,handles.filepath,filetype] = uigetfile({'*.mat;*.hdf5;*.nex'},'Choose file with timestamps...',handles.filepath);
    else
        [handles.pars.filename,handles.filepath,filetype] = uigetfile({'*.mat;*.hdf5;*.nex'},'Choose file with timestamps...');
    end
    [~,~,ext] = fileparts(handles.pars.filename);
catch
    disp(['User cancelled'])
    return;
end
prev_data = 0;
% checking data and loading files
if handles.pars.filename~=0
    %try
        switch ext
            case '.mat'
                %  ---------- TO DO: allow loading file previously
                %  analyzed --------------
                matfile = load([handles.filepath,handles.pars.filename]); % loading mat on matfile structure
                fnames = fieldnames(matfile); % fields of the structure matfile
                
                if isstruct(matfile.(fnames{1})) % previousaly analyzed data
                    handles.data = matfile.(fnames{1});
                    handles.data=rmfield(handles.data,'pars');
                    handles.pars = matfile.(fnames{1}).pars;
                    
                    % updating parameters
                    set(handles.npcs,'String',num2str(handles.pars.npcs));
                    set(handles.dc,'String',num2str(handles.pars.dc));
                    set(handles.minspk,'String',num2str(handles.pars.minspk));
                    set(handles.minsize,'String',num2str(handles.pars.minsize));
                    set(handles.cent_thr,'String',num2str(handles.pars.cent_thr));
                    set(handles.nsur,'String',num2str(handles.pars.nsur));
                    set(handles.prct,'String',num2str(handles.pars.prct));
                    set(handles.inner_corr,'String',num2str(handles.pars.inner_corr));
                    prev_data = 1;
                    
                else
                    data = matfile.(fnames{1}); % loading the first field
                    if isstruct(data) % if data is not cell array
                        fnames = fieldnames(data); % fields of the structure matfile
                        handles.data.spk = data.(fnames{structfun(@iscell, data)==1}); % spikes data in time points is on a cell array
                    else
                        handles.data.spk = data;
                    end
                    %if size(handles.data.spk,1) == 1
                    %    handles.data.spk = cellfun(@transpose, handles.data.spk, 'uniformoutput', 0);
                    %end
                end
                
            case '.hdf5'% hdf5 file
                handles.data.spk = loadTSfromHdf5(handles.filepath,handles.pars.filename,1);
                
            case '.nex' % nex file
                %handles.data.spk = readNexFile([handles.filepath,handles.pars.filename],0);
                [handles.data.spk] = nex2TS([handles.filepath,handles.pars.filename]);
                    
        end
        if prev_data==0
            handles.data.N = length(handles.data.spk);
            
            % dialog for sampling rate, bin size and interval
            prompt = {'Recording Sampling Rate [Hz]:','bin size [s]:','Interval Start [s]:','Interval End (0 for end of recording)[s]:'};
            def = {num2str(handles.pars.sr),num2str(handles.pars.bin),'0','100'};
            output = inputdlg(prompt,'',1,def);
            
            % assinging variables
            handles.pars.sr = str2double(output{1});
            handles.pars.bin = str2double(output{2});
            handles.pars.ini = str2double(output{3});
            handles.pars.end = str2double(output{4});

            if (handles.data.spk{1}(1))~=floor(handles.data.spk{1}(1)) % checking if timestamps are in sampling points or time basis
                handles.data.spk = cellfun(@(x) x.*handles.pars.sr,handles.data.spk,'uni',0);
            end
            
            if handles.pars.end ==0 % end of recording, assigns to full recording
                handles.pars.end = max(cellfun(@max, handles.data.spk))/handles.pars.sr;
            end
            % using [ini,end] segment of the recording
            handles.data.spk = cellfun(@(x) x(x>=handles.pars.ini*handles.pars.sr & x<=handles.pars.end*handles.pars.sr)...
                - handles.pars.ini*handles.pars.sr,handles.data.spk,'uni',0);
            % computing raster
            %[handles.data.raster] = ts2raster_poprate(handles.data.spk,handles.pars.sr,handles.pars.bin,0)>0; % logical raster
            [handles.data.raster] = ts2binaryraster(handles.data.spk,handles.pars.bin,handles.pars.sr);
            
            
            % update data info
            set(handles.filename_txt,'String', ['File: ',handles.pars.filename])
            set(handles.bin_txt,'String', ['Bin size: ',num2str(handles.pars.bin), ' s']);
            set(handles.sr_txt,'String', ['Sampling Rate: ',num2str(handles.pars.sr), ' Hz'] )
            set(handles.min_x_time,'String',num2str(0));
            set(handles.max_x_time,'String',num2str(handles.pars.end - handles.pars.ini));
            
            handles.data.ens_cols = [];
            handles.data.labels = [];
            handles.data.sel_labels = [];
            % Plot raster
            handles.data.tt = linspace(handles.pars.ini,handles.pars.end,size(handles.data.raster,2)) - handles.pars.ini;
            plot_raster(handles.raster_axes,handles.data.tt,handles.data.raster,handles.data.sel_labels,handles.data.ens_cols);
            
            handles.min_x = handles.data.tt(1);
            handles.max_x = handles.data.tt(end);
            xlim(handles.ens_ras_axes,[handles.min_x handles.max_x])
            xlabel(handles.ens_ras_axes,'Time (s)')
            
        else
            selcols = handles.data.ens_cols(handles.data.final_sel_ens,:);
            sellabs = 1:handles.data.Nens;
            sellabs = sellabs(handles.data.final_sel_ens);
            % ------------------- PLOTTING DATA ------------------------
            % Plotting raster
            plot_raster(handles.raster_axes,handles.data.tt,handles.data.raster,handles.data.sel_labels,selcols);
            if isfield(handles,'min_x')
                xlim(handles.raster_axes,[handles.min_x handles.max_x])
            else
                xlim(handles.raster_axes,[handles.data.tt(1) handles.data.tt(end)])
            end
            
            % Plotting ens_sequence
            plot_ens_seq(handles.ens_ras_axes,handles.data.tt,handles.data.sel_labels,selcols,sellabs);
            if isfield(handles,'min_x')
                xlim(handles.ens_ras_axes,[handles.min_x handles.max_x])
            else
                xlim(handles.ens_ras_axes,[handles.data.tt(1) handles.data.tt(end)])
            end
                        
            % plotting pca
            plot_pca(handles.pca_axes,handles.data.pcs,handles.data.labels,handles.data.ens_cols);
            
            % plotting delta vs rho
            plot_delta_rho(handles.delta_rho_axes,handles.data.rho,handles.data.delta,handles.data.cents,handles.data.predbounds,handles.data.ens_cols);
            
            % plotting corr(n,e)
            plot_core_cells(handles.corr_axes,handles.data.ens_cel_corr,[min(handles.data.ens_cel_corr(:)) max(handles.data.ens_cel_corr(:))]);
            title(handles.corr_axes,'corr(n,e)')
            
            % plotting core-cells
            plot_core_cells(handles.core_axes,handles.data.core_cells,[-1 1]);
            title(handles.core_axes,'Core-Cells')
            
            % plotting inner-ens correlations
            plot_ens_corr(handles.inner_corr_axes,handles.data.ens_corr,handles.data.corr_thr,handles.data.ens_cols);
            
            
        end
        disp('Data properly loaded')
        
        
    %catch
     %   errordlg('Error reading file');
    %end
end


guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over load_butt.
function load_butt_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to load_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ------------- ANALYZING DATA --------------------------------------------
%
%

% --- Executes on button press in analyze_butt.
function analyze_butt_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% gets old inputs parameters
if ~isnan(handles.pars.old_npcs)
    handles.pars.old_npcs = handles.pars.npcs;
    handles.pars.old_dc = handles.pars.dc;
    handles.pars.old_minspk = handles.pars.minspk;
    handles.pars.old_minsize = handles.pars.minsize;
    handles.pars.old_cent_thr = handles.pars.cent_thr;
    handles.pars.old_nsur = handles.pars.nsur;
    handles.pars.old_prct = handles.pars.prct;
    handles.pars.old_inner_corr = handles.pars.inner_corr;
end

% gets new parameters
handles.pars.npcs = str2double(get(handles.npcs,'String'));
handles.pars.dc = str2double(get(handles.dc,'String'));
handles.pars.minspk = str2double(get(handles.minspk,'String'));
handles.pars.minsize = str2double(get(handles.minsize,'String'));
handles.pars.cent_thr = str2double(get(handles.cent_thr,'String'));
handles.pars.nsur = str2double(get(handles.nsur,'String'));
handles.pars.prct = str2double(get(handles.prct,'String'));
handles.pars.inner_corr = str2double(get(handles.inner_corr,'String'));

% comparing parameters to avoid re-calculation of data
pc_op = handles.pars.old_npcs == handles.pars.npcs;
dc_op = handles.pars.old_dc == handles.pars.dc;
minspk_op = handles.pars.old_minspk == handles.pars.minspk;
minsize_op = handles.pars.old_minsize == handles.pars.minsize;
cent_op = handles.pars.old_cent_thr == handles.pars.cent_thr;
nsur_op = handles.pars.old_nsur == handles.pars.nsur;
prct_op = handles.pars.old_prct == handles.pars.prct;
corr_op = handles.pars.old_inner_corr == handles.pars.inner_corr;

% --------- Analyzing data
%if pc_op==0 && dc_op==0 && minspk_op==0 && minsize_op==0 && cent_op==0 && nsur_op==0 && prct_op==0 % compute everything
T = size(handles.data.raster,2);
wb = waitbar(0,'','Name','Analyzing...');
tot_steps = 6;

if minspk_op==0 % if the minimal population spike is change re-calculate everything
    % 1.- Remove low magnitude pop events, compute pca and dist mat
    waitbar(1/tot_steps,wb,'Removing low magnitude events')
    ras = handles.data.raster;
    handles.data.selbins = sum(ras)>handles.pars.minspk;
    ras=ras(:,handles.data.selbins)*1;
    
    % 2.- pca and distance matrix
    waitbar(2/tot_steps,wb,'Computing PCA')
    [~,pcs,~,~,handles.data.exp_var] = pca(ras'); % pca with npcs num components
    handles.data.pcs = pcs(:,1:handles.pars.npcs);
    handles.data.bincor = pdist2(handles.data.pcs,handles.data.pcs); % euclidean distance on principal component space
    
    % 3.- rho and delta computation
    waitbar(3/tot_steps,wb,'Clustering data')
    [~, handles.data.rho] = paraSetv2(handles.data.bincor, handles.pars.dc);
    handles.data.delta = delta_from_dist_mat(handles.data.bincor, handles.data.rho);
    [handles.data.Nens,handles.data.cents,handles.data.predbounds] = cluster_by_pow_fit(handles.data.delta,handles.data.rho,handles.pars.cent_thr);
    if handles.data.Nens==1
        handles.data.labels = ones(length(handles.data.delta),1);
    else
        dist2cent = handles.data.bincor(handles.data.cents>0,:); % distance from centroid to any other point
        [~,handles.data.labels] = min(dist2cent);
    end
    
    % 4.- ensemble raster
    handles.data.ensmat_out = zeros(handles.data.Nens,T);
    handles.data.ensmat_out(:,handles.data.selbins) = bsxfun(@eq,handles.data.labels',(1:handles.data.Nens))';
    
    % 5.- core-cells computation
    waitbar(4/tot_steps,wb,'Core-Cells computation')
    [handles.data.core_cells,~,handles.data.ens_cel_corr,handles.data.sur_cel_cor] = find_core_cells_by_correlation(handles.data.raster,handles.data.ensmat_out,handles.pars.nsur,handles.pars.prct);
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    % ----------------- END OF MINSPK MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==0
    % 2.- pca and distance matrix
    waitbar(2/tot_steps,wb,'Computing PCA')
    ras = handles.data.raster(:,handles.data.selbins)*1;
    [~,handles.data.pcs] = pca(ras','Numcomponents',handles.pars.npcs); % pca with npcs num components     
    handles.data.bincor = pdist2(handles.data.pcs,handles.data.pcs); % euclidean distance on principal component space
    
    % 3.- rho and delta computation
    waitbar(3/tot_steps,wb,'Clustering data')
    [~, handles.data.rho] = paraSetv2(handles.data.bincor, handles.pars.dc);
    handles.data.delta = delta_from_dist_mat(handles.data.bincor, handles.data.rho);
    [handles.data.Nens,handles.data.cents,handles.data.predbounds] = cluster_by_pow_fit(handles.data.delta,handles.data.rho,handles.pars.cent_thr);
    if handles.data.Nens==1
        handles.data.labels = ones(length(handles.data.delta),1);
    else
        dist2cent = handles.data.bincor(handles.data.cents>0,:); % distance from centroid to any other point
        [~,handles.data.labels] = min(dist2cent);
    end
    
    % 4.- ensemble raster
    handles.data.ensmat_out = zeros(handles.data.Nens,T);
    handles.data.ensmat_out(:,handles.data.selbins) = bsxfun(@eq,handles.data.labels',(1:handles.data.Nens))';
    
    % 5.- core-cells computation
    waitbar(4/tot_steps,wb,'Core-Cells computation')
    [handles.data.core_cells,~,handles.data.ens_cel_corr,handles.data.sur_cel_cor] = find_core_cells_by_correlation(handles.data.raster,handles.data.ensmat_out,handles.pars.nsur,handles.pars.prct);
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
    % ----------------- END OF MINSPK  AND PCA MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==0
    % 3.- rho and delta computation
    waitbar(3/tot_steps,wb,'Clustering data')
    [~, handles.data.rho] = paraSetv2(handles.data.bincor, handles.pars.dc);
    handles.data.delta = delta_from_dist_mat(handles.data.bincor, handles.data.rho);
    [handles.data.Nens,handles.data.cents,handles.data.predbounds] = cluster_by_pow_fit(handles.data.delta,handles.data.rho,handles.pars.cent_thr);
    if handles.data.Nens==1
        handles.data.labels = ones(length(handles.data.delta),1);
    else
        dist2cent = handles.data.bincor(handles.data.cents>0,:); % distance from centroid to any other point
        [~,handles.data.labels] = min(dist2cent);
    end
    
    % 4.- ensemble raster
    handles.data.ensmat_out = zeros(handles.data.Nens,T);
    handles.data.ensmat_out(:,handles.data.selbins) = bsxfun(@eq,handles.data.labels',(1:handles.data.Nens))';
    
    % 5.- core-cells computation
    waitbar(4/tot_steps,wb,'Core-Cells computation')
    [handles.data.core_cells,~,handles.data.ens_cel_corr,handles.data.sur_cel_cor] = find_core_cells_by_correlation(handles.data.raster,handles.data.ensmat_out,handles.pars.nsur,handles.pars.prct);
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
    % ----------------- END OF MINSPK, PCA, AND dc MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==1 && cent_op==0
    % 3.- pred bound computation
    [handles.data.Nens,handles.data.cents,handles.data.predbounds] = cluster_by_pow_fit(handles.data.delta,handles.data.rho,handles.pars.cent_thr);
    if handles.data.Nens==1
        handles.data.labels = ones(length(handles.data.delta),1);
    else
        dist2cent = handles.data.bincor(handles.data.cents>0,:); % distance from centroid to any other point
        [~,handles.data.labels] = min(dist2cent);
    end
    
    % 4.- ensemble raster
    handles.data.ensmat_out = zeros(handles.data.Nens,T);
    handles.data.ensmat_out(:,handles.data.selbins) = bsxfun(@eq,handles.data.labels',(1:handles.data.Nens))';
    
    % 5.- core-cells computation
    waitbar(4/tot_steps,wb,'Core-Cells computation')
    [handles.data.core_cells,~,handles.data.ens_cel_corr,handles.data.sur_cel_cor] = find_core_cells_by_correlation(handles.data.raster,handles.data.ensmat_out,handles.pars.nsur,handles.pars.prct);
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
    % ----------------- END OF MINSPK, PCA, dc AND PREDBOUNDS MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==1 && cent_op==1 && nsur_op==0
    % 5.- core-cells computation
    waitbar(4/tot_steps,wb,'Core-Cells computation')
    [handles.data.core_cells,~,handles.data.ens_cel_corr,handles.data.sur_cel_cor] = find_core_cells_by_correlation(handles.data.raster,handles.data.ensmat_out,handles.pars.nsur,handles.pars.prct);
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
    % ----------------- END OF MINSPK, PCA, dc, PREDBOUNDS AND NSUR MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==1 && cent_op==1 && nsur_op==1 && prct_op==0
    % 5.- only changing core-cell threshoild
    idthr = prctile(handles.data.sur_cel_cor,handles.pars.prct,3);
    handles.data.core_cells = handles.data.ens_cel_corr>idthr;
    handles.data.id_sel_core = sum(handles.data.core_cells,1)>handles.pars.minsize;
    
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
    % ----------------- END OF MINSPK, PCA, dc, PREDBOUNDS, NSUR AND PRCT MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==1 && cent_op==1 && nsur_op==1 && prct_op==1 && corr_op==0
    % 6.- filtering core cells
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    [handles.data.ens_corr,handles.data.corr_thr,handles.data.corr_selection] = filter_ens_by_inner_corr(handles.data.raster,handles.data.core_cells,handles.pars.inner_corr);
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    % ----------------- END OF MINSPK, PCA, dc, PREDBOUNDS, NSUR, PRCT AND ENS-CORR MODIFICATION -------------%%%%
    
elseif minspk_op==1 && pc_op==1 && dc_op==1 && cent_op==1 && nsur_op==1 && prct_op==1 && corr_op==1 && minsize_op==0
    % 6.- filtering core cells BY SIZE
    waitbar(5/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.final_sel_ens = handles.data.corr_selection & handles.data.id_sel_core;
    
    % 7.- final ensemble outputs
    waitbar(6/tot_steps,wb,'Filtering ensemble by core-cells')
    handles.data.sel_ensmat_out = handles.data.ensmat_out(handles.data.final_sel_ens,:); % filtering by magnitude &  inner cell correlation
    handles.data.sel_core_cells = handles.data.core_cells(:,handles.data.final_sel_ens);
    handles.data.Nens_final = size(handles.data.sel_ensmat_out,1);
    if handles.data.Nens_final>1
        handles.data.sel_labels = sum(bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)'));
    else
        handles.data.sel_labels = bsxfun(@times,handles.data.sel_ensmat_out,(1:handles.data.Nens_final)');
    end
    
    
end

delete(wb)

% colors
%handles.data.ens_cols = jet(handles.data.Nens);
handles.data.ens_cols = othercolor('Dark28',handles.data.Nens*2);
selcols = handles.data.ens_cols(handles.data.final_sel_ens,:);
sellabs = 1:handles.data.Nens;
sellabs = sellabs(handles.data.final_sel_ens);

% ------------------- PLOTTING DATA ------------------------
% Plotting raster
plot_raster(handles.raster_axes,handles.data.tt,handles.data.raster,handles.data.sel_labels,selcols);
xlim(handles.raster_axes,[handles.min_x handles.max_x])

% Plotting ens_sequence
plot_ens_seq(handles.ens_ras_axes,handles.data.tt,handles.data.sel_labels,selcols,sellabs);
xlim(handles.ens_ras_axes,[handles.min_x handles.max_x])

% plotting eigs
plot_eigs(handles.npcs_axes,handles.data.exp_var,handles.pars.npcs)

% plotting pca
plot_pca(handles.pca_axes,handles.data.pcs,handles.data.labels,handles.data.ens_cols);

% plotting delta vs rho
plot_delta_rho(handles.delta_rho_axes,handles.data.rho,handles.data.delta,handles.data.cents,handles.data.predbounds,handles.data.ens_cols);

% plotting corr(n,e)
plot_core_cells(handles.corr_axes,handles.data.ens_cel_corr,[min(handles.data.ens_cel_corr(:)) max(handles.data.ens_cel_corr(:))]);
title(handles.corr_axes,'corr(n,e)')

% plotting core-cells
plot_core_cells(handles.core_axes,handles.data.core_cells,[-1 1]);
title(handles.core_axes,'Core-Cells')

% plotting inner-ens correlations
plot_ens_corr(handles.inner_corr_axes,handles.data.ens_corr,handles.data.corr_thr,handles.data.ens_cols);

% ------------ UPDATING OLD PARAMETERS
handles.pars.old_npcs = handles.pars.npcs;
handles.pars.old_dc = handles.pars.dc;
handles.pars.old_minspk = handles.pars.minspk;
handles.pars.old_minsize = handles.pars.minsize;
handles.pars.old_cent_thr = handles.pars.cent_thr;
handles.pars.old_nsur = handles.pars.nsur;
handles.pars.old_prct = handles.pars.prct;
handles.pars.old_inner_corr = handles.pars.inner_corr;




guidata(hObject, handles);



function npcs_Callback(hObject, eventdata, handles)
% hObject    handle to npcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npcs as text
%        str2double(get(hObject,'String')) returns contents of npcs as a double


% --- Executes during object creation, after setting all properties.
function npcs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dc_Callback(hObject, eventdata, handles)
% hObject    handle to dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dc as text
%        str2double(get(hObject,'String')) returns contents of dc as a double


% --- Executes during object creation, after setting all properties.
function dc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minspk_Callback(hObject, eventdata, handles)
% hObject    handle to minspk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minspk as text
%        str2double(get(hObject,'String')) returns contents of minspk as a double


% --- Executes during object creation, after setting all properties.
function minspk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minspk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minsize_Callback(hObject, eventdata, handles)
% hObject    handle to minsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minsize as text
%        str2double(get(hObject,'String')) returns contents of minsize as a double


% --- Executes during object creation, after setting all properties.
function minsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cent_thr_Callback(hObject, eventdata, handles)
% hObject    handle to cent_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cent_thr as text
%        str2double(get(hObject,'String')) returns contents of cent_thr as a double


% --- Executes during object creation, after setting all properties.
function cent_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cent_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nsur_Callback(hObject, eventdata, handles)
% hObject    handle to nsur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsur as text
%        str2double(get(hObject,'String')) returns contents of nsur as a double


% --- Executes during object creation, after setting all properties.
function nsur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prct_Callback(hObject, eventdata, handles)
% hObject    handle to prct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prct as text
%        str2double(get(hObject,'String')) returns contents of prct as a double


% --- Executes during object creation, after setting all properties.
function prct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inner_corr_Callback(hObject, eventdata, handles)
% hObject    handle to inner_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inner_corr as text
%        str2double(get(hObject,'String')) returns contents of inner_corr as a double


% --- Executes during object creation, after setting all properties.
function inner_corr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inner_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_butt.
function save_butt_Callback(hObject, eventdata, handles)
% hObject    handle to save_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% saving data
try
    [saveFileName,savePathName] = uiputfile('*.mat','Save data to mat file',handles.filepath);
    ens_data = handles.data;
    ens_data = rmfield(ens_data,'sur_cel_cor');
    ens_data.pars = handles.pars;
    ens_data.filename = handles.pars.filename;
    ens_data.filepath = handles.filepath;
    save([savePathName,saveFileName,'.mat'],'ens_data','-v7.3');
    disp('Data properly saved')
catch
    disp('User cancelled saving')
end



% --- Executes on button press in print_butt.
function print_butt_Callback(hObject, eventdata, handles)
% hObject    handle to print_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes_handles = [handles.pca_axes,handles.delta_rho_axes,handles.raster_axes,...
    handles.ens_ras_axes,handles.corr_axes,handles.core_axes,handles.inner_corr_axes,handles.npcs_axes];
try
    SaveFigures(axes_handles)
catch
    disp('User Cancelled Saving')
end



function min_x_time_Callback(hObject, eventdata, handles)
% hObject    handle to min_x_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_x_time as text
%        str2double(get(hObject,'String')) returns contents of min_x_time as a double


% --- Executes during object creation, after setting all properties.
function min_x_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_x_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_x_time_Callback(hObject, eventdata, handles)
% hObject    handle to max_x_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_x_time as text
%        str2double(get(hObject,'String')) returns contents of max_x_time as a double


% --- Executes during object creation, after setting all properties.
function max_x_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_x_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_xlims.
function apply_xlims_Callback(hObject, eventdata, handles)
% hObject    handle to apply_xlims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.min_x = str2double(get(handles.min_x_time,'String'));
handles.max_x = str2double(get(handles.max_x_time,'String'));

set(handles.ens_ras_axes,'xlim',[handles.min_x handles.max_x]);
set(handles.raster_axes,'xlim',[handles.min_x handles.max_x]);

guidata(hObject, handles);
