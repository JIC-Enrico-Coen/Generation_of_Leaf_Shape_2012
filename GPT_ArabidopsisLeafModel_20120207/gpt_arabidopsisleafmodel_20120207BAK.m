function m = gpt_arabidopsisleafmodel_20120207( m )
%m = gpt_arabidopsisleafmodel_20120207( m )
%   Morphogen interaction function.
%   Written at 2012-02-08 13:37:14.
%   GFtbox revision 3898, 2011-12-05 10:25:20.692402.
%   Model last saved to SVN as revision 3912, 2011-12-14 16:56:38.561077.

% The user may edit any part of this function between delimiters
% of the form "USER CODE..." and "END OF USER CODE...".  The
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    try
        m = local_setproperties( m );
    catch
    end

    realtime = m.globalDynamicProps.currenttime;

%%% USER CODE: INITIALISATION
% model combinations explored in the paper: 
% Non-Deforming Model:                                             modelname='NONDEFORMING',     modelversion='DYNAMIC',   modeloperations='NONE';
% Organiser-Based Fixed Model:                                     modelname='PROXORG',              modelversion='FIXED',   modeloperations='NONE';
% Organiser-Based Dynamic Model:                                   modelname='PROXORG',              modelversion='DYNAMIC', modeloperations='NONE';
% Organiser-Based Dynamic Model with transient polarity removal:   modelname='PROXORG',              modelversion='DYNAMIC', modeloperations='TRANSIENTPOLREMOVAL';
% Organiser-Based Model with excision:                             modelname='PROXORG',              modelversion='FIXED',   modeloperations='EXCISION'; 


if Steps(m) == 0
    
    % Model Visualisation and Setup Parameters
    m = leaf_setproperty( m, 'mingradient', 0, 'timestep', 1, 'userpolarisation',false, 'bioAsplitcells', false, 'allowSplitBio',false);
    m = leaf_plotoptions( m, 'hiresdpi', 300,'FEthinlinesize', 1, 'layeroffset', 0.45 );
    m = leaf_fix_vertex( m, 'vertex', find(m.nodes(:,2)<=min(m.nodes(:,2)+0.001)), 'dfs', 'y' );
    
    % Model output times to plot images for the paper.
    m.userdata.outputtimes = [87 117 148 180 205];
    m.userdata.output = 0; % make high resultion images when time = m.userdata.outputtimes; 0 = no output.
    
    % Model types differing in their main polarity system
    m.userdata.ranges.modelname.range{1}='NONDEFORMING';
    m.userdata.ranges.modelname.range{2}='PROXORG';
    m.userdata.ranges.modelname.index =1;

    % Different polarity options for the organiser-based submodels
    m.userdata.ranges.modelversion.range{1} = 'FIXED'; 
    m.userdata.ranges.modelversion.range{2} = 'DYNAMIC';  
    m.userdata.ranges.modelversion.index =2;
    
    % modeloperation - ways to maniputale the canvas or polarity
    m.userdata.ranges.modeloperation.range{1} = 'NONE';
    m.userdata.ranges.modeloperation.range{2} = 'TRANSIENTPOLREMOVAL';
    m.userdata.ranges.modeloperation.range{3} = 'EXCISION';
    m.userdata.ranges.modeloperation.index =1;
    
    % second layer superimposed onto the canvas
    m.userdata.ranges.modeldeco.range{1} = 'NONE';
    m.userdata.ranges.modeldeco.range{2} = 'CELLSEARLY';
    m.userdata.ranges.modeldeco.range{3} = 'CELLSLATE';
    m.userdata.ranges.modeldeco.range{4} = 'CIRCLES';
    m.userdata.ranges.modeldeco.index =1;
    
    
    % parameters and values, with norm value and +/- 20%
    m.userdata.ranges.bpgrad.range = [0.039 0.156 0.195 0.234 0.351]; % for 20% and 80%
    m.userdata.ranges.bpgrad.index = 3;
    m.userdata.ranges.ppgrad.range = [0.0328 0.041 0.0492];
    m.userdata.ranges.ppgrad.index = 2;
    m.userdata.ranges.plam.range = [0.0188 0.0235 0.0282]; 
    m.userdata.ranges.plam.index = 2;
    m.userdata.ranges.glate.range = [0.0038 0.0048 0.0058]; 
    m.userdata.ranges.glate.index = 3;
    m.userdata.ranges.hlate.range = [1.76 2.2 2.64]; 
    m.userdata.ranges.hlate.index = 2;
    m.userdata.ranges.hmid.range = [0.4 0.5 0.6];
    m.userdata.ranges.hmid.index = 2;
    m.userdata.ranges.plate.range = [0.56 0.7 0.84]; 
    m.userdata.ranges.plate.index = 2;
    m.userdata.ranges.bpol.range = [0.08 0.1 0.12];
    m.userdata.ranges.bpol.index = 2;
    
end

modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index};
modelversion = m.userdata.ranges.modelversion.range{m.userdata.ranges.modelversion.index};
modeloperation = m.userdata.ranges.modeloperation.range{m.userdata.ranges.modeloperation.index};
modeldeco = m.userdata.ranges.modeldeco.range{m.userdata.ranges.modeldeco.index};
disp(sprintf('\n------> modelname %s modelversion %s modeloperation %s modeldeco %s \n',...
    modelname,modelversion,modeloperation,modeldeco));

bpol = m.userdata.ranges.bpol.range(m.userdata.ranges.bpol.index);

m.plotdefaults.perelement = [];
m.plotdefaults.perelementaxes = [];
m.plotdefaults.perelementcomponents = [];
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

    if isempty(m), return; end

    setGlobals();
    global gNEW_KA_PAR gNEW_KA_PER gNEW_KB_PAR gNEW_KB_PER
    global gNEW_K_NOR gNEW_POLARISER gNEW_STRAINRET gNEW_ARREST
    dt = m.globalProps.timestep;
    polariser_i = gNEW_POLARISER;
    P = m.morphogens(:,polariser_i);
    [kapar_i,kapar_p,kapar_a,kapar_l] = getMgenLevels( m, 'KAPAR' );
    [kaper_i,kaper_p,kaper_a,kaper_l] = getMgenLevels( m, 'KAPER' );
    [kbpar_i,kbpar_p,kbpar_a,kbpar_l] = getMgenLevels( m, 'KBPAR' );
    [kbper_i,kbper_p,kbper_a,kbper_l] = getMgenLevels( m, 'KBPER' );
    [knor_i,knor_p,knor_a,knor_l] = getMgenLevels( m, 'KNOR' );
    [strainret_i,strainret_p,strainret_a,strainret_l] = getMgenLevels( m, 'STRAINRET' );
    [arrest_i,arrest_p,arrest_a,arrest_l] = getMgenLevels( m, 'ARREST' );
    [v_leaf_i,v_leaf_p,v_leaf_a,v_leaf_l] = getMgenLevels( m, 'V_LEAF' );
    [id_proxorg_i,id_proxorg_p,id_proxorg_a,id_proxorg_l] = getMgenLevels( m, 'ID_PROXORG' );
    [id_inc_i,id_inc_p,id_inc_a,id_inc_l] = getMgenLevels( m, 'ID_INC' );
    [id_pgrad_i,id_pgrad_p,id_pgrad_a,id_pgrad_l] = getMgenLevels( m, 'ID_PGRAD' );
    [id_mid_i,id_mid_p,id_mid_a,id_mid_l] = getMgenLevels( m, 'ID_MID' );
    [id_lam_i,id_lam_p,id_lam_a,id_lam_l] = getMgenLevels( m, 'ID_LAM' );
    [karea_i,karea_p,karea_a,karea_l] = getMgenLevels( m, 'KAREA' );
    [id_late_i,id_late_p,id_late_a,id_late_l] = getMgenLevels( m, 'ID_LATE' );
    [id_distal_i,id_distal_p,id_distal_a,id_distal_l] = getMgenLevels( m, 'ID_DISTAL' );
    [id_cutedge_i,id_cutedge_p,id_cutedge_a,id_cutedge_l] = getMgenLevels( m, 'ID_CUTEDGE' );

% Mesh type: lobes
%            base: 0
%        cylinder: 0
%          height: 0.053
%           lobes: 1
%          radius: 0.042
%      randomness: 0
%           rings: 15
%          strips: 35
%         version: 1

%            Morphogen   Diffusion   Decay   Dilution   Mutant
%            -------------------------------------------------
%                KAPAR        ----    ----       ----     ----
%                KAPER        ----    ----       ----     ----
%                KBPAR        ----    ----       ----     ----
%                KBPER        ----    ----       ----     ----
%                 KNOR        ----    ----       ----     ----
%            POLARISER        ----    ----       ----     ----
%            STRAINRET        ----    ----       ----     ----
%               ARREST        ----    ----       ----     ----
%               V_LEAF        ----    ----       ----     ----
%           ID_PROXORG        ----    ----       ----     ----
%               ID_INC        ----    ----       ----     ----
%             ID_PGRAD        ----    ----       ----     ----
%               ID_MID        ----    ----       ----     ----
%               ID_LAM        ----    ----       ----     ----
%                KAREA        ----    ----       ----     ----
%              ID_LATE        ----    ----       ----     ----
%            ID_DISTAL        ----    ----       ----     ----
%           ID_CUTEDGE        ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

if Steps(m) == 1 % a bit later to give the cluster the chance to catch up.
    
    %% PRN setup
    
    switch modelname
        
        case 'NONDEFORMING'
            % setup polarity parallel to the external y-axis and canvas
            % midline
            P = max(m.nodes(:,2)) - m.nodes(:,2);
            P = P/max(P)*bpol;
            
            % POL levels remain constant throughout the simulation
            m.mgen_production(:,polariser_i) = 0;
            m = leaf_mgen_conductivity( m, 'polariser', 0 );
            m = leaf_mgen_absorption( m, 'polariser', 0 );
            
        case 'PROXORG'
                    % PROXORG and POLARISER (POL) setup 
                    proxorg_ind= m.nodes(:,2)<=min(m.nodes(:,2))+0.001;
                    id_proxorg_p(proxorg_ind) = 1;
                    P(proxorg_ind) = bpol;
                    m = leaf_fix_mgen( m, polariser_i,'vertex',find(proxorg_ind),'fix', 1);
                    m = leaf_mgen_absorption(m, polariser_i, 0.1);
                    
                    switch modelversion
                        % vary POL diffusion rate according to the
                        % model version to establish visible gradients
                        % through the simulation in both model versions.
                        case 'FIXED'
                    m = leaf_mgen_conductivity(m, polariser_i, 0.0001); 
                    
                        case 'DYNAMIC'
                    m = leaf_mgen_conductivity(m, polariser_i, 0.01); 
                    end
        otherwise
            % the model will grow isotropically;
    end

    %% GRN setup   
    
    switch modeloperation            
        case 'EXCISION'
            m.userdata.outputtimes = [148 205];
            
            % setting up DISTAL and CUTEDGE. Regions containing DISTAL will be removed at a later stage. 
            origin = [0 -0.15];
            radius = 0.14; 
            dist = sqrt((m.nodes(:,1)-origin(1)).^2 + (m.nodes(:,2)-origin(2)).^2);
            id_distal_p(dist > radius) = 1;
            id_cutedge_p(dist > (radius-radius/30)) = 1;    
    end
    
    
    %% KRN FACTORS
    
    % Visulalisation Factors setup
    v_leaf_p(:) =1;
    
    % Identity Factors setup
    
    % LAMINA (LAM)
    id_lam_p(m.nodes(:,2)>(-0.038)& m.nodes(:,2)<0.005) = 1.5;
    m = leaf_mgen_conductivity(m, id_lam_i, 0.00082);
    
    % MIDVEIN (MID)
    id_mid_p(:) = (max(m.nodes(:,2)) - m.nodes(:,2));
    id_mid_p = id_mid_p/max(id_mid_p);
    id_mid_p(abs(m.nodes(:,1))>0.0038) = 0; 
    
    % PGRAD
    % set up a linear gradient
    bpgrad = m.userdata.ranges.bpgrad.range(m.userdata.ranges.bpgrad.index);
    id_pgrad_p(:) = max(m.nodes(:,2)) - m.nodes(:,2);
    id_pgrad_p(:) = id_pgrad_p/max(id_pgrad_p)*(1-bpgrad);
    id_pgrad_p(:) = id_pgrad_p+bpgrad;
    
    
elseif realtime > 67 && realtime < 68+ dt
    
    % fixing LAM distribution and resetting LAM to a low value in
    % the petiole
    m = leaf_mgen_conductivity(m, id_lam_i, 0);
    id_lam_p(m.nodes(:,2)<(min(m.nodes(:,2))+0.011)) = 0.4;
    
    
elseif realtime > 85 && realtime < 86+ dt
    % Setup Parameters - fixing polarity directions
    
    switch modelversion
        case 'FIXED'
             m = leaf_setproperty( m, 'mingradient', 2);
             
             m = leaf_mgen_conductivity(m, polariser_i, 0);
             m = leaf_mgen_absorption(m, polariser_i, 0);
        otherwise 
            % the polarity will remain dynamic and directions will adjust
            % with the levels of POL. 
    end
        
    %% Biological Growth
elseif realtime >= 87
    
    % read in the parameters and simplify the parameter names
    ppgrad = m.userdata.ranges.ppgrad.range(m.userdata.ranges.ppgrad.index);
    plam = m.userdata.ranges.plam.range(m.userdata.ranges.plam.index);
    glate = m.userdata.ranges.glate.range(m.userdata.ranges.glate.index);
    plate = m.userdata.ranges.plate.range(m.userdata.ranges.plate.index);
    hlate = m.userdata.ranges.hlate.range(m.userdata.ranges.hlate.index);
    hmid = m.userdata.ranges.hmid.range(m.userdata.ranges.hmid.index);
    
    % model operations
    if realtime > 147 && realtime < 148+ dt
        switch modeloperation
            
            case 'TRANSIENTPOLREMOVAL'                
                ind = find(id_proxorg_l~=1);
                ind2 = randperm(length(ind));
                P(ind(ind2)) = P(ind)/10; 
                
            case 'EXCISION'
                id_inc_p(id_cutedge_l==1)=1;
        end
    end
    
    %% GRN
    % increase of LATE levels over time.
    if realtime >= 148 
        id_late_p = (realtime - 148)*glate;
    end
    
    
    %% KRN
        % growth interaction network for all models
            kapar_p(:) = ppgrad.*id_pgrad_l... Kpar promotion by PGRAD (parallel to POL gradient)
                        .*inh(hlate,id_late_l)... Kpar inhibition by LATE
                        .*inh(5,id_inc_l); % in case of distal excision, growth restriction at the cut.
            kbpar_p(:) = kapar_p;
            kaper_p(:) = plam*id_lam_l...Kper promotion by LAM (perpendicular to POL gradient)
                        .*pro(plate,id_late_l)...Kper promotion by LATE
                        .*inh(hmid, id_mid_l)...Kper inhibition along the midline by MID
                        .*inh(5,id_inc_l); % in case of distal excision, growth restriction at the cut.
            kbper_p(:) = kaper_p(:);

end

%% PRN continuous
    % in the non-deforming model, the POL gradient is reset after every iteration.
switch modelname
    case 'NONDEFORMING'
        P = max(m.nodes(:,2)) - m.nodes(:,2);
        P = P/max(P)*bpol;
    otherwise
        % in all other cases polarity interaction have been established during setup.
end

%% Model Decorations - optional

switch modeldeco
    case 'CIRCLES'  % make simple circles. 
         if Steps(m) == 1 
            m = leaf_makesecondlayer( m, ... 
                'mode', 'each', ...  
                'relarea', 1/80,...
                'probpervx', 'V_LEAF', ...
                'numcells',50,...
                'sides', 20, ...  
                'colors', [0.8 0.8 0.8],...
                'colorvariation', 0,...
                'add', false,...
                'allowoverlap', false);  
        end
        
    case 'CELLSEARLY' % load in cell outlines taken from a 3-day old leaf.
        if Steps(m) == 1 

            % rotation and rescaling parameters to project cells onto the
            % canvas
            ang = pi+0.1; scale = 1000;
            
            % load in the position of each vertex and rescale (from meters
            % to mm)
            path = [fileparts(which(m.globalProps.modelname)),filesep,'EarlyClones'];
            point = load([path,filesep,'allpoints.mat']);
            points = point.all_points(:,:,1)*scale;
            pointsC1 = points(:,1)-mean(points(:,1))+0.0023;
            pointsC2 = points(:,2)-mean(points(:,2))+0.0045;
            
            % rotate points
            R = [cos(ang) sin(ang); -sin(ang) cos(ang)];
            pointsC = [pointsC1 pointsC2]*R;
            points(:,3) = zeros(length(points),1); % to get z-axis component
            pointsC = [pointsC points(:,3)];
            
            % load in the cell vertices
            cells = load([path,filesep,'cells.mat']);
            cells = cells.cells;
            
            l = 0;
            cellcount = [86 76 83 73 68 72 0 89 62 60 56 53 5 50 13 19 22 91 96 95 24 28 94 82 35 41 8 4 122 115 109 125 116]; % from 122
            for i=cellcount+1
                l = l+1;
                cell{l} = cells(i).pts;
            end
            
            % position the cells in the canvas
            m = leaf_makesecondlayer( m, ...
                'vertexdata',pointsC,...
                'celldata',cell);
            
            % recolour superimposed cells
            cols = load([path,filesep,'1946Ycolours']);
            cols = cols.colours;
            for j=1:(size(cellcount,2)-size(cols,1))
            cols(end+1,:) = cols(j,:);
            end
            cols = cols./255;
            cols = cols(1:length(cellcount),:);
            m.secondlayer.cellcolor = cols;
            
        end
        
    case 'CELLSLATE' % load in cell outlines taken from a 6-day old leaf.
        if realtime > 156 && realtime < 157+ dt
            
            C2 = -0.31;
            ang = pi; scale = 1e3;
            
            path = [fileparts(which(m.globalProps.modelname)),filesep,'LateClones'];
            point = load([path,filesep,'allpoints.mat']);
            points = point.all_points(:,:,1)*scale;
            
            % reposition
            points(:,3) = zeros(length(points),1);
            pointsC1 = points(:,1)-mean(points(:,1));
            pointsC2 = points(:,2)-mean(points(:,2))+C2;
            
            % rotate
            R = [cos(ang) sin(ang); -sin(ang) cos(ang)];
            pointsC = [pointsC1 pointsC2]*R;
            pointsC = [pointsC points(:,3)];
            
            % cells
            cells = load([path,filesep,'cells.mat']);
            cells = cells.cells;
            
            for i=1:size(cells,2)
                cell{i} = cells(i).pts;
            end
            
            m = leaf_makesecondlayer( m, ...
                'vertexdata',pointsC,...
                'celldata',cell);
            
            % recolour superimposed cells
            cols = load([path,filesep,'1946Ycolours']);
            cols = cols.colours./255;
            m.secondlayer.cellcolor = cols;
            
        end
    otherwise 
        % the model will not have a second layer. 
end

%% Growth Rate Calculations

% Specified karea
karea_p(:) = kapar_p+kaper_p;


% m = calculateOutputs(m);
% [acA,frames] = tensorsToComponents( m.outputs.actualstrain.A, m.cellFrames );
% kmax = perFEtoperVertex( m, max(acA(:,[1 2]),[],2));
% kmin = perFEtoperVertex( m, min(acA(:,[1 2]),[],2));
% %aniso_p(:) = abs(kmax-kmin)./max(kmax,kmin);
% karea_p(:) = kmin+kmax;

% [inplane,outofplane] = splitVector( m.outputs.rotations,m.unitcellnormals );
% rot_p(:) = perFEtoperVertex( m, inplane);


output = m.userdata.outputtimes(2:end);

% calculate growth over 24 hour intervals
if sum(ismember(output-24,realtime))
    m.userdata.oldpos = m.prismnodes;
    m.userdata.starttime = realtime;
    
elseif sum(ismember(output-1,realtime))
    displacements = m.prismnodes - m.userdata.oldpos;
    [growth,gf] = leaf_computeGrowthFromDisplacements( m, displacements, ...
        realtime - m.userdata.starttime,'axisorder', 'maxminnor', ...
        'anisotropythreshold', 0.05);
    
     % plot resultant areal growth rates over 24-h intervals.    
     %m = leaf_plotoptions( m, 'pervertex',perFEtoperVertex(m,sum(growth(:,1:2),2)) ,'perelementaxes', gf(:,1,:), 'drawtensoraxes', true );
  
    
elseif sum(ismember(output,realtime)) && m.userdata.output ==1
    % Output - plot an image at high resolution   
    path = fileparts(which(m.globalProps.modelname));
    [m,ok] = leaf_snapshot( m,[path,filesep,'snapshots',filesep,modelname,'_',modeltype,'.png'], 'resolution',[]); %

end
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.
    m.morphogens(:,polariser_i) = P;
    m.morphogens(:,kapar_i) = kapar_p;
    m.morphogens(:,kaper_i) = kaper_p;
    m.morphogens(:,kbpar_i) = kbpar_p;
    m.morphogens(:,kbper_i) = kbper_p;
    m.morphogens(:,knor_i) = knor_p;
    m.morphogens(:,strainret_i) = strainret_p;
    m.morphogens(:,arrest_i) = arrest_p;
    m.morphogens(:,v_leaf_i) = v_leaf_p;
    m.morphogens(:,id_proxorg_i) = id_proxorg_p;
    m.morphogens(:,id_inc_i) = id_inc_p;
    m.morphogens(:,id_pgrad_i) = id_pgrad_p;
    m.morphogens(:,id_mid_i) = id_mid_p;
    m.morphogens(:,id_lam_i) = id_lam_p;
    m.morphogens(:,karea_i) = karea_p;
    m.morphogens(:,id_late_i) = id_late_p;
    m.morphogens(:,id_distal_i) = id_distal_p;
    m.morphogens(:,id_cutedge_i) = id_cutedge_p;

%%% USER CODE: FINALISATION
switch modeloperation
    case 'EXCISION' % take of part of the canvas (expressing DISTAL)
        if realtime > 147 && realtime < 148+ dt % cut
            m = leaf_deletenodes( m, id_distal_l==1);
        end
end

% In this section you may modify the mesh in any way whatsoever.
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS


function m = initproperties( m )
% This function is called at time zero in the INITIALISATION section of the
% interaction function.  It provides commands to set each of the properties
% that are contained in m.globalProps.  Uncomment whichever ones you would
% like to set yourself, and put in whatever value you want.
%
% Some of these properties are for internal use only and should never be
% set by the user.  At some point these will be moved into a different
% component of m, but for the present, just don't change anything unless
% you know what it is you're changing.

%    m = leaf_setproperty( m, 'trinodesvalid', true );
%    m = leaf_setproperty( m, 'prismnodesvalid', true );
%    m = leaf_setproperty( m, 'thicknessRelative', 0.010700 );
%    m = leaf_setproperty( m, 'thicknessArea', 1.000000 );
%    m = leaf_setproperty( m, 'activeGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedGrowth', 17.000000 );
%    m = leaf_setproperty( m, 'displayedMulti', [] );
%    m = leaf_setproperty( m, 'allowNegativeGrowth', true );
%    m = leaf_setproperty( m, 'usePrevDispAsEstimate', true );
%    m = leaf_setproperty( m, 'mingradient', 1.000000 );
%    m = leaf_setproperty( m, 'relativepolgrad', false );
%    m = leaf_setproperty( m, 'userpolarisation', false );
%    m = leaf_setproperty( m, 'thresholdsq', 0.000061 );
%    m = leaf_setproperty( m, 'splitmargin', 1.400000 );
%    m = leaf_setproperty( m, 'splitmorphogen', '' );
%    m = leaf_setproperty( m, 'thresholdmgen', 0.500000 );
%    m = leaf_setproperty( m, 'bulkmodulus', 1.000000 );
%    m = leaf_setproperty( m, 'unitbulkmodulus', true );
%    m = leaf_setproperty( m, 'poissonsRatio', 0.300000 );
%    m = leaf_setproperty( m, 'starttime', 68.000000 );
%    m = leaf_setproperty( m, 'timestep', 1.000000 );
%    m = leaf_setproperty( m, 'timeunitname', 'hours' );
%    m = leaf_setproperty( m, 'distunitname', 'mm' );
%    m = leaf_setproperty( m, 'scalebarvalue', 0.000000 );
%    m = leaf_setproperty( m, 'validateMesh', true );
%    m = leaf_setproperty( m, 'rectifyverticals', false );
%    m = leaf_setproperty( m, 'allowSplitLongFEM', false );
%    m = leaf_setproperty( m, 'longSplitThresholdPower', 0.000000 );
%    m = leaf_setproperty( m, 'allowSplitBentFEM', false );
%    m = leaf_setproperty( m, 'allowSplitBio', false );
%    m = leaf_setproperty( m, 'allowFlipEdges', false );
%    m = leaf_setproperty( m, 'allowElideEdges', false );
%    m = leaf_setproperty( m, 'mincellangle', 0.200000 );
%    m = leaf_setproperty( m, 'alwaysFlat', 1.000000 );
%    m = leaf_setproperty( m, 'flattenforceconvex', true );
%    m = leaf_setproperty( m, 'flatten', false );
%    m = leaf_setproperty( m, 'flattenratio', 1.000000 );
%    m = leaf_setproperty( m, 'useGrowthTensors', false );
%    m = leaf_setproperty( m, 'plasticGrowth', false );
%    m = leaf_setproperty( m, 'totalinternalrotation', 0.000000 );
%    m = leaf_setproperty( m, 'stepinternalrotation', 2.000000 );
%    m = leaf_setproperty( m, 'showinternalrotation', false );
%    m = leaf_setproperty( m, 'performinternalrotation', false );
%    m = leaf_setproperty( m, 'internallyrotated', false );
%    m = leaf_setproperty( m, 'maxFEcells', 0.000000 );
%    m = leaf_setproperty( m, 'inittotalcells', 0.000000 );
%    m = leaf_setproperty( m, 'maxBioAcells', 0.000000 );
%    m = leaf_setproperty( m, 'maxBioBcells', 0.000000 );
%    m = leaf_setproperty( m, 'colors', (3 values) );
%    m = leaf_setproperty( m, 'colorvariation', 0.000000 );
%    m = leaf_setproperty( m, 'colorparams', (6 values) );
%    m = leaf_setproperty( m, 'freezing', 0.000000 );
%    m = leaf_setproperty( m, 'canceldrift', 0.000000 );
%    m = leaf_setproperty( m, 'mgen_interaction', (unknown type ''function_handle'') );
%    m = leaf_setproperty( m, 'mgen_interactionName', 'arleaf_ver3187_ad0p4goodunfixrim_101021' );
%    m = leaf_setproperty( m, 'allowInteraction', true );
%    m = leaf_setproperty( m, 'interactionValid', true );
%    m = leaf_setproperty( m, 'gaussInfo', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'stitchDFs', [] );
%    m = leaf_setproperty( m, 'D', (36 values) );
%    m = leaf_setproperty( m, 'C', (36 values) );
%    m = leaf_setproperty( m, 'G', (6 values) );
%    m = leaf_setproperty( m, 'solver', 'cgs' );
%    m = leaf_setproperty( m, 'solvertolerance', 0.000100 );
%    m = leaf_setproperty( m, 'diffusiontolerance', 0.000001 );
%    m = leaf_setproperty( m, 'allowsparse', true );
%    m = leaf_setproperty( m, 'maxIters', 166.000000 );
%    m = leaf_setproperty( m, 'maxsolvetime', 1000.000000 );
%    m = leaf_setproperty( m, 'cgiters', 232.000000 );
%    m = leaf_setproperty( m, 'simsteps', 0.000000 );
%    m = leaf_setproperty( m, 'stepsperrender', 0.000000 );
%    m = leaf_setproperty( m, 'growthEnabled', true );
%    m = leaf_setproperty( m, 'diffusionEnabled', true );
%    m = leaf_setproperty( m, 'makemovie', 0.000000 );
%    m = leaf_setproperty( m, 'moviefile', '' );
%    m = leaf_setproperty( m, 'codec', 'None' );
%    m = leaf_setproperty( m, 'autonamemovie', true );
%    m = leaf_setproperty( m, 'overwritemovie', false );
%    m = leaf_setproperty( m, 'framesize', [] );
%    m = leaf_setproperty( m, 'mov', [] );
%    m = leaf_setproperty( m, 'jiggleProportion', 1.000000 );
%    m = leaf_setproperty( m, 'cvtperiter', 0.200000 );
%    m = leaf_setproperty( m, 'boingNeeded', false );
%    m = leaf_setproperty( m, 'initialArea', 0.008228 );
%    m = leaf_setproperty( m, 'bendunitlength', 0.090709 );
%    m = leaf_setproperty( m, 'targetRelArea', 1.000000 );
%    m = leaf_setproperty( m, 'defaultinterp', 'min' );
%    m = leaf_setproperty( m, 'readonly', false );
%    m = leaf_setproperty( m, 'projectdir', 'D:\DArT_Models\trunk\Growth Arrest\PaperWorthy\UnfixedBase' );
%    m = leaf_setproperty( m, 'modelname', 'ArLeaf_ver3187_Ad0p4GoodUnfixRim_101021' );
%    m = leaf_setproperty( m, 'allowsave', 0.000000 );
%    m = leaf_setproperty( m, 'addedToPath', true );
%    m = leaf_setproperty( m, 'bendsplit', 0.300000 );
%    m = leaf_setproperty( m, 'usepolfreezebc', false );
%    m = leaf_setproperty( m, 'dorsaltop', true );
%    m = leaf_setproperty( m, 'defaultazimuth', -45.000000 );
%    m = leaf_setproperty( m, 'defaultelevation', 33.750000 );
%    m = leaf_setproperty( m, 'defaultroll', 0.000000 );
%    m = leaf_setproperty( m, 'defaultViewParams', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'comment', '' );
%    m = leaf_setproperty( m, 'legendTemplate', '%T: %q\n%m' );
%    m = leaf_setproperty( m, 'bioAsplitcells', false );
%    m = leaf_setproperty( m, 'bioApullin', 0.142857 );
%    m = leaf_setproperty( m, 'bioAfakepull', 0.202073 );
%    m = leaf_setproperty( m, 'interactive', false );
%    m = leaf_setproperty( m, 'coderevision', 3187 );
%    m = leaf_setproperty( m, 'coderevisiondate', '2010-10-15 08:03:44.157047' );
%    m = leaf_setproperty( m, 'modelrevision', 3172.000000 );
%    m = leaf_setproperty( m, 'modelrevisiondate', '2010-10-12 08:38:13.898165' );
%    m = leaf_setproperty( m, 'vxgrad', (108 values) );
%    m = leaf_setproperty( m, 'lengthscale', 0.107000 );
%    m = leaf_setproperty( m, 'targetAbsArea', 0.007388 );
%    m = leaf_setproperty( m, 'RecordMeshes', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'usefrozengradient', true );
%    m = leaf_setproperty( m, 'solverprecision', 'double' );
%    m = leaf_setproperty( m, 'thicknessMode', 'physical' );
%    m = leaf_setproperty( m, 'savedrunname', '' );
%    m = leaf_setproperty( m, 'savedrundesc', '' );
%    m = leaf_setproperty( m, 'solvertolerancemethod', 'norm' );
%    m = leaf_setproperty( m, 'physicalThickness', 0.000000 );
end





% Here you may add any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.
% This comment can be deleted.of sink