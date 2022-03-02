function [model_units, model_w, model_y_hat, model_corr, model_n_train, model_n_valid, model_rmse, model_y_train] = ...
        get_subpop(subpop, model_units, model_w, model_y_hat, model_corr,model_n_train, model_n_valid, model_rmse, model_y_train)
%     function [model_units, model_w, model_y_hat, model_corr, model_n_train, model_n_valid, model_rmse, model_y_train] = ...
%         get_subpop(model_units, model_w, model_y_hat, model_corr,model_n_train, model_n_valid, model_rmse, model_y_train)
%     Gets subpopulation of dLGN cells, specifically units located in the top or the bottom of dLGN, as specified
%     by the input parameter 'subpop'.
%     TODO: avoid large variable list.


    %% Get list of all units (815 dLGN cells)
    % Either load unit_list or generate from previous scripts
    load('../data/dlgn_unit_list.mat');
    units_all = unitList;
    
    %% Get list of dLGN units for which we have information about RF location (divided top or bottom)
    
    % Load mat file
%     load('../data/top_bottom_units.mat');
    load('../data/top_1ch_bottom_units.mat');
    units_RF = myrec;
   
    %% Get subpopulation of units_RF depending on input par 'subpop'
    units_RF_sub = units_RF(strcmp({units_RF.classification}, subpop));

    %% Check whether subpopulation selection has already been done
    if size(model_units,1) == size(units_RF_sub,2)
        fprintf('Subpopulation of dLGN cells already selected: %s\n', subpop);
    else
        fprintf('Getting subpopulation of dLGN cells: %s\n', subpop); 
        
        %% Find units in units_all which are in that subpopulation units_RF_sub

        for i = 1:length(units_RF_sub)
            idx(i) = find([units_all.mouse_counter] == units_RF_sub(i).mouse_counter & ...
                        [units_all.series_num] == units_RF_sub(i).series_num &...
                        [units_all.unit_id] == units_RF_sub(i).unit_id);
        end

        %% Restrict model variables to subpopulation

        model_units = model_units(idx,:,:);
        model_w = model_w(idx,:,:);
        model_y_hat = model_y_hat(idx,:,:);
        model_corr = model_corr(idx,:);
        model_n_train = model_n_train(idx,:);
        model_n_valid = model_n_valid(idx,:);
        model_rmse = model_rmse(idx,:);
        model_y_train = model_y_train(idx,:,:);
        
        fprintf('Subpopulation of dLGN cells (%s): %i cells.\n', subpop, length(idx)); 

    end
end