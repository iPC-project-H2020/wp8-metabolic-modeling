function data_dir = constructDirectory(jobInfo)
    data_dir = '';
    if jobInfo.("is.cell_line")=="TRUE"
        data_dir = strcat(data_dir, '\cell_line\', jobInfo.model_name, '\', ...
            jobInfo.input_type_clustering);
        if strcmp(jobInfo.input_type_clustering, 'ic10')
            data_dir = strcat(data_dir, '\','ic_score_by_',...
                jobInfo.ic_score_by, '\');
        end
        data_dir = strcat(data_dir, '\cl_aggregation_by_', jobInfo.aggregation_by, '\');
    else
        data_dir = strcat(data_dir, '\pdx\', jobInfo.model_name, '\', ...
            'cl_resolution_', num2str(cell2mat(jobInfo.clustering_resolution)), '\', ...
            'ic_score_by_', jobInfo.ic_score_by, '\', ...
            'cl_aggregation_by_', jobInfo.aggregation_by, '\');
    end 
end