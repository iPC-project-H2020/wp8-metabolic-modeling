function [data, info] = addFilenameToData(data, info)
    [~, filename, ext] = fileparts(info.Filename);
    data.Filename(:, 1) = string([filename ext]);
end