
function browse_product()

    % Usage:
    
        % browse_product()

    % Nathan Weinstein, Weizmann Institute of Science, 2023
    %
    % Feel free to reach me out:
    % nathan.wainstein@weizmann.ac.il
    % https://linkedin.com/in/nathan-weinstein


    summary_params_names = select_bands();
    
    function summary_params_names = select_bands()

        list = {'SINDEX2', 'D2300', 'BD1900_2'};    % more in future...
        RGB = {'Red', 'Green', 'Blue'};
        summary_params_names = cell(size(RGB));

        for i = 1:length(RGB)
            [idx, tf] = listdlg('ListString', list,...
                'PromptString', [RGB(i), ' channel'],...
                'SelectionMode', 'single',...
                'InitialValue', 1,...
                'Name', 'Select RGB summary params',...
                'ListSize', [160, 300],...
                'OKString', 'Apply');

            if ~tf; disp('Canceled'); return; end

            summary_params_names(i) = list(idx);    % selected summary params

            % remove the selected from the list:
            list(1) = [];
        end

    end

    assert(length(summary_params_names) == 1 || length(summary_params_names) == 3, 'Invalid number of bands. Must be 1 or 3.')

    option_1 = 'Load SU/SR';
    option_2 = 'Calculate on IF cube';
    option_3 = 'Calculate Lab Reference';
    
    opt.Default = option_2;
    opt.Interpreter = 'tex';
    answer = questdlg('\fontname{Calibri} \fontsize{15} Load pre-calculated SU/SR summary parameters, or calculate them?', '', option_1, option_2, option_3, opt);

    if strcmpi(answer, option_1)

        [filename, path] = uigetfile({'*su*.hdr;*sr*.hdr'}, 'Select SR or SU Scene');
        filepath         = fullfile(path, filename);

        header  = enviinfo(filepath);
        nBands  = header.Bands;
        bands   = 1:nBands;
        hcube   = hypercube(filepath, bands);

        summary_params = get_bands(hcube, summary_params_names);

    elseif strcmpi(answer, option_2)

        hcube = load_scene();

        metadata = hcube.Metadata;

        answer = questdlg('Exact Reflectance or Median?', '', 'Exact', 'Median', 'Exact');

        if strcmpi(answer, 'Exact')
            cube_datatype = 'SU';
        elseif strcmpi(answer, 'Median')
            cube_datatype = 'SR';
        end

        metadata.cube_datatype = cube_datatype;         

        hcube = hypercube(hcube.DataCube, hcube.Wavelength, metadata);

        summary_params = calc_summary_parameters(summary_params_names, hcube);

    elseif strcmpi(answer, option_3)

        speclib = load_speclib();
        mineral = name_enter(speclib);

        wavelength  = column_data(speclib, 'Wavelength');
        reflectance = column_data(speclib, mineral);

        wavelength = nm2um_cube(wavelength);

        libData.Wavelength  = wavelength;
        libData.Reflectance = reflectance;

        summary_params = calc_summary_parameters(summary_params_names, libData);

        img = concatenate_bands(summary_params);

        figure();
        image(img);
        axis image

        return

    end

    summary_params = ignore_values(summary_params, hcube.Metadata.DataIgnoreValue);
    summary_params = bands2double(summary_params);
                
    stretch_params = input_stretch(summary_params_names);

    function stretch_params = input_stretch(summary_params_names)

        prompt = {[summary_params_names{1}, ' (Red; default)'],...
                    [summary_params_names{2}, ' (Green; default)'],...
                    [summary_params_names{3}, ' (Blue; default)']};

        definput = {'0.003, 0.020', '0.010, 0.038', '0.028, 0.058'};

        opts.Interpreter = 'none';
        opts.Resize = 'on';

        answer = inputdlg(prompt, 'Select strecth params for the RGB channels', [3, 50], definput, opts);

        if(isempty(answer) || isempty(answer{1})); disp('Canceled');   return; end

        stretch_params = cellfun(@(x) str2double(x), cellfun(@(x) strsplit(x, ', '), answer, 'UniformOutput', false), 'UniformOutput', false);
        stretch_params = stretch_params';
        
    end

    summary_params = stretch(summary_params, stretch_params);

    img = concatenate_bands(summary_params);

    figure();
    image(img);
    axis image

end


function summary_params = calc_summary_parameters(summary_params_names, hcube)

    summary_params = cellfun(@(x) calc_summary_parameter(x, hcube), summary_params_names, 'UniformOutput', false);

end

function summary_param = calc_summary_parameter(summary_params_names, hcube)

    if ~exist('hcube', 'var')
        hcube = load_scene();
    end

    if isobject(hcube)
        hcube = nm2um_cube(hcube);
    end

    switch summary_params_names
        
        case 'SINDEX2'
            summary_param = sindex2(hcube);

        case 'D2300'
            summary_param = bd2300(hcube);

        case 'BD1900_2'
            summary_param = bd1900_2(hcube);

    end


    function summary_param = sindex2(hcube)

        % Inverse Band Depth Parameter

        wl_S = 2.120;
        wl_C = 2.290;
        wl_L = 2.400;

        summary_param = band_depth(hcube, wl_S, wl_C, wl_L, 'Inverse');


        % checks:
%                 squeeze(hcube.DataCube(2, 32, gama_S_i-n:gama_S_i+n))
%                 squeeze(hcube.DataCube(2, 32, gama_S_i))
%                 hcube.Wavelength(gama_S_i)
%                 hcube.Wavelength(gama_C_i)
%                 hcube.Wavelength(gama_L_i-n:gama_L_i+n)
%                 check = hcube.Wavelength(index)

    end

    function summary_param = bd2300(hcube)

        W1815 = 1.815;
        W2120 = 2.120;
        W2170 = 2.170;
        W2210 = 2.210;
        W2290 = 2.290;
        W2320 = 2.320;
        W2330 = 2.330;
        W2530 = 2.530;

        % retrieve the CRISM wavelengths nearest the requested
        % values:
        [~, W1815_i] = get_wavelength(hcube.Wavelength, W1815);
        [~, W2120_i] = get_wavelength(hcube.Wavelength, W2120);
        [~, W2170_i] = get_wavelength(hcube.Wavelength, W2170);
        [~, W2210_i] = get_wavelength(hcube.Wavelength, W2210);
        [~, W2290_i] = get_wavelength(hcube.Wavelength, W2290);
        [~, W2320_i] = get_wavelength(hcube.Wavelength, W2320);
        [~, W2330_i] = get_wavelength(hcube.Wavelength, W2330);
        [~, W2530_i] = get_wavelength(hcube.Wavelength, W2530);

        % extract individual channels:
        if isobject(hcube)

            if strcmpi(hcube.Metadata.cube_datatype, 'SU')

                % for producing SU cube:
                R1815 = hcube.DataCube(:, :, W1815_i);
                R2120 = hcube.DataCube(:, :, W2120_i);
                R2170 = hcube.DataCube(:, :, W2170_i);
                R2210 = hcube.DataCube(:, :, W2210_i);
                R2290 = hcube.DataCube(:, :, W2290_i);
                R2320 = hcube.DataCube(:, :, W2320_i);
                R2330 = hcube.DataCube(:, :, W2330_i);
                R2530 = hcube.DataCube(:, :, W2530_i);

            elseif strcmpi(hcube.Metadata.cube_datatype, 'SR')

                % for SR cube:
                R1815 = median_reflectance(hcube, W1815_i, 5);
                R2120 = median_reflectance(hcube, W2120_i, 5);
                R2170 = median_reflectance(hcube, W2170_i, 5);
                R2210 = median_reflectance(hcube, W2210_i, 5);
                R2290 = median_reflectance(hcube, W2290_i, 3);
                R2320 = median_reflectance(hcube, W2320_i, 3);
                R2330 = median_reflectance(hcube, W2330_i, 3);
                R2530 = median_reflectance(hcube, W2530_i, 5);

            end

        else

            % for lab reference:
            R1815 = hcube.Reflectance(W1815_i);
            R2120 = hcube.Reflectance(W2120_i);
            R2170 = hcube.Reflectance(W2170_i);
            R2210 = hcube.Reflectance(W2210_i);
            R2290 = hcube.Reflectance(W2290_i);
            R2320 = hcube.Reflectance(W2320_i);
            R2330 = hcube.Reflectance(W2330_i);
            R2530 = hcube.Reflectance(W2530_i);

        end

        % compute the interpolated continuum values at selected wavelengths between 1815 and 2530
        slope  = ( R2530 - R1815 ) ./ ( W2530 - W1815 );
        CR2120 = R1815 + slope .* ( W2120 - W1815 );
        CR2170 = R1815 + slope .* ( W2170 - W1815 );
        CR2210 = R1815 + slope .* ( W2210 - W1815 );
        CR2290 = R1815 + slope .* ( W2290 - W1815 );
        CR2320 = R1815 + slope .* ( W2320 - W1815 );
        CR2330 = R1815 + slope .* ( W2330 - W1815 );

        summary_param = 1 - (((R2290 ./ CR2290) + (R2320 ./ CR2320) + (R2330 ./ CR2330)) ./ ...
                            ((R2120 ./ CR2120) + (R2170 ./ CR2170) + (R2210 ./ CR2210)));

    end
    
    function summary_param = bd1900_2(hcube)

        % Band Depth Parameter

        wl_S = 1.850;
        wl_C = 1.930;
        wl_L = 2.067;

        A = band_depth(hcube, wl_S, wl_C, wl_L);

        wl_S = 1.850;
        wl_C = 1.985;
        wl_L = 2.067;

        B = band_depth(hcube, wl_S, wl_C, wl_L);

        summary_param = 0.5 * A + 0.5 * B;

    end

    function BD = band_depth(hcube, wl_S, wl_C, wl_L, inverse)

        [~, wl_S_i] = get_wavelength(hcube.Wavelength, wl_S);
        [~, wl_C_i] = get_wavelength(hcube.Wavelength, wl_C);
        [~, wl_L_i] = get_wavelength(hcube.Wavelength, wl_L);

        if isobject(hcube)

            if strcmpi(hcube.Metadata.cube_datatype, 'SU')
                % for producing SU cube
                R_S = hcube.DataCube(:, :, wl_S_i);
                R_C = hcube.DataCube(:, :, wl_C_i);
                R_L = hcube.DataCube(:, :, wl_L_i);

            elseif strcmpi(hcube.Metadata.cube_datatype, 'SR')

                % for SR cube:
                R_S = median_reflectance(hcube, wl_S_i, 5);
                R_C = median_reflectance(hcube, wl_C_i, 7);
                R_L = median_reflectance(hcube, wl_L_i, 3);
            end

        else

            % for lab reference:
            R_S = hcube.Reflectance(wl_S_i);
            R_C = hcube.Reflectance(wl_C_i);
            R_L = hcube.Reflectance(wl_L_i);

        end

        % BD:
        b = (wl_C - wl_S) / (wl_L - wl_S);
        a = 1 - b;

        if exist('inverse', 'var') && strcmpi(inverse, 'Inverse')
            BD = 1 - (((a .* R_S) + (b .* R_L)) ./ R_C);
        else
            BD = 1 - (R_C ./ ( a * R_S + b * R_L));
        end
    end
end


function [val, idx] = get_wavelength(wavelengths, wl)
    if isobject(wavelengths)
        wavelengths = wavelengths.Wavelength;
    end
    [val, idx] = min(abs(wavelengths - wl));
end

function R = median_reflectance(hcube, R_idx, kernel_width)

%             minidx = (R_indx - half_width) > 0  % clamp to 0 if R-half_width is less than 0
%             maxidx = (R_indx + half_width) < (n_elements(wvt)-1)  % clamp to max index if R+half_width is past the end

    half_width  = fix(kernel_width / 2);
    minidx      = (R_idx - half_width);
    maxidx      = (R_idx + half_width);

    if isobject(hcube)
        subcube = hcube.DataCube(:, :, minidx:maxidx);
    else
        subcube = hcube.Reflectance(minidx:maxidx);
    end

% %             boxcar      = ones(1, 1, kernel_width)/kernel_width;
% %             R           = conv(subcube, boxcar, 'same');
% 
% %             R           = smoothdata(subcube, 3, 'movmean', kernel_width, 'omitnan');
% %             R           = R(:, :, ceil(size(R,3)/2));

    R           = median(subcube, 3);

%             R = subcube;

end

function bands = get_bands(hcube, bands_to_load)

%             bands = cell(1, length(bands_to_load));
%             for i = 1:length(bands_to_load)
%                 bands{i} = get_band(hcube, bands_to_load{i});
%             end
    bands = cellfun(@(x) get_band(hcube, x), bands_to_load, 'UniformOutput', false);

    
end

function band = get_band(hcube, band_name)
% get_DDR_slice
    %     textdata    = hcube.Metadata.BandNames;
    %     expression  = 'Latitude';

    % todo: generalize to LBL relevant filed with tab delimited to use on X =
    % multibandread()

%             a       = regexpi(hcube.Metadata.BandNames, band_name);
%             line    = find(~ cellfun(@isempty, a));
%             band_i  = line;

    band_i  = find(strcmp(hcube.Metadata.BandNames, band_name));
    band    = hcube.DataCube(:, :, band_i);

end


function A = ignore_values(A, ignore_value, whattobe)

    if ~exist("whattobe", 'var')
        whattobe = NaN;
    end
    A = cellfun(@(x) ignore_values_helper(x, ignore_value, whattobe), A, 'UniformOutput', false);

end

function A = ignore_values_helper(A, a, b)
    A(A == a) = b;
end

function bands = bands2double(bands)
    bands = cellfun(@im2double, bands, 'UniformOutput', false);
end

function bands = stretch(bands, stretch_params)
%             bands = arrayfun(@(i) imadjust(bands{i}, stretch_params{i}), 1:length(bands), 'UniformOutput', false);
    bands = cellfun(@(band, stretch) imadjust(band, stretch), bands, stretch_params, 'UniformOutput', false);

end

function img = concatenate_bands(bands)

%             img = cat(3, bands{1}, bands{2}, bands{3});
    img = cat(3, bands{:});

    % Or generaly:
%             function img = concatenate_bands(bands)
%                 img = bands{1};
%                 for i = 2:length(bands)
%                     img = cat(3, img, bands{i});
%                 end
%             end

end


function speclib = load_speclib(option)
    % load Spectral Library (should be resampled before to the hcube wavelengths):

    arguments

        option.speclib_path (1, :) char = ''

    end

    if ~ strcmpi(option.speclib_path, "")

        speclib         = importdata(option.speclib_path);
    else
        speclib_folder   = 'Spectral Libraries';
        filename         = 'Sulfates and Sulfite';
        [filename, path] = uigetfile('*.txt', 'Select Input Spectral Library - ASCII', [speclib_folder, '\', filename]);     if isequal(path, 0); speclib = 0; disp('Canceled or aborted'); return; end
        speclib          = importdata(fullfile(path, filename));
    end

end

function varargout = column_data(ASCII, column_name_or_i)
% function [data, i, name] = column_data(ASCII, column_name_or_i)

% get data from column

% [i, name] = findASCIIcolumn_i_name(ASCII, column_name_or_i);
[varargout{1}, varargout{2}] = findASCIIcolumn_i_name(ASCII, column_name_or_i);

a = varargout{1};

data      = ASCII.data(:, a);
varargout = [{data}, varargout];


end

function varargout = findASCIIcolumn_i_name(ASCII, column_name_or_i)
    % function [i, line_num, name] = findASCIIcolumn_i_name(ASCII, column_name_or_i)

    % find column #:

    textdata = ASCII.textdata;

    if ischar(column_name_or_i)

        expression  = [': *', column_name_or_i];             % later it takes the index of one before ':' char to extract the column #
        a           = regexpi(textdata, expression);

        line_num = find(~ cellfun(@isempty, a));
        a        = a{line_num};
        line     = textdata{line_num};

        i = str2double(line(a - 1));

        varargout{2} = line_num;

    elseif isa(column_name_or_i, 'double')
        line = textdata{column_name_or_i};
        i    = column_name_or_i;

        %         words    = split(line);
        %         b        = contains(words, ':');
        %         name     = words{find(b) + 1};

        idx  = strfind(line, ':');
        name = line(idx + 2:end);

        varargout{2} = name;
    end

    %     words    = split(line);
    %     b        = contains(words, ':');
    %     name     = words{find(b) + 1};

    varargout{1} = i;

end

function [mineralname, mineral_lab_name] = name_enter(speclib, option)

            %     answer = inputdlg('Mineral Name of spectrum:', 'Input');
            %     if(isempty(answer) || isempty(answer{1}));   disp('Canceled, aborted or empty string');   return; end
            %     mineralname = answer{1};

            arguments
                speclib (1, 1) struct

                option.results (1, :) string = ''
            end

            textdata = speclib.textdata;

            % remove non relevants:
            [~, line_num]        = findASCIIcolumn_i_name(speclib, 'Wavelength');
            textdata(1:line_num) = [];

            % find minerals names:
            expression = 'Column\s\d:\s?(?<mineral>[a-zA-Z]+)_?';       % find mineral name (substring) in a a string
            minerals   = regexpi(textdata, expression, 'names');
            minerals   = cellfun(@(minerals) char(minerals.mineral), minerals, 'UniformOutput', false);    % extract names

            if ~ strcmpi(option.results, "")
                results_name = option.results;
                A            = cellfun(@(minerals) strfind(lower(results_name), lower(minerals)), minerals, 'UniformOutput', false);    % extract names
                B            = cellfun(@(A) ~ isempty(A), A);
                mineralname  = minerals{B};

                line             = textdata{B};
                idx              = strfind(line, ':');
                mineral_lab_name = line(idx + 2:end);

                return
            end

            % list dialog:
            [indx, tf] = listdlg('PromptString', {'Select a mineral.', ...
                                 'Only one mineral can be selected at a time.', ''}, ...
                                 'SelectionMode', 'single', ...
                                 'ListString', minerals);

            if isequal(tf, 0); mineralname = 0; disp('Canceled or aborted'); return; end

            mineralname = minerals{indx};

            line             = textdata{indx};
            idx              = strfind(line, ':');
            mineral_lab_name = line(idx + 2:end);

%             mineral_lab_name = '';

        end

function hcube = load_scene(filepath)
    % load scene:
    if ~exist("filepath", 'var')
        [filename, path] = uigetfile({'*if*.img'}, 'Select Scene');     if isequal(path, 0); disp('Canceled or aborted'); return; end
        filepath = fullfile(path, filename);
    end
    hcube = hypercube(filepath);
end

function hcube_or_wave = nm2um_cube(hcube_or_wave)

    if isobject(hcube_or_wave) % wavelength == hcube...
        wavelength = hcube_or_wave.Wavelength;
    else
        wavelength = hcube_or_wave;
    end
 
    if min(wavelength) >= 100
        if isobject(hcube_or_wave) % wavelength == hcube...
            %                 if strcmp(hcube_or_wave.Metadata.WavelengthUnits, 'Nanometers')
            metadata                 = hcube_or_wave.Metadata;
            Wavelength_um            = wavelength ./ 1000;      % values convertion
            metadata.WavelengthUnits = 'Micrometers';                 % units convertion
            hcube_or_wave            = hypercube(hcube_or_wave.DataCube, Wavelength_um, metadata);
        else
            hcube_or_wave = wavelength / 1000;
        end
    end

end























