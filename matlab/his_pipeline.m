

% read in the data exported from R
base_dir = "~/NEMpipelineoutput/040_prepared"
old_dir = cd(base_dir)
listing = dir('lfc*.csv')
cd(old_dir)

for idx=1:length(listing)


    filename = listing(idx).name
    file = fullfile(base_dir, filename)

    opts = detectImportOptions(file)
    data = readtable(file)
    data = removevars(data, "Var1")


    % calculate the HIS on this data using the published his() function, with
    % GenomeRNAi Hs settings. Depending on your computer this may take a minute
    % or so. Increasing intMin and decreasing intSteps will make it go faster,
    % but find less interactions.


    intMin = 0.4;
    intMax = 2.5;
    intSteps = 150;
    strTail = 'both'; % default = both
    boolDrawAFigure = false; % switching on figure drawing will make it slower
    tic;matHIS = his(transpose(table2array(data)),intMin,intMax,intSteps,strTail,boolDrawAFigure);toc

    writematrix(full(matHIS), fullfile(base_dir, strcat("his.", filename)))

end


% Now we can calculate the enrichment p-values against a reference set of
% interactions using get_his_pvalue()
%boolDrawAFigure = true;
%matPs = get_his_pvalue(matHIS,reference_david,boolDrawAFigure);