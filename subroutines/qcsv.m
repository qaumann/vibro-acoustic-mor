function qcsv(fname,m,varargin)
%QCSV Write csv with precision and header
%   Detailed explanation goes here

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

p = inputParser();
p.addParameter('delimiter',',');
p.addParameter('precision',10);
p.addParameter('header','');
p.parse(varargin{:});
delim = p.Results.delimiter;
prec = p.Results.precision;
head = p.Results.header;

if isnumeric(m)
    if head
        if length(strsplit(head,delim)) ~= size(m,2)
            error('qcsv: header does not fit data')
        end
        fid = fopen(fname,'wt');
        fprintf(fid, [head '\n']);
        fclose(fid);
        dlmwrite(fname,m,'delimiter',delim,'precision',prec,'-append');
    else
        dlmwrite(fname,m,'delimiter',delim,'precision',prec);  
    end
elseif istable(m)
    writetable(m, fname, 'delimiter', delim);
else
    if head
        if length(strsplit(head,delim)) ~= size(m,2)
            error('qcsv: header does not fit data')
        end
        T = table(m{:}, 'VariableNames', strsplit(head, delim));
    else
        T = table(m{:});
    end
    writetable(T, fname, 'delimiter', delim);%, 'QuoteStrings', true);
end
end