function [names_clusters, names_groups] = get_rgc_names()
    % Gets RGC cluster and group names
    names_clusters = {'1 OFF local, OS    1'; '2 OFF DS    2'; '3 OFF step    3'; '4a OFF slow    4';...
        '4b OFF slow    5';'5a OFF alpha sust.    6';'5b OFF alpha sust.    7';...
        '5c OFF alpha sust.    8';'6 (ON-)OFF "JAM-B" mix    9';'7 OFF sust.   10';...
        '8a OFF alpha trans.   11'; '8b OFF alpha trans.   12'; '9 OFF "mini" alpha trans.   13';...
        '10 ON-OFF local-edge "W3"   14';'11a ON-OFF local   15'; '11b ON-OFF local   16';...
        '12a ON-OFF DS 1   17'; '12b ON-OFF DS 1   18'; '13 ON-OFF DS 2   19'; ...
        '14 (ON-)OFF local, OS   20'; '15 ON step   21'; '16 ON DS trans.   22';...
        '17a ON local trans., OS   23'; '17b ON local trans., OS   24'; ...
        '17c ON local trans., OS   25'; '18a ON trans.   26'; '18b ON trans.   27';...
        '19 ON trans., large   28'; '20 ON high freq.   29'; '21 ON low freq.   30';...
        '22a ON sust.   31'; '22b ON sust   32'; '23 ON "mini" alpha   33'; '24 ON alpha   34'; ...
        '25 ON DS sust. 1   35'; '26 ON DS sust. 2   36'; '27 ON slow   37';...
        '28a ON contrast suppr.   38'; '28b ON contrast suppr.   39'; '29 ON DS sust. 3   40';...
        '30 ON local sust., OS   41'; '31a OFF suppr. 1   42'; '31b OFF suppr. 1   43';...
        '31c OFF suppr. 1   44'; '31d OFF suppr. 1   45'; '31e OFF suppr. 1   46';...
        '32a OFF suppr. 2   47'; '32b OFF suppr. 2   48'; '32c OFF suppr. 2   49'};

    names_groups = {'OFF local, OS    1'; 'OFF DS    2'; 'OFF step    3'; ' OFF slow     4';...
        'OFF alpha sust.    5';'(ON-)OFF "JAM-B" mix    6';'OFF sust.   7';...
        'OFF alpha trans.   8'; 'OFF "mini" alpha trans.    9'; ' ON-OFF local-edge "W3"   10';...
        'ON-OFF local   11'; 'ON-OFF DS 1   12'; 'ON-OFF DS 2   13'; ...
        '(ON-)OFF local, OS   14'; 'ON step   15'; 'ON DS trans.   16';...
        'ON local trans., OS   17'; 'ON trans.   18'; 'ON trans., large   19';...
        'ON high freq.   20'; 'ON low freq.   21'; 'ON sust.   22'; 'ON "mini" alpha   23';...
        'ON alpha   24'; 'ON DS sust. 1   25'; 'ON DS sust. 2   26'; 'ON slow   27';...
        'ON contrast suppr.   28'; 'ON DS sust. 3   29'; 'ON local sust., OS   30';...
        'OFF suppr. 1   31'; 'OFF suppr. 2   32'};
end