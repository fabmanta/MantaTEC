function [ output ] = collapseData( data, gap_value )
%{
Removes gaps from time series data. Does not account for timing issues.
This is to be used in specific cases when you don't care about the timing
of the signal.
%}

collapsed_a = data(data~=gap_value);
collapsed_b = collapsed_a(~isnan(collapsed_a));
collapsed = collapsed_b;

output = collapsed;


end