function y=Q(x)
% co-error function: 1/sqrt(2*pi) * int_x^inf exp(-t^2/2) dt.£¨Îó²îº¯Êý£©
y=erfc(x/sqrt(2))/2;    %erfc»¥²¹Îó²îº¯Êý¡¢Qº¯Êý¡¢erfÎó²îº¯Êý