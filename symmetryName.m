function nam = symmetryName(symcode)
% turns the relatively arbitrary numeric system I made into
% standard Schoenflies notation.
K = num2str(symcode(2));
switch symcode(1)
    case 1
        nam = ['C_{' K '}'];
    case 2
        if symcode(2)==1
            nam = 'C_s';
        else
            nam = ['C_{' K 'v}'];
        end
    case 3
        if symcode(2)==1
            nam = 'C_s';
        else
            nam = ['C_{' K 'h}'];
        end
    case 4
        if symcode(2)==1
            nam = 'C_i';
        else
            nam = ['R_{' K '}'];
        end
    case 5
        nam = ['D_{' K '}'];
    case 6 
        nam = ['D_{' K 'd}'];
    case 7
        nam = ['D_{' K 'h}'];
    case 8
        nam = 'T';
    case 9
        nam = 'T_h';
    case 10
        nam = 'T_d';
    case 11
        nam = 'O';
    case 12
        nam = 'O_h';
    case 13
        nam = 'I';
    case 14
        nam = 'I_h';
end
end