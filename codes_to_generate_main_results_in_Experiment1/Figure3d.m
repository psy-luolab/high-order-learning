% get the averaged result across all EV (explained variance)



matrix = alldata;

[numRows, numCols] = size(matrix);

sumOddCols = zeros(numRows, 1);
sumEvenCols = zeros(numRows, 1);


for col = 1:numCols
    
    if rem(col, 2) == 1
        
        sumOddCols = sumOddCols + matrix(:, col);
    else
        
        sumEvenCols = sumEvenCols + matrix(:, col);
    end
end


avgOddCols = sumOddCols / (numCols / 2);
avgEvenCols = sumEvenCols / (numCols / 2);

[h,p,ci,stats]=ttest(avgEvenCols,avgOddCols)
