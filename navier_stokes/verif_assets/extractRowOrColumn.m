function extractedData = extractRowOrColumn(matrix, row_or_col, index)
    if (row_or_col=="row")
        extractedData = matrix(index, :); 
    else
        extractedData = matrix(:, index); 
    end 
    % Convert the extracted data to a vertical vector
    extractedData = extractedData(:);
end