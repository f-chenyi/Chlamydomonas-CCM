% given the barcodes of a set of engineering states, generates a matrix of
% distances between each engineering state (ie., how many steps it would
% take to get from a given state to another.

function dist = gen_distance(barcode)

N = size(barcode, 1);
dist = zeros([N,N]); 

for i = 1:N
    for j = 1:i
        
        dist(i,j) = dist(i,j) + sum(abs(barcode(i, 1:5) - barcode(j, 1:5)));
        dist(j,i) = dist(j,i) + sum(abs(barcode(i, 1:5) - barcode(j, 1:5)));
        
        % account for separate paths of starch vs. thylakoids
        if barcode(i, 6) == barcode(j, 6)

        elseif barcode(i, 6) == 1 && barcode(j, 6) == 2
            dist(i, j) = dist(i, j) + 2;
            dist(j, i) = dist(j, i) + 2;
        elseif barcode(j, 6) == 1 && barcode(i, 6) == 2
            dist(i, j) = dist(i, j) + 2;
            dist(j, i) = dist(j, i) + 2;
            
        else 
            dist(i, j) = dist(i, j) + 1;
            dist(j, i) = dist(j, i) + 1;
        end 
        
        % don't add a starch sheath before Rub localized
        if barcode(j, 5) == 1 && barcode(j, 6) == 2
            if i ~= j
                dist(i,j) = dist(i,j) + 1;
                dist(j,i) = dist(j,i) + 1;
            end 
        end 
        
        % LCIB cannot localize before Rubisco
        if barcode(j, 4) == 2 && barcode(j, 5) == 1
            if barcode(i, 4) == 1 && barcode(i, 5) == 1
                dist(i, j) = dist(i, j) + 1;
                dist(j, i) = dist(j, i) + 1;
            elseif barcode(i, 4) == 2 && barcode(i, 5) == 2
                dist(i, j) = dist(i, j) + 1;
                dist(j, i) = dist(j, i) + 1;
            end 
        end 
               
        % LCIB and Rubisco can localize / delocalize in the same step        
        if barcode(i, 4) == 1 && barcode(i, 5) == 1
            if barcode(j, 4) == 2 && barcode(j, 5) == 2
                dist(i, j) = dist(i, j) - 1;
                dist(j, i) = dist(j, i) - 1;
            end
        elseif barcode(i, 4) == 2 && barcode(i, 5) == 2
            if barcode(j, 4) == 1 && barcode(j, 5) == 1
                dist(i, j) = dist(i, j) - 1;
                dist(j, i) = dist(j, i) - 1;
            end
        end
        
    
    end
end 

            
    