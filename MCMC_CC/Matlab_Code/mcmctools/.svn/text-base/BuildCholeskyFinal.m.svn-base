function [L] = BuildCholeskyFinal(V)
%This function takes a vector V of sample parameters and constructs a
%Cholesky Factorization from them using a recursive algorithm. Futhermore,
%it will track if the factorization is valid (i.e. real) and alert you
%appropriately.

%Include configure options and static parameters:
global CONFIG;

persistent Q;

%Start Debugging or timing, if appropriate
if(CONFIG.OPTIONS.DEBUG.Debug || CONFIG.OPTIONS.DEBUG.DebugFactor)
    DEBUG = true;
else
    DEBUG = false;
end
if(CONFIG.OPTIONS.TIMING.Timing || CONFIG.OPTIONS.TIMING.TimeFactor)
    tic;
end
%Log Messages if necessary
if(CONFIG.OPTIONS.DEBUG.HandleMessages)
    HM = true;
else
    HM = false;
end

%Use a queue to store recent factorizations to improve effiency and
%eliminate double-work
if(isempty(Q))
    if(DEBUG)
        if(HM)
            HandleMessages(1, 'Building Queue');
        else
            disp('Building Queue');
        end
    end
    for i = 1:CONFIG.OPTIONS.LIKELIHOOD.Qsize
        Q(i).vector = [];
        Q(i).factor = [];
    end
else
    for q = 1:length(Q)
        if ~isempty(Q(q).vector) 
            if(V == Q(q).vector)
                if(DEBUG)
                    if(HM)
                        HandleMessages(1, 'Returning Factorization From Queue');
                    else
                        disp('Returning Factorization From Queue');
                    end
                end
                L = full(Q(q).factor);
                return;
            end
        end
    end
end

if(DEBUG)
    if(HM)
        HandleMessages(1,'Building New Factorization')
    else
        disp('Building New Factorization');
    end
end

%Check to see if the toeplitz matrix has cusps or not
if(CONFIG.OPTIONS.LIKELIHOOD.Cusping)
    if(DEBUG)
        if(HM)
            HandleMessages(1, 'Cusping turned on');
        else
            disp('Cusping turned on');
        end
    end
    datarowsize = CONFIG.STATIC.nl;
end

%N is total number of data points (dim of covariance matrix)
N = CONFIG.STATIC.nl*CONFIG.STATIC.nz;
V2 = V ==0;
V2 = +V2;
%If the initial vector has a zero, that entire diagonal of the covariance
%matrix is zero, so we don't have to calculate them
IsDefined = toeplitz(V2);
%... as long as L is initialized to zero.
L = zeros(N,N);

for row=1:N
    for col=1:row
        if(~IsDefined(row,col))
            if(col==1)
                %First column is the vector element
                L(row,col)=V(1,row);
            elseif(row==col)
                tempsum=0;
                for l=1:col-1
                    tempsum = tempsum+(L(row,l))^2;
                end
                sumordiff = V(1,1)^2-tempsum;
                if(sumordiff <= 0)  %not a valid factorization
                    if(HM)
                        HandleMessages(3, 'Complex Factorization. Check parameters.');
                    end
                    error('Complex Factorization. Check parameters.');
                end
                L(row,col) =(-1)^(floor(rand*2))*(sumordiff)^(.5);
                if(L(row,col) == 0) %just in case
                    if(HM)
                        HandlMessages(3, 'Zero variance error.');
                    end
                    error('Zero variance error.');
                end
                
            else
                if(CONFIG.OPTIONS.LIKELIHOOD.Cusping)
                    if((mod(row-1,datarowsize)==0) && (mod(col, datarowsize)==0) && ((row-1)/datarowsize == col/datarowsize))
                        %in a cusp
                        if(DEBUG)
                            if(HM)
                                HandleMessages(1, 'in a cusp');
                            else
                                disp('in a cusp');
                            end
                        end
                        L(row,col)=0;
                    else
                        tempsum=0;
                        for k=1:(col-1)
                            tempsum = tempsum + L(col,k)*L(row,k);
                        end;
                        
                        L(row,col) = (L(1,1)*L(row-col+1,1)-tempsum)/L(col,col);
                    end
                else
                    tempsum=0;
                    for k=1:(col-1)
                        tempsum = tempsum + L(col,k)*L(row,k);
                    end;
                    
                    L(row,col) = (L(1,1)*L(row-col+1,1)-tempsum)/L(col,col);
                    
                end
            end
        end
        if(DEBUG)
            temp = ['Storing ' num2str(L(row,col)) ' at point: ' num2str(row) ', ' num2str(col)];
            if(HM)
                HandleMessages(1, temp);
            else
                disp(temp);
            end
        end
    end
end
L = sparse(L);
%store the vector and factorization into the queue
Q(2:length(Q))= Q(1:length(Q)-1);
Q(1).vector = V;
Q(1).factor = L;

if(CONFIG.OPTIONS.TIMING.Timing || CONFIG.OPTIONS.TIMING.TimeFactor)
    temp = toc;
    disp(['Building Cholesky took ' num2str(temp) ' seconds.']);
end
