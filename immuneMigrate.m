function [state, DeathCounter] = immuneMigrate(state, i, j, k, type, prop, DeathCounter)
  %type = type of immune cell (1=innate, 2=Adap Tumor, 3=Adap Viral)  
  %prop = propensity for infected/tumor cell (depending on type)
  

% dim of state: nxnxnx3
n=size(state,1);

%First set with 3's, and then update only the ones in the domain
westNghbrs = [3 3 3 3]; 
eastNghbrs = [3 3 3 3]; 
northNghbrs = [3 3 3 3]; 
southNghbrs = [3 3 3 3]; 
northWestNghbrs = [3 3 3 3]; 
northEastNghbrs = [3 3 3 3]; 
southEastNghbrs = [3 3 3 3]; 
southWestNghbrs = [3 3 3 3]; 
% nbrs in layer above
upNghbrs = [3 3 3 3]; 
upwestNghbrs = [3 3 3 3]; 
upeastNghbrs = [3 3 3 3]; 
upnorthNghbrs = [3 3 3 3]; 
upsouthNghbrs = [3 3 3 3]; 
upnorthWestNghbrs = [3 3 3 3]; 
upnorthEastNghbrs = [3 3 3 3]; 
upsouthEastNghbrs = [3 3 3 3]; 
upsouthWestNghbrs = [3 3 3 3]; 
%nbrs in layer below
downNghbrs = [3 3 3 3]; 
downwestNghbrs = [3 3 3 3]; 
downeastNghbrs = [3 3 3 3]; 
downnorthNghbrs = [3 3 3 3]; 
downsouthNghbrs = [3 3 3 3]; 
downnorthWestNghbrs = [3 3 3 3]; 
downnorthEastNghbrs = [3 3 3 3]; 
downsouthEastNghbrs = [3 3 3 3]; 
downsouthWestNghbrs = [3 3 3 3]; 

% 
% % nbrs in same layer
% if i==1
%    southNghbrs = [3 3 3 3]; 
%    southEastNghbrs = [3 3 3 3];
%    southWestNghbrs = [3 3 3 3];
%    upsouthNghbrs = [3 3 3 3];
%    upsouthEastNghbrs = [3 3 3 3];
%    upsouthWestNghbrs = [3 3 3 3];
%    downsouthNghbrs = [3 3 3 3];
%    downsouthEastNghbrs = [3 3 3 3];
%    downsouthWestNghbrs = [3 3 3 3];
%    if j==1
%        eastNghbrs = [3 3 3 3];
%        northEastNghbrs = [3 3 3 3];
%        upeastNghbrs = [3 3 3 3];
%        upnorthEastNghbrs = [3 3 3 3];
%        downeastNghbrs = [3 3 3 3];
%        downnorthEastNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end   
%    elseif j==n
%        westNghbrs = [3 3 3 3];
%        northWestNghbrs = [3 3 3 3];
%        upwestNghbrs = [3 3 3 3];
%        upnorthWestNghbrs = [3 3 3 3];
%        downwestNghbrs = [3 3 3 3];
%        downnorthWestNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end
%    else %j not on boundary
%        westNghbrs = state(i,j+1,k,:);
%        eastNghbrs = state(i,j-1,k,:);
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end  
%    end
% elseif i==n
%    northNghbrs = [3 3 3 3];
%    northEastNghbrs = [3 3 3 3];
%    northWestNghbrs = [3 3 3 3];
%    upnorthNghbrs = [3 3 3 3];
%    upnorthEastNghbrs = [3 3 3 3];
%    upnorthWestNghbrs = [3 3 3 3];
%    downnorthNghbrs = [3 3 3 3];
%    downnorthEastNghbrs = [3 3 3 3];
%    downnorthWestNghbrs = [3 3 3 3];
%    if j==1
%        eastNghbrs = [3 3 3 3];
%        southEastNghbrs = [3 3 3 3];
%        upeastNghbrs = [3 3 3 3];
%        upsouthEastNghbrs = [3 3 3 3];
%        downeastNghbrs = [3 3 3 3];
%        downsouthEastNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end 
%     elseif j==n
%        westNghbrs = [3 3 3 3];
%        northWestNghbrs = [3 3 3 3];
%        upwestNghbrs = [3 3 3 3];
%        upnorthWestNghbrs = [3 3 3 3];
%        downwestNghbrs = [3 3 3 3];
%        downnorthWestNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end
%    else %j not on boundary
%        westNghbrs = state(i,j+1,k,:);
%        eastNghbrs = state(i,j-1,k,:);
%        if k==1
%           upNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%        end
%    end
%    
% else %i not on boundary
%     northNghbrs = state(i+1,j,k,:);
%     southNghbrs = state(i-1,j,k,:);
%     if j==1
%         southEastNghbrs = [3 3 3 3];
%         upsouthEastNghbrs = [3 3 3 3];
%         downsouthEastNghbrs = [3 3 3 3];
%         eastNghbrs = [3 3 3 3];
%         northEastNghbrs = [3 3 3 3];
%         upeastNghbrs = [3 3 3 3];
%         upnorthEastNghbrs = [3 3 3 3];
%         downeastNghbrs = [3 3 3 3];
%         downnorthEastNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%           upnorthNghbrs = [3 3 3 3];
%           upsouthNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%           downnorthNghbrs = [3 3 3 3];
%           downsouthNghbrs = [3 3 3 3];
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%           upsouthNghbrs = state(i-1,j,k+1,:);
%           upnorthNghbrs = state(i+1,j,k+1,:);
%        end 
%     elseif j==n    
%         southWestNghbrs = [3 3 3 3];
%         upsouthWestNghbrs = [3 3 3 3];
%         downsouthWestNghbrs = [3 3 3 3];
%         WestNghbrs = [3 3 3 3];
%         northWestNghbrs = [3 3 3 3];
%         upWestNghbrs = [3 3 3 3];
%         upnorthWestNghbrs = [3 3 3 3];
%         downWestNghbrs = [3 3 3 3];
%         downnorthWestNghbrs = [3 3 3 3];
%        if k==1
%           upNghbrs = [3 3 3 3];
%           upnorthNghbrs = [3 3 3 3];
%           upsouthNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%           downnorthNghbrs = [3 3 3 3];
%           downsouthNghbrs = [3 3 3 3];
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%           upsouthNghbrs = state(i-1,j,k+1,:);
%           upnorthNghbrs = state(i+1,j,k+1,:);
%        end
%     else   %j not on boundary
%        westNghbrs = state(i,j+1,k,:);
%        eastNghbrs = state(i,j-1,k,:);
%        if k==1
%           upNghbrs = [3 3 3 3];
%           upnorthNghbrs = [3 3 3 3];
%           upsouthNghbrs = [3 3 3 3];
%        elseif k==n
%           downNghbrs = [3 3 3 3]; 
%           downnorthNghbrs = [3 3 3 3];
%           downsouthNghbrs = [3 3 3 3];
%        else
%           upNghbrs = state(i,j,k+1,:);
%           downNghbrs = state(i,j,k-1,:); 
%           upsouthNghbrs = state(i-1,j,k+1,:);
%           upnorthNghbrs = state(i+1,j,k+1,:);
%        end
%        
%     end
%     
% end
   

if i~=1
    southNghbrs = state(i-1,j,k,:); 
    if j~=1
        eastNghbrs = state(i,j-1,k,:);
        southEastNghbrs = state(i-1,j-1,k,:);
        if k~=1 
            downNghbrs = state(i,j,k-1,:); 
            downeastNghbrs = state(i,j-1,k-1,:);
            downsouthNghbrs = state(i-1,j,k-1,:);
            downsouthEastNghbrs = state(i-1,j-1,k-1,:);
        end
        if k~=n
            upNghbrs = state(i,j,k+1,:);
            upeastNghbrs = state(i,j-1,k+1,:);
            upsouthNghbrs = state(i-1,j,k+1,:);
            upsouthEastNghbrs = state(i-1,j-1,k+1,:);
        end  
    end
    
    if j~=n
        westNghbrs = state(i,j+1,k,:);
        southWestNghbrs = state(i-1,j+1,k,:);
        if k~=1 
            downNghbrs = state(i,j,k-1,:); 
            downwestNghbrs = state(i,j+1,k-1,:);
            downsouthNghbrs = state(i-1,j,k-1,:);
            downsouthWestNghbrs = state(i-1,j+1,k-1,:);
        end
        if k~=n
            upNghbrs = state(i,j,k+1,:);
            upwestNghbrs = state(i,j+1,k+1,:);
            upsouthNghbrs = state(i-1,j,k+1,:);
            upsouthWestNghbrs = state(i-1,j+1,k+1,:);
        end  
    end
end 

if i~=n
    northNghbrs = state(i+1,j,k,:);
    if j~=1
        eastNghbrs = state(i,j-1,k,:);
        northEastNghbrs = state(i+1,j-1,k,:);
        if k~=1 
            downNghbrs = state(i,j,k-1,:); 
            downeastNghbrs = state(i,j-1,k-1,:);
            downnorthNghbrs = state(i+1,j,k-1,:);
            downnorthEastNghbrs = state(i+1,j-1,k-1,:);
        end
        if k~=n
            upNghbrs = state(i,j,k+1,:);
            upeastNghbrs = state(i,j-1,k+1,:);
            upnorthNghbrs = state(i+1,j,k+1,:);
            upnorthEastNghbrs = state(i+1,j-1,k+1,:);
        end  
    end
    
    if j~=n
        westNghbrs = state(i,j+1,k,:);
        northWestNghbrs = state(i+1,j+1,k,:);
        if k~=1 
            downNghbrs = state(i,j,k-1,:); 
            downwestNghbrs = state(i,j+1,k-1,:);
            downnorthNghbrs = state(i+1,j,k-1,:);
            downnorthWestNghbrs = state(i+1,j+1,k-1,:);
        end
        if k~=n
            upNghbrs = state(i,j,k+1,:);
            upwestNghbrs = state(i,j+1,k+1,:);
            upnorthNghbrs = state(i+1,j,k+1,:);
            upnorthWestNghbrs = state(i+1,j+1,k+1,:);
        end  
    end
end


% % nbrs in layer above
% if k~=n
%     upNghbrs = state(i,j,k+1,:);
%     if i~=1
%         upsouthNghbrs = state(i-1,j,k+1,:);
%         
%         if j~=1
%             upeastNghbrs = state(i,j-1,k+1,:);
%             upsouthEastNghbrs = state(i-1,j-1,k+1,:);
%         else
%             upeastNghbrs = [3 3 3 3];
%             upsouthEastNghbrs = [3 3 3 3];
%         end
% 
%         if j~=n
%             upwestNghbrs = state(i,j+1,k+1,:);
%             upsouthWestNghbrs = state(i-1,j+1,k+1,:);
%         else
%             upwestNghbrs = [3 3 3 3];
%             upsouthWestNghbrs = [3 3 3 3];
%         end
%     else
%         upsouthNghbrs = [3 3 3 3];
%     end    
%     if i~=n
%         upnorthNghbrs = state(i+1,j,k+1,:);
%         if j~=1
%             upnorthEastNghbrs = state(i+1,j-1,k+1,:);
%         else
%             upnorthEastNghbrs = [3 3 3 3];
%         end
%         if j~=n
%             upnorthWestNghbrs = state(i+1,j+1,k+1,:);
%         else
%             upnorthWestNghbrs = [3 3 3 3];
%         end
%     else
%         upnorthNghbrs = [3 3 3 3];
%     end 
% else
%     upNghbrs = [3 3 3 3];
% end
% %nbrs in layer below
% if k~=1
%     downNghbrs = state(i,j,k-1,:);
%     if i~=1
%         downsouthNghbrs = state(i-1,j,k-1,:);
%         
%         if j~=1
%             downeastNghbrs = state(i,j-1,k-1,:);
%             downsouthEastNghbrs = state(i-1,j-1,k-1,:);
%         else
%             downeastNghbrs = [3 3 3 3];
%             downsouthEastNghbrs = [3 3 3 3];
%         end
% 
%         if j~=n
%             downwestNghbrs = state(i,j+1,k-1,:);
%             downsouthWestNghbrs = state(i-1,j+1,k-1,:);
%         else
%             downwestNghbrs = [3 3 3 3];
%             downsouthWestNghbrs = [3 3 3 3];
%         end
%     else
%         downsouthNghbrs = [3 3 3 3];
%     end    
%     if i~=n
%         downnorthNghbrs = state(i+1,j,k-1,:);
%         if j~=1
%             downnorthEastNghbrs = state(i+1,j-1,k-1,:);
%         else
%             downnorthEastNghbrs = [3 3 3 3];
%         end
%         if j~=n
%             downnorthWestNghbrs = state(i+1,j+1,k-1,:);
%         else
%             downnorthWestNghbrs = [3 3 3 3];
%         end
%     else
%         downnorthNghbrs = [3 3 3 3];
%     end 
% else
%     downNghbrs = [3 3 3 3];
% end


% westNghbrs = state(i,j+1,k,:);
% eastNghbrs = state(i,j-1,k,:);
% northNghbrs = state(i+1,j,k,:);
% southNghbrs = state(i-1,j,k,:);
% northWestNghbrs = state(i+1,j+1,k,:);
% northEastNghbrs = state(i+1,j-1,k,:);
% southEastNghbrs = state(i-1,j-1,k,:);
% southWestNghbrs = state(i-1,j+1,k,:);
% % nbrs in layer above
% upNghbrs = state(i,j,k+1,:);
% upwestNghbrs = state(i,j+1,k+1,:);
% upeastNghbrs = state(i,j-1,k+1,:);
% upnorthNghbrs = state(i+1,j,k+1,:);
% upsouthNghbrs = state(i-1,j,k+1,:);
% upnorthWestNghbrs = state(i+1,j+1,k+1,:);
% upnorthEastNghbrs = state(i+1,j-1,k+1,:);
% upsouthEastNghbrs = state(i-1,j-1,k+1,:);
% upsouthWestNghbrs = state(i-1,j+1,k+1,:);
% %nbrs in layer below
% downNghbrs = state(i,j,k-1,:);
% downwestNghbrs = state(i,j+1,k-1,:);
% downeastNghbrs = state(i,j-1,k-1,:);
% downnorthNghbrs = state(i+1,j,k-1,:);
% downsouthNghbrs = state(i-1,j,k-1,:);
% downnorthWestNghbrs = state(i+1,j+1,k-1,:);
% downnorthEastNghbrs = state(i+1,j-1,k-1,:);
% downsouthEastNghbrs = state(i-1,j-1,k-1,:);
% downsouthWestNghbrs = state(i-1,j+1,k-1,:);



%find free immune cell spots
if type==1
    state_ind=2;    %index in state matrix for imm cell 
elseif type==2  %adap tumor
    state_ind=3;
else        %adap viral
    state_ind=3;
end

emptyLoc = []; %stores location of empty sites (no imm cells) without tumor cells
tumorLoc = [];  %stores location of empty sites (no imm cells) with (inf/any) tumor cells

% figure out which sites are empty (no imm cells), and which have tumor cells 
if westNghbrs(state_ind)==0
    if westNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 1 0];
    elseif (type==2 && westNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 1 0];
    else
        emptyLoc = [emptyLoc; 0 1 0];
    end
end
if eastNghbrs(state_ind)==0
    if eastNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 -1 0];
    elseif (type==2 && eastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 -1 0];
    else
        emptyLoc = [emptyLoc; 0 -1 0];
    end
end
if northNghbrs(state_ind)==0
    if northNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 0 0];
    elseif (type==2 && northNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 0 0];
    else
        emptyLoc = [emptyLoc; 1 0 0];
    end       
end
if southNghbrs(state_ind)==0
    if southNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 0 0];
    elseif (type==2 && southNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 0 0];
    else
        emptyLoc = [emptyLoc; -1 0 0];
    end
end
if northWestNghbrs(state_ind)==0
    if northWestNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 1 0];
    elseif (type==2 && northWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 1 0];
    else
        emptyLoc = [emptyLoc; 1 1 0];
    end
end
if northEastNghbrs(state_ind)==0
    if northEastNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 -1 0];
    elseif (type==2 && northEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 -1 0];
    else
        emptyLoc = [emptyLoc; 1 -1 0];
    end
end
if southEastNghbrs(state_ind)==0
    if southEastNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 -1 0];
    elseif (type==2 && southEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 -1 0];
    else
        emptyLoc = [emptyLoc; -1 -1 0];
    end
end
if southWestNghbrs(state_ind)==0
    if southWestNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 1 0];
    elseif (type==2 && southWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 1 0];
    else
        emptyLoc = [emptyLoc; -1 1 0];
    end
end
%same with nbrs above:
if upNghbrs(state_ind)==0
    if upNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 0 1];
    elseif (type==2 && upNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 0 1];
    else
        emptyLoc = [emptyLoc; 0 0 1];
    end
end
if upwestNghbrs(state_ind)==0
    if upwestNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 1 1];
    elseif (type==2 && upwestNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 1 1];
    else
        emptyLoc = [emptyLoc; 0 1 1];
    end
end
if upeastNghbrs(state_ind)==0
    if upeastNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 -1 1];
    elseif (type==2 && upeastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 -1 1];
    else
        emptyLoc = [emptyLoc; 0 -1 1];
    end
end
if upnorthNghbrs(state_ind)==0
    if upnorthNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 0 1];
    elseif (type==2 && upnorthNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 0 1];
    else
        emptyLoc = [emptyLoc; 1 0 1];
    end       
end
if upsouthNghbrs(state_ind)==0
    if upsouthNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 0 1];
    elseif (type==2 && upsouthNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 0 1];
    else
        emptyLoc = [emptyLoc; -1 0 1];
    end
end
if upnorthWestNghbrs(state_ind)==0
    if upnorthWestNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 1 1];
    elseif (type==2 && upnorthWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 1 1];
    else
        emptyLoc = [emptyLoc; 1 1 1];
    end
end
if upnorthEastNghbrs(state_ind)==0
    if upnorthEastNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 -1 1];
    elseif (type==2 && upnorthEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 -1 1];
    else
        emptyLoc = [emptyLoc; 1 -1 1];
    end
end
if upsouthEastNghbrs(state_ind)==0
    if upsouthEastNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 -1 1];
    elseif (type==2 && upsouthEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 -1 1];
    else
        emptyLoc = [emptyLoc; -1 -1 1];
    end
end
if upsouthWestNghbrs(state_ind)==0
    if upsouthWestNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 1 1];
    elseif (type==2 && upsouthWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 1 1];
    else
        emptyLoc = [emptyLoc; -1 1 1];
    end
end
%same with nbrs below:
if downNghbrs(state_ind)==0
    if downNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 0 -1];
    elseif (type==2 && downNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 0 -1];
    else
        emptyLoc = [emptyLoc; 0 0 -1];
    end
end
if downwestNghbrs(state_ind)==0
    if downwestNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 1 -1];
    elseif (type==2 && downwestNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 1 -1];
    else
        emptyLoc = [emptyLoc; 0 1 -1];
    end
end
if downeastNghbrs(state_ind)==0
    if downeastNghbrs(1)==2
        tumorLoc = [tumorLoc; 0 -1 -1];
    elseif (type==2 && downeastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 0 -1 -1];
    else
        emptyLoc = [emptyLoc; 0 -1 -1];
    end
end
if downnorthNghbrs(state_ind)==0
    if downnorthNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 0 -1];
    elseif (type==2 && downnorthNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 0 -1];
    else
        emptyLoc = [emptyLoc; 1 0 -1];
    end       
end
if downsouthNghbrs(state_ind)==0
    if downsouthNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 0 -1];
    elseif (type==2 && downsouthNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 0 -1];
    else
        emptyLoc = [emptyLoc; -1 0 -1];
    end
end
if downnorthWestNghbrs(state_ind)==0
    if downnorthWestNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 1 -1];
    elseif (type==2 && downnorthWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 1 -1];
    else
        emptyLoc = [emptyLoc; 1 1 -1];
    end
end
if downnorthEastNghbrs(state_ind)==0
    if downnorthEastNghbrs(1)==2
        tumorLoc = [tumorLoc; 1 -1 -1];
    elseif (type==2 && downnorthEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; 1 -1 -1];
    else
        emptyLoc = [emptyLoc; 1 -1 -1];
    end
end
if downsouthEastNghbrs(state_ind)==0
    if downsouthEastNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 -1 -1];
    elseif (type==2 && downsouthEastNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 -1 -1];
    else
        emptyLoc = [emptyLoc; -1 -1 -1];
    end
end
if downsouthWestNghbrs(state_ind)==0
    if downsouthWestNghbrs(1)==2
        tumorLoc = [tumorLoc; -1 1 -1];
    elseif (type==2 && downsouthWestNghbrs(1)==1)
        tumorLoc = [tumorLoc; -1 1 -1];
    else
        emptyLoc = [emptyLoc; -1 1 -1];
    end
end





numEmpty = size(emptyLoc,1);    %number of free spots with no tum cells
numTum = size(tumorLoc,1);      %number of free spots with (inf/any) tum cells


% Decide which site to move to (higher propensity for sites occupied by tumor cells 
rateSum = numEmpty+prop*numTum;
siteRand = rand*rateSum;    %multiply U(0,1) rand number by rateSum
if siteRand>0
    if siteRand <= numEmpty  %choose one of non-tum sites to move to
        siteInd=max(ceil(siteRand),1);
        newi=i+emptyLoc(siteInd,1);
        newj=j+emptyLoc(siteInd,2);
        newk=k+emptyLoc(siteInd,3);
        if type==3
            state(newi,newj,newk,state_ind)=2;
        else
            state(newi,newj,newk,state_ind)=1;
        end
        state(i,j,k,state_ind)=0;
        DeathCounter(newi,newj,newk,:)=DeathCounter(i,j,k,:);
        DeathCounter(i,j,k,:)=[0; 0];
    else    %choose tumor site to move to
        siteInd = max(ceil((rateSum-siteRand)/prop),1);
        newi=i+tumorLoc(siteInd,1);
        newj=j+tumorLoc(siteInd,2);
        newk=k+tumorLoc(siteInd,3);
        if type==3
            state(newi,newj,newk,state_ind)=2;
        else
            state(newi,newj,newk,state_ind)=1;
        end
        state(i,j,k,state_ind)=0;
        DeathCounter(newi,newj,newk,:)=DeathCounter(i,j,k,:);
        DeathCounter(i,j,k,:)=[0; 0];
    end
end


end