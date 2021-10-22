function emptyLoc = checkTumorNeighborsLoc(state,i,j,k)
%CHECKNEIGHBORSLOC 3D finds empty neighbors of a specific location (i,j,k)

state_ind=1;

% nbrs in same layer
westNghbrs = state(i,j+1,k,state_ind);
eastNghbrs = state(i,j-1,k,state_ind);
northNghbrs = state(i+1,j,k,state_ind);
southNghbrs = state(i-1,j,k,state_ind);
northWestNghbrs = state(i+1,j+1,k,state_ind);
northEastNghbrs = state(i+1,j-1,k,state_ind);
southEastNghbrs = state(i-1,j-1,k,state_ind);
southWestNghbrs = state(i-1,j+1,k,state_ind);
% nbrs in layer above
upNghbrs = state(i,j,k+1,state_ind);
upwestNghbrs = state(i,j+1,k+1,state_ind);
upeastNghbrs = state(i,j-1,k+1,state_ind);
upnorthNghbrs = state(i+1,j,k+1,state_ind);
upsouthNghbrs = state(i-1,j,k+1,state_ind);
upnorthWestNghbrs = state(i+1,j+1,k+1,state_ind);
upnorthEastNghbrs = state(i+1,j-1,k+1,state_ind);
upsouthEastNghbrs = state(i-1,j-1,k+1,state_ind);
upsouthWestNghbrs = state(i-1,j+1,k+1,state_ind);
%nbrs in layer below
downNghbrs = state(i,j,k-1,state_ind);
downwestNghbrs = state(i,j+1,k-1,state_ind);
downeastNghbrs = state(i,j-1,k-1,state_ind);
downnorthNghbrs = state(i+1,j,k-1,state_ind);
downsouthNghbrs = state(i-1,j,k-1,state_ind);
downnorthWestNghbrs = state(i+1,j+1,k-1,state_ind);
downnorthEastNghbrs = state(i+1,j-1,k-1,state_ind);
downsouthEastNghbrs = state(i-1,j-1,k-1,state_ind);
downsouthWestNghbrs = state(i-1,j+1,k-1,state_ind);


% % nbrs in same layer
% westNghbrs = circshift(state(i,j,k,state_ind),[0,1,0,0]);
% eastNghbrs = circshift(state(i,j,k,state_ind),[0,-1,0,0]);
% northNghbrs = circshift(state(i,j,k,state_ind),[1,0,0,0]);
% southNghbrs = circshift(state(i,j,k,state_ind),[-1,0,0,0]);
% northWestNghbrs = circshift(state(i,j,k,state_ind),[1,1,0,0]);
% northEastNghbrs = circshift(state(i,j,k,state_ind),[1,-1,0,0]);
% southEastNghbrs = circshift(state(i,j,k,state_ind),[-1,-1,0,0]);
% southWestNghbrs = circshift(state(i,j,k,state_ind),[-1,1,0,0]);
% % nbrs in layer above
% upNghbrs = circshift(state(i,j,k,state_ind),[0,0,1,0]);
% upwestNghbrs = circshift(state(i,j,k,state_ind),[0,1,1,0]);
% upeastNghbrs = circshift(state(i,j,k,state_ind),[0,-1,1,0]);
% upnorthNghbrs = circshift(state(i,j,k,state_ind),[1,0,1,0]);
% upsouthNghbrs = circshift(state(i,j,k,state_ind),[-1,0,1,0]);
% upnorthWestNghbrs = circshift(state(i,j,k,state_ind),[1,1,1,0]);
% upnorthEastNghbrs = circshift(state(i,j,k,state_ind),[1,-1,1,0]);
% upsouthEastNghbrs = circshift(state(i,j,k,state_ind),[-1,-1,1,0]);
% upsouthWestNghbrs = circshift(state(i,j,k,state_ind),[-1,1,1,0]);
% %nbrs in layer below
% downNghbrs = circshift(state(i,j,k,state_ind),[0,0,-1,0]);
% downwestNghbrs = circshift(state(i,j,k,state_ind),[0,1,-1,0]);
% downeastNghbrs = circshift(state(i,j,k,state_ind),[0,-1,-1,0]);
% downnorthNghbrs = circshift(state(i,j,k,state_ind),[1,0,-1,0]);
% downsouthNghbrs = circshift(state(i,j,k,state_ind),[-1,0,-1,0]);
% downnorthWestNghbrs = circshift(state(i,j,k,state_ind),[1,1,-1,0]);
% downnorthEastNghbrs = circshift(state(i,j,k,state_ind),[1,-1,-1,0]);
% downsouthEastNghbrs = circshift(state(i,j,k,state_ind),[-1,-1,-1,0]);
% downsouthWestNghbrs = circshift(state(i,j,k,state_ind),[-1,1,-1,0]);


% westNghbrs = circshift(state,[0,1,0,0]);
% eastNghbrs = circshift(state,[0,-1,0,0]);
% northNghbrs = circshift(state,[1,0,0,0]);
% southNghbrs = circshift(state,[-1,0,0,0]);
% northWestNghbrs = circshift(state,[1,1,0,0]);
% northEastNghbrs = circshift(state,[1,-1,0,0]);
% southEastNghbrs = circshift(state,[-1,-1,0,0]);
% southWestNghbrs = circshift(state,[-1,1,0,0]);
% % nbrs in layer above
% upNghbrs = circshift(state,[0,0,1,0]);
% upwestNghbrs = circshift(state,[0,1,1,0]);
% upeastNghbrs = circshift(state,[0,-1,1,0]);
% upnorthNghbrs = circshift(state,[1,0,1,0]);
% upsouthNghbrs = circshift(state,[-1,0,1,0]);
% upnorthWestNghbrs = circshift(state,[1,1,1,0]);
% upnorthEastNghbrs = circshift(state,[1,-1,1,0]);
% upsouthEastNghbrs = circshift(state,[-1,-1,1,0]);
% upsouthWestNghbrs = circshift(state,[-1,1,1,0]);
% %nbrs in layer below
% downNghbrs = circshift(state,[0,0,-1,0]);
% downwestNghbrs = circshift(state,[0,1,-1,0]);
% downeastNghbrs = circshift(state,[0,-1,-1,0]);
% downnorthNghbrs = circshift(state,[1,0,-1,0]);
% downsouthNghbrs = circshift(state,[-1,0,-1,0]);
% downnorthWestNghbrs = circshift(state,[1,1,-1,0]);
% downnorthEastNghbrs = circshift(state,[1,-1,-1,0]);
% downsouthEastNghbrs = circshift(state,[-1,-1,-1,0]);
% downsouthWestNghbrs = circshift(state,[-1,1,-1,0]);

% % nbrs of the specific site (i,j,k) for its immune type
% westNghbrs = westNghbrs(i,j,k,state_ind);
% eastNghbrs = eastNghbrs(i,j,k,state_ind);
% northNghbrs = northNghbrs(i,j,k,state_ind);
% southNghbrs = southNghbrs(i,j,k,state_ind);
% northWestNghbrs = northWestNghbrs(i,j,k,state_ind);
% northEastNghbrs = northEastNghbrs(i,j,k,state_ind);
% southEastNghbrs = southEastNghbrs(i,j,k,state_ind);
% southWestNghbrs = southWestNghbrs(i,j,k,state_ind);
% 
% upNghbrs = upNghbrs(i,j,k,state_ind);
% upwestNghbrs = upwestNghbrs(i,j,k,state_ind);
% upeastNghbrs = upeastNghbrs(i,j,k,state_ind);
% upnorthNghbrs = upnorthNghbrs(i,j,k,state_ind);
% upsouthNghbrs = upsouthNghbrs(i,j,k,state_ind);
% upnorthWestNghbrs = upnorthWestNghbrs(i,j,k,state_ind);
% upnorthEastNghbrs = upnorthEastNghbrs(i,j,k,state_ind);
% upsouthEastNghbrs = upsouthEastNghbrs(i,j,k,state_ind);
% upsouthWestNghbrs = upsouthWestNghbrs(i,j,k,state_ind);
% 
% downNghbrs = downNghbrs(i,j,k,state_ind);
% downwestNghbrs = downwestNghbrs(i,j,k,state_ind);
% downeastNghbrs = downeastNghbrs(i,j,k,state_ind);
% downnorthNghbrs = downnorthNghbrs(i,j,k,state_ind);
% downsouthNghbrs = downsouthNghbrs(i,j,k,state_ind);
% downnorthWestNghbrs = downnorthWestNghbrs(i,j,k,state_ind);
% downnorthEastNghbrs = downnorthEastNghbrs(i,j,k,state_ind);
% downsouthEastNghbrs = downsouthEastNghbrs(i,j,k,state_ind);
% downsouthWestNghbrs = downsouthWestNghbrs(i,j,k,state_ind);








% Check how many empty cells around me
% westFree = zeros(n,n,n); eastFree = zeros(n,n,n); northFree = zeros(n,n,n); southFree = zeros(n,n,n);
% northWestFree = zeros(n,n,n); northEastFree = zeros(n,n,n); southWestFree = zeros(n,n,n); southEastFree = zeros(n,n,n);
% upFree = zeros(n,n,n); upwestFree = zeros(n,n,n); upeastFree = zeros(n,n,n); upnorthFree = zeros(n,n,n); upsouthFree = zeros(n,n,n);
% upnorthWestFree = zeros(n,n,n); upnorthEastFree = zeros(n,n,n); upsouthWestFree = zeros(n,n,n); upsouthEastFree = zeros(n,n,n);
% downFree = zeros(n,n,n); downwestFree = zeros(n,n,n); downeastFree = zeros(n,n,n); downnorthFree = zeros(n,n,n); downsouthFree = zeros(n,n,n);
% downnorthWestFree = zeros(n,n,n); downnorthEastFree = zeros(n,n,n); downsouthWestFree = zeros(n,n,n); downsouthEastFree = zeros(n,n,n);

% westFree(westNghbrs(i,j,k,1)==0)=1;
% eastFree(eastNghbrs(i,j,k,1)==0)=1;
% northFree(northNghbrs(i,j,k,1)==0)=1;
% southFree(southNghbrs(i,j,k,1)==0)=1;
% northWestFree(northWestNghbrs(i,j,k,1)==0)=1;
% northEastFree(northEastNghbrs(i,j,k,1)==0)=1;
% southWestFree(southWestNghbrs(i,j,k,1)==0)=1;
% southEastFree(southEastNghbrs(i,j,k,1)==0)=1;
% %above nbrs:
% upFree(upNghbrs(i,j,k,1)==0)=1;
% upwestFree(upwestNghbrs(i,j,k,1)==0)=1;
% upeastFree(upeastNghbrs(i,j,k,1)==0)=1;
% upnorthFree(upnorthNghbrs(i,j,k,1)==0)=1;
% upsouthFree(upsouthNghbrs(i,j,k,1)==0)=1;
% upnorthWestFree(upnorthWestNghbrs(i,j,k,1)==0)=1;
% upnorthEastFree(upnorthEastNghbrs(i,j,k,1)==0)=1;
% upsouthWestFree(upsouthWestNghbrs(i,j,k,1)==0)=1;
% upsouthEastFree(upsouthEastNghbrs(i,j,k,1)==0)=1;
% %below nbrs:
% downFree(downNghbrs(i,j,k,1)==0)=1;
% downwestFree(downwestNghbrs(i,j,k,1)==0)=1;
% downeastFree(downeastNghbrs(i,j,k,1)==0)=1;
% downnorthFree(downnorthNghbrs(i,j,k,1)==0)=1;
% downsouthFree(downsouthNghbrs(i,j,k,1)==0)=1;
% downnorthWestFree(downnorthWestNghbrs(i,j,k,1)==0)=1;
% downnorthEastFree(downnorthEastNghbrs(i,j,k,1)==0)=1;
% downsouthWestFree(downsouthWestNghbrs(i,j,k,1)==0)=1;
% downsouthEastFree(downsouthEastNghbrs(i,j,k,1)==0)=1;


%totalFree = westFree+eastFree+northFree+southFree+northEastFree+northWestFree+southEastFree+southWestFree...
% +upFree+upwestFree+upeastFree+upnorthFree+upsouthFree+upnorthEastFree+upnorthWestFree...
% +upsouthEastFree+upsouthWestFree+downFree+downwestFree+downeastFree+downnorthFree...
% +downsouthFree+downnorthEastFree+downnorthWestFree+downsouthEastFree+downsouthWestFree;






emptyLoc = [];
if westNghbrs==0
    emptyLoc = [emptyLoc; 0 1 0];
end
if eastNghbrs==0
    emptyLoc = [emptyLoc; 0 -1 0];
end
if northNghbrs==0
    emptyLoc = [emptyLoc; 1 0 0];
end
if southNghbrs==0
    emptyLoc = [emptyLoc; -1 0 0];
end
if northWestNghbrs==0
    emptyLoc = [emptyLoc; 1 1 0];
end
if northEastNghbrs==0
    emptyLoc = [emptyLoc; 1 -1 0];
end
if southEastNghbrs==0
    emptyLoc = [emptyLoc; -1 -1 0];
end
if southWestNghbrs==0
    emptyLoc = [emptyLoc; -1 1 0];
end
%nbrs above:
if upNghbrs==0
    emptyLoc = [emptyLoc; 0 0 1];
end
if upwestNghbrs==0
    emptyLoc = [emptyLoc; 0 1 1];
end
if upeastNghbrs==0
    emptyLoc = [emptyLoc; 0 -1 1];
end
if upnorthNghbrs==0
    emptyLoc = [emptyLoc; 1 0 1];
end
if upsouthNghbrs==0
    emptyLoc = [emptyLoc; -1 0 1];
end
if upnorthWestNghbrs==0
    emptyLoc = [emptyLoc; 1 1 1];
end
if upnorthEastNghbrs==0
    emptyLoc = [emptyLoc; 1 -1 1];
end
if upsouthEastNghbrs==0
    emptyLoc = [emptyLoc; -1 -1 1];
end
if upsouthWestNghbrs==0
    emptyLoc = [emptyLoc; -1 1 1];
end
%nbrs below:
if downNghbrs==0
    emptyLoc = [emptyLoc; 0 0 -1];
end
if downwestNghbrs==0
    emptyLoc = [emptyLoc; 0 1 -1];
end
if downeastNghbrs==0
    emptyLoc = [emptyLoc; 0 -1 -1];
end
if downnorthNghbrs==0
    emptyLoc = [emptyLoc; 1 0 -1];
end
if downsouthNghbrs==0
    emptyLoc = [emptyLoc; -1 0 -1];
end
if downnorthWestNghbrs==0
    emptyLoc = [emptyLoc; 1 1 -1];
end
if downnorthEastNghbrs==0
    emptyLoc = [emptyLoc; 1 -1 -1];
end
if downsouthEastNghbrs==0
    emptyLoc = [emptyLoc; -1 -1 -1];
end
if downsouthWestNghbrs==0
    emptyLoc = [emptyLoc; -1 1 -1];
end

% freeSpots=[westFree; eastFree; northFree; southFree; northEastFree; northWestFree; southEastFree; southWestFree;...
%  upFree; upwestFree; upeastFree; upnorthFree; upsouthFree; upnorthEastFree; upnorthWestFree;...
%  upsouthEastFree; upsouthWestFree; downFree; downwestFree; downeastFree; downnorthFree;...
%  downsouthFree; downnorthEastFree; downnorthWestFree; downsouthEastFree; downsouthWestFree];


end

