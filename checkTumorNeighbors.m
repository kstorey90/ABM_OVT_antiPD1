function totalFree = checkTumorNeighbors(state,n)
%CHECKNEIGHBORS 3D Computes the number of empty cells in the first
%order Moore neighborhood


% nbrs in same layer
westNghbrs = circshift(state,[0,1,0,0]);
eastNghbrs = circshift(state,[0,-1,0,0]);
northNghbrs = circshift(state,[1,0,0,0]);
southNghbrs = circshift(state,[-1,0,0,0]);
northWestNghbrs = circshift(state,[1,1,0,0]);
northEastNghbrs = circshift(state,[1,-1,0,0]);
southEastNghbrs = circshift(state,[-1,-1,0,0]);
southWestNghbrs = circshift(state,[-1,1,0,0]);
% nbrs in layer above
upNghbrs = circshift(state,[0,0,1,0]);
upwestNghbrs = circshift(state,[0,1,1,0]);
upeastNghbrs = circshift(state,[0,-1,1,0]);
upnorthNghbrs = circshift(state,[1,0,1,0]);
upsouthNghbrs = circshift(state,[-1,0,1,0]);
upnorthWestNghbrs = circshift(state,[1,1,1,0]);
upnorthEastNghbrs = circshift(state,[1,-1,1,0]);
upsouthEastNghbrs = circshift(state,[-1,-1,1,0]);
upsouthWestNghbrs = circshift(state,[-1,1,1,0]);
%nbrs in layer below
downNghbrs = circshift(state,[0,0,-1,0]);
downwestNghbrs = circshift(state,[0,1,-1,0]);
downeastNghbrs = circshift(state,[0,-1,-1,0]);
downnorthNghbrs = circshift(state,[1,0,-1,0]);
downsouthNghbrs = circshift(state,[-1,0,-1,0]);
downnorthWestNghbrs = circshift(state,[1,1,-1,0]);
downnorthEastNghbrs = circshift(state,[1,-1,-1,0]);
downsouthEastNghbrs = circshift(state,[-1,-1,-1,0]);
downsouthWestNghbrs = circshift(state,[-1,1,-1,0]);

% Check how many empty cells around me
westFree = zeros(n,n,n); eastFree = zeros(n,n,n); northFree = zeros(n,n,n); southFree = zeros(n,n,n);
northWestFree = zeros(n,n,n); northEastFree = zeros(n,n,n); southWestFree = zeros(n,n,n); southEastFree = zeros(n,n,n);
upFree = zeros(n,n,n); upwestFree = zeros(n,n,n); upeastFree = zeros(n,n,n); upnorthFree = zeros(n,n,n); upsouthFree = zeros(n,n,n);
upnorthWestFree = zeros(n,n,n); upnorthEastFree = zeros(n,n,n); upsouthWestFree = zeros(n,n,n); upsouthEastFree = zeros(n,n,n);
downFree = zeros(n,n,n); downwestFree = zeros(n,n,n); downeastFree = zeros(n,n,n); downnorthFree = zeros(n,n,n); downsouthFree = zeros(n,n,n);
downnorthWestFree = zeros(n,n,n); downnorthEastFree = zeros(n,n,n); downsouthWestFree = zeros(n,n,n); downsouthEastFree = zeros(n,n,n);

westFree(westNghbrs(:,:,:,1)==0)=1;
eastFree(eastNghbrs(:,:,:,1)==0)=1;
northFree(northNghbrs(:,:,:,1)==0)=1;
southFree(southNghbrs(:,:,:,1)==0)=1;
northWestFree(northWestNghbrs(:,:,:,1)==0)=1;
northEastFree(northEastNghbrs(:,:,:,1)==0)=1;
southWestFree(southWestNghbrs(:,:,:,1)==0)=1;
southEastFree(southEastNghbrs(:,:,:,1)==0)=1;
%above nbrs:
upFree(upNghbrs(:,:,:,1)==0)=1;
upwestFree(upwestNghbrs(:,:,:,1)==0)=1;
upeastFree(upeastNghbrs(:,:,:,1)==0)=1;
upnorthFree(upnorthNghbrs(:,:,:,1)==0)=1;
upsouthFree(upsouthNghbrs(:,:,:,1)==0)=1;
upnorthWestFree(upnorthWestNghbrs(:,:,:,1)==0)=1;
upnorthEastFree(upnorthEastNghbrs(:,:,:,1)==0)=1;
upsouthWestFree(upsouthWestNghbrs(:,:,:,1)==0)=1;
upsouthEastFree(upsouthEastNghbrs(:,:,:,1)==0)=1;
%below nbrs:
downFree(downNghbrs(:,:,:,1)==0)=1;
downwestFree(downwestNghbrs(:,:,:,1)==0)=1;
downeastFree(downeastNghbrs(:,:,:,1)==0)=1;
downnorthFree(downnorthNghbrs(:,:,:,1)==0)=1;
downsouthFree(downsouthNghbrs(:,:,:,1)==0)=1;
downnorthWestFree(downnorthWestNghbrs(:,:,:,1)==0)=1;
downnorthEastFree(downnorthEastNghbrs(:,:,:,1)==0)=1;
downsouthWestFree(downsouthWestNghbrs(:,:,:,1)==0)=1;
downsouthEastFree(downsouthEastNghbrs(:,:,:,1)==0)=1;


totalFree = westFree+eastFree+northFree+southFree+northEastFree+northWestFree+southEastFree+southWestFree...
+upFree+upwestFree+upeastFree+upnorthFree+upsouthFree+upnorthEastFree+upnorthWestFree...
+upsouthEastFree+upsouthWestFree+downFree+downwestFree+downeastFree+downnorthFree...
+downsouthFree+downnorthEastFree+downnorthWestFree+downsouthEastFree+downsouthWestFree;




end

