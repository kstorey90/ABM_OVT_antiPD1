function totals = checkAllNeighbors(state,n)
%CHECKNEIGHBORS Computes the number of all neighbor types in the first
%order Moore neighborhood, 
% third index values: 1-susc, 2-inf, 3-innate, 4-antitumor adap,
% 5-antiviral adap

numtypes=5;

totals=zeros(n,n,n,numtypes);


%nbrs in same layer
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



% westFree = zeros(n,n,n); eastFree = zeros(n,n,n); northFree = zeros(n,n,n); southFree = zeros(n,n,n);
% northWestFree = zeros(n,n,n); northEastFree = zeros(n,n,n); southWestFree = zeros(n,n,n); southEastFree = zeros(n,n,n);
% Check how many empty cells around me:
% westFree(westNghbrs(:,:,:,1)==0)=1;
% eastFree(eastNghbrs(:,:,:,1)==0)=1;
% northFree(northNghbrs(:,:,:,1)==0)=1;
% southFree(southNghbrs(:,:,:,1)==0)=1;
% northWestFree(northWestNghbrs(:,:,:,1)==0)=1;
% northEastFree(northEastNghbrs(:,:,:,1)==0)=1;
% southWestFree(southWestNghbrs(:,:,:,1)==0)=1;
% southEastFree(southEastNghbrs(:,:,:,1)==0)=1;
% totalFree = westFree+eastFree+northFree+southFree+northEastFree+northWestFree...
%     +southEastFree+southWestFree;

% %Check how many of each type around me:
% westTypes = zeros(n,n,n,5); eastTypes = zeros(n,n,n,5); northTypes = zeros(n,n,n,5); southTypes = zeros(n,n,n,5);
% northWestTypes = zeros(n,n,n,5); northEastTypes = zeros(n,n,n,5); southWestTypes = zeros(n,n,n,5); southEastTypes = zeros(n,n,n,5);
% upTypes = zeros(n,n,n,5); upwestTypes = zeros(n,n,n,5); upeastTypes = zeros(n,n,n,5); upnorthTypes = zeros(n,n,n,5); upsouthTypes = zeros(n,n,n,5);
% upnorthWestTypes = zeros(n,n,n,5); upnorthEastTypes = zeros(n,n,n,5); upsouthWestTypes = zeros(n,n,n,5); upsouthEastTypes = zeros(n,n,n,5);
% downTypes = zeros(n,n,n,5); downwestTypes = zeros(n,n,n,5); downeastTypes = zeros(n,n,n,5); downnorthTypes = zeros(n,n,n,5); downsouthTypes = zeros(n,n,n,5);
% downnorthWestTypes = zeros(n,n,n,5); downnorthEastTypes = zeros(n,n,n,5); downsouthWestTypes = zeros(n,n,n,5); downsouthEastTypes = zeros(n,n,n,5);

for ii=1:numtypes
    westTypes = zeros(n,n,n); eastTypes = zeros(n,n,n); northTypes = zeros(n,n,n); southTypes = zeros(n,n,n);
    northWestTypes = zeros(n,n,n); northEastTypes = zeros(n,n,n); southWestTypes = zeros(n,n,n); southEastTypes = zeros(n,n,n);
    upTypes = zeros(n,n,n); upwestTypes = zeros(n,n,n); upeastTypes = zeros(n,n,n); upnorthTypes = zeros(n,n,n); upsouthTypes = zeros(n,n,n);
    upnorthWestTypes = zeros(n,n,n); upnorthEastTypes = zeros(n,n,n); upsouthWestTypes = zeros(n,n,n); upsouthEastTypes = zeros(n,n,n);
    downTypes = zeros(n,n,n); downwestTypes = zeros(n,n,n); downeastTypes = zeros(n,n,n); downnorthTypes = zeros(n,n,n); downsouthTypes = zeros(n,n,n);
    downnorthWestTypes = zeros(n,n,n); downnorthEastTypes = zeros(n,n,n); downsouthWestTypes = zeros(n,n,n); downsouthEastTypes = zeros(n,n,n);

    if (ii==1||ii==2)
        stateind=1;
        typeind=ii;
    elseif ii==3
        stateind=2;
        typeind=1;
    else
        stateind=3;
        typeind=ii-3;
    end
    
    %check many of type ii around me
    westTypes(westNghbrs(:,:,:,stateind)==typeind)=1;
    eastTypes(eastNghbrs(:,:,:,stateind)==typeind)=1;
    northTypes(northNghbrs(:,:,:,stateind)==typeind)=1;
    southTypes(southNghbrs(:,:,:,stateind)==typeind)=1;
    northWestTypes(northWestNghbrs(:,:,:,stateind)==typeind)=1;
    northEastTypes(northEastNghbrs(:,:,:,stateind)==typeind)=1;
    southWestTypes(southWestNghbrs(:,:,:,stateind)==typeind)=1;
    southEastTypes(southEastNghbrs(:,:,:,stateind)==typeind)=1;
    %in above nbrs:
    upTypes(upNghbrs(:,:,:,stateind)==typeind)=1;
    upwestTypes(upwestNghbrs(:,:,:,stateind)==typeind)=1;
    upeastTypes(upeastNghbrs(:,:,:,stateind)==typeind)=1;
    upnorthTypes(upnorthNghbrs(:,:,:,stateind)==typeind)=1;
    upsouthTypes(upsouthNghbrs(:,:,:,stateind)==typeind)=1;
    upnorthWestTypes(upnorthWestNghbrs(:,:,:,stateind)==typeind)=1;
    upnorthEastTypes(upnorthEastNghbrs(:,:,:,stateind)==typeind)=1;
    upsouthWestTypes(upsouthWestNghbrs(:,:,:,stateind)==typeind)=1;
    upsouthEastTypes(upsouthEastNghbrs(:,:,:,stateind)==typeind)=1;
    %in below nbrs:
    downTypes(downNghbrs(:,:,:,stateind)==typeind)=1;
    downwestTypes(downwestNghbrs(:,:,:,stateind)==typeind)=1;
    downeastTypes(downeastNghbrs(:,:,:,stateind)==typeind)=1;
    downnorthTypes(downnorthNghbrs(:,:,:,stateind)==typeind)=1;
    downsouthTypes(downsouthNghbrs(:,:,:,stateind)==typeind)=1;
    downnorthWestTypes(downnorthWestNghbrs(:,:,:,stateind)==typeind)=1;
    downnorthEastTypes(downnorthEastNghbrs(:,:,:,stateind)==typeind)=1;
    downsouthWestTypes(downsouthWestNghbrs(:,:,:,stateind)==typeind)=1;
    downsouthEastTypes(downsouthEastNghbrs(:,:,:,stateind)==typeind)=1;
    
    totals(:,:,:,ii)=westTypes+eastTypes+northTypes+southTypes...
    +northEastTypes+northWestTypes+southEastTypes+southWestTypes...
    +upTypes+upwestTypes+upeastTypes+upnorthTypes+upsouthTypes...
    +upnorthEastTypes+upnorthWestTypes+upsouthEastTypes+upsouthWestTypes...
    +downTypes+downwestTypes+downeastTypes+downnorthTypes+downsouthTypes...
    +downnorthEastTypes+downnorthWestTypes+downsouthEastTypes+downsouthWestTypes;
end

%totalSusc = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
%     +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
%     +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
%     +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
%     +downTypes(:,:,:,curind)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
%     +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);

    

%Check how many susc around me
% westSusc = zeros(n,n); eastSusc = zeros(n,n); northSusc = zeros(n,n); southSusc = zeros(n,n);
% northWestSusc = zeros(n,n); northEastSusc = zeros(n,n); southWestSusc = zeros(n,n); southEastSusc = zeros(n,n);
%curind=1;   %susc
% westTypes(westNghbrs(:,:,:,1)==1,curind)=1;
% eastTypes(eastNghbrs(:,:,:,1)==1,curind)=1;
% northTypes(northNghbrs(:,:,:,1)==1,curind)=1;
% southTypes(southNghbrs(:,:,:,1)==1,curind)=1;
% northWestTypes(northWestNghbrs(:,:,:,1)==1,curind)=1;
% northEastTypes(northEastNghbrs(:,:,:,1)==1,curind)=1;
% southWestTypes(southWestNghbrs(:,:,:,1)==1,curind)=1;
% southEastTypes(southEastNghbrs(:,:,:,1)==1,curind)=1;
% %in above nbrs:
% upTypes(upNghbrs(:,:,:,1)==1,curind)=1;
% upwestTypes(upwestNghbrs(:,:,:,1)==1,curind)=1;
% upeastTypes(upeastNghbrs(:,:,:,1)==1,curind)=1;
% upnorthTypes(upnorthNghbrs(:,:,:,1)==1,curind)=1;
% upsouthTypes(upsouthNghbrs(:,:,:,1)==1,curind)=1;
% upnorthWestTypes(upnorthWestNghbrs(:,:,:,1)==1,curind)=1;
% upnorthEastTypes(upnorthEastNghbrs(:,:,:,1)==1,curind)=1;
% upsouthWestTypes(upsouthWestNghbrs(:,:,:,1)==1,curind)=1;
% upsouthEastTypes(upsouthEastNghbrs(:,:,:,1)==1,curind)=1;
% %in below nbrs:
% downTypes(downNghbrs(:,:,:,1)==1,curind)=1;
% downwestTypes(downwestNghbrs(:,:,:,1)==1,curind)=1;
% downeastTypes(downeastNghbrs(:,:,:,1)==1,curind)=1;
% downnorthTypes(downnorthNghbrs(:,:,:,1)==1,curind)=1;
% downsouthTypes(downsouthNghbrs(:,:,:,1)==1,curind)=1;
% downnorthWestTypes(downnorthWestNghbrs(:,:,:,1)==1,curind)=1;
% downnorthEastTypes(downnorthEastNghbrs(:,:,:,1)==1,curind)=1;
% downsouthWestTypes(downsouthWestNghbrs(:,:,:,1)==1,curind)=1;
% downsouthEastTypes(downsouthEastNghbrs(:,:,:,1)==1,curind)=1;
% 
% 
% totalSusc = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
% +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
% +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
% +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
% +downTypes(:,:,:,1)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
% +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);

%totalInf = 8-totalFree-totalSusc;

% %check how many inf around me
% curind=2;   %infected
% westTypes(westNghbrs(:,:,:,1)==2,curind)=1;
% eastTypes(eastNghbrs(:,:,:,1)==2,curind)=1;
% northTypes(northNghbrs(:,:,:,1)==2,curind)=1;
% southTypes(southNghbrs(:,:,:,1)==2,curind)=1;
% northWestTypes(northWestNghbrs(:,:,:,1)==2,curind)=1;
% northEastTypes(northEastNghbrs(:,:,:,1)==2,curind)=1;
% southWestTypes(southWestNghbrs(:,:,:,1)==2,curind)=1;
% southEastTypes(southEastNghbrs(:,:,:,1)==2,curind)=1;
% %in above nbrs:
% upTypes(upNghbrs(:,:,:,1)==2,curind)=1;
% upwestTypes(upwestNghbrs(:,:,:,1)==2,curind)=1;
% upeastTypes(upeastNghbrs(:,:,:,1)==2,curind)=1;
% upnorthTypes(upnorthNghbrs(:,:,:,1)==2,curind)=1;
% upsouthTypes(upsouthNghbrs(:,:,:,1)==2,curind)=1;
% upnorthWestTypes(upnorthWestNghbrs(:,:,:,1)==2,curind)=1;
% upnorthEastTypes(upnorthEastNghbrs(:,:,:,1)==2,curind)=1;
% upsouthWestTypes(upsouthWestNghbrs(:,:,:,1)==2,curind)=1;
% upsouthEastTypes(upsouthEastNghbrs(:,:,:,1)==2,curind)=1;
% %in below nbrs:
% downTypes(downNghbrs(:,:,:,1)==2,curind)=1;
% downwestTypes(downwestNghbrs(:,:,:,1)==2,curind)=1;
% downeastTypes(downeastNghbrs(:,:,:,1)==2,curind)=1;
% downnorthTypes(downnorthNghbrs(:,:,:,1)==2,curind)=1;
% downsouthTypes(downsouthNghbrs(:,:,:,1)==2,curind)=1;
% downnorthWestTypes(downnorthWestNghbrs(:,:,:,1)==2,curind)=1;
% downnorthEastTypes(downnorthEastNghbrs(:,:,:,1)==2,curind)=1;
% downsouthWestTypes(downsouthWestNghbrs(:,:,:,1)==2,curind)=1;
% downsouthEastTypes(downsouthEastNghbrs(:,:,:,1)==2,curind)=1;
% 
% 
% totalInf = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
% +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
% +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
% +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
% +downTypes(:,:,:,1)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
% +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);
% 
% 
% %Check innate around me
% 
% curind=3;   %innate
% westTypes(westNghbrs(:,:,:,2)==1,curind)=1;
% eastTypes(eastNghbrs(:,:,:,2)==1,curind)=1;
% northTypes(northNghbrs(:,:,:,2)==1,curind)=1;
% southTypes(southNghbrs(:,:,:,2)==1,curind)=1;
% northWestTypes(northWestNghbrs(:,:,:,2)==1,curind)=1;
% northEastTypes(northEastNghbrs(:,:,:,2)==1,curind)=1;
% southWestTypes(southWestNghbrs(:,:,:,2)==1,curind)=1;
% southEastTypes(southEastNghbrs(:,:,:,2)==1,curind)=1;
% %in above nbrs:
% upTypes(upNghbrs(:,:,:,2)==1,curind)=1;
% upwestTypes(upwestNghbrs(:,:,:,2)==1,curind)=1;
% upeastTypes(upeastNghbrs(:,:,:,2)==1,curind)=1;
% upnorthTypes(upnorthNghbrs(:,:,:,2)==1,curind)=1;
% upsouthTypes(upsouthNghbrs(:,:,:,2)==1,curind)=1;
% upnorthWestTypes(upnorthWestNghbrs(:,:,:,2)==1,curind)=1;
% upnorthEastTypes(upnorthEastNghbrs(:,:,:,2)==1,curind)=1;
% upsouthWestTypes(upsouthWestNghbrs(:,:,:,2)==1,curind)=1;
% upsouthEastTypes(upsouthEastNghbrs(:,:,:,2)==1,curind)=1;
% %in below nbrs:
% downTypes(downNghbrs(:,:,:,2)==1,curind)=1;
% downwestTypes(downwestNghbrs(:,:,:,2)==1,curind)=1;
% downeastTypes(downeastNghbrs(:,:,:,2)==1,curind)=1;
% downnorthTypes(downnorthNghbrs(:,:,:,2)==1,curind)=1;
% downsouthTypes(downsouthNghbrs(:,:,:,2)==1,curind)=1;
% downnorthWestTypes(downnorthWestNghbrs(:,:,:,2)==1,curind)=1;
% downnorthEastTypes(downnorthEastNghbrs(:,:,:,2)==1,curind)=1;
% downsouthWestTypes(downsouthWestNghbrs(:,:,:,2)==1,curind)=1;
% downsouthEastTypes(downsouthEastNghbrs(:,:,:,2)==1,curind)=1;
% 
% 
% totalInn = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
% +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
% +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
% +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
% +downTypes(:,:,:,1)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
% +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);
% 
% 
% %Check antitumor adapative around me
% 
% curind=4;   %adap tumor
% westTypes(westNghbrs(:,:,:,3)==1,curind)=1;
% eastTypes(eastNghbrs(:,:,:,3)==1,curind)=1;
% northTypes(northNghbrs(:,:,:,3)==1,curind)=1;
% southTypes(southNghbrs(:,:,:,3)==1,curind)=1;
% northWestTypes(northWestNghbrs(:,:,:,3)==1,curind)=1;
% northEastTypes(northEastNghbrs(:,:,:,3)==1,curind)=1;
% southWestTypes(southWestNghbrs(:,:,:,3)==1,curind)=1;
% southEastTypes(southEastNghbrs(:,:,:,3)==1,curind)=1;
% %in above nbrs:
% upTypes(upNghbrs(:,:,:,3)==1,curind)=1;
% upwestTypes(upwestNghbrs(:,:,:,3)==1,curind)=1;
% upeastTypes(upeastNghbrs(:,:,:,3)==1,curind)=1;
% upnorthTypes(upnorthNghbrs(:,:,:,3)==1,curind)=1;
% upsouthTypes(upsouthNghbrs(:,:,:,3)==1,curind)=1;
% upnorthWestTypes(upnorthWestNghbrs(:,:,:,3)==1,curind)=1;
% upnorthEastTypes(upnorthEastNghbrs(:,:,:,3)==1,curind)=1;
% upsouthWestTypes(upsouthWestNghbrs(:,:,:,3)==1,curind)=1;
% upsouthEastTypes(upsouthEastNghbrs(:,:,:,3)==1,curind)=1;
% %in below nbrs:
% downTypes(downNghbrs(:,:,:,3)==1,curind)=1;
% downwestTypes(downwestNghbrs(:,:,:,3)==1,curind)=1;
% downeastTypes(downeastNghbrs(:,:,:,3)==1,curind)=1;
% downnorthTypes(downnorthNghbrs(:,:,:,3)==1,curind)=1;
% downsouthTypes(downsouthNghbrs(:,:,:,3)==1,curind)=1;
% downnorthWestTypes(downnorthWestNghbrs(:,:,:,3)==1,curind)=1;
% downnorthEastTypes(downnorthEastNghbrs(:,:,:,3)==1,curind)=1;
% downsouthWestTypes(downsouthWestNghbrs(:,:,:,3)==1,curind)=1;
% downsouthEastTypes(downsouthEastNghbrs(:,:,:,3)==1,curind)=1;
% 
% 
% totalAT = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
% +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
% +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
% +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
% +downTypes(:,:,:,1)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
% +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);
% 
% 
% %Check antiviral adapative around me
% 
% curind=5;   %adap viral
% westTypes(westNghbrs(:,:,:,3)==2,curind)=1;
% eastTypes(eastNghbrs(:,:,:,3)==2,curind)=1;
% northTypes(northNghbrs(:,:,:,3)==2,curind)=1;
% southTypes(southNghbrs(:,:,:,3)==2,curind)=1;
% northWestTypes(northWestNghbrs(:,:,:,3)==2,curind)=1;
% northEastTypes(northEastNghbrs(:,:,:,3)==2,curind)=1;
% southWestTypes(southWestNghbrs(:,:,:,3)==2,curind)=1;
% southEastTypes(southEastNghbrs(:,:,:,3)==2,curind)=1;
% %in above nbrs:
% upTypes(upNghbrs(:,:,:,3)==2,curind)=1;
% upwestTypes(upwestNghbrs(:,:,:,3)==2,curind)=1;
% upeastTypes(upeastNghbrs(:,:,:,3)==2,curind)=1;
% upnorthTypes(upnorthNghbrs(:,:,:,3)==2,curind)=1;
% upsouthTypes(upsouthNghbrs(:,:,:,3)==2,curind)=1;
% upnorthWestTypes(upnorthWestNghbrs(:,:,:,3)==2,curind)=1;
% upnorthEastTypes(upnorthEastNghbrs(:,:,:,3)==2,curind)=1;
% upsouthWestTypes(upsouthWestNghbrs(:,:,:,3)==2,curind)=1;
% upsouthEastTypes(upsouthEastNghbrs(:,:,:,3)==2,curind)=1;
% %in below nbrs:
% downTypes(downNghbrs(:,:,:,3)==2,curind)=1;
% downwestTypes(downwestNghbrs(:,:,:,3)==2,curind)=1;
% downeastTypes(downeastNghbrs(:,:,:,3)==2,curind)=1;
% downnorthTypes(downnorthNghbrs(:,:,:,3)==2,curind)=1;
% downsouthTypes(downsouthNghbrs(:,:,:,3)==2,curind)=1;
% downnorthWestTypes(downnorthWestNghbrs(:,:,:,3)==2,curind)=1;
% downnorthEastTypes(downnorthEastNghbrs(:,:,:,3)==2,curind)=1;
% downsouthWestTypes(downsouthWestNghbrs(:,:,:,3)==2,curind)=1;
% downsouthEastTypes(downsouthEastNghbrs(:,:,:,3)==2,curind)=1;
% 
% 
% totalAV = westTypes(:,:,:,curind)+eastTypes(:,:,:,curind)+northTypes(:,:,:,curind)+southTypes(:,:,:,curind)...
% +northEastTypes(:,:,:,curind)+northWestTypes(:,:,:,curind)+southEastTypes(:,:,:,curind)+southWestTypes(:,:,:,curind)...
% +upTypes(:,:,:,curind)+upwestTypes(:,:,:,curind)+upeastTypes(:,:,:,curind)+upnorthTypes(:,:,:,curind)+upsouthTypes(:,:,:,curind)...
% +upnorthEastTypes(:,:,:,curind)+upnorthWestTypes(:,:,:,curind)+upsouthEastTypes(:,:,:,curind)+upsouthWestTypes(:,:,:,curind)...
% +downTypes(:,:,:,1)+downwestTypes(:,:,:,curind)+downeastTypes(:,:,:,curind)+downnorthTypes(:,:,:,curind)+downsouthTypes(:,:,:,curind)...
% +downnorthEastTypes(:,:,:,curind)+downnorthWestTypes(:,:,:,curind)+downsouthEastTypes(:,:,:,curind)+downsouthWestTypes(:,:,:,curind);
% 
% 
% 
% totals(:,:,:,1)=totalSusc;
% totals(:,:,:,2)=totalInf;
% totals(:,:,:,3)=totalInn;
% totals(:,:,:,4)=totalAT;
% totals(:,:,:,5)=totalAV;

end

