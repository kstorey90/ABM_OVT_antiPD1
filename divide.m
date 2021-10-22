function [state,divisionCounter,divisionFlag] = divide(state,divisionCounter,divisionFlag,row,col,ht,n,Sm,Ssd)
%DIVIDE NO RESTRICTION 3D Produces a daughter cell if space is available. 
%If no space available, randomly choose neighboring space to push outward
%and continue shifting until you reach open site

currState=state(row,col,ht,1); 
currDivisionCounter=initializeCycleDur(Sm,Ssd); 
currDivisionFlag=0;

emptysite=0;
curi=row;
curj=col;
curk=ht;


nbrChoices=[0 1 0; 0 -1 0; 1 0 0; -1 0 0; 1 1 0; 1 -1 0; -1 -1 0; -1 1 0;...
    0 0 1; 0 1 1; 0 -1 1; 1 0 1; -1 0 1; 1 1 1; 1 -1 1; -1 -1 1; -1 1 1;...
    0 0 -1; 0 1 -1; 0 -1 -1; 1 0 -1; -1 0 -1; 1 1 -1; 1 -1 -1; -1 -1 -1; -1 1 -1];
numNbrs=26;


while ~emptysite
    emptyLoc = checkTumorNeighborsLoc(state,curi,curj,curk);
    numFree=size(emptyLoc,1);
    if numFree>0
        % randomly choose one of the free spots to place daughter cell
        newSite=ceil(rand*numFree);
        %newSite=randsample(numFree,1);
        curi=curi+emptyLoc(newSite,1);
        curj=curj+emptyLoc(newSite,2);
        curk=curk+emptyLoc(newSite,3);
        state(curi,curj,curk,1)=currState;
        divisionCounter(curi,curj,curk)=currDivisionCounter;
        divisionFlag(curi,curj,curk)=currDivisionFlag;
        emptysite=1;
    else
        %randomly choose a neighboring spot, and then keep going
        
        %check whether site is on border:
        if (curi==1 || curi==n || curj==1 || curj==n || curk==1 || curk==n)
            emptysite=1;
        else    
            newSite=ceil(rand*numNbrs);

            curi=curi+nbrChoices(newSite,1);
            curj=curj+nbrChoices(newSite,2);
            curk=curk+nbrChoices(newSite,3);
            
            %continue pushing in the same direction
            nbrChoices=nbrChoices(newSite,:);
            numNbrs=1;

            newState=state(curi,curj,curk,1);
            newDivCtr=divisionCounter(curi,curj,curk);
            newDivFlag=divisionFlag(curi,curj,curk);

            state(curi,curj,curk,1)=currState;
            divisionCounter(curi,curj,curk)=currDivisionCounter;
            divisionFlag(curi,curj,curk)=currDivisionFlag;

            currState=newState;
            currDivisionCounter=newDivCtr;
            currDivisionFlag=newDivFlag;
        
    end
    
end

% while (state(row-1,col,1)~=0 && state(row-1,col-1,1)~=0 && state(row-1,col+1,1)~=0 &&...
%         state(row,col-1,1)~=0 && state(row,col+1,1)~=0 &&...
%         state(row+1,col-1,1)~=0 && state(row+1,col,1)~=0 && state(row+1,col+1,1)~=0)
%     
%     % Look in all 8 directions and find the direction with lowest distance
%     % to the spheroid surface.
%     % North...
%     i=row; j=col; northDistance=0;
%     while state(i-1,j,1)~=0
%         i=i-1; northDistance=northDistance+1;
%     end
%     % East...
%     i=row; j=col; eastDistance=0;
%     while state(i,j+1,1)~=0
%         j=j+1; eastDistance=eastDistance+1;
%     end
%     % South...
%     i=row; j=col; southDistance=0;
%     while state(i+1,j,1)~=0
%         i=i+1; southDistance=southDistance+1;
%     end
%     % West...
%     i=row; j=col; westDistance=0;
%     while state(i,j-1,1)~=0
%         j=j-1; westDistance=westDistance+1;
%     end
%     % North-east...
%     i=row; j=col; northEastDistance=0;
%     while state(i-1,j+1,1)~=0
%         i=i-1; j=j+1; northEastDistance=northEastDistance+1;
%     end
%     % North-west...
%     i=row; j=col; northWestDistance=0;
%     while state(i-1,j-1,1)~=0
%         i=i-1; j=j-1; northWestDistance=northWestDistance+1;
%     end
%     % South-east
%     i=row; j=col; southEastDistance=0;
%     while state(i+1,j+1,1)~=0
%         i=i+1; j=j+1; southEastDistance=southEastDistance+1;
%     end
%     % South-west...
%     i=row; j=col; southWestDistance=0;
%     while state(i+1,j-1,1)~=0
%         i=i+1; j=j-1; southWestDistance=southWestDistance+1;
%     end
%     
%     % Push cell in the biased random direction
%     distances=[northDistance eastDistance southDistance westDistance...
%         northEastDistance northWestDistance southEastDistance southWestDistance];
%     tmp=zeros(1,8);
%     sortedDistances=sort(distances);
%     
%     tmp(distances==sortedDistances(8))=rand(size(tmp(distances==sortedDistances(8))));
%     tmp(distances==sortedDistances(7))=rand(size(tmp(distances==sortedDistances(7))));
%     tmp(distances==sortedDistances(6))=rand(size(tmp(distances==sortedDistances(6))));
%     tmp(distances==sortedDistances(5))=rand(size(tmp(distances==sortedDistances(5))))+0.6;
%     tmp(distances==sortedDistances(4))=rand(size(tmp(distances==sortedDistances(4))))+0.7;
%     tmp(distances==sortedDistances(3))=rand(size(tmp(distances==sortedDistances(3))))+0.8;
%     tmp(distances==sortedDistances(2))=rand(size(tmp(distances==sortedDistances(2))))+0.9;
%     tmp(distances==sortedDistances(1))=rand(size(tmp(distances==sortedDistances(1))))+1;
%     
%     northDistance=tmp(1);
%     eastDistance=tmp(2);
%     southDistance=tmp(3);
%     westDistance=tmp(4);
%     northEastDistance=tmp(5);
%     northWestDistance=tmp(6);
%     southEastDistance=tmp(7);
%     southWestDistance=tmp(8);
%     
%     %Shift cell along the lowest distance to nearest empty cell
%     if northDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                   = state(row-1,col,1);
%         tmpDivisionCounter         = divisionCounter(row-1,col);
%         tmpDivisionFlag            = divisionFlag(row-1,col);
%         state(row-1,col,1)           = currState;
%         divisionCounter(row-1,col) = currDivisionCounter;
%         divisionFlag(row-1,col)    = currDivisionFlag;
%         currState                  = tmpState;
%         currDivisionCounter        = tmpDivisionCounter;
%         currDivisionFlag           = tmpDivisionFlag;
%         row                        = row-1;
%     elseif eastDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                   = state(row,col+1,1);
%         tmpDivisionCounter         = divisionCounter(row,col+1);
%         tmpDivisionFlag            = divisionFlag(row,col+1);
%         state(row,col+1,1)           = currState;
%         divisionCounter(row,col+1) = currDivisionCounter;
%         divisionFlag(row,col+1)    = currDivisionFlag;
%         currState                  = tmpState;
%         currDivisionCounter        = tmpDivisionCounter;
%         currDivisionFlag           = tmpDivisionFlag;
%         col                        = col+1;
%     elseif southDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                   = state(row+1,col,1);
%         tmpDivisionCounter         = divisionCounter(row+1,col);
%         tmpDivisionFlag            = divisionFlag(row+1,col);
%         state(row+1,col,1)           = currState;
%         divisionCounter(row+1,col) = currDivisionCounter;
%         divisionFlag(row+1,col)    = currDivisionFlag;
%         currState                  = tmpState;
%         currDivisionCounter        = tmpDivisionCounter;
%         currDivisionFlag           = tmpDivisionFlag;
%         row                        = row+1;
%     elseif westDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                   = state(row,col-1,1);
%         tmpDivisionCounter         = divisionCounter(row,col-1);
%         tmpDivisionFlag            = divisionFlag(row,col-1);
%         state(row,col-1,1)           = currState;
%         divisionCounter(row,col-1) = currDivisionCounter;
%         divisionFlag(row,col-1)    = currDivisionFlag;
%         currState                  = tmpState;
%         currDivisionCounter        = tmpDivisionCounter;
%         currDivisionFlag           = tmpDivisionFlag;
%         col                        = col-1;
%     elseif northEastDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                     = state(row-1,col+1,1);
%         tmpDivisionCounter           = divisionCounter(row-1,col+1);
%         tmpDivisionFlag              = divisionFlag(row-1,col+1);
%         state(row-1,col+1,1)           = currState;
%         divisionCounter(row-1,col+1) = currDivisionCounter;
%         divisionFlag(row-1,col+1)    = currDivisionFlag;
%         currState                    = tmpState;
%         currDivisionCounter          = tmpDivisionCounter;
%         currDivisionFlag             = tmpDivisionFlag;
%         row                          = row-1;
%         col                          = col+1;
%     elseif northWestDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                     = state(row-1,col-1,1);
%         tmpDivisionCounter           = divisionCounter(row-1,col-1);
%         tmpDivisionFlag              = divisionFlag(row-1,col-1);
%         state(row-1,col-1,1)           = currState;
%         divisionCounter(row-1,col-1) = currDivisionCounter;
%         divisionFlag(row-1,col-1)    = currDivisionFlag;
%         currState                    = tmpState;
%         currDivisionCounter          = tmpDivisionCounter;
%         currDivisionFlag             = tmpDivisionFlag;
%         row                          = row-1;
%         col                          = col-1;
%     elseif southEastDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                     = state(row+1,col+1,1);
%         tmpDivisionCounter           = divisionCounter(row+1,col+1);
%         tmpDivisionFlag              = divisionFlag(row+1,col+1);
%         state(row+1,col+1,1)           = currState;
%         divisionCounter(row+1,col+1) = currDivisionCounter;
%         divisionFlag(row+1,col+1)    = currDivisionFlag;
%         currState                    = tmpState;
%         currDivisionCounter          = tmpDivisionCounter;
%         currDivisionFlag             = tmpDivisionFlag;
%         row                          = row+1;
%         col                          = col+1;
%     elseif southWestDistance==max([northDistance eastDistance southDistance westDistance...
%             northEastDistance northWestDistance southEastDistance southWestDistance])
%         tmpState                     = state(row+1,col-1,1);
%         tmpDivisionCounter           = divisionCounter(row+1,col-1);
%         tmpDivisionFlag              = divisionFlag(row+1,col-1);
%         state(row+1,col-1,1)           = currState;
%         divisionCounter(row+1,col-1) = currDivisionCounter;
%         divisionFlag(row+1,col-1)    = currDivisionFlag;
%         currState                    = tmpState;
%         currDivisionCounter          = tmpDivisionCounter;
%         currDivisionFlag             = tmpDivisionFlag;
%         row                          = row+1;
%         col                          = col-1;
%     end
% end
% 
% % Now that we know there is at least one free space in the Moore neighborhood
% % we assign random numbers to the empty ones
% nRnd=0; sRnd=0; wRnd=0; eRnd=0; nwRnd=0; neRnd=0; seRnd=0; swRnd=0;
% totalFree = checkTumorNeighbors(state,n);
% %nbrs in cur layer
% n=state(row-1,col,ht,1); s=state(row+1,col,ht,1); w=state(row,col-1,ht,1); e=state(row,col+1,ht,1);
% nw=state(row-1,col-1,ht,1); ne=state(row-1,col+1,ht,1); se=state(row+1,col+1,ht,1); sw=state(row+1,col-1,ht,1);
% %%above nbrs:
% u=state(row,col,ht+1,1); un=state(row-1,col,ht+1,1); us=state(row+1,col,ht+1,1); uw=state(row,col-1,ht+1,1); ue=state(row,col+1,ht+1,1);
% unw=state(row-1,col-1,ht+1,1); une=state(row-1,col+1,ht+1,1); use=state(row+1,col+1,ht+1,1); usw=state(row+1,col-1,ht+1,1);
% %%below nbrs:
% d=state(row,col,ht-1,1); dn=state(row-1,col,ht-1,1); ds=state(row+1,col,ht-1,1); dw=state(row,col-1,ht-1,1); de=state(row,col+1,ht-1,1);
% dnw=state(row-1,col-1,ht-1,1); dne=state(row-1,col+1,ht-1,1); dse=state(row+1,col+1,ht-1,1); dsw=state(row+1,col-1,ht-1,1);
% 
% nRnd(n==0)   = rand+(26-totalFree(row-1,col,ht));
% sRnd(s==0)   = rand+(26-totalFree(row+1,col,ht));
% wRnd(w==0)   = rand+(26-totalFree(row,col-1,ht));
% eRnd(e==0)   = rand+(26-totalFree(row,col+1,ht));
% nwRnd(nw==0) = rand+(26-totalFree(row-1,col-1,ht));
% neRnd(ne==0) = rand+(26-totalFree(row-1,col+1,ht));
% seRnd(se==0) = rand+(26-totalFree(row+1,col+1,ht));
% swRnd(sw==0) = rand+(26-totalFree(row+1,col-1,ht));
% %FINISH ADDING HERE (FIGURE OUT WHAT THESE ARE USED FOR)
% 
% % nRnd(n=='E')   = rand+1;
% % sRnd(s=='E')   = rand+1;
% % wRnd(w=='E')   = rand+1;
% % eRnd(e=='E')   = rand+1;
% % nwRnd(nw=='E') = rand;
% % neRnd(ne=='E') = rand;
% % seRnd(se=='E') = rand;
% % swRnd(sw=='E') = rand;
% 
% if nRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row-1,col,1)=currState;
%     divisionCounter(row-1,col)=currDivisionCounter;
%     divisionFlag(row-1,col)=currDivisionFlag;
% elseif sRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row+1,col,1)=currState;
%     divisionCounter(row+1,col)=currDivisionCounter;
%     divisionFlag(row+1,col)=currDivisionFlag;
% elseif wRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row,col-1,1)=currState;
%     divisionCounter(row,col-1)=currDivisionCounter;
%     divisionFlag(row,col-1)=currDivisionFlag;
% elseif eRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row,col+1,1)=currState;
%     divisionCounter(row,col+1)=currDivisionCounter;
%     divisionFlag(row,col+1)=currDivisionFlag;
% elseif nwRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row-1,col-1,1)=currState;
%     divisionCounter(row-1,col-1)=currDivisionCounter;
%     divisionFlag(row-1,col-1)=currDivisionFlag;
% elseif neRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row-1,col+1,1)=currState;
%     divisionCounter(row-1,col+1)=currDivisionCounter;
%     divisionFlag(row-1,col+1)=currDivisionFlag;
% elseif seRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row+1,col+1,1)=currState;
%     divisionCounter(row+1,col+1)=currDivisionCounter;
%     divisionFlag(row+1,col+1)=currDivisionFlag;
% elseif swRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
%     state(row+1,col-1,1)=currState;
%     divisionCounter(row+1,col-1)=currDivisionCounter;
%     divisionFlag(row+1,col-1)=currDivisionFlag;
% end

end


