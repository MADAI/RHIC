#ifndef __ACTION_PERFORM_REGEN_CC__
#define __ACTION_PERFORM_REGEN_CC__
#include "b3d.h"

void CAction::PerformMuTCalc(){
	int ix,iy,ieta;
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			for(ieta=0;ieta<2*b3d->NETA;ieta++)
				b3d->muTinfo[ix][iy]->UpdateNMPE(b3d->cell[ix][iy][ieta]);
		}
	}
	// only uses x, y coordinates of cell for particle placement
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			if(b3d->regen->CheckForRegeneration(b3d->cell[ix][iy][0],b3d->muTinfo[ix][iy]))
				b3d->nregenerate+=1;	
		}
	}
}

#endif
