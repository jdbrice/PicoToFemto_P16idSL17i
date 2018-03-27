#ifndef FEMTO_TRACK_PROXY_H
#define FEMTO_TRACK_PROXY_H

#include "FemtoTrack.h"
#include "FemtoMcTrack.h"
#include "FemtoMtdPidTraits.h"
#include "TClonesArrayReader.h"



class FemtoTrackProxy {
public:
	FemtoTrack        * _track   = nullptr;
	FemtoMcTrack      * _mcTrack = nullptr;
	FemtoMtdPidTraits * _mtdPid  = nullptr;

	// extras
	bool _fMcTrack = false;
	bool _fTrack   = false;
	bool _fMtdPid  = false;
	float _mass    = 0.0;	// mass assumption

	void setMtdPidTraits( TClonesArrayReader<FemtoMtdPidTraits> &_rMtd ){
		if ( nullptr != this->_track && this->_track->mMtdPidTraitsIndex >= 0 )
			this->_mtdPid = _rMtd.get( this->_track->mMtdPidTraitsIndex );
		else 
			this->_mtdPid = nullptr;
	}

	void assemble( uint i,
					TClonesArrayReader<FemtoTrack> &_rTracks,
					TClonesArrayReader<FemtoMtdPidTraits> &_rMtd ) {
		this->_mcTrack = nullptr;
		this->_track = _rTracks.get( i );
		setMtdPidTraits( _rMtd );

	}

	void assemble( uint i, 
					TClonesArrayReader<FemtoMcTrack> &_rMc, 
					TClonesArrayReader<FemtoTrack> &_rTracks, 
					TClonesArrayReader<FemtoMtdPidTraits> &_rMtd ){

		this->_mcTrack   = _rMc.get( i );
		if ( nullptr == this->_mcTrack ) return;
		
		if ( this->_mcTrack->mAssociatedIndex >= 0 )
			this->_track = _rTracks.get( this->_mcTrack->mAssociatedIndex );
		else 
			this->_track = nullptr;

		setMtdPidTraits( _rMtd );
	}
};

#endif