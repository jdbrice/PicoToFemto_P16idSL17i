#ifndef FEMTO_DST_WRITER_H
#define FEMTO_DST_WRITER_H

#include "PicoDstSkimmer/PicoDstSkimmer.h"

#include "ProductionUtils/RunMapFactory.h"

#include "PicoDstP16id/StPicoMtdPidTraits.h"

#include "FemtoDstFormat/BranchWriter.h"
#include "FemtoDstFormat/TClonesArrayWriter.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoBTofPidTraits.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"

#include "TTree.h"
#include "TAxis.h"
#include "Compression.h"

#include <set>

#include "vendor/loguru.h"

class FemtoDstWriter : public PicoDstSkimmer
{
protected:
	FemtoEvent _event;
	RunMapFactory rmf;

	TTree * tree;
	BranchWriter<FemtoEvent> _wEvent;
	TClonesArrayWriter<FemtoTrack> _wTrack;
	TClonesArrayWriter<FemtoMtdPidTraits> _wMtdPid;
	TClonesArrayWriter<FemtoBTofPidTraits> _wBTofPid;
	TClonesArrayWriter<FemtoTrackHelix>    _wHelix;
	FemtoTrack _track;
	FemtoMtdPidTraits _mtdPid;
	FemtoBTofPidTraits _btofPid;
	FemtoTrackHelix _helix;

	bool require_mtdPid = false;
	bool require_btofPid = false;

	double max_pT = 1000;

	std::set<unsigned int> triggerIds;

public:
	virtual const char* classname() const { return "FemtoDstWriter"; }
	FemtoDstWriter() {}
	~FemtoDstWriter() {}

	virtual void initialize(){
		PicoDstSkimmer::initialize();

		TFile * f = book->getOutputFile();
		f->SetCompressionLevel( 9 );

		book->cd();
		tree = new TTree( "FemtoDst", "FemtoDst" );
		_wEvent.createBranch( tree, "Event" );
		_wTrack.createBranch( tree, "Tracks" );
		_wMtdPid.createBranch( tree, "MtdPidTraits" );
		_wBTofPid.createBranch( tree, "BTofPidTraits" );
		_wHelix.createBranch( tree, "Helices" );


		require_mtdPid = config.getBool( nodePath + ".Require:mtd", false );
		require_btofPid = config.getBool( nodePath + ".Require:btof", false );
		max_pT = config.getDouble( nodePath + ".Require:max_pT", max_pT );

	}
protected:	

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

	}
	virtual bool keepEvent( StPicoEvent *event ){

		book->fill( "events", "All" );
		if ( makeTriggerWord( event ) == 0)
			return false;

		book->fill( "events", "Trigger" );

		// TODO more bad run rejection
		// LOG_F( INFO, "runId: %d --> day: %d", event->runId(), RunMapFactory::day( event->runId() ) );
		if ( RunMapFactory::day( event->runId() ) <= 50 )
			return false;
		book->fill( "events", "bad" );

		auto vtx = event->primaryVertex();
		float deltaVz = event->vzVpd() - vtx.z();

		book->fill( "vtx_z", vtx.z() );
		book->fill( "delta_z", deltaVz );
		book->fill( "ranking", event->ranking() * 1e-6 );
		book->fill( "mtd_mult", _rMtdPid.N() );

		if ( fabs( vtx.z() ) > 100 )
			return false;
		book->fill( "events", "vtx" );

		if ( fabs(deltaVz) > 6.0 )
			return false;
		book->fill( "events", "vtx_delta" );

		
		// Ranking cut
		if ( event->ranking() < 0 )
			return false;
		book->fill( "events", "ranking" );


		// require at least 1 MTD match
		if ( _rMtdPid.N() < 1 )
			return false;
		book->fill( "events", "mtd_match" );

		book->fill( "pass_mtd_mult", _rMtdPid.N() );
		book->fill( "pass_vtx_z", vtx.z() );
		book->fill( "pass_delta_z", deltaVz );
		book->fill( "pass_ranking", event->ranking() * 1e-6 );

		return true;
	}

	virtual void analyzeEvent(){
		StPicoEvent *event = _rEvent.get( 0 );
		if ( nullptr == event ){
			LOG_F( ERROR, "NULL EVENT" );
		}

		if ( false == keepEvent( event ) )
			return;

		
		auto vtx = event->primaryVertex();
		float deltaVz = event->vzVpd() - vtx.z();

		_event.mPrimaryVertex_mX1 = vtx.x();
		_event.mPrimaryVertex_mX2 = vtx.y();
		_event.mPrimaryVertex_mX3 = vtx.z();

		_event.mRunId       = event->runId();
		_event.mRunIndex    = rmf.indexForRun( _event.mRunId );
		_event.mEventId     = event->eventId();
		_event.mGRefMult    = event->grefMult();
		_event.mTriggerWord = makeTriggerWord( event );
		_event.mWeight      = deltaVz;
		LOG_F( 9, "RunId: %d --> %d", _event.mRunId, rmf.indexForRun( _event.mRunId ) );

		_wEvent.set( _event );

		_wTrack.reset();
		_wMtdPid.reset();
		_wBTofPid.reset();
		_wHelix.reset();
		size_t nTracks = _rTrack.N();
		// LOG_F( INFO, "N TRACKS = %lu", nTracks );;
		size_t nMtdTracks = 0;
		size_t nBTofTracks = 0;
		size_t nTracksKeep = 0;
		for ( size_t i = 0; i < nTracks; i++ ){
			
			StPicoTrack * track = _rTrack.get( i );

			fillTrack( nTracksKeep, track, event );
			
			if ( fabs(_track.mPt) > 25 ){
				LOG_F( INFO, "pT=%f", _track.mPt );
			}
			
			if ( fabs(_track.mPt) > 0.01 ){
				nTracksKeep++;

				_wTrack.add( _track );
				_wHelix.add( _helix );

				if ( _track.mMtdPidTraitsIndex >= 0 ){
					_wMtdPid.add( _mtdPid );
					nMtdTracks++;
				}

				if ( _track.mBTofPidTraitsIndex >= 0 ){
					nBTofTracks++;
					_wBTofPid.add( _btofPid );
				}
			
			}// mPt > 0.01
		} // loop on tracks

		if ( nMtdTracks >= 1 ){
			book->fill( "events", "tree" );
			tree->Fill();
		}
	}

	virtual void fillTrack( size_t idx, StPicoTrack *track, StPicoEvent *event ){
		// LOG_F( INFO, "Filling %lu", i );
		_track.reset();

		_track.mId  = track->id();
		auto pMom   = track->pMom();
		_track.mPt  = pMom.perp();
		_track.mEta = pMom.pseudoRapidity();
		_track.mPhi = pMom.phi();
		
		_track.mNHitsFit  = track->nHitsFit() * track->charge();
		_track.mNHitsMax  = track->nHitsMax();
		_track.mNHitsDedx = track->nHitsDedx();

		_track.nSigmaElectron ( track->nSigmaElectron ( ) );
		_track.nSigmaPion     ( track->nSigmaPion     ( ) );
		_track.nSigmaKaon     ( track->nSigmaKaon     ( ) );
		_track.nSigmaProton   ( track->nSigmaProton   ( ) );

		double dca_full_prec = track->globalDCA( event->bField(), event->primaryVertex() );
		_track.gDCA( dca_full_prec );

		_track.dEdx( track->dEdx() );

		fillTrackHelix( idx, track, dca_full_prec );

		if ( track->mtdPidTraitsIndex() >= 0 ){
			fillMtdPid( idx, track );
		}

		if ( track->bTofPidTraitsIndex() >= 0 ){
			fillBTofPid( idx, track );
		}
		
	}


	virtual void fillMtdPid( size_t idx, StPicoTrack *track ){
		_mtdPid.reset();

		size_t index = _wMtdPid.N();
		StPicoMtdPidTraits *mtdPid = _rMtdPid.get( track->mtdPidTraitsIndex() );
		_mtdPid.mMatchFlag = mtdPid->matchFlag();
		_mtdPid.mDeltaY = mtdPid->deltaY();
		_mtdPid.mDeltaZ = mtdPid->deltaZ();
		_mtdPid.mDeltaTimeOfFlight = mtdPid->deltaTimeOfFlight();
		_mtdPid.mMtdHitChan = mtdPid->mMtdHitChan;


		size_t nMtdHits = _rMtdHit.N();
		for ( size_t iMtdHit = 0; iMtdHit < nMtdHits; iMtdHit++ ){
			StPicoMtdHit *mtdHit = _rMtdHit.get( iMtdHit );
			
			if ( nullptr != mtdHit && mtdHit->gChannel() == _mtdPid.mMtdHitChan ){
				_mtdPid.mTriggerFlag = mtdHit->triggerFlag();
				break;
			}
		} // loop iMtdHit

		_track.mMtdPidTraitsIndex = index;

	}


	virtual void fillBTofPid( size_t idx, StPicoTrack *track ){
		_btofPid.reset();
		StPicoBTofPidTraits *btofPid = _rBTofPid.get( track->bTofPidTraitsIndex() );

		_btofPid.beta( btofPid->btofBeta() );
		_btofPid.matchFlag( btofPid->btofMatchFlag() );
		_btofPid.yLocal( btofPid->btofYLocal() );
		_btofPid.zLocal( btofPid->btofZLocal() );

		size_t index = _wBTofPid.N();
		_track.mBTofPidTraitsIndex = index;

	}

	virtual void fillTrackHelix( size_t i, StPicoTrack *track, double dca ){
		_helix.reset();

		_helix.mPar[0] = track->gMom().x();
		_helix.mPar[1] = track->gMom().y();
		_helix.mPar[2] = track->gMom().z();
		_helix.mPar[3] = track->origin().x();
		_helix.mPar[4] = track->origin().y();
		_helix.mPar[5] = track->origin().z();
		_helix.mMap0 = track->topologyMap(0);
		_helix.mMap1 = track->topologyMap(1);
		_helix.mDCA = dca;

		size_t index = _wHelix.N();
		_track.mHelixIndex = index;
	}


	virtual void postEventLoop(){
		
		

	}


	unsigned int makeTriggerWord( StPicoEvent * event ){
		unsigned int word = 0;

		// dimuon
		if ( event->isTrigger( 470602 ) || event->isTrigger( 480602 ) || event->isTrigger(490602) ){
			word = 1 << 0;
		}

		// // singlemuon-5
		// if ( event->isTrigger( 470600 ) || event->isTrigger( 480600 ) || event->isTrigger(490600) ){
		// 	word = 1 << 1;
		// }

		// // emuon-30
		// if ( event->isTrigger( 470601 ) || event->isTrigger( 480601 ) || event->isTrigger(490601) ){
		// 	word = 1 << 2;
		// }

		return word;
	}

};


#endif