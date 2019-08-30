// -*- C++ -*-
//
// Package:    Producer/DiamondTimingProducer
// Class:      DiamondTimingProducer
//
/**\class DiamondTimingProducer DiamondTimingProducer.cc Producer/DiamondTimingProducer/plugins/DiamondTimingProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piotr Maciej Cwiklicki
//         Created:  Thu, 22 Aug 2019 13:48:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include <vector> // ##
#include "DataFormats/Math/interface/Point3D.h" // ##

enum Station_id_
{
	STATION_210_M_ID,
	STATION_TIMING_ID,
	STATION_220_M_ID
};

enum Sector_id_
{
  SECTOR_45_ID,
  SECTOR_56_ID
};
//
// class declaration
//

class DiamondTimingProducer : public edm::stream::EDProducer<> {
   public:
      explicit DiamondTimingProducer(const edm::ParameterSet&);
      ~DiamondTimingProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      typedef std::vector<bool> TBA;
      typedef std::map< std::pair< int , int >, int> Pmux_map;
			//typedef edm::DetSetVector<CTPPSDiamondRecHit> RecHit;
			typedef edm::DetSetVector<CTPPSPixelLocalTrack> PixelLocal;
			//typedef edm::DetSetVector<CTPPSDiamondDigi> DiamondDigi;
			//typedef edm::DetSetVector<CTPPSDiamondLocalTrack> DiamondLocalTracks;

			/*edm::EDGetTokenT< DiamondDigi > tokenDigi_;
      edm::EDGetTokenT< RecHit > tokenRecHit_;
			edm::EDGetTokenT< DiamondLocalTracks > tokenLocalTrack_;*/
      edm::EDGetTokenT< PixelLocal > tokenPixelLocalTrack_;

      std::map< std::pair< int , int >, std::pair< int , int > > Ntracks_cuts_map_;
      Pmux_map Pixel_Mux_map_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
DiamondTimingProducer:: DiamondTimingProducer(const edm::ParameterSet& iConfig)
  :
	//tokenDigi_            ( consumes< edm::DetSetVector<CTPPSDiamondDigi> >      ( iConfig.getParameter<edm::InputTag>( "tagDigi" ) ) ),
  //tokenRecHit_          ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >      ( iConfig.getParameter<edm::InputTag>( "tagRecHit" ))),
  //tokenLocalTrack_      ( consumes< edm::DetSetVector<CTPPSDiamondLocalTrack> >      ( iConfig.getParameter<edm::InputTag>( "tagLocalTrack" ))),
  tokenPixelLocalTrack_      ( consumes< edm::DetSetVector<CTPPSPixelLocalTrack> >      ( iConfig.getParameter<edm::InputTag>( "tagPixelLocalTrack" )))
 {

   /*for (int i=0; i<MAX_SECTOR_NUMBER; i++){
     for (int j=0; j<2; j++){
       if (j==1) j=2;
       Ntracks_cuts_map_[std::make_pair((Sector_id_)i,(Station_id_)j)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[i+j],
                                                iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[i+j]);
     }
   }*/

   Ntracks_cuts_map_[std::make_pair(SECTOR_45_ID,STATION_210_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[0],
                                            iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[0]);
   Ntracks_cuts_map_[std::make_pair(SECTOR_45_ID,STATION_220_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[1],
                                            iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[1]);
   Ntracks_cuts_map_[std::make_pair(SECTOR_56_ID,STATION_210_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[2],
                                            iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[2]);
   Ntracks_cuts_map_[std::make_pair(SECTOR_56_ID,STATION_220_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[3],
                                            iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[3]);

    produces<Pmux_map>( "PixelMuxmap" ).setBranchAlias( "pixel_Mux_map_");
    produces<TBA>( "SectorTBA" ).setBranchAlias( "sector_TBA");
		/*produces<RecHit>( "ctppsDiamondRecHits" ).setBranchAlias("ctpps_Diamond_RecHits");
		produces<PixelLocal>( "ctppsPixelLocalTracks" ).setBranchAlias("ctpps_Pixel_LocalTracks");
		produces<DiamondDigi>( "ctppsDiamondDigi" ).setBranchAlias("ctpps_Diamond_Digi");
		produces<DiamondLocalTracks>( "ctppsDiamondLocalTracks" ).setBranchAlias("ctpps_Diamond_Local_Track");*/

}


DiamondTimingProducer::~DiamondTimingProducer()
{

   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DiamondTimingProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //edm::Handle< edm::DetSetVector<CTPPSDiamondDigi> > timingDigi;
  //edm::Handle<RecHit> timingRecHitt;
  edm::Handle<PixelLocal> pixelLocalTrack;
	//edm::Handle<DiamondDigi> diamondDigit;
	//edm::Handle<DiamondLocalTracks> diamondLocalTrackt;
  //iEvent.getByToken( tokenDigi_, timingDigi );
  //iEvent.getByToken( tokenRecHit_, timingRecHitt );
  iEvent.getByToken( tokenPixelLocalTrack_, pixelLocalTrack );
	//iEvent.getByToken( tokenDigi_, diamondDigit);
	//iEvent.getByToken( tokenLocalTrack_, diamondLocalTrackt);

	/*RecHit timingRecHit=*timingRecHitt;
	PixelLocal pixelLocalTrackt=*pixelLocalTrack;
	DiamondDigi diamondDigi=*diamondDigit;
	DiamondLocalTracks diamondLocalTrack=*diamondLocalTrackt;*/
//*run_ = iEvent.id().run();
//*ev_ = iEvent.id().event();
//*lumiblock_ = iEvent.luminosityBlock();

 ////////////////////////////////////////////////////////////////
//
//		EXTRACT PIXELS TRACK NUMBER
//
/////////////////////////////////////////////////////////////////

  Pixel_Mux_map_.clear();

  TBA Sector_TBA(2,true);

   for (const auto& RP_trks : *pixelLocalTrack) //array of tracks
  {
      const CTPPSDetId detId( RP_trks.detId() );

 //	std::cout << "Tracks in arm " << detId.arm() << ", station " << detId.station() << ", rp " << detId.rp() << std::endl;

   for ( const auto& trk : RP_trks )
   {
   if ( !trk.isValid() ) continue;
   //pixelsppslist.push_back( std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>(&trk, det_id) );
   //std::cout << "track found in " << trk.getX0() << " , "<< trk.getY0() << std::endl;

   //Pixel_Hit_Hmap_[ std::make_pair(detId.arm(), detId.station()) ] -> Fill(trk.getX0(),trk.getY0());
   //if (included("Pixel_Hit_Hmap_")) Histos_TH2F[MapKey("Pixel_Hit_Hmap_",detId.station(),detId.arm(),-1,-1)] -> Fill(trk.getX0(),trk.getY0());
   Pixel_Mux_map_[ std::make_pair(detId.arm(), detId.station()) ]++;
      }

 }

  for (const auto& Ntracks_cuts_iter_ :  Ntracks_cuts_map_)
  {
   if ( (Ntracks_cuts_iter_.second.first < 0) || (Ntracks_cuts_iter_.second.second < 0) ) continue; // don't care condition
   if ( (Pixel_Mux_map_[Ntracks_cuts_iter_.first] < Ntracks_cuts_iter_.second.first) ||
      (Pixel_Mux_map_[Ntracks_cuts_iter_.first] > Ntracks_cuts_iter_.second.second))  //condition violated
   {
     Sector_TBA[Ntracks_cuts_iter_.first.first] = false;
   }
  }

	std::unique_ptr<Pmux_map> Pmm(new Pmux_map(Pixel_Mux_map_));
  std::unique_ptr<TBA> ST(new TBA(Sector_TBA));
	/*std::unique_ptr<RecHit> RH(new RecHit(timingRecHit));
	std::unique_ptr<PixelLocal> PL(new PixelLocal(pixelLocalTrackt));
	std::unique_ptr<DiamondDigi> DD(new DiamondDigi(diamondDigi));
	std::unique_ptr<DiamondLocalTracks> DLT(new DiamondLocalTracks(diamondLocalTrack));*/

  iEvent.put(std::move(Pmm),"PixelMuxmap");
  iEvent.put(std::move(ST),"SectorTBA");
	/*iEvent.put(std::move(RH),"ctppsDiamondRecHits");
	iEvent.put(std::move(PL),"ctppsPixelLocalTracks");
	iEvent.put(std::move(DD),"ctppsDiamondDigi");
	iEvent.put(std::move(DLT),"ctppsDiamondLocalTracks");*/
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DiamondTimingProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DiamondTimingProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DiamondTimingProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
DiamondTimingProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DiamondTimingProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DiamondTimingProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiamondTimingProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingProducer);
