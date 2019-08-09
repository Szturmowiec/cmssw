// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"


#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DiamondDetectorClass.h"

#include <map>
#include <cmath>
#include <string>
#include <utility>

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"

enum Station_id_
{
	STATION_210_M_ID,
	STATION_TIMING_ID,
	STATION_220_M_ID

	};

struct MapKey
	{
    std::string name;
		int station, sector, plane, channel;

		bool operator < (const MapKey &rhs) const {
			return std::tie(this->name, this->station, this->sector, this->plane, this->channel) < std::tie(rhs.name, rhs.station, rhs.sector, rhs.plane, rhs.channel);
		}

		MapKey(const std::string name="noname", const int station=-1, const int sector=-1, const int plane=-1, const int channel=-1)
		: name(name), station(station), sector(sector), plane(plane), channel(channel)
		{};
	};

	std::vector<std::string> excluded_histos; //list of histograms that are not to be created, its content depends on the debug_ variable value

Double_t FermiDirac(Double_t *x, Double_t *par) {

  Double_t value=0.;
  value=(par[0]/(exp((x[0]-par[1])/par[2])+1))+par[3];

  return value;
}

class DiamondTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

	typedef std::map< std::pair<int,int> , TH2F* >  TH2F_map;
	typedef std::map< std::pair<int,int> , TH1F* >  TH1F_map;


   public:
      explicit DiamondTimingAnalyzer(const edm::ParameterSet&);
      ~ DiamondTimingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:



      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      virtual void endJob() override;

      void initHistograms( const CTPPSDiamondDetId&);

      // ---------- constants ---------------------------
      static const int CHANNELS_X_PLANE=12;
      static const int PLANES_X_DETECTOR=4;
      static const int MAX_SECTOR_NUMBER=2;
			static const int MAX_STATION_NUMBER=2;
	  edm::Service<TFileService> fs_;

	enum Sector_id_
	{
		SECTOR_45_ID,
		SECTOR_56_ID
    };

	enum Plane_id_
	{
		PLANE_0_ID,
		PLANE_1_ID,
		PLANE_2_ID,
		PLANE_3_ID

    };

      // ---------- objects to retrieve ---------------------------
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondDigi> > tokenDigi_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSPixelLocalTrack> > tokenPixelLocalTrack_;

      // ---------- directories ---------------------------
      std::map< CTPPSDiamondDetId, TFileDirectory > mainDir_map_;
      std::map< CTPPSDiamondDetId, TFileDirectory > planeDir_map_;
      std::map< ChannelKey, TFileDirectory > channelDir_map_;

	  std::vector<TFileDirectory> dir_cyl_V ;
	  std::vector<TFileDirectory> dir_cyl_DDTI_V ;
	  std::vector < std::map<std::pair<int,int>, TFileDirectory > > dir_cyl_TPTI_V ;
	  std::vector < std::map<int, TFileDirectory > > dir_cyl_L2_3p_Res_V ;
	  std::vector < std::map<int, TFileDirectory > > dir_cyl_L2Res_V ;
	  std::vector < std::vector<TFileDirectory> > dir_plane_VV;


      // ---------- global histograms ---------------------------

			TH1F* RunNumber_Histo_;
      std::map<MapKey,TH1F*> Histos_TH1F;
			std::map<MapKey,TH2F*> Histos_TH2F;
      std::map< ChannelKey, TProfile*> TOTvsT_Profile_map_;


      // ---------- correlation histograms ---------------------------

      std::map< std::pair<ChannelKey,ChannelKey> , TH1F* > Raw_DT_TrackedPlane_Hmap_;  //<<ChannelKey,ChannelKey>,histogram>
      std::map< std::pair<ChannelKey,ChannelKey> , TH1F* > SPC_DT_TrackedPlane_Hmap_;  //<<ChannelKey,ChannelKey>,histogram>
      std::map< std::pair<ChannelKey,ChannelKey> , TH2F* > ToT_vs_ToT_TrackedPlane_Hmap_; //<<ChannelKey,ChannelKey>,histogram>

      std::vector < TH2F_map > TimingCorrelation_mux_Hmap_V_; //plane 1, plane 2

      // ---------- pixel mux map ---------------------------

      std::map< std::pair< int , int >, int> Pixel_Mux_map_; //arm, station

      // ---------- channel graphs ---------------------------

      std::vector <TGraph*> Raw_resolution_DD_V_;
      std::vector <TGraph*> SPC_resolution_DD_V_;
      std::vector <TGraph*> Saturated_percentage_V_;

      std::vector <TH2F_map> Raw_resolution_TP_V_; // sector, <pl 1, pl 2>
      std::vector <TH2F_map> SPC_resolution_TP_V_;
      std::vector <TH2F_map> Correlation_TP_V_;


	  std::vector<TH1F*>  Tracks_time_RAW_histo_V_;
	  std::vector<TH1F*>  Tracks_time_SPC_histo_V_;
	  std::vector<TH1F*>  Tracks_resolution_histo_V_;

	  std::map<ChannelKey, double> Resolution_L1_map_;
	  std::map<ChannelKey, double> Resolution_L2_map_;
	  std::map<ChannelKey, double> Resolution_L2_3p_map_;

	  std::map<std::pair<int,int>, TGraph*>  Resolution_L2_graph_map_;
	  std::map<std::pair<int,int>, TGraph*>  Resolution_L2_3p_graph_map_;
	  std::map<std::pair<int,int>, TGraph*>  ImprouvedRes_L2_graph_map_;
	  std::map<std::pair<int,int>, TGraph*>  ImprouvedRes_L2_3p_graph_map_;

	  TTree  *TCalibration_;
	  int cal_channel_;
	  int cal_plane_;
	  int cal_sector_;
	  double cal_par_0_;
	  double cal_par_1_;
	  double cal_par_2_;
	  double cal_par_3_;
	  double Chi_2_;
	  double resolution_L1_;
	  double resolution_L2_;
	  double resolution_L2_3p_;

	  //external
	  DiamondDetectorClass DiamondDet;
	  int valid_OOT_;
		int debug_;

	  // selection parameters


     std::map< std::pair< int , int >, std::pair< int , int > > Ntracks_cuts_map_; //arm, station ,, Lcut,Ucut

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
DiamondTimingAnalyzer:: DiamondTimingAnalyzer(const edm::ParameterSet& iConfig)
 :
 tokenDigi_            ( consumes< edm::DetSetVector<CTPPSDiamondDigi> >      ( iConfig.getParameter<edm::InputTag>( "tagDigi" ) ) ),
 tokenRecHit_          ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >      ( iConfig.getParameter<edm::InputTag>( "tagRecHit" ))),
 tokenLocalTrack_          ( consumes< edm::DetSetVector<CTPPSDiamondLocalTrack> >      ( iConfig.getParameter<edm::InputTag>( "tagLocalTrack" ))),
 tokenPixelLocalTrack_      ( consumes< edm::DetSetVector<CTPPSPixelLocalTrack> >      ( iConfig.getParameter<edm::InputTag>( "tagPixelLocalTrack" ))),
 TimingCorrelation_mux_Hmap_V_(2),
 Raw_resolution_DD_V_(2),
 SPC_resolution_DD_V_(2),
 Saturated_percentage_V_(2),
 Raw_resolution_TP_V_(2),
 SPC_resolution_TP_V_(2),
 Correlation_TP_V_(2),
 Tracks_time_RAW_histo_V_(2),
 Tracks_time_SPC_histo_V_(2),
 Tracks_resolution_histo_V_(2),
 DiamondDet(iConfig,tokenRecHit_,tokenLocalTrack_),
valid_OOT_ (iConfig.getParameter< int >( "tagValidOOT" )),
debug_ (iConfig.getParameter< int >( "debug" ))
{
  usesResource("TFileService");

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

}


DiamondTimingAnalyzer::~ DiamondTimingAnalyzer()
{}

bool included(std::string name){
	return !(std::find(excluded_histos.begin(),excluded_histos.end(),name)!=excluded_histos.end());
}

//std::map<ChannelKey,TH1F*>
// member functions
//

void
DiamondTimingAnalyzer::initHistograms(const CTPPSDiamondDetId& detId)
{
   ChannelKey recHitKey(detId.arm(),detId.plane(),detId.channel());
  //now do what ever initialization is needed
  if (mainDir_map_.find(detId.getRPId()) == mainDir_map_.end())
  {
	std::cout << "RP is new, creating directory and histos" << std::endl;
    std::string dirName;
    detId.rpName(dirName, CTPPSDiamondDetId::nPath);
    std::string rpName;
    detId.rpName(rpName, CTPPSDiamondDetId::nFull);

    // create directory for the detector, if not already done
    mainDir_map_[ detId.getRPId() ] = fs_->mkdir( dirName );



   }

   if (planeDir_map_.find( detId.getPlaneId() ) == planeDir_map_.end())
   {
	  std::cout << "Plane is new, creating directory and histos" << std::endl;
      std::string dirName;
      detId.planeName(dirName, CTPPSDiamondDetId::nPath);
      std::string plName;
      detId.planeName(plName, CTPPSDiamondDetId::nFull);

      // create directory for the plane, if not already done
      planeDir_map_[ detId.getPlaneId() ] = fs_->mkdir( dirName );
   }
    //create channels histograms

   if (channelDir_map_.find(recHitKey) == channelDir_map_.end())
   {
	  std::cout << "Channels is new, creating directory and histos" << std::endl;
      std::string dirName;
      detId.channelName(dirName, CTPPSDiamondDetId::nPath);
      std::string chName;
      detId.channelName(chName, CTPPSDiamondDetId::nFull);

      // create directory for the channel, if not already done
      channelDir_map_[recHitKey] = fs_->mkdir( dirName );

	 // create all histograms

      std::vector<std::array<int,3>> sizes_1;
      sizes_1.push_back({{1200,-60,60}});
      sizes_1.push_back({{1200,-60,60}});
      sizes_1.push_back({{1200,-60,60}});
      sizes_1.push_back({{10,-5,5}});
      sizes_1.push_back({{10,-5,5}});
      sizes_1.push_back({{100,-20,20}});
      sizes_1.push_back({{100,-20,20}});

      std::vector<std::array<int,6>> sizes_2;
      sizes_2.push_back({{240, 0, 60, 450, -20, 25}});
      sizes_2.push_back({{240, 0, 60, 12, 0, 12}});
      sizes_2.push_back({{240, 0, 60, 2, 0, 2}});
      sizes_2.push_back({{240, 0, 60, 20, 0, 20}});
      sizes_2.push_back({{240, 0, 60, 12, 0, 12}});
      sizes_2.push_back({{1000, -100, 100, 1000, -100, 100}});
      sizes_2.push_back({{1000, -100, 100, 1000, -100, 100}});
      sizes_2.push_back({{240, 0, 60, 450, -20, 25}});

			std::string histo_names_1[]={"RAW_T_Distribution_","VALID_RAW_T_Distribution_","SinglePadCorr_T_Distribution_","OOT_Distribution_","ValidOOT_Distribution_",
      "TOT_Distribution_","ValidTOT_Distribution_"};

      std::string histo_names_2[]={"TvsTOT_Distribution_","TOTvsPlaneMux_Distribution_","TOTvsMH_Distribution_","TOTvsPT_Distribution_","TOTvsTT_Distribution_",
      "Pad tomography 220","Pad tomography 210","TOTvsSPCT_Distribution_"};

			if (debug_==5){
				for (std::string x : histo_names_1){
					//excluded_histos.push_back(x);
				}
				excluded_histos.push_back("RAW_T_Distribution_");
			}

      for (unsigned int i=0; i<std::size(histo_names_1); i++){
				if (included(histo_names_1[i])){
	        TH1F* m;
	        m=channelDir_map_[recHitKey].make<TH1F>((histo_names_1[i]+chName).c_str(), (histo_names_1[i]+chName).c_str(), sizes_1[i][0], sizes_1[i][1], sizes_1[i][2] );
	        MapKey mk(histo_names_1[i],-1,recHitKey.sector,recHitKey.plane,recHitKey.channel);
	        Histos_TH1F.insert(std::pair<MapKey,TH1F*>(mk,m));
				}
      }

      for (unsigned int i=0; i<std::size(histo_names_2); i++){
				if (included(histo_names_2[i])){
	        TH2F* m;
	        m=channelDir_map_[recHitKey].make<TH2F>((histo_names_2[i]+chName).c_str(), (histo_names_2[i]+chName).c_str(), sizes_2[i][0], sizes_2[i][1], sizes_2[i][2], sizes_2[i][3], sizes_2[i][4], sizes_2[i][5] );
					MapKey mk(histo_names_2[i],-1,recHitKey.sector,recHitKey.plane,recHitKey.channel);
	        Histos_TH2F.insert(std::pair<MapKey,TH2F*>(mk,m));
				}
      }
   }
}

void planes_res(int planes,DiamondDetectorClass DiamondDet,std::vector<bool> Sector_TBA,std::vector<std::map<int,TFileDirectory>> &dir,std::map<MapKey,TH1F*> &Histos_TH1F,
	int planes_x_detector,std::vector<TH1F*>  &Tracks_time_RAW_histo_V_,std::vector<TH1F*>  &Tracks_time_SPC_histo_V_,std::vector<TH1F*>  &Tracks_resolution_histo_V_){

		enum Sector_id_
		{
			SECTOR_45_ID,
			SECTOR_56_ID
	    };

		enum Plane_id_
		{
			PLANE_0_ID,
			PLANE_1_ID,
			PLANE_2_ID,
			PLANE_3_ID

	    };

	for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map()) // loop on predigested tracks
	{
	  int sec_number = LocalTrack_mapIter.first.getZ0() > 0.0 ? SECTOR_45_ID : SECTOR_56_ID;
	  if (!(Sector_TBA[sec_number])) continue;

		std::map<int,int> active_pl_map;
		if (planes==3){
		  if (sec_number==0)
		  {
		    active_pl_map[0]=0;
		    active_pl_map[1]=1;
		    active_pl_map[2]=1;
		    active_pl_map[3]=1;
		  }
		  else
		  {
		    active_pl_map[0]=1;
		    active_pl_map[1]=1;
		    active_pl_map[2]=1;
		    active_pl_map[3]=0;
		  }
		}

	  bool mark_tag=false;

	//std::cout << "check for mark" << std::endl;
	if (planes==4){
	  if (DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)==1 &&
	    DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)==1 &&
	    DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)==1 &&
	    DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)==1 &&
	    DiamondDet.GetTrackMuxInSector(sec_number)==1) mark_tag=true;
	}

	if (planes==3){
		if (!((DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)==1 || active_pl_map[PLANE_0_ID]==0 )&&
	    (DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)==1 || active_pl_map[PLANE_1_ID]==0 )&&
	    (DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)==1 || active_pl_map[PLANE_2_ID]==0 )&&
	    (DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)==1 || active_pl_map[PLANE_3_ID]==0 )&&
	    DiamondDet.GetTrackMuxInSector(sec_number)==1)) continue;
	}

		std::map<int,ChannelKey>  hit_selected;
	  //std::vector<ChannelKey>  hit_selected(4);
	  std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>::const_iterator  hit_iter;

	  //std::cout << "check for size" << std::endl;
	  if	(LocalTrack_mapIter.second.size() ==0) continue;

	  double Track_time_SPC = 100.0;
	  double Track_time_RAW = 100.0;
	  double Track_precision_SPC = 100.0;
	  double Track_precision_RAW = 1.0; // fixed weight

	  //std::cout << "start loop on hit" << std::endl;
	  for (hit_iter=LocalTrack_mapIter.second.begin();hit_iter < LocalTrack_mapIter.second.end(); hit_iter++ ) // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> >
	  {
	  //std::cout << "looping on hits" << std::endl;
	    double hit_time_SPC= DiamondDet.GetTime_SPC((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
	    double hit_time_RAW= DiamondDet.GetTime((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
	    double hit_prec_SPC= DiamondDet.GetPadPrecision((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
	    double hit_weig_SPC= DiamondDet.GetPadWeight((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);

			if (planes==3){
				if (active_pl_map[(*hit_iter).first.plane] ==1) hit_selected[(*hit_iter).first.plane] = (*hit_iter).first;
			}

			if (planes==4 && mark_tag){
		        hit_selected[(*hit_iter).first.plane] = (*hit_iter).first;  // save for resolution reco
	    }

	    if (hit_iter == LocalTrack_mapIter.second.begin())
	    {

	      //	std::cout << "setting starting timing" << std::endl;
	      Track_time_RAW = hit_time_RAW;
	      Track_precision_RAW =  1.0;
	      Track_time_SPC = hit_time_SPC;
	      Track_precision_SPC =  hit_prec_SPC;
	    }
	    else
	    {
	        //std::cout << "updating  timing" << std::endl;
	      Track_time_RAW = (Track_time_RAW*pow(Track_precision_RAW,-2) + hit_time_RAW*1.0)/(pow(Track_precision_RAW,-2)+1.0);
	      Track_precision_RAW = pow((pow(Track_precision_RAW,-2)+1.0),-0.5);
	      Track_time_SPC = (Track_time_SPC*pow(Track_precision_SPC,-2) + hit_time_SPC*hit_weig_SPC)/(pow(Track_precision_SPC,-2)+hit_weig_SPC);
	      Track_precision_SPC = pow((pow(Track_precision_SPC,-2)+hit_weig_SPC),-0.5);
	    }

			}

	        //std::cout << "saving timing" << std::endl;
		if (planes==4){
		  Tracks_time_SPC_histo_V_[ sec_number ]->Fill(Track_time_SPC);
		  Tracks_resolution_histo_V_[ sec_number ]->Fill(Track_precision_SPC);
		  Tracks_time_RAW_histo_V_[ sec_number ]->Fill(Track_time_RAW);
		}

	        //std::cout << "timing saved" << std::endl;

	  if (planes==3 || (planes==4 && mark_tag))
	  {
	    for (int pl_mark = 0 ; pl_mark < planes_x_detector; pl_mark++)
	    {
				if (planes==3 && active_pl_map[pl_mark] ==0) continue;

	      double Marked_track_time = 12.5;
	      double Marked_track_precision = 25.0;
	      double Marked_hit_time=0.0;
	      int Marked_hit_channel=-1;

	      for (int pl_loop = 0 ; pl_loop < planes_x_detector; pl_loop++)
	      {
					if (planes==3 && active_pl_map[pl_loop] ==0) continue;

	        if (pl_loop == pl_mark)
	        {
	          Marked_hit_time= DiamondDet.GetTime_SPC(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
	          Marked_hit_channel = hit_selected[pl_loop].channel;
	          continue;
	        }
	        double Others_hit_time= DiamondDet.GetTime_SPC(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
	        double Others_hit_prec= DiamondDet.GetPadPrecision(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
	        double Others_hit_weig= DiamondDet.GetPadWeight(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);

	        if (hit_iter == LocalTrack_mapIter.second.begin())
	        {

	          Marked_track_time = Others_hit_time;
	          Marked_track_precision = Others_hit_prec;
	        }
	        else
	        {
	          Marked_track_time = (Marked_track_time*pow(Marked_track_precision,-2) + Others_hit_time*Others_hit_weig)/(pow(Marked_track_precision,-2)+Others_hit_weig);
	          Marked_track_precision = pow((pow(Marked_track_precision,-2)+Others_hit_weig),-0.5);
	        }
	      }

	      double Marked_hit_difference = Marked_hit_time-Marked_track_time;

				std::pair<std::string,std::string> names_histos_L2_1[]={std::make_pair("L2_resolution_DT_sec_","Resolution_L2_histo_map_"),std::make_pair("L2_partialTrk_Time_sec_","TrkTime_L2_histo_map_"),
	      std::make_pair("L2_partialTrk_ExpectedRes_sec_","TrkExpectedRes_L2_histo_map_")};

				if (planes==3){
					for (unsigned int i=0; i<size(names_histos_L2_1); i++){
						std::get<0>(names_histos_L2_1[i]).insert(std::get<0>(names_histos_L2_1[i]).find("L2_")+3,"3p_");
						std::get<1>(names_histos_L2_1[i]).insert(std::get<1>(names_histos_L2_1[i]).find("L2_")+3,"3p_");
					}
				}

				std::vector<std::array<int,3>> sizes_L2_1;
	      sizes_L2_1.push_back({{1200,-60,60}});
	      sizes_L2_1.push_back({{1200,-60,60}});
	      sizes_L2_1.push_back({{1000,0,2}});

	      if (Histos_TH1F.find(MapKey(std::get<1>(names_histos_L2_1[0]),-1,sec_number,pl_mark,Marked_hit_channel))==Histos_TH1F.end()) // creta pair histos if do not exist
	      {
					for (unsigned int i=0; i<size(names_histos_L2_1); i++){
						std::string name=std::get<0>(names_histos_L2_1[i])+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+std::to_string(Marked_hit_channel);
						MapKey mk=MapKey(std::get<1>(names_histos_L2_1[i]),-1,sec_number,pl_mark,Marked_hit_channel);
						Histos_TH1F[mk]=dir[sec_number][pl_mark].make<TH1F>(name.c_str(), name.c_str(), sizes_L2_1[i][0], sizes_L2_1[i][1], sizes_L2_1[i][2] );
					}
	      }

	      Histos_TH1F[MapKey(std::get<1>(names_histos_L2_1[0]),-1,sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_hit_difference);
	      Histos_TH1F[MapKey(std::get<1>(names_histos_L2_1[1]),-1,sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_time);
	      Histos_TH1F[MapKey(std::get<1>(names_histos_L2_1[2]),-1,sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_precision);

	    	}
	  	}
		}
	}

// ------------ method called for each event  ------------
void
DiamondTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //edm::Handle< edm::DetSetVector<CTPPSDiamondDigi> > timingDigi;
  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > timingRecHit;
  edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelLocalTrack;
  //iEvent.getByToken( tokenDigi_, timingDigi );
  iEvent.getByToken( tokenRecHit_, timingRecHit );
  iEvent.getByToken( tokenPixelLocalTrack_, pixelLocalTrack );





//*run_ = iEvent.id().run();
//*ev_ = iEvent.id().event();
//*lumiblock_ = iEvent.luminosityBlock();


 ////////////////////////////////////////////////////////////////
//
//		EXTRACT PIXELS TRACK NUMBER
//
/////////////////////////////////////////////////////////////////


  Pixel_Mux_map_.clear();

  std::vector<bool> Sector_TBA(2,true);

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
		if (included("Pixel_Hit_Hmap_")) Histos_TH2F[MapKey("Pixel_Hit_Hmap_",detId.station(),detId.arm(),-1,-1)] -> Fill(trk.getX0(),trk.getY0());
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

if (!(Sector_TBA[0] || Sector_TBA[1])) return;


  RunNumber_Histo_->Fill(iEvent.id().run());


 ////////////////////////////////////////////////////////////////
//
//		EXTRACT Dimoand detector info
//
/////////////////////////////////////////////////////////////////

  DiamondDet.ExtractData(iEvent);

 ////////////////////////////////////////////////////////////////
//
//		Tommography
//
/////////////////////////////////////////////////////////////////

 for (const auto& RP_trks : *pixelLocalTrack) //array of tracks
 {
      const CTPPSDetId detId( RP_trks.detId() );

	//	std::cout << "Tracks in arm " << detId.arm() << ", station " << detId.station() << ", rp " << detId.rp() << std::endl;
	  for ( const auto& trk : RP_trks )
	  {
		if ( !trk.isValid() ) continue;
		if ((Pixel_Mux_map_[ std::make_pair(detId.arm(), detId.station()) ]==1) && (DiamondDet.GetTrackMuxInSector(detId.arm()) > 0) && (Sector_TBA[detId.arm()]))
				if (included("Pixel_Tomography_Hmap_")) Histos_TH2F[MapKey("Pixel_Tomography_Hmap_",detId.station(),detId.arm(),-1,-1)] -> Fill(trk.getX0(),trk.getY0());
      }

}



////////////////////////////////////////////////////////////////
//
//		Single pad analisys
//
/////////////////////////////////////////////////////////////////




  for (const auto& recHits : *timingRecHit) //rechits = array of hits in one channel
  {
    const CTPPSDiamondDetId detId( recHits.detId() );
	if (!(Sector_TBA[detId.arm()])) continue;

	ChannelKey  recHitKey(detId.arm(),detId.plane(),detId.channel());

    if (channelDir_map_.find(recHitKey) == channelDir_map_.end() && recHits.size() > 0)
    {
		std::cout << "Found new channel with data! Creating Histograms..." << std::endl;
		initHistograms( detId );
	}

	// Perform channel histogram

    for (const auto& recHit : recHits ) //rechit
    {
		//std::cout << "Hits in channel " << detId.channel() << ", plane " << detId.plane() << ", rp " << detId.rp()<< ", station " << detId.station()<<
	//", arm " << detId.arm() << std::endl;

		if (((recHit.getOOTIndex() != (int)((DiamondDet.GetSPCMap())[recHitKey].offset/25) ) &&  valid_OOT_!=-1) ||  recHit.getMultipleHits()) continue;

    if (included("RAW_T_Distribution_")) Histos_TH1F[MapKey("RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->Fill( recHit.getT() );
    if (included("TOT_Distribution_")) Histos_TH1F[MapKey("TOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->Fill( recHit.getToT() );
    if (included("OOT_Distribution_")) Histos_TH1F[MapKey("OOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->Fill( recHit.getOOTIndex() );

		if (DiamondDet.PadActive(detId.arm(), detId.plane(),detId.channel())) // T and Tot present
		{

			if (included("ValidTOT_Distribution_")) Histos_TH1F[MapKey("ValidTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill( DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()) );
			if (included("SinglePadCorr_T_Distribution_")) Histos_TH1F[MapKey("SinglePadCorr_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill( DiamondDet.GetTime_SPC(detId.arm(), detId.plane(),detId.channel()));
			if (included("ValidOOT_Distribution_")) Histos_TH1F[MapKey("ValidOOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill( DiamondDet.GetPadOOT(detId.arm(), detId.plane(),detId.channel()));
			if (included("TOTvsSPCT_Distribution_")) Histos_TH2F[MapKey("TOTvsSPCT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetTime_SPC(detId.arm(),detId.plane(),detId.channel()));
			if (included("TOTvsPlaneMux_Distribution_")) Histos_TH2F[MapKey("TOTvsPlaneMux_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetMux(detId.arm(), detId.plane()));
			if (included("TOTvsMH_Distribution_")) Histos_TH2F[MapKey("TOTvsMH_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetMH(detId.arm(), detId.plane(),detId.channel()));
			if (included("TOTvsPT_Distribution_")) Histos_TH2F[MapKey("TOTvsPT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), Pixel_Mux_map_[std::make_pair(detId.arm(),STATION_220_M_ID)]);
			if (included("TOTvsTT_Distribution_")) Histos_TH2F[MapKey("TOTvsTT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()),  DiamondDet.GetTrackMuxInSector(detId.arm()));

			if (included("TvsTOT_Distribution_")) Histos_TH2F[MapKey("TvsTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()));

			if (included("VALID_RAW_T_Distribution_")) Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]-> Fill( DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()) );
		}
    }
  }



////////////////////////////////////////////////////////////////
//
//		Multiplicity correlation for all events
//
/////////////////////////////////////////////////////////////////

   // plane multiplicity correlation
for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
{

	if (!(Sector_TBA[sec_number])) continue;

  	for  (int pl_number=0; pl_number < PLANES_X_DETECTOR; pl_number++)
	{
		for  (int pl_number_2 = pl_number+1; pl_number_2 < PLANES_X_DETECTOR; pl_number_2++)

		{
			//multiplicity corralation between planes (all hits)
			TimingCorrelation_mux_Hmap_V_[sec_number][ std::make_pair(pl_number,pl_number_2) ]->Fill(DiamondDet.GetMux(sec_number,pl_number),DiamondDet.GetMux(sec_number,pl_number_2));
		}

		// multiplicity distribution (all hits)
		Histos_TH1F[ MapKey("PlaneMux_sSector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMux(sec_number,pl_number));
		// multiplicity distribution (valid time hits)
		Histos_TH1F[ MapKey("PlaneMuxValidT_sector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number));
		// hits multiplicity (valid T) Vs Spread (valid T)
		Histos_TH2F[ MapKey("Mux-Spread correlation sector ",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number),DiamondDet.GetSpread(sec_number,pl_number));
		// hits multiplicty (valid T) Vs 220 pixel tracks multiplicity
		Histos_TH2F[MapKey("Diamond-Pixel 220 mux correlation sector ",-1,sec_number,pl_number,-1) ]-> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_220_M_ID)]);
		// hits multiplicty (valid T) Vs 210 pixel tracks multiplicity
		Histos_TH2F[MapKey("Diamond-Pixel 210 mux correlation sector ",-1,sec_number,pl_number,-1) ]-> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_210_M_ID)]);
        // hits multiplicty (in tracks) Vs Spread (in tracks)
		Histos_TH2F[ MapKey("Tracks Mux-Spread correlation sector",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxInTrack(sec_number,pl_number),DiamondDet.GetTrackSpread(sec_number,pl_number));
        // multiplicity distribution (hits contained in tracks)
		Histos_TH1F[ MapKey("Track_PlaneMuxValidT_sector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxInTrack(sec_number,pl_number));

	}


	Histos_TH1F[MapKey("Timing tracks mux sector",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number));
	Histos_TH2F[MapKey("DiamondTracks-Pixel 220 mux correlation sector ",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_220_M_ID)]);
	Histos_TH2F[MapKey("DiamondTracks-Pixel 210 mux correlation sector ",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_210_M_ID)]);



	if (DiamondDet.SecIsSaturated(sec_number))
	{
		for  (int pl_number=0; pl_number < PLANES_X_DETECTOR; pl_number++)
			{

				Histos_TH1F[ MapKey("PlaneMuxSaturated_sector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMux(sec_number,pl_number));
				Histos_TH1F[ MapKey("PlaneMuxSaturatedValidT_sector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number));
				Histos_TH1F[ MapKey("Track_PlaneMuxSaturatedValidT_sector_",-1,sec_number,pl_number,-1) ] -> Fill(DiamondDet.GetMuxInTrack(sec_number,pl_number));
				Histos_TH2F[MapKey("Diamond-Pixel 220 Saturated mux correlation sector ",-1,sec_number,pl_number,-1) ]-> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_220_M_ID)]);
				Histos_TH2F[MapKey("Diamond-Pixel 210 Saturated mux correlation sector ",-1,sec_number,pl_number,-1) ]-> Fill(DiamondDet.GetMuxValidT(sec_number,pl_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_210_M_ID)]);

			}

		Histos_TH1F[MapKey("Timing tracks saturated mux sector",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number));
		Histos_TH2F[MapKey("DiamondTracks-Pixel 220 Saturated mux correlation sector ",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_220_M_ID)]);
		Histos_TH2F[MapKey("DiamondTracks-Pixel 210 Saturated mux correlation sector ",-1,sec_number,-1,-1)]-> Fill(DiamondDet.GetTrackMuxInSector(sec_number),Pixel_Mux_map_[std::make_pair(sec_number,STATION_210_M_ID)]);
	}
}

for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map())
{
	int sec_number = LocalTrack_mapIter.first.getZ0() > 0.0 ? SECTOR_45_ID : SECTOR_56_ID;


	if (!(Sector_TBA[sec_number])) continue;

	Histos_TH1F[MapKey(" Reconstructed Constituent Pad Mux sector",-1,sec_number,-1,-1)]->Fill(LocalTrack_mapIter.second.size());
	Histos_TH2F[ MapKey("Timing track hit map sector ",-1,sec_number,-1,-1) ]-> Fill (LocalTrack_mapIter.first.getX0(),LocalTrack_mapIter.first.getY0());
	// start TBR
	/*if ((LocalTrack_mapIter.first.getX0()<6.0) && (std::abs(LocalTrack_mapIter.first.getY0())>0.2))
	{
		std::cout << "2017 anoumaluos track" << std::endl;
		DiamondDet.DumpEvent();
	}*/
}

////////////////////////////////////////////////////////////////
//
//		TIMING STUDIES
//
/////////////////////////////////////////////////////////////////


	// double diamond analysis method

for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
{

	if (!(Sector_TBA[sec_number])) continue;

	for  (int ch_number=0; ch_number < CHANNELS_X_PLANE; ch_number++)
	{
		if (DiamondDet.PairActive(sec_number, PLANE_2_ID,ch_number,PLANE_3_ID,ch_number))
		{
			Histos_TH1F[MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_number)]->Fill(DiamondDet.GetTime(sec_number,PLANE_2_ID,ch_number)-DiamondDet.GetTime(sec_number,PLANE_3_ID,ch_number));



			double TimeA_corr=DiamondDet.GetTime_SPC(sec_number, PLANE_2_ID,ch_number);
			double TimeB_corr=DiamondDet.GetTime_SPC(sec_number, PLANE_3_ID,ch_number);

			Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_number)]->Fill(TimeA_corr-TimeB_corr);
			Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,sec_number,-1,ch_number)]->Fill(DiamondDet.GetToT(sec_number, PLANE_2_ID,ch_number),DiamondDet.GetToT(sec_number,PLANE_3_ID,ch_number));

		}
	}
}

//typedef std::map< CTPPSDiamondLocalTrack , std::pair<ChannelKey,CTPPSDiamondRecHit> >   LocalTrack_map;
//std::vector<TH1F_map> Raw_DT_TrackedPlane_Hmap_;  //<<ChannelKey1,ChannelKey>,histogram>

// loop on all tracks --> generation of time resolution timing track based
for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map()) // loop on predigested tracks
{
	int sec_number = LocalTrack_mapIter.first.getZ0() > 0.0 ? SECTOR_45_ID : SECTOR_56_ID;


	if (!(Sector_TBA[sec_number])) continue;


	std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>::const_iterator  hitA_iter;
	std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>::const_iterator  hitB_iter;



//std::cout << "Before hit loop: #hit " << LocalTrack_mapIter.second.size()  << std::endl;
	if	(LocalTrack_mapIter.second.size() < 2) continue;

	for (hitA_iter=LocalTrack_mapIter.second.begin();hitA_iter < LocalTrack_mapIter.second.end()-1; hitA_iter++ ) // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> >
	{

		for (hitB_iter=hitA_iter+1; hitB_iter < LocalTrack_mapIter.second.end(); hitB_iter++ ) // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> >
		{
			if ((*hitA_iter).first.plane == (*hitB_iter).first.plane) continue;
			if((*hitA_iter).first.sector != (*hitB_iter).first.sector)  std::cout << "WARNING: hits from different sectors in same tracks!!" << std::endl;

			if (Raw_DT_TrackedPlane_Hmap_.find(std::make_pair((*hitA_iter).first, (*hitB_iter).first))==Raw_DT_TrackedPlane_Hmap_.end()) // creta pair histos if do not exist
			{
				std::string Histo_name = "Raw_resolution_TP_sec_"+std::to_string(sec_number)+"_planes_"+std::to_string((*hitA_iter).first.plane)+"_"+std::to_string((*hitB_iter).first.plane)+
					"_channels_"+std::to_string((*hitA_iter).first.channel)+"_"+std::to_string((*hitB_iter).first.channel);
				Raw_DT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)] =
					dir_cyl_TPTI_V[sec_number][std::make_pair((*hitA_iter).first.plane,(*hitB_iter).first.plane)].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );

				Histo_name = "SPC_resolution_TP_sec_"+std::to_string(sec_number)+"_planes_"+std::to_string((*hitA_iter).first.plane)+"_"+std::to_string((*hitB_iter).first.plane)+
					"_channels_"+std::to_string((*hitA_iter).first.channel)+"_"+std::to_string((*hitB_iter).first.channel);
				SPC_DT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)] =
					dir_cyl_TPTI_V[sec_number][std::make_pair((*hitA_iter).first.plane,(*hitB_iter).first.plane)].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );

				Histo_name = "TOT_vs_TOT_TP_sector_"+std::to_string(sec_number)+"_planes_"+std::to_string((*hitA_iter).first.plane)+"_"+std::to_string((*hitB_iter).first.plane)+
					"_channels_"+std::to_string((*hitA_iter).first.channel)+"_"+std::to_string((*hitB_iter).first.channel);
				ToT_vs_ToT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)] =
						dir_cyl_TPTI_V[sec_number][std::make_pair((*hitA_iter).first.plane,(*hitB_iter).first.plane)].make<TH2F>(Histo_name.c_str(), Histo_name.c_str(), 100, 5, 25, 100, 5, 25 );

			}


			Raw_DT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)]->
				Fill(DiamondDet.GetTime((*hitA_iter).first.sector,(*hitA_iter).first.plane,(*hitA_iter).first.channel)-DiamondDet.GetTime((*hitB_iter).first.sector,(*hitB_iter).first.plane,(*hitB_iter).first.channel));

			double TimeA_corr=DiamondDet.GetTime_SPC((*hitA_iter).first.sector,(*hitA_iter).first.plane,(*hitA_iter).first.channel);
			double TimeB_corr=DiamondDet.GetTime_SPC((*hitB_iter).first.sector,(*hitB_iter).first.plane,(*hitB_iter).first.channel);

			SPC_DT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)]->Fill(TimeA_corr-TimeB_corr);
			ToT_vs_ToT_TrackedPlane_Hmap_[std::make_pair((*hitA_iter).first, (*hitB_iter).first)]->
				Fill(DiamondDet.GetToT((*hitA_iter).first.sector,(*hitA_iter).first.plane,(*hitA_iter).first.channel),DiamondDet.GetToT((*hitB_iter).first.sector,(*hitB_iter).first.plane,(*hitB_iter).first.channel));
			Correlation_TP_V_[sec_number][std::make_pair((*hitA_iter).first.plane,(*hitB_iter).first.plane)] ->Fill((*hitA_iter).first.channel, (*hitB_iter).first.channel);

		}

	}

}



////////////////////////////////////////////////////////////////
//
//		TRACK TIMING STUDIES
//
/////////////////////////////////////////////////////////////////

//std::cout << "start new block" << std::endl;


////////////////////////////////////////////////////////////////
//
//		4 planes!!!!!
//
/////////////////////////////////////////////////////////////////

//planes_4(DiamondDet,Sector_TBA,dir_cyl_L2Res_V,Histos_L2_1,PLANES_X_DETECTOR,Tracks_time_RAW_histo_V_,Tracks_time_SPC_histo_V_,Tracks_resolution_histo_V_);
planes_res(4,DiamondDet,Sector_TBA,dir_cyl_L2Res_V,Histos_TH1F,PLANES_X_DETECTOR,Tracks_time_RAW_histo_V_,Tracks_time_SPC_histo_V_,Tracks_resolution_histo_V_);

////////////////////////////////////////////////////////////////
//
//		3 planes!!!!!
//
/////////////////////////////////////////////////////////////////

//planes_3(DiamondDet,Sector_TBA,dir_cyl_L2Res_V,Histos_L2_1,PLANES_X_DETECTOR,Tracks_time_RAW_histo_V_,Tracks_time_SPC_histo_V_,Tracks_resolution_histo_V_,dir_cyl_L2_3p_Res_V);
planes_res(3,DiamondDet,Sector_TBA,dir_cyl_L2_3p_Res_V,Histos_TH1F,PLANES_X_DETECTOR,Tracks_time_RAW_histo_V_,Tracks_time_SPC_histo_V_,Tracks_resolution_histo_V_);
//std::cout << "end new block" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void
DiamondTimingAnalyzer::beginJob()
{
   // edm::Service<TFileService> fs_;

	dir_cyl_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr" ));
	dir_cyl_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr" ));

	dir_cyl_DDTI_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/DDTimeInfo" ));
	dir_cyl_DDTI_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/DDTimeInfo" ));

	dir_cyl_TPTI_V.resize(2);

  std::string sector_names[]={"sector 45","sector 56"};
  std::string plane_names_1[]={"plane0_1","plane0_2","plane0_3","plane1_2","plane1_3","plane2_3"};

  for (unsigned int i=0; i<std::size(sector_names); i++){
    for (unsigned int j=0; j<std::size(plane_names_1); j++){
      std::string p=plane_names_1[j].substr(5,1);
      std::string q=plane_names_1[j].substr(7,1);
      int a=std::atoi(p.c_str());
      int b=std::atoi(q.c_str());
      dir_cyl_TPTI_V[i][std::make_pair(a,b)]=fs_->mkdir( "CTPPS/TimingDiamond/"+sector_names[i]+"/station 220cyl/cyl_hr/TPTimeInfo/"+plane_names_1[j] );
    }
  }

	dir_cyl_L2Res_V.resize(2);
	dir_cyl_L2_3p_Res_V.resize(2);

  int pl_num=4;

  for (unsigned int i=0; i<std::size(sector_names); i++){
    for (int j=0; j<pl_num; j++){
      dir_cyl_L2Res_V[i][j]=fs_->mkdir( "CTPPS/TimingDiamond/"+sector_names[i]+"/station 220cyl/cyl_hr/L2Res/plane"+std::to_string(j) );
      dir_cyl_L2_3p_Res_V[i][j]=fs_->mkdir( "CTPPS/TimingDiamond/"+sector_names[i]+"/station 220cyl/cyl_hr/L2_3p_Res/plane"+std::to_string(j) );
    }
  }



	dir_plane_VV.push_back(std::vector<TFileDirectory>());
	dir_plane_VV.push_back(std::vector<TFileDirectory>());

  for (unsigned int i=0; i<std::size(sector_names); i++){
    for (int j=0; j<pl_num; j++){
      dir_plane_VV[i].push_back(fs_->mkdir( "CTPPS/TimingDiamond/"+sector_names[i]+"/station 220cyl/cyl_hr/plane "+std::to_string(j) ));
    }
  }

	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
	{
		for  (int ch_number=0; ch_number < CHANNELS_X_PLANE; ch_number++)
		{


			std::string Raw_DT_name("RAW_DeltaT_DD_sector_"+std::to_string(sec_number)+"_channel_"+std::to_string(ch_number));
			std::string SPC_DT_name("SPC_DeltaT_DD_sector_"+std::to_string(sec_number)+"_channel_"+std::to_string(ch_number));
			std::string ToT_vs_ToT_name("TOT_vs_TOT_DD_sector_"+std::to_string(sec_number)+"_channel_"+std::to_string(ch_number));

			Histos_TH1F[ MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_number) ] = dir_cyl_DDTI_V[sec_number].make<TH1F>(Raw_DT_name.c_str(), Raw_DT_name.c_str(), 1200, -60, 60 );
			Histos_TH1F[ MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_number) ] = dir_cyl_DDTI_V[sec_number].make<TH1F>(SPC_DT_name.c_str(), SPC_DT_name.c_str(), 1200, -60, 60 );
			Histos_TH2F[ MapKey("ToT_vs_ToT_Hmap_",-1,sec_number,-1,ch_number) ] = dir_cyl_DDTI_V[sec_number].make<TH2F>(ToT_vs_ToT_name.c_str(), ToT_vs_ToT_name.c_str(), 100, 5, 25, 100, 5, 25 );

		}
	}

	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
	{

		for  (int pl_number=0; pl_number < PLANES_X_DETECTOR; pl_number++)
		{

      std::vector<std::array<int,3>> sizes_mux_1;
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{16,0,16}});
      sizes_mux_1.push_back({{12,0,12}});

      std::vector<std::array<int,6>> sizes_mux_2;
      sizes_mux_2.push_back({{12, 0, 12, 12, 0, 12}});
			sizes_mux_2.push_back({{12, 0, 12, 13, 0, 12}});
      sizes_mux_2.push_back({{12, 0, 12, 20, 0, 20}});
      sizes_mux_2.push_back({{12, 0, 12, 20, 0, 20}});
      sizes_mux_2.push_back({{12, 0, 12, 20, 0, 20}});
      sizes_mux_2.push_back({{12, 0, 12, 20, 0, 20}});

			std::string histo_mux_names_1[]={"PlaneMux_sSector_","PlaneMuxValidT_sector_","PlaneMuxSaturated_sector_","PlaneMuxSaturatedValidT_sector_","Track_PlaneMuxValidT_sector_",
      "Track_PlaneMuxSaturatedValidT_sector_","Valid TOT mean sector "};

      std::string histo_mux_names_2[]={"Mux-Spread correlation sector ","Tracks Mux-Spread correlation sector","Diamond-Pixel 220 mux correlation sector ",
      "Diamond-Pixel 210 mux correlation sector ","Diamond-Pixel 220 Saturated mux correlation sector ","Diamond-Pixel 210 Saturated mux correlation sector "};

      for (unsigned int i=0; i<std::size(histo_mux_names_1); i++){
				std::string p="_plane_";
				if (histo_mux_names_1[i]=="Valid TOT mean sector ") p=" plane ";
        std::string Plane_Mux_Hname(histo_mux_names_1[i]+std::to_string(sec_number)+p+std::to_string(pl_number));
        Histos_TH1F[ MapKey(histo_mux_names_1[i],-1,sec_number,pl_number,-1) ] = dir_plane_VV[sec_number][pl_number].make<TH1F>(Plane_Mux_Hname.c_str(), Plane_Mux_Hname.c_str(), sizes_mux_1[i][0], sizes_mux_1[i][1], sizes_mux_1[i][2] );
      }

      for (unsigned int i=0; i<std::size(histo_mux_names_2); i++){
        std::string TP_Mux_Hname(histo_mux_names_2[i]+std::to_string(sec_number)+" plane "+std::to_string(pl_number));
        Histos_TH2F[ MapKey(histo_mux_names_2[i],-1,sec_number,pl_number,-1) ] = dir_plane_VV[sec_number][pl_number].make<TH2F>(TP_Mux_Hname.c_str(), TP_Mux_Hname.c_str(), sizes_mux_2[i][0], sizes_mux_2[i][1], sizes_mux_2[i][2], sizes_mux_2[i][3], sizes_mux_2[i][4], sizes_mux_2[i][5] );
      }

			std::string Hname = "Valid TOT mean sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);


			Resolution_L2_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="L2 resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			Resolution_L2_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());

			Resolution_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2_3p_Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="L2_3p resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			Resolution_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());


			ImprouvedRes_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2_3p_Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="Improuved L2_3p resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			ImprouvedRes_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());


			ImprouvedRes_L2_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="Improuved L2 resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			ImprouvedRes_L2_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());


			for (int pl_number_2 = pl_number+1; pl_number_2 < PLANES_X_DETECTOR; pl_number_2++)
			{

				Hname = "Diamond plane mux correlation sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number)+ "_" + std::to_string(pl_number_2);
				TimingCorrelation_mux_Hmap_V_[sec_number][ std::make_pair(pl_number,pl_number_2) ] = dir_cyl_V[sec_number].make<TH2F>(Hname.c_str(), Hname.c_str(), 12, 0, 12, 12, 0, 12 );


				Hname = "Raw_resolution_TP_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_number)+ "_" + std::to_string(pl_number_2);
				Raw_resolution_TP_V_[sec_number][std::make_pair(pl_number,pl_number_2)] = dir_cyl_TPTI_V[sec_number][std::make_pair(pl_number,pl_number_2)].make<TH2F>(Hname.c_str(), Hname.c_str(), 12, 0, 12,12, 0, 12);

				Hname = "SPC_resolution_TP_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_number)+ "_" + std::to_string(pl_number_2);
				SPC_resolution_TP_V_[sec_number][std::make_pair(pl_number,pl_number_2)] = dir_cyl_TPTI_V[sec_number][std::make_pair(pl_number,pl_number_2)].make<TH2F>(Hname.c_str(), Hname.c_str(), 12, 0, 12,12, 0, 12);

				Hname = "Occupancy_TP_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_number)+ "_" + std::to_string(pl_number_2);
				Correlation_TP_V_[sec_number][std::make_pair(pl_number,pl_number_2)] = dir_cyl_TPTI_V[sec_number][std::make_pair(pl_number,pl_number_2)].make<TH2F>(Hname.c_str(), Hname.c_str(), 12, 0, 12,12, 0, 12);

			}

		}
	}

	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
	{

    std::vector<std::array<int,3>> sizes_timing_1;
    sizes_timing_1.push_back({{12, 0, 12}});
    sizes_timing_1.push_back({{12, 0, 12}});
    sizes_timing_1.push_back({{12, 0, 12}});

    std::vector<std::array<int,6>> sizes_timing_2;
    sizes_timing_2.push_back({{12, 0, 12, 20, 0, 20}});
    sizes_timing_2.push_back({{12, 0, 12, 20, 0, 20}});
    sizes_timing_2.push_back({{12, 0, 12, 20, 0, 20}});
    sizes_timing_2.push_back({{12, 0, 12, 20, 0, 20}});
    sizes_timing_2.push_back({{100, 0, 20, 10, -1, 1}});

    std::string timing_histos_names_1[]={" Reconstructed Constituent Pad Mux sector","Timing tracks mux sector","Timing tracks saturated mux sector"};
    std::string timing_histos_names_2[]={"DiamondTracks-Pixel 220 mux correlation sector ","DiamondTracks-Pixel 210 mux correlation sector ",
    "DiamondTracks-Pixel 220 Saturated mux correlation sector ","DiamondTracks-Pixel 210 Saturated mux correlation sector ","Timing track hit map sector "};

    for (unsigned int i=0; i<std::size(timing_histos_names_1); i++){
      std::string name=timing_histos_names_1[i]+ std::to_string(sec_number);
      Histos_TH1F[MapKey(timing_histos_names_1[i],-1,sec_number,-1,-1)] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), sizes_timing_1[i][0], sizes_timing_1[i][1], sizes_timing_1[i][2]);
    }

    for (unsigned int i=0; i<std::size(timing_histos_names_2); i++){
      std::string name=timing_histos_names_2[i]+ std::to_string(sec_number);
      Histos_TH2F[MapKey(timing_histos_names_2[i],-1,sec_number,-1,-1)] = dir_cyl_V[sec_number].make<TH2F>(name.c_str(), name.c_str(), sizes_timing_2[i][0], sizes_timing_2[i][1], sizes_timing_2[i][2], sizes_timing_2[i][3], sizes_timing_2[i][4], sizes_timing_2[i][5]);
    }

		std::pair<std::string,std::string> names_pixel_hitmap[]={std::make_pair("Pixel hit map sector ","Pixel_Hit_Hmap_"),std::make_pair("Tomography sector ","Pixel_Tomography_Hmap_")};

		std::vector<std::array<int,6>> sizes_hitmap;
		sizes_hitmap.push_back({{1000, -100, 100, 1000, -100, 100}});
		sizes_hitmap.push_back({{1000, -100, 100, 1000, -100, 100}});
		sizes_hitmap.push_back({{1000, -100, 100, 1000, -100, 100}});
		sizes_hitmap.push_back({{1000, -100, 100, 1000, -100, 100}});

			for(unsigned int i=0; i<size(names_pixel_hitmap); i++){
				for(int j=0; j<MAX_STATION_NUMBER; j++){
					Station_id_ id=STATION_220_M_ID;
					std::string name=std::get<0>(names_pixel_hitmap[i])+std::to_string(sec_number);
					if (j==0){
						name+=" station 210";
						id=STATION_210_M_ID;
					}
					else name+=" station 220";
					MapKey mk=MapKey(std::get<1>(names_pixel_hitmap[i]),id,sec_number,-1,-1);
					Histos_TH2F[mk]=dir_cyl_V[sec_number].make<TH2F>(name.c_str(), name.c_str(), sizes_hitmap[i][0], sizes_hitmap[i][1], sizes_hitmap[i][2], sizes_hitmap[i][3], sizes_hitmap[i][4], sizes_hitmap[i][5]);
			}
		}

		std::string name;
		Raw_resolution_DD_V_[sec_number] = dir_cyl_DDTI_V[sec_number].make<TGraph>(12);
		name="Raw_resolution_DD_sec_"+std::to_string(sec_number);
		Raw_resolution_DD_V_[sec_number] -> SetNameTitle(name.c_str(),name.c_str());

		SPC_resolution_DD_V_[sec_number] = dir_cyl_DDTI_V[sec_number].make<TGraph>(12);
		name="SPC_resolution_DD_sec_"+std::to_string(sec_number);
		SPC_resolution_DD_V_[sec_number] -> SetNameTitle(name.c_str(),name.c_str());

		Saturated_percentage_V_[sec_number] = dir_cyl_V[sec_number].make<TGraph>(12);
		name="Saturated_percentage_sec_"+std::to_string(sec_number);
		Saturated_percentage_V_[sec_number] -> SetNameTitle(name.c_str(),name.c_str());

	    name="Timing track time RAW sector "+std::to_string(sec_number);
		Tracks_time_RAW_histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), 1200, -60, 60);
	    name="Timing track time SPC sector "+std::to_string(sec_number);
		Tracks_time_SPC_histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), 1200, -60, 60);
	    name="Timing track resolution sector "+std::to_string(sec_number);
		Tracks_resolution_histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), 1000, 0, 1);
	}



	TFileDirectory  BaseDir = fs_->mkdir( "Calib/" );

	TCalibration_ = BaseDir.make<TTree>("CalibTree","rootuple");

  std::map<std::string,double> branches={{"Sector",this->cal_sector_},{"Plane",this->cal_plane_},{"Channel",this->cal_channel_},{"Par_0",this->cal_par_0_},
  {"Par_1",this->cal_par_1_},{"Par_2",this->cal_par_2_},{"Par_3",this->cal_par_3_},{"Chi_2",this->Chi_2_},{"Res_L1",this->resolution_L1_},
  {"Res_L2",this->resolution_L2_},{"Res_L2_3p",this->resolution_L2_3p_}};

  /*for (auto& b : branches){
    std::string str=b.first;
    if(str=="Sector" || str=="Plane" || str=="Channel"){
      str.append("/I");
      int var=(int) b.second;
      TCalibration_ -> Branch(str.c_str(),&(var),str.c_str());
    }
    else str.append("/D");
    TCalibration_ -> Branch(str.c_str(),&(b.second),str.c_str());
  }*/

	TCalibration_ -> Branch("Sector"  	, &(this->cal_sector_), "Sector/I");
	TCalibration_ -> Branch("Plane"		, &(this->cal_plane_), "Plane/I");
	TCalibration_ -> Branch("Channel"  	, &(this->cal_channel_), "Channel/I");
	TCalibration_ -> Branch("Par_0"		, &(this->cal_par_0_), "Par_0/D" );
	TCalibration_ -> Branch("Par_1"		, &(this->cal_par_1_), "Par_1/D" );
	TCalibration_ -> Branch("Par_2"		, &(this->cal_par_2_), "Par_2/D" );
	TCalibration_ -> Branch("Par_3"		, &(this->cal_par_3_), "Par_3/D" );
	TCalibration_ -> Branch("Chi_2"		, &(this->Chi_2_), "Chi_2/D" );
	TCalibration_ -> Branch("Res_L1"		, &(this->resolution_L1_), "Res_L1/D" );
	TCalibration_ -> Branch("Res_L2"		, &(this->resolution_L2_), "Res_L2/D" );
	TCalibration_ -> Branch("Res_L2_3p"		, &(this->resolution_L2_3p_), "Res_L2_3p/D" );


	TFileDirectory  GlobalInfoDir = fs_->mkdir( "GlobalInfo/" );
	RunNumber_Histo_ = GlobalInfoDir.make<TH1F>("RunNumber", "RunNumber", 300000, 300000, 330000);

}

// ------------ method called once each job just after ending the event loop  ------------
void
DiamondTimingAnalyzer::endJob()
{


	std::cout << "Starting end job task" << std::endl;
	int X_sat_start = Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,0,-1,0)]-> GetXaxis()->FindBin(15);
	int X_sat_stop = Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,0,-1,0)]-> GetNbinsX();
	int Y_sat_start = Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,0,-1,0)]-> GetYaxis()->FindBin(15);
	int Y_sat_stop = Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,0,-1,0)]-> GetNbinsY();

	///////////////////////////////////////////
    // deriving simple DD channels resolution
    ///////////////////////////////////////////

	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
	{
		for  (int ch_idx=0; ch_idx < CHANNELS_X_PLANE; ch_idx++)
		{
			if (Histos_TH1F[MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetEntries() >100)
			{
				Saturated_percentage_V_[sec_number] ->SetPoint(ch_idx+1, ch_idx, (Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,sec_number,-1,ch_idx)]->
				Integral(X_sat_start,X_sat_stop,Y_sat_start,Y_sat_stop))/(Histos_TH2F[MapKey("ToT_vs_ToT_Hmap_",-1,sec_number,-1,ch_idx)]->Integral()));

				Histos_TH1F[MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_idx)]->Fit("gaus","+","",-3,3);
				if (Histos_TH1F[MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")!= NULL)
					Raw_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, Histos_TH1F[MapKey("Raw_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")->GetParameter(2)/sqrt(2));
				else
				{
					std::cout << "WARNING: RAW DT fit unseccessfull for DD channel pair " << ch_idx << std::endl;
					Raw_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, 0);
				}

				Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->Fit("gaus","+Q","",-10,10);
				if (Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")!= NULL)
				{
					double SPC_mean = Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")->GetParameter(1);
					double SPC_sigma = Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")->GetParameter(2);
					Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->Fit("gaus","","",SPC_mean-(2.2*SPC_sigma),SPC_mean+(2.2*SPC_sigma));
					SPC_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, Histos_TH1F[MapKey("SPC_DT_Hmap_",-1,sec_number,-1,ch_idx)]->GetFunction("gaus")->GetParameter(2)/sqrt(2));
				}
				else
				{
					std::cout << "WARNING: SPC DT fit unseccessfull for DD channel pair " << ch_idx << std::endl;
					SPC_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, 0);
				}

			}
			else
			{
				Raw_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, 0);
				SPC_resolution_DD_V_[sec_number]->SetPoint(ch_idx+1, ch_idx, 0);
			}
		}
	}

	std::cout << "Analyzing channel couplings " << std::endl;

	////////////////////////////////////////////
    // computing track based dual L2 resolution
    ////////////////////////////////////////////

	for (const auto& Raw_histo_iter	: Raw_DT_TrackedPlane_Hmap_) // <<ChannelKey,ChannelKey>,histogram>
	{
		std::pair<ChannelKey,ChannelKey> histo_key = Raw_histo_iter.first;


		if(histo_key.first.sector != histo_key.second.sector)  std::cout << "WARNING: hits from different sectors in same tracks!!" << std::endl;


		int plane_1 = histo_key.first.plane;
		int plane_2 = histo_key.second.plane;
		int chann_1 = histo_key.first.channel;
		int chann_2 = histo_key.second.channel;

		int sector = histo_key.first.sector;

		if (Raw_DT_TrackedPlane_Hmap_[histo_key]->GetEntries()>100) // controlla se l'histo ha abbastanya  hit
		{
			Raw_DT_TrackedPlane_Hmap_[histo_key]->Fit("gaus","+","",-3,3);
			if (Raw_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")!= NULL)
				Raw_resolution_TP_V_[sector][std::make_pair(plane_1,plane_2)]->SetBinContent(chann_1+1, chann_2+1, Raw_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")->GetParameter(2)/sqrt(2));
			else
				std::cout << "WARNING: RAW DT fit unseccessfull for sector  " << sector <<" (plane " << plane_1 <<" channel " << chann_1 << " VS plane " << plane_2 <<" channel " << chann_2 << ")" << std::endl;


			SPC_DT_TrackedPlane_Hmap_[histo_key]->Fit("gaus","+Q","",-10,10);

			if (SPC_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")!= NULL)
			{
				double SPC_mean = SPC_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")->GetParameter(1);
				double SPC_sigma = SPC_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")->GetParameter(2);
				SPC_DT_TrackedPlane_Hmap_[histo_key]->Fit("gaus","","",SPC_mean-(2.2*SPC_sigma),SPC_mean+(2.2*SPC_sigma));
				SPC_resolution_TP_V_[sector][std::make_pair(plane_1,plane_2)]->SetBinContent(chann_1+1, chann_2+1, SPC_DT_TrackedPlane_Hmap_[histo_key]->GetFunction("gaus")->GetParameter(2)/sqrt(2));
			}
			else
				std::cout << "WARNING: SPC DT fit unseccessfull for sector  " << sector <<" (plane " << plane_1 <<" channel " << chann_1 << " VS plane " << plane_2 <<" channel " << chann_2 << ")" << std::endl;


		}
	}

    //////////////////////////////////////////////////////
    // extracting L2 precision (most populated criteria)
    //////////////////////////////////////////////////////

	for (int sec_id=0; sec_id < MAX_SECTOR_NUMBER; sec_id ++)
	{
		for (int plA_id=0; plA_id < PLANES_X_DETECTOR; plA_id +=2)
		{
			int plB_id = plA_id+1;
			for (int ch_id = 0; ch_id<CHANNELS_X_PLANE; ch_id++)
			{
				//select the bin in X
				Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetXaxis()->SetRangeUser(ch_id,ch_id+1);
				int max_bin = Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetMaximumBin();
				int max_value = Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetBinContent(max_bin);

				if (max_value > 100)
				{
					double res = SPC_resolution_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetBinContent(max_bin);
					Resolution_L1_map_[ChannelKey(sec_id,plA_id,ch_id)] = res;
				}
				else
				{
					Resolution_L1_map_[ChannelKey(sec_id,plA_id,ch_id)] = 0.400;
				}
				Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetXaxis()->SetRange(0,0);
				std::cout << "L1 resolution ("<<max_value<<" events) sec " << sec_id << " plane " << plA_id << " ch " << ch_id << ": " << Resolution_L1_map_[ChannelKey(sec_id,plA_id,ch_id)] << std::endl;

				//select the bin in Y
				Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetYaxis()->SetRangeUser(ch_id,ch_id+1);
				max_bin = Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetMaximumBin();
				max_value = Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetBinContent(max_bin);

				if (max_value > 100)
				{
					double res = SPC_resolution_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetBinContent(max_bin);
					Resolution_L1_map_[ChannelKey(sec_id,plB_id,ch_id)] = res;
				}
				else
				{
					Resolution_L1_map_[ChannelKey(sec_id,plB_id,ch_id)] = 0.400;
				}
				Correlation_TP_V_[sec_id][std::make_pair(plA_id,plB_id)]->GetYaxis()->SetRange(0,0);

				std::cout << "L1 resolution ("<<max_value<<" events) sec " << sec_id << " plane " << plB_id << " ch " << ch_id << ": " << Resolution_L1_map_[ChannelKey(sec_id,plB_id,ch_id)] << std::endl;

			}
		}
	}


    ///////////////////////////////////////////
    // deriving full track based L2 resolution
    ///////////////////////////////////////////

	for (int sec_id=0; sec_id < MAX_SECTOR_NUMBER; sec_id ++)
	{
		for (int pl_id=0; pl_id < PLANES_X_DETECTOR; pl_id ++)
		{
			for (int ch_id=0; ch_id < CHANNELS_X_PLANE; ch_id ++)
			{

				ChannelKey chk(sec_id,pl_id,ch_id);
        MapKey histo_key2=MapKey("Resolution_L2_histo_map_",-1,sec_id,pl_id,ch_id);

				if (Histos_TH1F.find(histo_key2)!=Histos_TH1F.end())
				{
					if (Histos_TH1F[histo_key2]->GetEntries() > 100)
					{
						Histos_TH1F[histo_key2]->Fit("gaus","+Q","",-10,10);

						if (Histos_TH1F[histo_key2]->GetFunction("gaus")!= NULL)
						{

							double ResL2_mean = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(1);
							double ResL2_sigma = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(2);
							Histos_TH1F[histo_key2]->Fit("gaus","","",ResL2_mean-(2.2*ResL2_sigma),ResL2_mean+(2.2*ResL2_sigma));
							ResL2_sigma = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(2);
							double Exp_sigma =  Histos_TH1F[MapKey("TrkExpectedRes_L2_histo_map_",-1,sec_id,pl_id,ch_id)]-> GetMean();
							Resolution_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, ResL2_sigma);

							if (ResL2_sigma > Exp_sigma)
							{

								ImprouvedRes_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2),0.5));
								Resolution_L2_map_[chk] = pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2),0.5);
							}
							else
							{
								ImprouvedRes_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0.05);
								Resolution_L2_map_[chk] = 0.05;
							}
						}
						else
						{
							Resolution_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0);
							Resolution_L2_map_[chk] = 0.400;
						}
					}
				}

				chk=ChannelKey(sec_id,pl_id,ch_id);
        histo_key2=MapKey("Resolution_L2_3p_histo_map_",-1,sec_id,pl_id,ch_id);

				if (Histos_TH1F.find(histo_key2)!=Histos_TH1F.end())
				{
					if (Histos_TH1F[histo_key2]->GetEntries() > 100)
					{
						Histos_TH1F[histo_key2]->Fit("gaus","+Q","",-10,10);

						if (Histos_TH1F[histo_key2]->GetFunction("gaus")!= NULL)
						{
							double ResL2_3p_mean = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(1);
							double ResL2_3p_sigma = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(2);
							Histos_TH1F[histo_key2]->Fit("gaus","","",ResL2_3p_mean-(2.2*ResL2_3p_sigma),ResL2_3p_mean+(2.2*ResL2_3p_sigma));
							ResL2_3p_sigma = Histos_TH1F[histo_key2]->GetFunction("gaus")->GetParameter(2);
							double Exp_sigma =  Histos_TH1F[MapKey("TrkExpectedRes_L2_3p_histo_map_",-1,sec_id,pl_id,ch_id)]-> GetMean();
							Resolution_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, ResL2_3p_sigma);
							if (ResL2_3p_sigma > Exp_sigma)
							{
								ImprouvedRes_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, pow(pow(ResL2_3p_sigma,2)-pow(Exp_sigma,2),0.5));
								Resolution_L2_3p_map_[chk] = pow(pow(ResL2_3p_sigma,2)-pow(Exp_sigma,2),0.5);
							}
							else
							{
								ImprouvedRes_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0.05);
								Resolution_L2_3p_map_[chk] = 0.05;
							}
						}
						else
						{
							std::cout << "WARNING: L2_3p resolution fit unseccessfull for sector  " << sec_id <<" plane " << pl_id <<" channel " << ch_id << std::endl;
							Resolution_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0);
							Resolution_L2_3p_map_[chk] = 0.400;
						}
					}
				}
			}

		}
	}

	///////////////////////////
    // deriving SPC corrections
    ///////////////////////////

    for(std::map<MapKey,TH2F*>::iterator iter=Histos_TH2F.begin(); iter!=Histos_TH2F.end(); ++iter)
    {

      if (iter->first.name!="TvsTOT_Distribution_") continue;
      ChannelKey recHitKey=ChannelKey(iter->first.sector,iter->first.plane,iter->first.channel);

		//TOT average
		Histos_TH1F[MapKey("Valid TOT mean sector ",-1,recHitKey.sector,recHitKey.plane,-1)]->SetBinContent(recHitKey.channel+1,
			Histos_TH1F[MapKey("ValidTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetMean());

		TOTvsT_Profile_map_[ recHitKey ] = channelDir_map_[ recHitKey ].make<TProfile>(*iter->second->ProfileX("_pfx",1,-1));
		if (iter->second->GetEntries() > 100 )
		{
			double Fit_UpRange = Histos_TH1F[MapKey("ValidTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetMean()+2.5;


			TF1 *myfermi = new TF1("MYfit","[0]/(exp((x-[1])/[2])+1)+[3]",10.5,25);
			//myfermi-> SetParameters(1.5,12.5,1,4);
			myfermi-> SetParameters(3*Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetRMS(),
									Histos_TH1F[MapKey("ValidTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetMean(),
									0.8,
									Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]
									->GetMean()-Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetRMS());

			myfermi->SetParLimits(1,9,15);
			myfermi->SetParLimits(2,0.2,2.5);

			std::cout << "starting parameters: " << 3*Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetRMS()
										<< " ; " << Histos_TH1F[MapKey("ValidTOT_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetMean()
										<< " ; " << 0.8
										<< " ; " << Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]
										->GetMean()-Histos_TH1F[MapKey("VALID_RAW_T_Distribution_",-1,recHitKey.sector,recHitKey.plane,recHitKey.channel)]->GetRMS()
										<< std::endl;

			TOTvsT_Profile_map_[ recHitKey ]->Fit(myfermi,"B+","",10.4,Fit_UpRange);

			//secondary fit
			/*TF1 *mygaus = new TF1("MYgaus","[3]-gaus",9.5,25);
			mygaus-> SetParameters(3*ValidT_Histo_map_[recHitKey]->GetRMS(),
									ValidTOT_Histo_map_[recHitKey]->GetMean()+ ValidTOT_Histo_map_[recHitKey]->GetRMS(),
									1.3,
									ValidT_Histo_map_[recHitKey]->GetMean()+ValidT_Histo_map_[recHitKey]->GetRMS());
			TOTvsT_Profile_map_[ recHitKey ]->Fit(mygaus,"B+","",9.5,Fit_UpRange);

			std::cout << "starting parameters gaus: " << 3*ValidT_Histo_map_[recHitKey]->GetRMS()
										<< " ; " << ValidTOT_Histo_map_[recHitKey]->GetMean()+ ValidTOT_Histo_map_[recHitKey]->GetRMS()
										<< " ; " << 1.3
										<< " ; " << ValidT_Histo_map_[recHitKey]->GetMean()+ValidT_Histo_map_[recHitKey]->GetRMS()
										<< std::endl;*/

			cal_sector_    = recHitKey.sector;
			cal_plane_  = recHitKey.plane;
			cal_channel_= recHitKey.channel;
			cal_par_0_  = myfermi-> GetParameter(0);
			cal_par_1_  = myfermi-> GetParameter(1);
			cal_par_2_  = myfermi-> GetParameter(2);
			cal_par_3_  = myfermi-> GetParameter(3);

			Chi_2_  = myfermi-> GetChisquare();

			if (Resolution_L1_map_.find(recHitKey ) != Resolution_L1_map_.end() )
				resolution_L1_ = Resolution_L1_map_[ recHitKey ];
			else
				resolution_L1_= 0.401;

			if (Resolution_L2_map_.find(recHitKey ) != Resolution_L2_map_.end() )
				resolution_L2_ = Resolution_L2_map_[ recHitKey ];
			else
				resolution_L2_= 0.999;

			if (Resolution_L2_3p_map_.find(recHitKey ) != Resolution_L2_3p_map_.end() )
				resolution_L2_3p_ = Resolution_L2_3p_map_[ recHitKey ];
			else
				resolution_L2_3p_= 0.999;


			TCalibration_ ->Fill();
		}
	}

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiamondTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingAnalyzer);
