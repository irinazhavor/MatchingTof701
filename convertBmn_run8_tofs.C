
#define BMNHYPNUCLPAIR_H
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;
using VecF2D = vector<vector<float>>;
//using RVecMap = RVec<TMap<int,double>>;
//vector<float> -> RVecF
//vector Int

TF1 *f1_m2_400_211 = new TF1( "m2_400_211", "pol1", 0, 10 );
TF1 *f1_sigma_400_211 = new TF1( "sigma_400_211", "pol2", 0, 10 );
TF1 *f1_m2_400_321 = new TF1( "m2_400_321", "pol1", 0, 10 );
TF1 *f1_sigma_400_321 = new TF1( "sigma_400_321", "pol2", 0, 10 );
TF1 *f1_m2_400_2212 = new TF1( "m2_400_2212", "pol1", 0, 10 );
TF1 *f1_sigma_400_2212 = new TF1( "sigma_400_2212", "pol2", 0, 10 );
TF1 *f1_m2_400_1000010020 = new TF1( "m2_400_1000010020", "pol1", 0, 10 );
TF1 *f1_sigma_400_1000010020 = new TF1( "sigma_400_1000010020", "pol2", 0, 10 );
TF1 *f1_m2_700_211 = new TF1( "m2_700_211", "pol1", 0, 10 );
TF1 *f1_sigma_700_211 = new TF1( "sigma_700_211", "pol2", 0, 10 );
TF1 *f1_m2_700_2212 = new TF1( "m2_700_2212", "pol1", 0, 10 );
TF1 *f1_sigma_700_2212 = new TF1( "sigma_700_2212", "pol2", 0, 10 );
TF1 *f1_m2_700_1000010020 = new TF1( "m2_700_1000010020", "pol1", 0, 10 );
TF1 *f1_sigma_700_1000010020 = new TF1( "sigma_700_1000010020", "pol2", 0, 10 );


BmnFieldMap* magField{nullptr};

TChain* makeChain(string& filename, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (filename.rfind(".root") < filename.size())
    chain->Add(filename.data());
  else {
    TFileCollection fc("fc", "", filename.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

bool isGoodSimTrack(const CbmMCTrack &track)
{
  if (track.GetStartZ() > 900) return false;
  return true;
}

XYZVector ExtrapolateStraightLine(const FairTrackParam *par, float z)
{
  float dz = z - par->GetZ();
  float x = par->GetX() + par->GetTx() * dz;
  float y = par->GetY() + par->GetTy() * dz;
  return {x, y, z};
}

double getMass(int pdg)
{
  if(pdg<1000000) 
    return TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  else 
    return 0.931*(pdg/10%1000);
}

int getCharge(int pdg)
{
  if(pdg<1000000)
    return TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
  else 
    return pdg/10000%1000;
}

vector<fourVector> simMomentum(const RVec<CbmMCTrack> tracks)
{
  vector<fourVector> momenta;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    TVector3 mom;
    track.GetMomentum(mom);
    double mass=getMass(track.GetPdgCode());
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),sqrt(mass*mass+mom.Mag2())});
  }
  return momenta;
}
 
vector<XYZTVector> simPosStart(const RVec<CbmMCTrack> tracks)
{
  vector<XYZTVector> pos;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pos.push_back({track.GetStartX(),track.GetStartY(),track.GetStartZ(),track.GetStartT()});
  }
  return pos;
}

RVec<int> simMotherId(const RVec<CbmMCTrack> tracks)
{
  vector<int> mothId;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    mothId.push_back(track.GetMotherId());
  }
  return mothId;
}

RVec<int> simPdg(const RVec<CbmMCTrack> tracks)
{
  vector<int> pdg;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pdg.push_back(track.GetPdgCode());
  }
  return pdg;
}

RVec<short> simCharge (const RVec<int> pdg)
{
  vector<short> ch;
  for (auto &p:pdg)
    ch.push_back(getCharge(p));
  return ch;
}

RVec<float> trackP(const RVec<BmnGlobalTrack> tracks)
try {
  vector<float> momenta;
  for (auto track:tracks)
    momenta.push_back(1./track.GetParamFirst()->GetQp());
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<fourVector> trackMomentum(const RVec<BmnGlobalTrack> tracks)
try {
  vector<fourVector> momenta;
  for (auto track:tracks) {
    auto *par = track.GetParamFirst();   
    TVector3 mom;
    par->Momentum(mom);
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),0});
  }   
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<short> recCharge(const RVec<BmnGlobalTrack> tracks)
try {
  vector<short> charge;
  for (auto track:tracks) {
    int q = track.GetParamFirst()->GetQp() > 0 ? 1 : -1;
    charge.push_back(q);
  }
  return charge;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> recDca(const RVec<BmnGlobalTrack> tracks, const CbmVertex vtx)
try {
  vector<XYZVector> dca;
  for (auto track:tracks) {
    auto par = track.GetParamFirst();
    dca.push_back({par->GetX()-vtx.GetX(),par->GetY()-vtx.GetY(),par->GetZ()-vtx.GetZ()});
  }
  return dca;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D covMatrix(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  VecF2D covariance_matrix;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);
    auto* par = track.GetParamFirst();
    covariance_matrix.emplace_back();
    for( int i=0; i<5; ++i ){
      for( int j=0; j<=i; ++j ){
        covariance_matrix.back().push_back( par->GetCovariance(i, j) );
      }
    }
    // Lower triangle of the symmetric covariance matrix
    // C[x, y, tx, ty, Qp]
    // { c_00, c1[0..1], c2[0..2], ... c4[0..4] }
  }
  return covariance_matrix;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D globalTrackCovMatrix(RVec<BmnGlobalTrack> global_tracks)
try {
  VecF2D covariance_matrix;
  for (auto& global_track : global_tracks) {
    auto* par = global_track.GetParamFirst();
    covariance_matrix.emplace_back();
    for( int i=0; i<5; ++i ){
      for( int j=0; j<=i; ++j ){
        covariance_matrix.back().push_back( par->GetCovariance(i, j) );
      }
    }
    // Lower triangle of the symmetric covariance matrix
    // C[x, y, tx, ty, Qp]
    // { c_00, c1[0..1], c2[0..2], ... c4[0..4] }
  }
  return covariance_matrix;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

float determinant3x3( const std::array<std::array<float, 3>, 3>& matrix ) try {
  auto x_0 = matrix[0][0] * ( matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]  );
  auto x_1 = matrix[0][1] * ( matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]  );
  auto x_2 = matrix[0][2] * ( matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]  );

  return x_0 - x_1 + x_2;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

std::array<float, 3> cramerFieldSolver3x3( std::array<float, 3> field, std::array<float, 3> coordinate ) try {
  // Solving the system of equation to extract parameters of quadratic extrapolation of the magnetic field
  // Ax = B
  // xi = detAi / detA
  std::array<std::array<float, 3>, 3> A;
  A[0] = {1.0f, 1.0f, 1.0f };
  A[1] = { coordinate[0], coordinate[1], coordinate[2] };
  A[2] = { coordinate[0]*coordinate[0], coordinate[1]*coordinate[1], coordinate[2]*coordinate[2] };

  auto A0 = A;
  A0[0] = field;
  auto A1 = A;
  A1[1] = field;
  auto A2 = A;
  A2[2] = field;

  auto detA = determinant3x3( A );
  auto detA0 = determinant3x3( A0 );
  auto detA1 = determinant3x3( A1 );
  auto detA2 = determinant3x3( A2 );

  auto p0 = detA0 / detA;
  auto p1 = detA1 / detA;
  auto p2 = detA2 / detA;

  return {p0, p1, p2};
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}


VecF2D magneticField(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks, RVec<CbmStsHit> sts_hits)
try {
  VecF2D magnetic_field;
  for (auto& global_track : global_tracks ) {

    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);

    std::array<float, 3> hit_z;
    std::array<float, 3> hit_bx;
    std::array<float, 3> hit_by;
    std::array<float, 3> hit_bz;

    for( int i=0; i<3; ++i ){
      // It seems size of the hitmap cannot be less than 4, but just to be safe
      if( i > track.GetStsHits()->GetSize() )
        magnetic_field.push_back( std::vector<float>(10, 0.0f) );

      auto sts_idx = track.GetStsHits()->At(i);
      auto x = sts_hits.at(sts_idx).GetX();
      auto y = sts_hits.at(sts_idx).GetY();
      auto z = sts_hits.at(sts_idx).GetZ();

      hit_z.at(i) = z;
      hit_bx.at(i) = magField->GetBx( x, y, z ); // kGs
      hit_by.at(i) = magField->GetBy( x, y, z ); // kGs
      hit_bz.at(i) = magField->GetBz( x, y, z ); // kGs
    }

    auto parameters_bx = cramerFieldSolver3x3( hit_bx, hit_z );
    auto parameters_by = cramerFieldSolver3x3( hit_by, hit_z );
    auto parameters_bz = cramerFieldSolver3x3( hit_bz, hit_z );

    magnetic_field.emplace_back();
    for( const auto& c : parameters_bx )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_by )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_bz )
      magnetic_field.back().push_back( c );
    magnetic_field.back().push_back( 0.0 ); // z0
  }
  return magnetic_field;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D stsTrackParameters(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  VecF2D parameters;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);

    auto* par = track.GetParamFirst();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D globalTrackParameters(RVec<BmnGlobalTrack> global_tracks)
try {
  VecF2D parameters;
  for (auto& global_track : global_tracks) {
    auto* par = global_track.GetParamFirst();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D trParamFirst(RVec<BmnGlobalTrack> global_tracks)
try {
  VecF2D parameters;
  for (auto& global_track : global_tracks) {
    auto* par = global_track.GetParamFirst();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
VecF2D trParamLast(RVec<BmnGlobalTrack> global_tracks)
try {
  vector<vector<float>> parameters;
  for (auto& global_track : global_tracks) {
    auto* par = global_track.GetParamLast();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D BeamTrackParameters(RVec<BmnTrack> beam_tracks)
try {
  VecF2D parameters;
  for (auto& beam_track : beam_tracks) {
    auto *par = beam_track.GetParamLast();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<fourVector> stsTrackMomentum(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  vector<fourVector> momenta;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);
    auto *par = track.GetParamFirst();
    TVector3 mom;
    par->Momentum(mom);
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),0});
  }
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> stsTrackChi2Ndf(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  RVec<float> vec_chi2;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);

    auto chi2 = track.GetChi2();
    auto ndf = track.GetNDF();

    vec_chi2.push_back( chi2/ndf );
  }
  return vec_chi2;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> stsTrackNdf(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  vector<int> vec_ndf;
  for (auto& global_track : global_tracks) {
        auto idx = global_track.GetGemTrackIndex();
        auto track = tracks.at(idx);

    auto ndf = track.GetNDF();

    vec_ndf.push_back( ndf );
  }
  return vec_ndf;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> stsTrackNhits(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
try {
  vector<int> vec_ndf;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);

    auto ndf = track.GetNStsHits();

    vec_ndf.push_back( ndf );
  }
  return vec_ndf;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> stsTrackSimPdg(ROOT::VecOps::RVec<int> sim_index, ROOT::VecOps::RVec<int> sim_pdg) try {
  vector<int> pdg;
  for( auto idx : sim_index ) {
    if( idx < 0 ) {
      pdg.push_back(-1);
      continue;
    }
    if( idx > sim_pdg.size() ) {
      pdg.push_back(-1);
      continue;
    }
    pdg.push_back(sim_pdg.at(idx));
  }
  return pdg;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

/// BeamHit
vector<XYZVector> beamHitXYZ(const RVec<BmnSiBTHit> tracks)
try {
  vector<XYZVector> pos;
  for (auto track:tracks){
    pos.push_back({track.GetX(), track.GetY(), track.GetZ()});
  }
  return pos;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> beamHitStation(const RVec<BmnSiBTHit> tracks)
try {
  vector<float> _station;
  for (auto track:tracks)
    _station.push_back(track.GetStation());
  return _station;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> beamHitIndex(const RVec<BmnSiBTHit> tracks)
try {
  vector<float> _index;
  for (auto track:tracks)
    _index.push_back(track.GetIndex());
  return _index;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> recSimIndex(const RVec<BmnGlobalTrack> recTracks, const RVec<CbmMCTrack> simTracks)
try {
  vector<int> newIndex;
  int shift=0;
  int nSimTracks = simTracks.size();
  for (int i=0;i<nSimTracks;i++) {
    if (!isGoodSimTrack(simTracks.at(i)))
    {
      shift++;
      newIndex.push_back(-1);
    }
    else
      newIndex.push_back(i-shift);
  }
  vector<int> simIndex;
  for (auto track:recTracks) {
    int oldIndex=track.GetRefIndex();
    if (oldIndex<0 || oldIndex>=nSimTracks)
      simIndex.push_back(-1);
    else
      simIndex.push_back(newIndex.at(oldIndex));
  }
  return simIndex;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> recPosLast(const RVec<BmnGlobalTrack> tracks)
try {
  vector<XYZVector> pos;
  for (auto track:tracks){
    auto par=track.GetParamLast();
    pos.push_back({par->GetX(), par->GetY(), par->GetZ()});
  }
  return pos;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> recPos450(const RVec<BmnGlobalTrack> tracks)
try {
  vector<XYZVector> pos;
  for (auto track:tracks) 
    pos.push_back(ExtrapolateStraightLine(track.GetParamLast(), 450));
  return pos;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> tofHitPosition(const TClonesArray hits)
try {
  vector<XYZVector> pos;
  for (const auto& hitObj:hits){
    auto hit=(BmnTofHit*)hitObj;
    pos.push_back({hit->GetX(),hit->GetY(),hit->GetZ()});
  }
  return pos;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

auto tofTrackPosition = [](const int coordIdx) { // o - x coordinate, 1 - y coordinate
  return [coordIdx](ROOT::VecOps::RVec<double> hitsRes, vector<vector<double>> hitsCoordinates)
  {
      std::vector<double> trCoordinate;
      trCoordinate.reserve(hitsCoordinates.size());

      // std::cout << "coordIdx: " << coordIdx << std::endl;
      // std::cout << "hitsCoordinates.size(): " << hitsCoordinates.size() << endl;
      // std::cout << "hitsRes.size(): " << hitsRes.size() << endl;
      // getchar();

      for (int i = 0; i < hitsCoordinates.size(); i++)
      {
          const auto hit_x_y = hitsCoordinates.at(i);
          if (fabs(hitsCoordinates.at(i).at(coordIdx)) > 300 || fabs(hitsRes.at(i)) > 100)
          {
              trCoordinate.push_back(-1000);
              continue;
          }
          trCoordinate.push_back(hit_x_y.at(coordIdx) + hitsRes.at(i));
      }
      // std::cout << "trCoordinate.size(): " << trCoordinate.size() << std::endl;
      return trCoordinate;
  };
};

vector<XYZVector> tofRes(const RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
try {
  vector<XYZVector> res;
  auto testHit=(BmnTofHit*)hits.At(0);
  if (testHit)
  {
    bool tof400 = testHit->GetZ()<550 ? true : false;
    for (auto track:tracks){
      int hitIndex = tof400 ? track.GetTof1HitIndex() : track.GetTof701HitIndex();
      if (hitIndex<0 || !hits.At(hitIndex)) 
      {
        res.push_back({-999,-999,-999});
        continue;
      }
      auto hit=(BmnTofHit*)hits.At(hitIndex);
      auto par=track.GetParamLast();
      TVector3 pos;
      par->Position(pos);
      auto posAtHitZ_ = ExtrapolateStraightLine(par, hit->GetZ());
      TVector3 posAtHitZ(posAtHitZ_.x(), posAtHitZ_.y(), posAtHitZ_.z());
      TVector3 hitPos;
      hit->Position(hitPos);
      auto posDiff = hitPos - posAtHitZ;
      res.push_back({posDiff.X(), posDiff.Y(), posDiff.Z()});
    }
  }
  return res;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> moduleId (const vector<XYZVector> modulePos)
try {
  vector <int> moduleIds;
  for (int i=0;i<modulePos.size();i++)
    moduleIds.push_back(i+1);
  return moduleIds;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> modulePos (const char *geoFile, const char *detectorTag)
try {
  bool verbose=false;
  map <int,XYZVector> modulePosMap;
  printf("Reading %s geometry from geometry file\n", detectorTag);
  TGeoManager* geoMan = TGeoManager::Import(geoFile, "FAIRGeom");
  if( !geoMan )
    throw runtime_error(Form("ERROR: No TGeoManager in file %s", geoFile));
  TGeoNode* caveNode = geoMan->GetTopNode();
  if( !caveNode )
    throw runtime_error(Form("ERROR: No cave node in file %s", geoFile));
  TGeoNode* detectorNode = nullptr;
  TString nodeName;
  
  bool nodeFound=false;
  for (int i = 0; i < caveNode->GetNdaughters(); i++) {
    detectorNode = caveNode->GetDaughter(i);
    nodeName = detectorNode->GetName();
    nodeName.ToLower();
    if (nodeName.Contains(detectorTag))
    {	
      nodeFound=true;
      break;
    }
  }
  if( !nodeFound )
    throw runtime_error(Form("ERROR: No detector node %s in cave", detectorTag));
  detectorNode = detectorNode->GetDaughter(0);

  auto geoMatrix = detectorNode->GetMatrix();
  auto geoBox = (TGeoBBox*) detectorNode->GetVolume()->GetShape();
  TVector3 frontFaceLocal(0, 0, -geoBox->GetDZ());
  TVector3 frontFaceGlobal;
  geoMatrix->LocalToMaster(&frontFaceLocal[0], &frontFaceGlobal[0]);

  nodeName=detectorNode->GetName();
  if (nodeName.Contains("box"))
    detectorNode = detectorNode->GetDaughter(detectorNode->GetNdaughters()-1);
  printf("%s node name: %s\n", detectorTag, detectorNode->GetName());

  int nModules = detectorNode->GetNdaughters();
  for (int i = 0; i < nModules; ++i) {
    auto* daughter = detectorNode->GetDaughter(i);
    auto geoMatrix = daughter->GetMatrix();
    TVector3 translation(geoMatrix->GetTranslation());

    int modId = daughter->GetNumber();
    double x  = translation.X();
    double y  = translation.Y();
    translation.SetZ(frontFaceGlobal.Z());
    double z  = translation.Z();
    modulePosMap.insert({modId, {x,y,z}});
  }

  geoMan->GetListOfVolumes()->Delete();
  geoMan->GetListOfShapes()->Delete();
  delete geoMan;
  nModules=modulePosMap.rbegin()->first;
  vector <XYZVector> modulePosVector(nModules,{0.,0.,0.});
  for(auto &modulePos:modulePosMap)
    modulePosVector.at(modulePos.first-1)=modulePos.second;
  if (verbose)
  {
    printf("%d module positions:\n", nModules);
    for(int i=0;i<nModules;i++)
      printf("%d: (%f, %f, %f)\n", i, modulePosVector.at(i).x(), modulePosVector.at(i).y(), modulePosVector.at(i).z());
  }
  return modulePosVector;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> fhcalModE(BmnFHCalEvent event)
try {
  vector<float> fhcalModEnergy_;
  for (int i = 0; i < 54; i++)
    fhcalModEnergy_.push_back(event.GetModule(i+1)->GetEnergy());
  return fhcalModEnergy_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

VecF2D fhcalSectionE(BmnFHCalEvent event)
try {
  VecF2D fhcalSecEnergy_;
  for (Int_t iModule = 1; iModule <= event.GetTotalModules(); iModule++)
  {
    BmnFHCalModule* module = event.GetModule(iModule);
    fhcalSecEnergy_.emplace_back();
    for (Int_t iSect = 1; iSect <= module->GetNsections(); iSect++)
    {
      fhcalSecEnergy_.back().push_back(module->GetSectionEnergy(iSect));
    }
  }
  return fhcalSecEnergy_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> scwallModQ(BmnScWallEvent event)
try {
  vector<float> scwallModCharge_;
  for (int i = 0; i < 174; i++)
    scwallModCharge_.push_back(event.GetCell(i+1)->GetSignal());
  return scwallModCharge_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> hodoModQ(BmnHodoEvent event)
try {
  vector<float> hodoModCharge_;
  for (int i = 0; i < 16; i++)
    hodoModCharge_.push_back(event.GetStrip(i+1)->GetSignal());
  return hodoModCharge_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> getEloss(const RVec<int> moduleId, const RVec<double> eLossDigis, const vector<XYZVector> modulePos)
try {
  int nModules=modulePos.size();
  vector<float> eLossModules(nModules,0);
  int nDigis=eLossDigis.size();
  for(int i=0;i<nDigis;i++)
    if(moduleId.at(i) <= nModules)
      eLossModules.at(moduleId.at(i)-1)=eLossDigis.at(i);
  return eLossModules;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> fdEloss (const RVec<BmnFDPoint> points)
{
  vector<float> eLoss;
  for (auto p:points)
    if (getCharge(p.GetPdgId()) > 0)
      eLoss.push_back(p.GetEnergyLoss());
  return eLoss;
}

RVec<bool> hasHitFhcal (RVec<CbmMCTrack> particles)
{
  vector<bool> hasHit;
  for(auto &part:particles)
    hasHit.push_back(part.GetNPoints(kFHCAL)>0);
  return hasHit;
}

RVec<short> modNhits (RVec<short> digiModIds, RVec<short> pointModIds)
{
  vector<short> nHits(digiModIds.size(),0);
  for (auto &pointModId:pointModIds)
    for (short i=0;i<digiModIds.size();i++)
      if(pointModId==digiModIds.at(i))
        nHits.at(i)++;
  return nHits;
}

RVec<float> mcPointEloss(const TClonesArray points)
{
  vector<float>el; 
  for (auto p:points)
    el.push_back(((FairMCPoint*)p)->GetEnergyLoss());
  return el;
}

int trigNSamples(const RVec<BmnTrigWaveDigit> trigger)
try {
  return (trigger.at(0)).GetNSamples();
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

float trigIntegral(const RVec<BmnTrigWaveDigit> trigger)
try {
  return (trigger.at(0)).GetIntegral(); // 
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

float trigAmp(const RVec<BmnTrigWaveDigit> trigger)
try {
  return (trigger.at(0)).GetPeak();
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<double> trigTdcTimes(const TClonesArray trigger)
try {
  return ((BmnTrigWaveDigit*)trigger.At(0))->TdcVector();
} catch( const std::exception& e ){
  std::cout << e.what() << std::endl;
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<short> trigValues(const RVec<BmnTrigWaveDigit> trigger)
try {
  vector<short>values;
  short* triggerValues = (trigger.at(0)).GetShortValue();
  for (int i = 0; i < 450; i++)
    values.push_back(triggerValues[i]);
  return values;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}


// m2
RVec<float> trM2(const vector<fourVector> mTr, const RVec<double> fBeta)
try {
  RVec<float> trM2_;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = mTr.at(itr).P();
    auto p2 = p*p;
    auto beta = fBeta.at(itr);
    auto beta2 = beta*beta;
    auto gamma2 = 1 - beta2;
    auto m2 = beta > -990. ? p2 / beta2 * gamma2 : -999.0;
    trM2_.push_back(m2);
  }
  return trM2_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}


// m2
RVec<double> trBetaTof701_(RVec<BmnGlobalTrack> global_tracks, const TClonesArray hits)
try {
  RVec<double> trBeta_;
  trBeta_.reserve(global_tracks.size());

  for (auto& global_track : global_tracks) {
    
    auto hitIndex = global_track.GetTof701HitIndex();
    
    if (hitIndex<0 || !hits.At(hitIndex)) {
      trBeta_.push_back(-999.);
      continue;
    }

    BmnTofHit* tof_hit = (BmnTofHit*)hits.UncheckedAt(hitIndex);
    Int_t mod = ((tof_hit->GetDetectorID() & 0x0000FF00) >> 8) - 1;

    if(mod==30){
      trBeta_.push_back(-999.);
      continue;
    }
    
    trBeta_.push_back(global_track.GetBeta(3));//beta from tof-701
  }
  return trBeta_;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> ModHitTof701(RVec<BmnGlobalTrack> global_tracks, const TClonesArray hits)
try {
  RVec<int> parameters;
  parameters.reserve(global_tracks.size());
  for (auto& global_track : global_tracks) {
    auto hitIndex = global_track.GetTof701HitIndex();
    if (hitIndex<0 || !hits.At(hitIndex)) {
      parameters.push_back(-999);
      continue;
    }
    BmnTofHit* hit = (BmnTofHit*)hits.At(hitIndex);
    Int_t mod = ((hit->GetDetectorID() & 0x0000FF00) >> 8) - 1;
    parameters.push_back(mod);

  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> ModHitTof400(RVec<BmnGlobalTrack> global_tracks, const TClonesArray hits)
try {
  RVec<int> parameters;
  parameters.reserve(global_tracks.size());
  for (auto& global_track : global_tracks) {
    auto hitIndex = global_track.GetTof1HitIndex();
    if (hitIndex<0 || !hits.At(hitIndex)) {
      parameters.push_back(-999);
      continue;
    }
    BmnTofHit* hit = (BmnTofHit*)hits.At(hitIndex);
    Int_t mod = ((hit->GetDetectorID() & 0x0000FF00) >> 8) - 1;
    parameters.push_back(mod);

  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> TofHitCoord701(RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
try
{
  vector<XYZVector> coords;
  coords.reserve(tracks.size());

    for (auto &track : tracks)
    {
      int hitIndex = track.GetTof701HitIndex();

      if (hitIndex < 0 || !hits.At(hitIndex))
      {
        coords.push_back({-999, -999, -999});
        continue;
      }

      BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
      coords.push_back({hit->GetX(), hit->GetY(), 0});
  }
  return coords;
}
catch (const std::exception &e)
{
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> TofHitCoord400(RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
try
{
  vector<XYZVector> coords;
  coords.reserve(tracks.size());

    for (auto &track : tracks)
    {
      int hitIndex = track.GetTof1HitIndex();

      if (hitIndex < 0 || !hits.At(hitIndex))
      {
        coords.push_back({-999, -999, -999});
        continue;
      }

      BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
      coords.push_back({hit->GetX(), hit->GetY(), 0});
  }
  return coords;
}
catch (const std::exception &e)
{
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> tofResFromHit701(const RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
try
{
  vector<XYZVector> res;
  res.reserve(tracks.size());
    for (auto &track : tracks)
    {
      int hitIndex = track.GetTof701HitIndex();
      if (hitIndex < 0 || !hits.At(hitIndex))
      {
        res.push_back({-999, -999, -999});
        continue;
      }
      BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
      res.push_back({hit->GetResX(), hit->GetResY(), 0});
    }
  return res;
}
catch (const std::exception &e)
{
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> tofResFromHit400(const RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
try
{
  vector<XYZVector> res;
  res.reserve(tracks.size());
    for (auto &track : tracks)
    {
      int hitIndex = track.GetTof1HitIndex();
      if (hitIndex < 0 || !hits.At(hitIndex))
      {
        res.push_back({-999, -999, -999});
        continue;
      }
      BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
      res.push_back({hit->GetResX(), hit->GetResY(), 0});
    }
  return res;
}
catch (const std::exception &e)
{
  std::cout << __func__ << std::endl;
  throw e;
}

// nSigma 211 TOF-400
RVec<float> nSigmaM2Tof400_211(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_400_211->Eval(p);
    auto sigma = f1_sigma_400_211->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 2212 TOF-400
RVec<float> nSigmaM2Tof400_2212(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_400_2212->Eval(p);
    auto sigma = f1_sigma_400_2212->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 321 TOF-400
RVec<float> nSigmaM2Tof400_321(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_400_321->Eval(p);
    auto sigma = f1_sigma_400_321->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 1000010020 TOF-400
RVec<float> nSigmaM2Tof400_1000010020(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_400_1000010020->Eval(p);
    auto sigma = f1_sigma_400_1000010020->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 211 TOF-700
RVec<float> nSigmaM2Tof700_211(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_700_211->Eval(p);
    auto sigma = f1_sigma_700_211->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 2212 TOF-700
RVec<float> nSigmaM2Tof700_2212(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_700_2212->Eval(p);
    auto sigma = f1_sigma_700_2212->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}
// nSigma 1000010020 TOF-700
RVec<float> nSigmaM2Tof700_1000010020(const vector<fourVector> mTr, const RVec<float> mTrM2)
try {
  RVec<float> nSigma;
  for (int itr=0; itr<mTr.size(); itr++) {
    auto p = abs(mTr.at(itr).P());
    auto _m2 = mTrM2.at(itr);
    auto mean = f1_m2_700_1000010020->Eval(p);
    auto sigma = f1_sigma_700_1000010020->Eval(p);
    auto _nsigma = _m2!=-999.0 ? (_m2 - mean) / sigma : -999.0;
    nSigma.push_back(_nsigma);
  }
  return nSigma;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

// pileup
int CentralHitIndexBC1S(BmnBC1hitInfo hit)
try {
  return hit.GetCentralHitIndexBC1S();
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

int BC1hitClasses(BmnBC1hitInfo hit, int centralHitIndex)
try {
  if (centralHitIndex < 0) return -999;
  return hit.GetBC1hitClasses().at(centralHitIndex)==BmnEventClass::k1;
} catch( const std::exception& e ){
  std::cout << e.what() << std::endl;
  std::cout << __func__ << std::endl;
  throw e;
}

std::vector<float> ClosestBC1hitsDt_k0(BmnBC1hitInfo hit)
try {
  return hit.GetClosestBC1hitsDt(BmnEventClass::k0);
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

std::vector<float> ClosestBC1hitsDt_kV0(BmnBC1hitInfo hit)
try {
  return hit.GetClosestBC1hitsDt(BmnEventClass::kV0);
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

std::vector<float> ClosestBC1hitsDt_k1(BmnBC1hitInfo hit)
try {
  return hit.GetClosestBC1hitsDt(BmnEventClass::k1);
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

std::vector<float> ClosestBC1hitsDt_kV1(BmnBC1hitInfo hit)
try {
  return hit.GetClosestBC1hitsDt(BmnEventClass::kV1);
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

// inTofTracks - independent matching file from the Matching_Tofs.C macro
// inReco      - dst file from eos
// inDigi      - digi file from eos

void convertBmn_run8_tofs(string inTofTracks = "" , string inReco="data/run8/rec.root", string inDigi="data/run8/digi.root",  std::string fileOut = "out.tree.root")
{
  
  TChain *chainRec=makeChain(inReco, "bmndata");
  TChain *chainDigi=makeChain(inDigi, "bmndata");
  TChain *chainTofGlobalTracks=makeChain(inTofTracks, "TofMatch");
  chainRec->AddFriend(chainDigi);
  chainRec->AddFriend(chainTofGlobalTracks);
  ROOT::RDataFrame d(*chainRec);

  int nEvents = chainRec->GetEntries();

  DstRunHeader* run_header = dynamic_cast<DstRunHeader*>( chainRec->GetCurrentFile()->Get("DstRunHeader") );
  if (run_header) {
    cout << "\n|||||||||||||||| RUN SUMMARY |||||||||||||||" << endl;
    cout << "||\t\t\t\t\t  ||" << endl;
    cout << "||   Period:        " << run_header->GetPeriodNumber() << "\t\t\t  ||" << endl;
    cout << "||   Number:        " << run_header->GetRunNumber() << "\t\t  ||" << endl;
    cout << "||   Start Time:    " << run_header->GetStartTime().AsString("s") << "\t  ||" << endl;
    cout << "||   End Time:      " << run_header->GetFinishTime().AsString("s") << "\t  ||" << endl;
    cout << "||   Beam:          A = " << run_header->GetBeamA() << ", Z = " << run_header->GetBeamA() << "\t  ||" << endl;
    cout << "||   Beam energy:   " << run_header->GetBeamEnergy() << " GeV\t\t  ||" << endl;
    cout << "||   Target:        A = " << run_header->GetTargetA() << ", Z = " << run_header->GetTargetZ() << "\t  ||" << endl;
    //cout << "||   Field voltage: " << setprecision(4) << run_header->GetMagneticField() << " mV\t\t  ||" << endl;
    cout << "||\t\t\t\t\t  ||" << endl;
    cout << "||||||||||||||||||||||||||||||||||||||||||||\n" << endl;
  }
  
  auto run_id = run_header->GetRunNumber();

  gRandom->SetSeed(0);
  TString geoFileName = Form("current_geo_file_%d.root", UInt_t(gRandom->Integer(UINT32_MAX)));
  Int_t res_code = UniRun::ReadGeometryFile(run_header->GetPeriodNumber(), run_header->GetRunNumber(), (char*)geoFileName.Data());
  if (res_code != 0) {
    cout << "ERROR: could not read geometry file from the database" << endl;
    exit(-3);
  }

  // get gGeoManager from ROOT file (if required)
  TFile* geoFile = new TFile(geoFileName, "READ");
  if (!geoFile->IsOpen()) {
    cout << "ERROR: could not open ROOT file with geometry: " + geoFileName << endl;
    exit(-4);
  }
  
  UniRun* pCurrentRun = UniRun::GetRun(run_header->GetPeriodNumber(), run_header->GetRunNumber());
  if (pCurrentRun == 0)
    exit(-6);
  
  Double_t* field_voltage = pCurrentRun->GetFieldVoltage();
  if (field_voltage == NULL) {
    cout << "ERROR: no field voltage was found for run " << run_header->GetPeriodNumber() << ":" <<  run_header->GetRunNumber() << endl;
    exit(-7);
  }

  Double_t map_current = 112.0;// run7 = 55.87; run8 =112.0
  Double_t fieldScale = 0.;
  if (*field_voltage < 10) {
    fieldScale = 0;
  } else
    fieldScale = (*field_voltage) / map_current;

  //magField = new BmnNewFieldMap("field_sp41v5_ascii_Extrap.root");
  //FieldMap_1900_extrap_noPed
  magField = new BmnNewFieldMap("FieldMap_1900_extrap_noPed.root");
  magField->SetScale(fieldScale);
  magField->Init();

  auto scwallModPos=modulePos(geoFileName.Data(),"scwall");
  auto hodoModPos=modulePos(geoFileName.Data(),"hodo");
  auto fhcalModPos=modulePos(geoFileName.Data(),"fhcal");

  //pion TOF-400
  //TF1 *f1_m2_400_211 = new TF1( "m2_400_211", "pol1", 0, 10 );
  f1_m2_400_211->SetParameter( 0, 0.03911 );
  f1_m2_400_211->SetParameter( 1, -0.02269 );
  //TF1 *f1_sigma_400_211 = new TF1( "sigma_400_211", "pol2", 0, 10 );
  f1_sigma_400_211->SetParameter( 0, 0.04135 );
  f1_sigma_400_211->SetParameter( 1, -0.03773 );
  f1_sigma_400_211->SetParameter( 2, 0.02134 );
  //kaon TOF-400
  //TF1 *f1_m2_400_321 = new TF1( "m2_400_321", "pol1", 0, 10 );
  f1_m2_400_321->SetParameter( 0, 0.2181 );
  f1_m2_400_321->SetParameter( 1, 0.001171 );
  //TF1 *f1_sigma_400_321 = new TF1( "sigma_400_321", "pol2", 0, 10 );
  f1_sigma_400_321->SetParameter( 0, 0.2596 );
  f1_sigma_400_321->SetParameter( 1, -0.2067 );
  f1_sigma_400_321->SetParameter( 2, 0.06465 );
  //proton TOF-400
  //TF1 *f1_m2_400_2212 = new TF1( "m2_400_2212", "pol1", 0, 10 );
  f1_m2_400_2212->SetParameter( 0, 0.9461 );
  f1_m2_400_2212->SetParameter( 1, -0.02658 );
  //TF1 *f1_sigma_400_2212 = new TF1( "sigma_400_2212", "pol2", 0, 10 );
  f1_sigma_400_2212->SetParameter( 0, 0.08764 );
  f1_sigma_400_2212->SetParameter( 1, -0.02868 );
  f1_sigma_400_2212->SetParameter( 2, 0.01334 );
  //d/He4 TOF-400
  //TF1 *f1_m2_400_1000010020 = new TF1( "m2_400_1000010020", "pol1", 0, 10 );
  f1_m2_400_1000010020->SetParameter( 0, 3.684 );
  f1_m2_400_1000010020->SetParameter( 1, -0.02706 );
  //TF1 *f1_sigma_400_1000010020 = new TF1( "sigma_400_1000010020", "pol2", 0, 10 );
  f1_sigma_400_1000010020->SetParameter( 0, 0.3756 );
  f1_sigma_400_1000010020->SetParameter( 1, -0.1047 );
  f1_sigma_400_1000010020->SetParameter( 2, 0.02066 );
  // TOF-700
  //pion TOF-700
  //TF1 *f1_m2_700_211 = new TF1( "m2_700_211", "pol1", 0, 10 );
  f1_m2_700_211->SetParameter( 0, -0.07376 );
  f1_m2_700_211->SetParameter( 1, 0.08399 );
  //TF1 *f1_sigma_700_211 = new TF1( "sigma_700_211", "pol2", 0, 10 );
  f1_sigma_700_211->SetParameter( 0, 0.01582 );
  f1_sigma_700_211->SetParameter( 1, -0.01641 );
  f1_sigma_700_211->SetParameter( 2, 0.03549 );
  //proton TOF-700
  //TF1 *f1_m2_700_2212 = new TF1( "m2_700_2212", "pol1", 0, 10 );
  f1_m2_700_2212->SetParameter( 0, 0.919 );
  f1_m2_700_2212->SetParameter( 1, 0.0007007 );
  //TF1 *f1_sigma_700_2212 = new TF1( "sigma_700_2212", "pol2", 0, 10 );
  f1_sigma_700_2212->SetParameter( 0, 0.01756 );
  f1_sigma_700_2212->SetParameter( 1, 0.02298 );
  f1_sigma_700_2212->SetParameter( 2, 0.01183 );
  //d/He4 TOF-700
  //TF1 *f1_m2_700_1000010020 = new TF1( "m2_700_1000010020", "pol1", 0, 10 );
  f1_m2_700_1000010020->SetParameter( 0, 3.692 );
  f1_m2_700_1000010020->SetParameter( 1, -0.0268 );
  //TF1 *f1_sigma_700_1000010020 = new TF1( "sigma_700_1000010020", "pol2", 0, 10 );
  f1_sigma_700_1000010020->SetParameter( 0, 0.2818 );
  f1_sigma_700_1000010020->SetParameter( 1, -0.111 );
  f1_sigma_700_1000010020->SetParameter( 2, 0.03891 );

  auto dd=d
    // pileup
    .Define("centralHitIndex",CentralHitIndexBC1S,{"BmnBC1hitInfo."})
    //.Filter("centralHitIndex>=0")
    .Define("k1central",BC1hitClasses,{"BmnBC1hitInfo.","centralHitIndex"})
    .Define("k0closest",  ClosestBC1hitsDt_k0, {"BmnBC1hitInfo."})
    .Define("kV0closest", ClosestBC1hitsDt_kV0,{"BmnBC1hitInfo."})
    .Define("k1closest",  ClosestBC1hitsDt_k1, {"BmnBC1hitInfo."})
    .Define("kV1closest", ClosestBC1hitsDt_kV1,{"BmnBC1hitInfo."})
    .Define("k0before", "return (k0closest.at(0)>-200 || kV0closest.at(0)>-200)")
    .Define("k1before", "return (k1closest.at(0)>-3000 || kV1closest.at(0)>-3000)")
    .Define("k0after", "return (k0closest.at(1)<200 || kV0closest.at(1)<200)")
    .Define("k1after", "return (k1closest.at(1)<700 || kV1closest.at(1)<700)")
    .Define("noPileup", "if(k1central && !k1before && !k1after && !k0before && !k0after) return true; else return false;")
    //
    .Define("vtxNtracks","PrimaryVertex.fNTracks")
    //.Filter("vtxNtracks > 1")
     //.Range(1000,1010)
    .Define("runId",[run_id](){ return run_id; }, {} )
    .Define("evtId","DstEventHeader.fEventId")
    .Define("triggerMapBR","BmnTrigInfo.fInputsBR")
    .Define("triggerMapAR","BmnTrigInfo.fInputsAR")
    .Define("bc1sNSamples",trigNSamples,{"TQDC_BC1S"})
    .Define("bc1sIntegral",trigIntegral,{"TQDC_BC1S"})
    .Define("bc1sAmplitude",trigAmp,{"TQDC_BC1S"})
    .Define("bc1sTdcValues",trigValues,{"TQDC_BC1S"})
    .Define("bc1sTdcTimes",trigTdcTimes,{"TQDC_BC1S"})
    .Define("bc1tNSamples",trigNSamples,{"TQDC_BC1T"})
    .Define("bc1tIntegral",trigIntegral,{"TQDC_BC1T"})
    .Define("bc1tAmplitude",trigAmp,{"TQDC_BC1T"})
    .Define("bc1tTdcValues",trigValues,{"TQDC_BC1T"})
    .Define("bc1tTdcTimes",trigTdcTimes,{"TQDC_BC1T"})
    .Define("bc1bNSamples",trigNSamples,{"TQDC_BC1B"})
    .Define("bc1bIntegral",trigIntegral,{"TQDC_BC1B"})
    .Define("bc1bAmplitude",trigAmp,{"TQDC_BC1B"})
    .Define("bc1bTdcValues",trigValues,{"TQDC_BC1B"})
    .Define("bc1bTdcTimes",trigTdcTimes,{"TQDC_BC1B"})
    .Define("vcsNSamples",trigNSamples,{"TQDC_VCS"})
    .Define("vcsIntegral",trigIntegral,{"TQDC_VCS"})
    .Define("vcsAmplitude",trigAmp,{"TQDC_VCS"})
    .Define("vcsTdcValues",trigValues,{"TQDC_VCS"})
    .Define("vcsTdcTimes",trigTdcTimes,{"TQDC_VCS"})
    .Define("vctNSamples",trigNSamples,{"TQDC_VCT"})
    .Define("vctIntegral",trigIntegral,{"TQDC_VCT"})
    .Define("vctAmplitude",trigAmp,{"TQDC_VCT"})
    .Define("vctTdcValues",trigValues,{"TQDC_VCT"})
    .Define("vctTdcTimes",trigTdcTimes,{"TQDC_VCT"})
    .Define("vcbNSamples",trigNSamples,{"TQDC_VCB"})
    .Define("vcbIntegral",trigIntegral,{"TQDC_VCB"})
    .Define("vcbAmplitude",trigAmp,{"TQDC_VCB"})
    .Define("vcbTdcValues",trigValues,{"TQDC_VCB"})
    .Define("vcbTdcTimes",trigTdcTimes,{"TQDC_VCB"})
    .Define("fdNSamples",trigNSamples,{"TQDC_FD"})
    .Define("fdIntegral",trigIntegral,{"TQDC_FD"})
    .Define("fdAmplitude",trigAmp,{"TQDC_FD"})
    .Define("fdTdcValues",trigValues,{"TQDC_FD"})
    .Define("fdTdcTimes",trigTdcTimes,{"TQDC_FD"})   
    .Define("bdMult","BmnTrigInfo.fBDMult")
    .Define("bdModId", "BD.fMod")
    .Define("bdModAmp", "BD.fAmp") 
    .Define("simdMult","BmnTrigInfo.fSiMDMult")
    .Define("vtxX","PrimaryVertex.fX")
    .Define("vtxY","PrimaryVertex.fY")
    .Define("vtxZ","PrimaryVertex.fZ")
    // .Define("vtxNtracks","PrimaryVertex.fNTracks")
    .Define("vtxChi2","PrimaryVertex.fChi2")
    .Define("vtxNdf","PrimaryVertex.fNDF")
    .Define("trMom",trackMomentum,{"BmnGlobalTrack"})
    .Define("trNhits","BmnGlobalTrack.fNhits")
    .Define("trNdf","BmnGlobalTrack.fNDF")
    .Define("trChi2","BmnGlobalTrack.fChi2")
    .Define("trP",trackP,{"BmnGlobalTrack"})
    .Define("trChi2vtx","BmnGlobalTrack.fChi2InVertex")
    .Define("trLength","BmnGlobalTrack.fLength")
    .Define("trCharge",recCharge,{"BmnGlobalTrack"})
    .Define("trDca",recDca,{"BmnGlobalTrack","PrimaryVertex."})
    .Define("trTof400hit","BmnGlobalTrack.fTof1Hit")
    .Define("trTof700hit","BmnGlobalTrack.fTof2Hit")
    .Define("trTof701hit","BmnGlobalTrack.fTof701Hit")
    .Define("trBeta400","BmnGlobalTrack.fBeta400")
    .Define("trBeta700","BmnGlobalTrack.fBeta700")
    .Define("trBeta701","BmnGlobalTrack.fBeta701")
    .Define("trBeta701_2", trBetaTof701_, {"BmnGlobalTrack", "BmnTof701Hit"})

    .Define("trBeta400_mod",ModHitTof400,{"BmnGlobalTrack","BmnTof400Hit"})
    .Define("trBeta701_mod",ModHitTof701,{"BmnGlobalTrack","BmnTof701Hit"})

    .Define("trBeta400Hit_coors", TofHitCoord400, {"BmnGlobalTrack", "BmnTof400Hit"})
    .Define("trBeta701Hit_coors", TofHitCoord701, {"BmnGlobalTrack", "BmnTof701Hit"})
    
    // separate procedure of sts tracks extrapolation and matching to the tofs


    .Define("vtxX_m", "Vertex.fX")
    .Define("vtxY_m", "Vertex.fY")
    .Define("vtxZ_m", "Vertex.fZ")
    .Define("vtxNtracks_m","Vertex.fNTracks")
    .Define("vtxChi2_m","Vertex.fChi2")
    .Define("vtxNdf_m","Vertex.fNDF")
    .Define("trMom_m",trackMomentum,{"Track"})
    .Define("trNhits_m","Track.fNhits")
    .Define("trP_m",trackP,{"Track"})
    .Define("trLength_m","Track.fLength")
    .Define("trCharge_m",recCharge,{"Track"})
    .Define("trPosLast_m",recPosLast,{"Track"})

    .Define("trParamFirst_m", trParamFirst, { "Track" })
    .Define("trParamLast_m", trParamLast, { "Track" })
    .Define("globalTrackParameters_m", globalTrackParameters, { "Track" })
    .Define("trDca_m",recDca,{"Track","Vertex"})

    .Define("trTof400hit_m","Track.fTof1Hit")
    .Define("trTof700hit_m","Track.fTof2Hit")
    .Define("trTof701hit_m","Track.fTof701Hit")
    .Define("trBeta400_m","Track.fBeta400")
    .Define("trBeta700_m","Track.fBeta700")
    .Define("trBeta701_m","Track.fBeta701")
    .Define("trBeta701_2_m", trBetaTof701_, {"Track", "Tof701"})

    .Define("trBeta701_mod_m",ModHitTof701,{"Track","Tof701"})
    .Define("trBeta400_mod_m",ModHitTof400,{"Track","Tof400"})

    .Define("trBeta400Hit_coors_m", TofHitCoord400, {"Track", "Tof400"})
    .Define("trBeta701Hit_coors_m", TofHitCoord701, {"Track", "Tof701"})

    .Define("trM2Tof400_m",trM2,{"trMom_m","trBeta400_m"})
    .Define("trM2Tof700_m",trM2,{"trMom_m","trBeta700_m"})
    .Define("trM2Tof701_m",trM2,{"trMom_m","trBeta701_m"})
    .Define("trM2Tof701_2_m",trM2,{"trMom_m","trBeta701_2_m"})

    //
    
    .Define("trPosLast",recPosLast,{"BmnGlobalTrack"})
    .Define("trPos450",recPos450,{"BmnGlobalTrack"})
    .Define("gemDigits","GEM.fUniqueID")
    .Define("stsDigits","SILICON.fUniqueID")
    .Define("stsTrackCovMatrix", covMatrix, { "BmnGlobalTrack", "StsVector" })
    .Define("stsTrackMagField", magneticField, { "BmnGlobalTrack", "StsVector", "StsHit" })
    .Define("stsTrackParameters", stsTrackParameters, { "BmnGlobalTrack", "StsVector" })
    .Define("trParamFirst", trParamFirst, { "BmnGlobalTrack" })
    .Define("trParamLast", trParamLast, { "BmnGlobalTrack" })
    .Define("globalTrackParameters", globalTrackParameters, { "BmnGlobalTrack" })
    .Define("globalTrackCovMatrix", globalTrackCovMatrix, { "BmnGlobalTrack" })
    .Define("stsTrackMomentum", stsTrackMomentum, { "BmnGlobalTrack", "StsVector" })
    .Define("stsTrackChi2Ndf", stsTrackChi2Ndf, { "BmnGlobalTrack", "StsVector" })
    .Define("stsTrackNdf", stsTrackNdf, { "BmnGlobalTrack", "StsVector" })
    .Define("stsTrackNhits", stsTrackNhits, { "BmnGlobalTrack", "StsVector" })
    //nSigma from TOF-400 and TOF-700
    .Define("trM2Tof400",trM2,{"trMom","trBeta400"})
    .Define("trM2Tof700",trM2,{"trMom","trBeta700"})
    .Define("trM2Tof701",trM2,{"trMom","trBeta701"})
    .Define("trM2Tof701_2",trM2,{"trMom","trBeta701_2"})

    // .Define("nSigmaM2Tof400_211", nSigmaM2Tof400_211, {"trMom","trM2Tof400"})
    // .Define("nSigmaM2Tof400_321", nSigmaM2Tof400_321, {"trMom","trM2Tof400"})
    // .Define("nSigmaM2Tof400_2212",nSigmaM2Tof400_2212,{"trMom","trM2Tof400"})
    // .Define("nSigmaM2Tof400_1000010020",nSigmaM2Tof400_1000010020,{"trMom","trM2Tof400"})
    // .Define("nSigmaM2Tof700_211", nSigmaM2Tof700_211, {"trMom","trM2Tof700"})
    // .Define("nSigmaM2Tof700_2212",nSigmaM2Tof700_2212,{"trMom","trM2Tof700"})
    // .Define("nSigmaM2Tof700_1000010020",nSigmaM2Tof700_1000010020,{"trMom","trM2Tof700"})

    .Define("beamHitXYZ", beamHitXYZ, { "BmnSiBTHit" })
    .Define("beamHitStation", beamHitStation, { "BmnSiBTHit" })
    .Define("beamHitIndex",   beamHitIndex,   { "BmnSiBTHit" })
    .Define("beamTrackChi2","BmnBeamTrack.fChi2")
    .Define("beamTrackNDF", "BmnBeamTrack.fNDF")
    .Define("beamTrackB",   "BmnBeamTrack.fB")
    .Define("beamTrackParameters", BeamTrackParameters, { "BmnBeamTrack" })

    .Define("tof400Digits","TOF400.fUniqueID")
    .Define("tof400Plane","TOF400.fPlane")
    .Define("tof400Strip","TOF400.fStrip")
    //.Define("tof400Side", "TOF400.fSide")
    .Define("tof400hitPos",tofHitPosition,{"BmnTof400Hit"})
    .Define("tof400hitT","BmnTof400Hit.fTimeStamp")
    .Define("tof400hitL","BmnTof400Hit.fLength")
    .Define("tof400hitResX","BmnTof400Hit.fResX")
    .Define("tof400hitResY","BmnTof400Hit.fResY")
    .Define("tof400hitRefIndex","BmnTof400Hit.fRefIndex")
    .Define("trTof400hitRes2D", tofResFromHit400, {"BmnGlobalTrack", "BmnTof400Hit"})


    .Define("tof700Digits","TOF700.fUniqueID")
    .Define("tof700Plane","TOF700.fPlane")
    .Define("tof700Strip","TOF700.fStrip")
    //.Define("tof700Side", "TOF700.fSide")
    .Define("tof700hitPos",tofHitPosition,{"BmnTof700Hit"})
    .Define("tof700hitT","BmnTof700Hit.fTimeStamp")
    .Define("tof700hitL","BmnTof700Hit.fLength")
    .Define("tof700hitResX","BmnTof700Hit.fResX")
    .Define("tof700hitResY","BmnTof700Hit.fResY")
    .Define("tof700hitRefIndex","BmnTof700Hit.fRefIndex")
 
    .Define("tof701Digits","TOF701.fUniqueID")
    .Define("tof701Plane","TOF701.fPlane")
    .Define("tof701Strip","TOF701.fStrip")
    //.Define("tof701Side", "TOF701.fSide")
    .Define("tof701hitPos",tofHitPosition,{"BmnTof701Hit"})
    .Define("tof701hitT","BmnTof701Hit.fTimeStamp")
    .Define("tof701hitL","BmnTof701Hit.fLength")
    .Define("tof701hitResX","BmnTof701Hit.fResX")
    .Define("tof701hitResY","BmnTof701Hit.fResY")
    .Define("tof701hitRefIndex","BmnTof701Hit.fRefIndex")
    .Define("trTof701hitRes2D", tofResFromHit701, {"BmnGlobalTrack", "BmnTof701Hit"})

  // fields with independent matching results "_m"
    
    .Define("tof701hitPos_m", tofHitPosition, {"Tof701"})
    .Define("tof701hitT_m", "Tof701.fTimeStamp")
    .Define("tof701hitL_m", "Tof701.fLength")
    .Define("tof701hitResX_m", "Tof701.fResX")
    .Define("tof701hitResY_m", "Tof701.fResY")
    .Define("tof701hitModule_m", "Tof701.fModule")
    .Define("trTof701hitRes2D_m", tofResFromHit701, {"Track", "Tof701"})

    .Define("tof400hitPos_m", tofHitPosition, {"Tof400"})
    .Define("tof400hitT_m", "Tof400.fTimeStamp")
    .Define("tof400hitL_m", "Tof400.fLength")
    .Define("tof400hitResX_m", "Tof400.fResX")
    .Define("tof400hitResY_m", "Tof400.fResY")
    .Define("tof400hitRefIndex_m", "Tof400.fRefIndex")
    .Define("tof400hitModule_m", "Tof400.fModule")
    .Define("trTof400hitRes2D_m", tofResFromHit400, {"Track", "Tof400"})
//

    .Define("scwallModPos",[scwallModPos](){return scwallModPos;})
    .Define("scwallModId",moduleId, {"scwallModPos"})
    .Define("scwallModQ", scwallModQ, {"ScWallEvent"})
    
    .Define("hodoModPos",[hodoModPos](){return hodoModPos;})
    .Define("hodoModId",moduleId, {"hodoModPos"})
    .Define("hodoModQ", hodoModQ, {"HodoEvent"})
    
    .Define("fhcalModPos",[fhcalModPos](){return fhcalModPos;})
    .Define("fhcalModId", moduleId, {"fhcalModPos"})
    .Define("fhcalModE",fhcalModE,{"FHCalEvent"})
    .Define("fhcalSecE",fhcalSectionE,{"FHCalEvent"}) 

//    .Define("fdQ","Sum(FDPoint.fCharge*FDPoint.fCharge)")
//    .Define("fdLight","Sum(FDPoint.fLightYield)")
//    .Define("fdEloss", fdEloss, {"FDPoint"})

//    .Define("bdPointEloss", mcPointEloss, {"BdPoint"})
//    .Define("bdPointModId", "BdPoint.nCopy")
//    .Define("bdPointPdg", "BdPoint.fPdgId")
//    .Define("bdPointIsPrimary", "BdPoint.fIsPrimary")
  ;
  // dd.Foreach([](uint evtId){if (evtId % 1000 == 0) cout << "\n" << evtId;}, {"evtId"}); // progress display 
  std::cout << std::endl;

  vector<string> definedNames;
  vector<string> toExclude={/*"scwallModPos","fhcalModPos","hodoModPos"*/};
  for (auto& definedName:dd.GetDefinedColumnNames())
  {
    bool exclude=false;
    for (auto &nameToExclude:toExclude)
      if (definedName==nameToExclude)
        exclude=true;
    if (!exclude)
      definedNames.push_back(definedName);
  }
  dd.Snapshot("t", fileOut, definedNames);

  std::cout<<"Convert_done"<<std::endl;
}
