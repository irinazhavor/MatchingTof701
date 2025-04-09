#include <cmath>
#include <vector>
#include <functional>
using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;

auto CalcTrackCoord = [](const int coordIdx) { // o - x coordinate, 1 - y coordinate
    return [coordIdx](std::vector<double> hitsRes, std::vector<double> hitsCoordinates)
    {
        std::vector<double> trCoordinate;
        trCoordinate.reserve(hitsCoordinates.size());

        // std::cout << "\ncoordIdx: " << coordIdx << std::endl;
        // std::cout << "hitsCoordinates.size(): " << hitsCoordinates.size() << endl;
        // std::cout << "hitsRes.size(): " << hitsRes.size() << endl;

        for (int i = 0; i < hitsCoordinates.size(); i++)
        {
            const auto hit_x_y = hitsCoordinates.at(i);
            if (fabs(hitsCoordinates.at(i)) > 300 || fabs(hitsRes.at(i)) > 100)
            {
                trCoordinate.push_back(-1000);
                continue;
            }
            trCoordinate.push_back(hit_x_y + hitsRes.at(i));
        }
        // std::cout << "trCoordinate.size(): " << trCoordinate.size() << " (coordIdx: " << coordIdx << ");" << std::endl;
        // getchar();

        return trCoordinate;
    };
};

auto BetaPerMod_generator = [](const auto mod_id)
{
    return [mod_id](ROOT::VecOps::RVec<double> tr_beta_vec, ROOT::VecOps::RVec<int> hit_idx_vec, ROOT::VecOps::RVec<int> mod_vec)
    {
        std::vector<double> beta_mod;
        beta_mod.reserve(tr_beta_vec.size());
        for (int i = 0; i < tr_beta_vec.size(); i++)
        {
            int hit_idx = hit_idx_vec.at(i);

            if (hit_idx < 0 || !hit_idx_vec.at(hit_idx))
            {
                beta_mod.push_back(-999.);
                continue;
            }

            beta_mod.push_back(mod_vec.at(hit_idx) == mod_id ? tr_beta_vec : -999);
        }
        return beta_mod;
    };
};

vector<vector<double>> tofHitPosition(const TClonesArray hits)
try
{
    vector<vector<double>> pos;
    for (const auto &hitObj : hits)
    {
        auto hit = (BmnTofHit *)hitObj;
        pos.push_back({hit->GetX(), hit->GetY(), hit->GetZ()});
    }
    return pos;
}
catch (const std::exception &e)
{
    std::cout << __func__ << std::endl;
    throw e;
}

TChain *makeChain(string &filename, const char *treename)
{
    cout << "Adding files to chain:" << endl;
    TChain *chain = new TChain(treename);
    if (filename.rfind(".root") < filename.size())
        chain->Add(filename.data());
    else
    {
        TFileCollection fc("fc", "", filename.c_str());
        chain->AddFileInfoList((TCollection *)fc.GetList());
    }
    chain->ls();
    return chain;
}

vector<int> ModHitTof701(RVec<BmnGlobalTrack> global_tracks, const TClonesArray hits)
try
{
    vector<int> parameters;
    for (auto &global_track : global_tracks)
    {
        auto hitIndex = global_track.GetTof701HitIndex();
        if (hitIndex < 0 || !hits.At(hitIndex))
        {
            parameters.push_back(-1000);
            continue;
        }
        BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
        parameters.push_back(hit->GetModule());
    }
    return parameters;
}
catch (const std::exception &e)
{
    std::cout << __func__ << std::endl;
    throw e;
}

vector<int> ModHitTof400(RVec<BmnGlobalTrack> global_tracks, const TClonesArray hits)
try
{
    vector<int> parameters;
    for (auto &global_track : global_tracks)
    {
        auto hitIndex = global_track.GetTof1HitIndex();
        if (hitIndex < 0 || !hits.At(hitIndex))
        {
            parameters.push_back(-1000);
            continue;
        }
        BmnTofHit *hit = (BmnTofHit *)hits.At(hitIndex);
        Int_t mod = ((hit->GetDetectorID() & 0x0000FF00) >> 8) - 1;
        parameters.push_back(mod);
    }
    return parameters;
}
catch (const std::exception &e)
{
    std::cout << __func__ << std::endl;
    throw e;
}

vector<int> BetaCorrelationHitMod(vector<int> mod_globtracks, vector<int> mod_matching)
try
{
    vector<int> isMod;
    for (int i = 0; i < mod_globtracks.size(); i++)
    {
        if (mod_globtracks.at(i) == 1 && mod_matching.at(i) == 1)
            isMod.push_back(1);
        else
            isMod.push_back(0);
    }
    return isMod;
}
catch (const std::exception &e)
{
    std::cout << e.what() << std::endl;
    std::cout << __func__ << std::endl;
    throw e;
}

auto tofTrackPosition = [](const int coordIdx) { // o - x coordinate, 1 - y coordinate
    return [coordIdx](ROOT::VecOps::RVec<double> hitsRes, vector<vector<double>> hitsCoordinates)
    {
        std::vector<double> trCoordinate;
        trCoordinate.reserve(hitsCoordinates.size());

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
        return trCoordinate;
    };
};

void qa_converter_tofs(std::string str_in_converted = "",
                       std::string out_file_name = "qa.root")
{

    TFileCollection collection("collection", "", str_in_converted.c_str());
    auto *chain = new TChain("t");
    auto *chain1 = new TChain("t");
    chain->AddFileInfoList(collection.GetList());
    chain1->AddFileInfoList(collection.GetList());
    ROOT::RDataFrame d1(*chain);
    ROOT::RDataFrame d2(*chain1);

    // ROOT::RDataFrame d("t", str_in_converted.c_str());

    std::cout << "Preparing the RDF" << endl;

    // RDF for the independent matching - "_m"
    auto ddd = d2
                   .Define("vtxR_m", "return sqrt(vtxX_m*vtxX_m + vtxY_m*vtxY_m);")
                   .Filter("fabs(vtxR_m) < 1")
                   .Filter("vtxZ_m < 0.2")
                   .Filter("vtxNtracks_m > 1")

                   .Define("tof400HitX_", "std::vector<double> Xhit400mcoord; for ( int i = 0; i < tof400hitPos_m.size(); i++ ) { Xhit400mcoord.push_back( tof400hitPos_m.at(i).X() ); } return Xhit400mcoord;")
                   .Define("tof400HitY_", "std::vector<double> Yhit400mcoord; for ( int i = 0; i < tof400hitPos_m.size(); i++ ) { Yhit400mcoord.push_back( tof400hitPos_m.at(i).Y() ); } return Yhit400mcoord;")
                   .Define("tof701HitX_", "std::vector<double> Xhit701mcoord; for ( int i = 0; i < tof701hitPos_m.size(); i++ ) { Xhit701mcoord.push_back( tof701hitPos_m.at(i).X() ); } return Xhit701mcoord;")
                   .Define("tof701HitY_", "std::vector<double> Yhit701mcoord; for ( int i = 0; i < tof701hitPos_m.size(); i++ ) { Yhit701mcoord.push_back( tof701hitPos_m.at(i).Y() ); } return Yhit701mcoord;")

                   .Define("trTof400ResX_", "std::vector<double> res; for ( int i = 0; i < trTof400hitRes2D_m.size(); i++ ){ res.push_back( trTof400hitRes2D_m.at(i).X() ); } return res;")
                   .Define("trTof400ResY_", "std::vector<double> res; for ( int i = 0; i < trTof400hitRes2D_m.size(); i++ ){ res.push_back( trTof400hitRes2D_m.at(i).Y() ); } return res;")
                   .Define("trTof701ResX_", "std::vector<double> res; for ( int i = 0; i < trTof701hitRes2D_m.size(); i++ ){ res.push_back( trTof701hitRes2D_m.at(i).X() ); } return res;")
                   .Define("trTof701ResY_", "std::vector<double> res; for ( int i = 0; i < trTof701hitRes2D_m.size(); i++ ){ res.push_back( trTof701hitRes2D_m.at(i).Y() ); } return res;")

                   .Define("trTof400HitX_", "std::vector<double> Xtrmcoord; for ( int i = 0; i < trBeta400Hit_coors_m.size(); i++ ) { Xtrmcoord.push_back( trBeta400Hit_coors_m.at(i).X() ); } return Xtrmcoord;")
                   .Define("trTof400HitY_", "std::vector<double> Ytrmcoord; for ( int i = 0; i < trBeta400Hit_coors_m.size(); i++ ) { Ytrmcoord.push_back( trBeta400Hit_coors_m.at(i).Y() ); } return Ytrmcoord;")
                   .Define("trTof701HitX_", "std::vector<double> Xtrmcoord; for ( int i = 0; i < trBeta701Hit_coors_m.size(); i++ ) { Xtrmcoord.push_back( trBeta701Hit_coors_m.at(i).X() ); } return Xtrmcoord;")
                   .Define("trTof701HitY_", "std::vector<double> Ytrmcoord; for ( int i = 0; i < trBeta701Hit_coors_m.size(); i++ ) { Ytrmcoord.push_back( trBeta701Hit_coors_m.at(i).Y() ); } return Ytrmcoord;")

                   .Define("trXfromResTof400_", CalcTrackCoord(0), {"trTof400ResX_", "trTof400HitX_"})
                   .Define("trYfromResTof400_", CalcTrackCoord(1), {"trTof400ResY_", "trTof400HitY_"})
                   .Define("trXfromResTof701_", CalcTrackCoord(0), {"trTof701ResX_", "trTof701HitX_"})
                   .Define("trYfromResTof701_", CalcTrackCoord(1), {"trTof701ResY_", "trTof701HitY_"});

    // RDF for the Global Tracking "
    auto dd = d1
                  .Define("vtxR", "return sqrt(vtxX*vtxX + vtxY*vtxY);")
                  .Filter("fabs(vtxR) < 1")
                  .Filter("vtxZ < 0.2")
                  .Filter("vtxNtracks > 1")

                  .Define("tof400HitX", "std::vector<double> Xhit400coord; for ( int i = 0; i < tof400hitPos.size(); i++ ) { Xhit400coord.push_back( tof400hitPos.at(i).X() ); } return Xhit400coord;")
                  .Define("tof400HitY", "std::vector<double> Yhit400coord; for ( int i = 0; i < tof400hitPos.size(); i++ ) { Yhit400coord.push_back( tof400hitPos.at(i).Y() ); } return Yhit400coord;")
                  .Define("tof701HitX", "std::vector<double> Xhit701coord; for ( int i = 0; i < tof701hitPos.size(); i++ ) { Xhit701coord.push_back( tof701hitPos.at(i).X() ); } return Xhit701coord;")
                  .Define("tof701HitY", "std::vector<double> Yhit701coord; for ( int i = 0; i < tof701hitPos.size(); i++ ) { Yhit701coord.push_back( tof701hitPos.at(i).Y() ); } return Yhit701coord;")

                  .Define("trTof400ResX", "std::vector<double> res; for ( int i = 0; i < trTof400hitRes2D.size(); i++ ){ res.push_back( trTof400hitRes2D.at(i).X() ); } return res;")
                  .Define("trTof400ResY", "std::vector<double> res; for ( int i = 0; i < trTof400hitRes2D.size(); i++ ){ res.push_back( trTof400hitRes2D.at(i).Y() ); } return res;")
                  .Define("trTof701ResX", "std::vector<double> res; for ( int i = 0; i < trTof701hitRes2D.size(); i++ ){ res.push_back( trTof701hitRes2D.at(i).X() ); } return res;")
                  .Define("trTof701ResY", "std::vector<double> res; for ( int i = 0; i < trTof701hitRes2D.size(); i++ ){ res.push_back( trTof701hitRes2D.at(i).Y() ); } return res;")

                  .Define("trTof400HitX", "std::vector<double> Xtrcoord; for ( int i = 0; i < trBeta400Hit_coors.size(); i++ ) { Xtrcoord.push_back( trBeta400Hit_coors.at(i).X() ); } return Xtrcoord;")
                  .Define("trTof400HitY", "std::vector<double> Ytrcoord; for ( int i = 0; i < trBeta400Hit_coors.size(); i++ ) { Ytrcoord.push_back( trBeta400Hit_coors.at(i).Y() ); } return Ytrcoord;")
                  .Define("trTof701HitX", "std::vector<double> Xtrcoord; for ( int i = 0; i < trBeta701Hit_coors.size(); i++ ) { Xtrcoord.push_back( trBeta701Hit_coors.at(i).X() ); } return Xtrcoord;")
                  .Define("trTof701HitY", "std::vector<double> Ytrcoord; for ( int i = 0; i < trBeta701Hit_coors.size(); i++ ) { Ytrcoord.push_back( trBeta701Hit_coors.at(i).Y() ); } return Ytrcoord;")

                  .Define("trXfromResTof400", CalcTrackCoord(0), {"trTof400ResX", "trTof400HitX"})
                  .Define("trYfromResTof400", CalcTrackCoord(1), {"trTof400ResY", "trTof400HitY"})
                  .Define("trXfromResTof701", CalcTrackCoord(0), {"trTof701ResX", "trTof701HitX"})
                  .Define("trYfromResTof701", CalcTrackCoord(1), {"trTof701ResY", "trTof701HitY"})
                  ;

    TString FieldName, Instruction;

    for (int iPl = 0; iPl < 59; iPl++)
    {
        FieldName = Form("trTof701mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(trBeta701_mod_m.size(), 0); for( int i = 0; i < trBeta701_mod_m.size(); i++ ){ if (trBeta701_mod_m.at(i) == %d) isWork.at(i) = 1; } return isWork;", iPl);
        ddd = ddd.Define(FieldName.Data(), Instruction.Data());

        FieldName = Form("Tof701mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(tof701hitModule_m.size(), 0); for( int i = 0; i < tof701hitModule_m.size(); i++ ){ if (tof701hitModule_m.at(i) == %d) isWork.at(i) = 1; }  return isWork;", iPl);
        ddd = ddd.Define(FieldName.Data(), Instruction.Data());

        FieldName = Form("gltrTof701mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(trBeta701_mod.size(), 0); for( int i = 0; i < trBeta701_mod.size(); i++ ){ if (trBeta701_mod.at(i) == %d) isWork.at(i) = 1; } return isWork;", iPl);
        dd = dd.Define(FieldName.Data(), Instruction.Data());
    }

    for (int iPl = 0; iPl < 20; iPl++)
    {
        FieldName = Form("trTof400mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(trBeta400_mod_m.size(), 0); for( int i = 0; i < trBeta400_mod_m.size(); i++ ){ if (trBeta400_mod_m.at(i) == %d) isWork.at(i) = 1; }  return isWork;", iPl);
        ddd = ddd.Define(FieldName.Data(), Instruction.Data());

        FieldName = Form("Tof400mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(tof400hitModule_m.size(), 0); for( int i = 0; i < tof400hitModule_m.size(); i++ ){ if (tof400hitModule_m.at(i) == %d) isWork.at(i) = 1; }  return isWork;", iPl);
        ddd = ddd.Define(FieldName.Data(), Instruction.Data());

        FieldName = Form("gltrTof400mod_%d", iPl);
        Instruction = Form("std::vector<int> isWork(trBeta400_mod.size(), 0); for( int i = 0; i < trBeta400_mod.size(); i++ ){ if (trBeta400_mod.at(i) == %d) isWork.at(i) = 1; }  return isWork;", iPl);
        dd = dd.Define(FieldName.Data(), Instruction.Data());
    }

    std::vector<ROOT::RDF::RResultPtr<::TH2D>> histo2d701_;
    std::vector<ROOT::RDF::RResultPtr<::TH2D>> histo2d400_;

    std::vector<ROOT::RDF::RResultPtr<::TH2D>> histo2d701;
    std::vector<ROOT::RDF::RResultPtr<::TH2D>> histo2d400;

    histo2d400_.push_back(ddd.Histo2D({"h2_hitY_hitX_tof400_", "hits (x, y) from tracks matched; x_{hit}, cm; y_{hit}, cm", 230, -200, 260, 200, -100, 100}, "tof400HitX_", "tof400HitY_"));
    histo2d400_.push_back(ddd.Histo2D({"h2_trY_trX_tof400_", "hits (x, y) from tracks + (ResX, ResY) ; x_{tr}, cm; y_{tr}, cm", 230, -200, 260, 200, -100, 100}, "trXfromResTof400_", "trYfromResTof400_"));
    histo2d400_.push_back(ddd.Histo2D({"h2_beta_pq_tof400_", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP_m", "trBeta400_m"));
    histo2d400_.push_back(ddd.Histo2D({"h2_m2_pq_tof400_", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP_m", "trM2Tof400_m"));
    histo2d400_.push_back(ddd.Histo2D({"h2_m2_ResX_tof400_", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof400ResX_", "trM2Tof400_m"));
    histo2d400_.push_back(ddd.Histo2D({"h2_m2_ResY_tof400_", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof400ResY_", "trM2Tof400_m"));

    histo2d400.push_back(dd.Histo2D({"h2_hitY_hitX_tof400", "hits (x, y) from tracks matched; x_{hit}, cm; y_{hit}, cm", 230, -200, 260, 200, -100, 100}, "tof400HitX", "tof400HitY"));
    histo2d400.push_back(dd.Histo2D({"h2_trY_trX_tof400", "hits (x, y) from tracks + (ResX, ResY) ; x_{tr}, cm; y_{tr}, cm", 230, -200, 260, 200, -100, 100}, "trXfromResTof400", "trYfromResTof400"));
    histo2d400.push_back(dd.Histo2D({"h2_beta_pq_tof400", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta400"));
    histo2d400.push_back(dd.Histo2D({"h2_m2_pq_tof400", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof400"));
    histo2d400.push_back(dd.Histo2D({"h2_m2_ResX_tof400_", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof400ResX", "trM2Tof400"));
    histo2d400.push_back(dd.Histo2D({"h2_m2_ResY_tof400_", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof400ResY", "trM2Tof400"));

    histo2d701_.push_back(ddd.Histo2D({"h2_hitY_hitX_tof701_", "hits (x, y) from tracks matched; x_{hit}, cm; y_{hit}, cm", 370, -150, 220, 100, -100, 100}, "tof701HitX_", "tof701HitY_"));
    histo2d701_.push_back(ddd.Histo2D({"h2_trY_trX_tof701_", "hits (x, y) from tracks + (ResX, ResY); x_{tr}, cm; y_{tr}, cm", 370, -150, 220, 100, -100, 100}, "trXfromResTof701_", "trYfromResTof701_"));
    histo2d701_.push_back(ddd.Histo2D({"h2_beta_pq_tof701_", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP_m", "trBeta701_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_beta_pq_tof701_wo30_", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP_m", "trBeta701_2_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_pq_tof701_", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP_m", "trM2Tof701_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_pq_tof701_wo30_", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP_m", "trM2Tof701_2_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_ResX_tof701_", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResX_", "trM2Tof701_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_ResY_tof701_", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResY_", "trM2Tof701_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_ResX_tof701_wo30_", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResX_", "trM2Tof701_2_m"));
    histo2d701_.push_back(ddd.Histo2D({"h2_m2_ResY_tof701_wo30_", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResY_", "trM2Tof701_2_m"));

    histo2d701.push_back(dd.Histo2D({"h2_hitY_hitX_tof701", "hits (x, y) from tracks matched; x_{hit}, cm; y_{hit}, cm", 370, -150, 220, 100, -100, 100}, "tof701HitX", "tof701HitY"));
    histo2d701.push_back(dd.Histo2D({"h2_trY_trX_tof701", "hits (x, y) from tracks + (ResX, ResY); x_{tr}, cm; y_{tr}; cm", 370, -150, 220, 100, -100, 100}, "trXfromResTof701", "trYfromResTof701"));
    histo2d701.push_back(dd.Histo2D({"h2_beta_pq_tof701", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta701"));
    histo2d701.push_back(dd.Histo2D({"h2_beta_pq_tof701_wo30", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta701_2"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_pq_tof701", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof701"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_pq_tof701_wo30", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof701_2"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_ResX_tof701", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResX", "trM2Tof701"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_ResY_tof701", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResY", "trM2Tof701"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_ResX_tof701_wo30", "; x_{tr} - x_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResX", "trM2Tof701_2"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_ResY_tof701_wo30", "; y_{tr} - y_{hit}, cm; m^{2} (GeV^{2}/c^{4})", 80, -20, 20, 320, -2., 6.}, "trTof701ResY", "trM2Tof701_2"));

    histo2d701.push_back(dd.Histo2D({"h2_beta_pq_tof700", "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta700"));
    histo2d701.push_back(dd.Histo2D({"h2_m2_pq_tof700", "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof700"));

    TString HistName, ModWeightName, HistTitle;

    for (int iPl = 0; iPl < 59; iPl++)
    {
        ModWeightName = Form("Tof701mod_%d", iPl);

        HistName = Form("h2_hitY_hitX_tof701_%d", iPl);
        histo2d701_.push_back(ddd.Histo2D({HistName.Data(), "; x_{hit}, cm; y_{hit}, cm", 370, -150, 220, 100, -100, 100}, "tof701HitX_", "tof701HitY_", ModWeightName.Data()));
        //
        ModWeightName = Form("trTof701mod_%d", iPl);

        HistName = Form("h2_trY_trX_tof701_%d", iPl);
        histo2d701_.push_back(ddd.Histo2D({HistName.Data(), "; x_{tr}, cm; y_{tr}, cm", 370, -150, 220, 100, -100, 100}, "trXfromResTof701_", "trYfromResTof701_", ModWeightName.Data()));

        HistName = Form("h2_beta_pq_tof701_%d", iPl);
        histo2d701_.push_back(ddd.Histo2D({HistName.Data(), "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP_m", "trBeta701_m", ModWeightName.Data()));

        HistName = Form("h2_m2_pq_tof701_%d", iPl);
        histo2d701_.push_back(ddd.Histo2D({HistName.Data(), "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP_m", "trM2Tof701_m", ModWeightName.Data()));
        //
        ModWeightName = Form("gltrTof701mod_%d", iPl);

        HistName = Form("h2gltr_beta_pq_tof701_%d", iPl);
        histo2d701.push_back(dd.Histo2D({HistName.Data(), "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta701", ModWeightName.Data()));

        HistName = Form("h2_m2_pq_tof701_%d", iPl);
        histo2d701.push_back(dd.Histo2D({HistName.Data(), "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof701", ModWeightName.Data()));

        HistName = Form("h2gltr_trY_trX_tof701_%d", iPl);
        histo2d701.push_back(dd.Histo2D({HistName.Data(), "; x_{tr}, cm; y_{tr}, cm", 370, -150, 220, 100, -100, 100}, "trXfromResTof701", "trYfromResTof701", ModWeightName.Data()));
    }

    for (int iPl = 0; iPl < 20; iPl++)
    {
        ModWeightName = Form("Tof400mod_%d", iPl);

        HistName = Form("h2_hitY_hitX_tof400_%d", iPl);
        histo2d400_.push_back(ddd.Histo2D({HistName.Data(), "; x_{hit}, cm; y_{hit}, cm", 230, -200, 260, 200, -100, 100}, "tof400HitX_", "tof400HitY_", ModWeightName.Data()));
        //
        ModWeightName = Form("trTof400mod_%d", iPl);

        HistName = Form("h2_trY_trX_tof400_%d", iPl);
        histo2d400_.push_back(ddd.Histo2D({HistName.Data(), "; x_{tr}, cm; y_{tr}, cm", 230, -200, 260, 200, -100, 100}, "trXfromResTof400_", "trYfromResTof400_", ModWeightName.Data()));

        HistName = Form("h2_beta_pq_tof400_%d", iPl);
        histo2d400_.push_back(ddd.Histo2D({HistName.Data(), "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP_m", "trBeta400_m", ModWeightName.Data()));

        HistName = Form("h2_m2_pq_tof400_%d", iPl);
        histo2d400_.push_back(ddd.Histo2D({HistName.Data(), "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP_m", "trM2Tof400_m", ModWeightName.Data()));
        //
        ModWeightName = Form("gltrTof400mod_%d", iPl);

        HistName = Form("h2gltr_beta_pq_tof400_%d", iPl);
        histo2d400.push_back(dd.Histo2D({HistName.Data(), "; p/q (GeV/c); #beta", 400, -0.5, 9.5, 400, 0.2, 1.2}, "trP", "trBeta400", ModWeightName.Data()));

        HistName = Form("h2gltr_m2_pq_tof400_%d", iPl);
        histo2d400.push_back(dd.Histo2D({HistName.Data(), "; p/q (GeV/c); m^{2} (GeV^{2}/c^{4})", 400, -0.5, 9.5, 320, -2., 6.}, "trP", "trM2Tof400", ModWeightName.Data()));

        HistName = Form("h2gltr_trY_trX_tof400_%d", iPl);
        histo2d400.push_back(dd.Histo2D({HistName.Data(), "; x_{hit}, cm; y_{hit}, cm", 230, -200, 260, 200, -100, 100}, "trXfromResTof400", "trYfromResTof400", ModWeightName.Data()));
    }

    std::cout << "Writing" << endl;

    auto file_out = TFile::Open(out_file_name.c_str(), "RECREATE");
    file_out->cd();
    file_out->mkdir("tof701");
    file_out->cd("tof701");

    std::for_each(histo2d701_.begin(), histo2d701_.end(), [](auto h1)
                  { h1->Write(); });

    file_out->mkdir("tof400");
    file_out->cd("tof400");

    std::for_each(histo2d400_.begin(), histo2d400_.end(), [](auto h1)
                  { h1->Write(); });

    file_out->mkdir("tof701_gltr");
    file_out->cd("tof701_gltr");

    std::for_each(histo2d701.begin(), histo2d701.end(), [](auto h1)
                  { h1->Write(); });

    file_out->mkdir("tof400_gltr");
    file_out->cd("tof400_gltr");

    std::for_each(histo2d400.begin(), histo2d400.end(), [](auto h1)
                  { h1->Write(); });

    file_out->Close();
}