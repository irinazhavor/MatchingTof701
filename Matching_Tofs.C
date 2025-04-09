/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CheckAlignment.C
 * Author: mikhailr
 *
 * Created on August 9, 2021, 1:24 PM
 */

#include <FairHit.h>
#include <cstdlib>
#include <map>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <stdlib.h>
// #include <BmnTOF1Point.h>

using namespace std;
using namespace TMath;

const Int_t run_period = 8;
const Double_t dX_small_WinMatch = 4.4, dX_big_WinMatch = 8.1, dY_small_WinMatch = 3.25, dY_big_WinMatch = 5.0;
const Double_t C_mPerNs = TMath::C() * TMath::Power(10, -9);

const Double_t TofYcorr[59] = {0};
const Double_t TofXcorr[59] = {0};
const Double_t TofZcorr[59] = {0};

Double_t TofTimecorr[59][48] = {0.};

Double_t TofX = 0, TofY = 0, TofZ = 0;
const Int_t NDet = 59;
const Int_t NStr = 48;
Double_t dX_win, dY_win;

int Matching_tofs(TString FileDstName = "", TString FileOutName = "", Long64_t nEvForRead = 0)
{

    TStopwatch timer;
    timer.Start();

    //--------------------------------------------------------------------------
    // --------- Read dst file -------------------------------------------------
    if (FileDstName == "")
    {
        cout << "File not specified!\n\n"
             << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
             << endl;
        return -1;
    }

    TFile *RootFileInDst = new TFile(FileDstName.Data(), "READ");
    if (!RootFileInDst->IsOpen())
    {
        cout << "File " << FileDstName.Data() << " are not find" << endl;
        return -4;
    }

    printf("Open %s \n", FileDstName.Data());

    TTree *eveTreeDst = (TTree *)RootFileInDst->Get("bmndata");

    // read first Run Header if present
    DstRunHeader *run_header = (DstRunHeader *)eveTreeDst->GetCurrentFile()->Get("DstRunHeader");

    DstEventHeader *evHeader = new DstEventHeader();
    eveTreeDst->SetBranchAddress("DstEventHeader.", &evHeader);

    TClonesArray *globTracks = nullptr;
    eveTreeDst->SetBranchAddress("BmnGlobalTrack", &globTracks);

    TClonesArray *StsTracks = nullptr;
    eveTreeDst->SetBranchAddress("StsVector", &StsTracks);

    CbmVertex *vertices = nullptr;
    eveTreeDst->SetBranchAddress("PrimaryVertex.", &vertices);

    TClonesArray *tof701Hits = nullptr;
    eveTreeDst->SetBranchAddress("BmnTof701Hit", &tof701Hits);

    TClonesArray *tof400Hits = nullptr;
    eveTreeDst->SetBranchAddress("BmnTof400Hit", &tof400Hits);

    if (run_header)
    {
        cout << "\n|||||||||||||||| RUN SUMMARY |||||||||||||||" << endl;
        cout << "||\t\t\t\t\t  ||" << endl;
        cout << "||   Period:        " << run_header->GetPeriodNumber() << "\t\t\t  ||" << endl;
        cout << "||   Number:        " << run_header->GetRunNumber() << "\t\t  ||" << endl;
        cout << "||   Start Time:    " << run_header->GetStartTime().AsString("s") << "\t  ||" << endl;
        cout << "||   End Time:      " << run_header->GetFinishTime().AsString("s") << "\t  ||" << endl;
        cout << "||   Beam:          A = " << run_header->GetBeamA() << ", Z = " << run_header->GetBeamA() << "\t  ||" << endl;
        cout << "||   Beam energy:   " << run_header->GetBeamEnergy() << " GeV\t\t  ||" << endl;
        cout << "||   Target:        A = " << run_header->GetTargetA() << ", Z = " << run_header->GetTargetZ() << "\t  ||" << endl;
        // cout << "||   Field voltage: " << setprecision(4) << run_header->GetMagneticField() << " mV\t\t  ||" << endl;
        cout << "||\t\t\t\t\t  ||" << endl;
        cout << "||||||||||||||||||||||||||||||||||||||||||||\n"
             << endl;
    }

    Double_t VertexPos[4] = {0, 0, 0, 0};
    VertexPos[0] = 0;
    VertexPos[1] = 0;
    VertexPos[2] = 0;
    Double_t dVertexPos[3];
    Double_t tempL = 0;
    BmnKalmanFilter *kalman = new BmnKalmanFilter();
    BmnStatus resultPropagate;
    Int_t run_period = run_header->GetPeriodNumber();
    Int_t run_number = run_header->GetRunNumber();

    //--------------------------------------------------------------------------
    // -----  Get geometry and field scale from db  ----------------------------
    FairRunAna *fRunAna = new FairRunAna();
    Double_t fieldScale = 0.;
    Bool_t isField = kFALSE, isTarget = kFALSE;

    // gRandom->SetSeed(0);
    TString geoFileStr = "full_geometry_run8.root";
    //    TString dir = Form("%s%s%s", getenv("VMCWORKDIR"), "/input/", NameFile.Data());

    TString geoFileName = Form("%s%s%s", "/scratch1/iazhavor/bmnroot_sample/bmnroot", "/macro/run/geometry_run/", geoFileStr.Data());
    // TString geoFileName = Form("current_geo_file_%d.root", UInt_t(gRandom->Integer(UINT32_MAX)));
    // Int_t res_code = UniRun::ReadGeometryFile(run_period, run_number, (char*) geoFileName.Data());
    // if (res_code != 0) {
    //     cout << "ERROR: could not read geometry file from the database" << endl;
    //     exit(-3);
    // }

    TFile *geoFile = new TFile(geoFileName, "READ");
    if (!geoFile->IsOpen())
    {
        cout << "ERROR: could not open ROOT file with geometry: " + geoFileName << endl;
        exit(-4);
    }
    TList *keyList = geoFile->GetListOfKeys();
    TIter next(keyList);
    TKey *key = (TKey *)next();
    TString className(key->GetClassName());
    if (className.BeginsWith("TGeoManager"))
        key->ReadObj();
    else
    {
        cout << "ERROR: TGeoManager is not top element in geometry file " + geoFileName << endl;
        exit(-5);
    }

    //--------------------------------------------------------------------------/*
    // set magnet field with factor corresponding to the given run
    UniRun *pCurrentRun = UniRun::GetRun(run_period, run_number);
    if (pCurrentRun == 0)
        exit(-6);
    Double_t *field_voltage = pCurrentRun->GetFieldVoltage();
    if (field_voltage == NULL)
    {
        cout << "ERROR: no field voltage was found for run " << run_period << ":" << run_number << endl;
        exit(-7);
    }
    Double_t map_current = (run_period == 8) ? 112.0 : (run_period == 7) ? 55.87
                                                                         : 0.0;
    if (*field_voltage < 10)
    {
        fieldScale = 0;
    }
    else
        fieldScale = (*field_voltage) / map_current;

    BmnFieldMap *magField = new BmnNewFieldMap("FieldMap_1900_extrap_noPed.root");
    magField->SetScale(fieldScale);
    magField->Init();
    fRunAna->SetField(magField);
    TString targ;
    if (pCurrentRun->GetTargetParticle() == NULL)
    {
        targ = "-";
        isTarget = kFALSE;
    }
    else
    {
        targ = (pCurrentRun->GetTargetParticle())[0];
        isTarget = kTRUE;
    }
    TString beam = pCurrentRun->GetBeamParticle();

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cout << "\n\n|||||||||||||||| EXPERIMENTAL RUN SUMMARY ||||||||||||||||" << endl;
    cout << "||\t\t\t\t\t\t\t||" << endl;
    cout << "||\t\tPeriod:\t\t" << run_period << "\t\t\t||" << endl;
    cout << "||\t\tNumber:\t\t" << run_number << "\t\t\t||" << endl;
    cout << "||\t\tBeam:\t\t" << beam << "\t\t\t||" << endl;
    cout << "||\t\tTarget:\t\t" << targ << "\t\t\t||" << endl;
    cout << "||\t\tField scale:\t" << setprecision(4) << fieldScale << "\t\t\t||" << endl;
    cout << "||\t\t\t\t\t\t\t||" << endl;
    cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n"
         << endl;

    //--------------------------------------------------------------------------
    //----- output file --------------------------------------------------------

    TFile *RootFileOut = new TFile(FileOutName.Data(), "RECREATE");

    TTree *t_out = new TTree("TofMatch", "TofMatch");

    Int_t oRunID = run_header->GetRunNumber();
    t_out->Branch("RunId", &oRunID);

    Int_t oEvID = -1;
    t_out->Branch("EvId", &oEvID);

    TClonesArray *oTrack = new TClonesArray("BmnGlobalTrack");
    t_out->Branch("Track", &oTrack);

    TClonesArray *oTof400 = new TClonesArray("BmnHit");
    t_out->Branch("Tof400", &oTof400);

    TClonesArray *oTof701 = new TClonesArray("BmnHit");
    t_out->Branch("Tof701", &oTof701);

    CbmVertex *oVertices = nullptr;
    t_out->Branch("Vertex", &oVertices);
    
    //--------------------------------------------------------------------------
    // event loop
    Int_t iEvAna;
    Long64_t nEvents = eveTreeDst->GetEntries();
    if (nEvForRead == 0 || nEvForRead > nEvents)
        nEvForRead = nEvents;
    cout << "Will be read " << nEvForRead << " events from " << nEvents << endl;

    for (Int_t iEv = 0; iEv < nEvForRead; iEv++)
    {

        if (iEv % 1000 == 0)
            cout << "Event# " << iEv << endl;

        oTrack->Delete();
        oTof400->Delete();
        oTof701->Delete();

        globTracks->Delete();
        StsTracks->Delete();
        tof701Hits->Delete();
        tof400Hits->Delete();

        eveTreeDst->GetEntry(iEv);
        oEvID = iEv;

        memset(dVertexPos, 0, sizeof(dVertexPos));
        memset(VertexPos, 0, sizeof(VertexPos));
        VertexPos[0] = vertices->GetX();
        VertexPos[1] = vertices->GetY();
        VertexPos[2] = vertices->GetZ();
        VertexPos[3] = vertices->GetNTracks();

        for (Int_t iTrack = 0; iTrack < StsTracks->GetEntriesFast(); iTrack++)
        {
            CbmStsTrack *gemTr = (CbmStsTrack *)StsTracks->UncheckedAt(iTrack);

            FairTrackParam *parLastGlob = gemTr->GetParamLast();
            FairTrackParam *parFirstGlob = gemTr->GetParamFirst();

            Int_t fPDG, NHits = 0;
            if (1. / parFirstGlob->GetQp() > 0)
                fPDG = 2212;
            else
                fPDG = -211;

            // copy parameters for propagation
            FairTrackParam pLastTemp = *parLastGlob;
            FairTrackParam pLastTemp2 = *parLastGlob;
            FairTrackParam pFirstTemp = *parFirstGlob;
            FairTrackParam pFirstTemp2 = *parFirstGlob;

            //------------------------------------------------------------------

            Int_t mindist2 = 1000, dist = 1100, dist2 = 1100, indexHit4 = -1, indexHit7 = -1;
            Double_t LtoTarget = 0, LpFirstToLast = 0, LtoTof400 = 0, LtoTof700 = 0;
            Double_t beta400 = -1000, beta701 = -1000;

            resultPropagate = kalman->TGeoTrackPropagate(&pFirstTemp2, VertexPos[2], fPDG, NULL, &LtoTarget, isField);
            if (resultPropagate == kBMNERROR)
                continue;

            dVertexPos[0] = VertexPos[0] - pFirstTemp2.GetX();
            dVertexPos[1] = VertexPos[1] - pFirstTemp2.GetY();

            Double_t dca = Sqrt(Power(dVertexPos[0], 2) + Power(dVertexPos[1], 2));

            BmnGlobalTrack ogemTr;
            ogemTr.SetParamLast(pLastTemp);
            ogemTr.SetParamFirst(pFirstTemp2);
            ogemTr.SetGemTrackIndex(iTrack);
            ogemTr.SetNHits(gemTr->GetNStsHits());
            ogemTr.SetDCAInVertex(dca);

            LpFirstToLast = 0;
            pFirstTemp = *parFirstGlob;
            resultPropagate = kalman->TGeoTrackPropagate(&pFirstTemp, parLastGlob->GetZ(), fPDG, NULL, &LpFirstToLast, isField);
            if (resultPropagate == kBMNERROR)
                continue;

            // matching tof400

            pLastTemp = *parLastGlob;
            resultPropagate = kalman->TGeoTrackPropagate(&pLastTemp, 424, fPDG, NULL, &LtoTof400, isField);
            if (resultPropagate == kBMNERROR)
                continue;

            for (Int_t iHit = 0; iHit < tof400Hits->GetEntriesFast(); iHit++)
            {
                BmnTofHit *hitTemp = (BmnTofHit *)tof400Hits->At(iHit); 
                BmnTofHit hit;
                hit = *hitTemp; 
                TofX = hit.GetX();
                TofY = hit.GetY();
                TofZ = hit.GetZ();

                if (Abs(pLastTemp.GetX() - TofX) > 20 || Abs(pLastTemp.GetY() - TofY) > 20)
                    continue; // skip the far hits

                pLastTemp2 = pLastTemp;
                resultPropagate = kalman->TGeoTrackPropagate(&pLastTemp2, TofZ, fPDG, NULL, &LtoTof400, isField);
                if( Abs( pLastTemp2.GetX() - TofX ) > 6 || Abs( pLastTemp2.GetY() - TofY ) > 6 ) continue;
                dist2 = Sqrt(Power(pLastTemp2.GetX() - TofX, 2) + Power(pLastTemp2.GetY() - TofY, 2));

                if (dist2 < mindist2)
                {
                    mindist2 = dist2;
                    indexHit4 = iHit;
                }
            }

            // get the nearest TofHit
            if (indexHit4 != -1)
            {
                BmnTofHit *hit4save = new BmnTofHit();

                BmnTofHit *hit = (BmnTofHit *)tof400Hits->At(indexHit4);
                hit4save->SetX(hit->GetX());
                hit4save->SetY(hit->GetY());
                hit4save->SetZ(hit->GetZ());
                hit4save->SetTimeStamp(hit->GetTimeStamp());

                // Calculate length from ParamLast to Hit
                LtoTof400 = 0;
                pLastTemp = *parLastGlob;
                resultPropagate = kalman->TGeoTrackPropagate(&pLastTemp, hit4save->GetZ(), fPDG, NULL, &LtoTof400, isField);
                if (resultPropagate == kBMNERROR)
                    continue; // correction

                Double_t Ltrack = LtoTarget + LpFirstToLast + LtoTof400;
                Double_t dx = pLastTemp.GetX() - hit4save->GetX();
                Double_t dy = pLastTemp.GetY() - hit4save->GetY();

                hit4save->SetResXY(dx, dy);
                hit4save->SetModule(((hit->GetDetectorID() & 0x0000FF00) >> 8) - 1);
                hit4save->SetStation((hit->GetDetectorID() & 0x000000FF) - 1);
                hit4save->SetLength(Ltrack);
                hit4save->SetDetectorID(hit->GetDetectorID());
                Double_t Vel_mPerNs = (Ltrack / 100.) / (hit4save->GetTimeStamp()); // m/ns
                beta400 = Vel_mPerNs / C_mPerNs;

                new ((*oTof400)[oTof400->GetEntriesFast()]) BmnTofHit(*hit4save);

                ogemTr.SetLength(Ltrack);
                ogemTr.SetTof1HitIndex(oTof400->GetEntriesFast() - 1);
                ogemTr.SetBeta(beta400, 1);
            }
            // tof700 matching
            mindist2 = 1000;
            pLastTemp = *parLastGlob;
            resultPropagate = kalman->TGeoTrackPropagate(&pLastTemp, 584, fPDG, NULL, &LtoTof700, isField); 
            if (resultPropagate == kBMNERROR)
                continue;

            // find the nearest Tof hit

            for (Int_t iHit = 0; iHit < tof701Hits->GetEntriesFast(); iHit++)
            {
                BmnTofHit *hitTemp = (BmnTofHit *)tof701Hits->At(iHit); // correction
                BmnTofHit hit;
                hit = *hitTemp;
                TofX = hit.GetX() + TofXcorr[BmnTOF1Point::GetModule(hit.GetDetectorID()) - 1];
                TofY = hit.GetY() + TofYcorr[BmnTOF1Point::GetModule(hit.GetDetectorID()) - 1];
                TofZ = hit.GetZ() + TofZcorr[BmnTOF1Point::GetModule(hit.GetDetectorID()) - 1];

                if (Abs(pLastTemp.GetX() - TofX) > 40 || Abs(pLastTemp.GetY() - TofY) > 20)
                    continue; // skip the far hits

                pLastTemp2 = pLastTemp;
                resultPropagate = kalman->TGeoTrackPropagate(&pLastTemp2, TofZ, fPDG, NULL, &LtoTof700, isField);
                if(Abs( pLastTemp2.GetX() - TofX ) > 4.2 || Abs( pLastTemp2.GetY() - TofY ) > 3. ) continue;
                dist2 = Sqrt(Power(pLastTemp2.GetX() - TofX, 2) + Power(pLastTemp2.GetY() - TofY, 2));

                if (dist2 < mindist2)
                {
                    mindist2 = dist2;
                    indexHit7 = iHit;
                }

            } 

            // get the nearest TofHit and save
            if (indexHit7 != -1)
            {
                BmnTofHit *hit4save = new BmnTofHit();

                BmnTofHit *hit = (BmnTofHit *)tof701Hits->At(indexHit7);
                hit4save->SetX(hit->GetX());
                hit4save->SetY(hit->GetY());
                hit4save->SetZ(hit->GetZ());
                hit4save->SetTimeStamp(hit->GetTimeStamp());

                // Calculate length from ParamLast to Hit
                LtoTof700 = 0;
                pLastTemp = *parLastGlob;
                resultPropagate =
                    kalman->TGeoTrackPropagate(&pLastTemp, hit4save->GetZ(), fPDG, NULL, &LtoTof700, isField);
                if (resultPropagate == kBMNERROR)
                    continue; // correction

                Double_t Ltrack = LtoTarget + LpFirstToLast + LtoTof700;
                Double_t dx = pLastTemp.GetX() - hit4save->GetX();
                Double_t dy = pLastTemp.GetY() - hit4save->GetY();

                hit4save->SetResXY(dx, dy);
                hit4save->SetModule((((hit->GetDetectorID() & 0x0000FF00) >> 8) - 1));
                hit4save->SetStation(((hit->GetDetectorID() & 0x000000FF) - 1));
                hit4save->SetDetectorID(hit->GetDetectorID());
                hit4save->SetLength(Ltrack);
                Double_t Vel_mPerNs = (Ltrack / 100.) / (hit4save->GetTimeStamp()); // m/ns
                beta701 = Vel_mPerNs / C_mPerNs;

                new ((*oTof701)[oTof701->GetEntriesFast()]) BmnTofHit(*hit4save);

                ogemTr.SetTof701HitIndex(oTof701->GetEntriesFast() - 1);
                ogemTr.SetBeta(beta701, 3);
            }
            new ((*oTrack)[oTrack->GetEntriesFast()]) BmnGlobalTrack(ogemTr);
        }

        oVertices = vertices;
        t_out->Fill();
    }

    RootFileInDst->Close();

    timer.Stop();
    cout << "Readed " << nEvForRead << " events" << endl;
    cout << "Time          = " << timer.RealTime() << " s" << endl;
    cout << "Time/Event    = " << timer.RealTime() / (Double_t)nEvForRead * 1000. << " ms/Event" << endl;
    cout << "TimeCPU       = " << timer.CpuTime() << " s" << endl;
    cout << "TimeCPU/Event = " << timer.CpuTime() / (Double_t)nEvForRead * 1000. << " ms/Event" << endl;

    timer.Reset();
    timer.Start();

    RootFileOut->cd();
    t_out->Write();
    RootFileOut->Close();

    timer.Stop();
    cout << "Time for write root = " << timer.RealTime() << " s" << endl;

    //*/
    return 0;
}
