/********************************************************************* *
*
 * Configfemtoanalysis.C - configuration macro for the femtoscopic       *
 * analysis, meant as a QA process for two-particle effects                              *
 *
*
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                                                                     *
 *
*
 *********************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODMultSelection.h"
 #include "AliFemtoEventReaderNanoAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoMJTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoCorrFctnDEtaDPhiCorrections.h"
#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections.h"
#include "AliFemtoCorrFctnDYDPhiSimple.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoModelCorrFctn3DSpherical.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoPairCutPt.h"
#include "AliFemtoCorrFctnDPhiStarDEta.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoCorrFctnDEtaDPhiStar.h"
#include "AliFemtoKKTrackCutFull.h"
#include "AliFemtoCorrFctnMezonPhi.h"

#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#endif

#include <stdio.h>
#include <string.h>

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

        double PionMass = 0.13956995;
        double KaonMass = 0.493677;
        double ProtonMass = 0.938272013;
        double LambdaMass = 1.115683;
        double XiMass = 1.32171;


        const int numOfMultBins = 5;
        const int numOfChTypes = 3;
        const int numOfkTbins = 5;

        bool performSharedDaughterCut = false;
        bool enablePairMonitors = true;

        char *par = new char[strlen(params)+1];
        strcpy(par,params);
        char *parameter[21];
        if(strlen(params)!=0)
          {
                        parameter[0] = strtok(par, ","); // Splits spaces between words in params
            cout<<"Parameter [0] (filterbit):"<<parameter[0]<<endl; // Writes first parameter
            parameter[1] = strtok(NULL, ",");
            cout<<"Parameter [1] (ktdep):"<<parameter[1]<<" "<<endl;
            parameter[2] = strtok(NULL, ",");
            cout<<"Parameter [2] (multdep):"<<parameter[2]<<" "<<endl;
            parameter[3] = strtok(NULL, ",");
            cout<<"Parameter [3]: (MinPlpContribSPD)"<<parameter[3]<<" "<<endl;
            parameter[4] = strtok(NULL, ",");
            cout<<"Parameter [4]: (multbino)"<<parameter[4]<<" "<<endl;
            parameter[5] = strtok(NULL, ",");
            cout<<"Parameter [5]: (zvertbino)"<<parameter[5]<<" "<<endl;
            parameter[6] = strtok(NULL, ",");
            cout<<"Parameter [6]: (ifGlobalTracks=true/false)"<<parameter[6]<<" "<<endl;
            parameter[7] = strtok(NULL, ",");
            cout<<"Parameter [7]: (shareQuality)"<<parameter[7]<<" "<<endl;
            parameter[8] = strtok(NULL, ",");
            cout<<"Parameter [8]: (shareFraction)"<<parameter[8]<<" "<<endl;
            parameter[9] = strtok(NULL, ",");
            cout<<"Parameter [9]: (ifElectronRejection)"<<parameter[9]<<" "<<endl;
            parameter[10] = strtok(NULL, ",");
            cout<<"Parameter [10]: (nSigma)"<<parameter[10]<<" "<<endl;
            parameter[11] = strtok(NULL, ",");
            cout<<"Parameter [12]: (etaMin)"<<parameter[11]<<" "<<endl;
            parameter[12] = strtok(NULL, ",");
            cout<<"Parameter [12]: (etaMax)"<<parameter[12]<<" "<<endl;
            parameter[13] = strtok(NULL, ",");
            cout<<"Parameter [13]: (ispileup)"<<parameter[13]<<" "<<endl;
            parameter[14] = strtok(NULL, ",");
            cout<<"Parameter [14]: (max pT)"<<parameter[14]<<" "<<endl;
            parameter[15] = strtok(NULL, ",");
            cout<<"Parameter [15]: (SetMostProbable 1)"<<parameter[15]<<" "<<endl;
            parameter[16] = strtok(NULL, ",");
            cout<<"Parameter [16]: (SetMostProbable 2)"<<parameter[16]<<" "<<endl;
            parameter[17] = strtok(NULL, ",");
            cout<<"Parameter [17]: (SetMostProbable 3)"<<parameter[17]<<" "<<endl;
            parameter[18] = strtok(NULL, ",");
            cout<<"Parameter [18]: (FILE no)"<<parameter[18]<<" "<<endl;
            parameter[19] = strtok(NULL, ",");
            cout<<"Parameter [19]: (monitors)"<<parameter[19]<<" "<<endl;
            parameter[20] = strtok(NULL, ",");
            cout<<"Parameter [20]: (nSigma2)"<<parameter[20]<<" "<<endl;
          }
        int filterbit = atoi(parameter[0]); //96 / 768 / 128
        int runktdep = atoi(parameter[1]); //0
        int runmultdep = atoi(parameter[2]); //0
        int minPlpContribSPD = atoi(parameter[3]); //3
        int multbino = atoi(parameter[4]); //30
        int zvertbino = atoi(parameter[5]); //10
        Bool_t ifGlobalTracks=kFALSE; if(atoi(parameter[6]))ifGlobalTracks=kTRUE;//kTRUE
        double shareQuality = atof(parameter[7]); //0.00
        double shareFraction = atof(parameter[8]); //0.05
        bool ifElectronRejection = atoi(parameter[9]); //true
        double nSigmaVal = atof(parameter[10]); //3.0
        double nEtaMin = atof(parameter[11]); //-0.8
        double nEtaMax = atof(parameter[12]);  //0.8
        bool ifIsPileUp = 0;//atoi(parameter[13]); //true
        double maxPt = atof(parameter[14]);  //4.0

        int setMostProb1 = atoi(parameter[15]);
        int setMostProb2 = atoi(parameter[16]);
        int setMostProb3 = atoi(parameter[17]);

        Bool_t ifMonitors=kFALSE; if(atoi(parameter[19]))ifMonitors=kTRUE;//kTRUE
        double nSigmaVal2 = atof(parameter[20]); //2.0 or 3.0
        /*
        printf("*** Connect to AliEn ***\n");
        TGrid::Connect("alien://");
        */
        int runmults[numOfMultBins] = {1, 0, 0, 0, 0};
        if(runmultdep)    {runmults[0]=1; runmults[1]=1; runmults[2]=1;   }
        int multbins[numOfMultBins+1] = {0, 5000, 1000, 0, 0, 0};

        int runch[numOfChTypes] = {/* kaons */ 0, 0, 1};
        const char *chrgs[numOfChTypes] = { "KpKp", "KmKm", "KpKm"};


        double ktrng[numOfkTbins+1] = {0.0, 0, 0, 0, 0, 0};
        double ktrngAll[numOfkTbins+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 100.0};
        double ktrngPion[numOfkTbins+1] = {0.0, 0.8, 1.2, 1.4, 2.5, 100.0};
        double ktrngKaon[numOfkTbins+1] = {0.0, 2.5, 3.5, 100.0, 0, 0};
        double ktrngProton[numOfkTbins+1] = {0.0, 2.75, 100, 0, 0, 0};

        int runqinv = 1;
        int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

        int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner       //global tracks ->mfit ITS+TPC
        int owncuts = 0;
        int owndca = 0;

        int gammacut = 1;       // cut na ee z gamma

        double shqmax = 1.0;
        int nbinssh = 100;

        //AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
        //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


        AliFemtoEventReaderAODMultSelection *Reader = new AliFemtoEventReaderAODMultSelection();
        Reader->SetFilterMask(filterbit);
        Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
        Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
        //Reader->SetMinPlpContribSPD(minPlpContribSPD);
        //Reader->SetIsPileUpEvent(ifIsPileUp);

        Reader->SetReadCascade(kFALSE);
        Reader->SetUseOutOfBunchPlpSelection(kFALSE);
        Reader->SetUseMVPlpSelection(ifIsPileUp);
        Reader->SetTrackPileUpRemoval(ifIsPileUp);

        /*
        AliFemtoEventReaderNanoAODChain *Reader = new AliFemtoEventReaderNanoAODChain();
        Reader->SetFilterMask(filterbit);
        Reader->SetCovMatPresent(false);
        Reader->SetDCAglobalTrack(1); //false for FB7, true for the rest //we do not use DCA at all
        //Reader->SetUseMultiplicity("MultSelection.RefMult08.Value");
        Reader->SetUseMultiplicity("V0M");
  Reader->SetReadV0(kTRUE);
        */

        AliFemtoManager* Manager = new AliFemtoManager();
        Manager->SetEventReader(Reader);

        AliFemtoVertexMultAnalysis      *anetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoBasicEventCut           *mecetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorEventMult     *cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorEventMult     *cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorEventVertex   *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorEventVertex   *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorCollections   *cutPassColletaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorCollections   *cutFailColletaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoKKTrackCutFull          *dtc1etaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoKKTrackCutFull          *dtc2etaphitpc[numOfMultBins*numOfChTypes];
        //AliFemtoMJTrackCut            *dtc1etaphitpc[numOfMultBins*numOfChTypes];
        //AliFemtoMJTrackCut            *dtc2etaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoMJTrackCut              *dtc3etaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0TrackCut              *dtc4etaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0TrackCut              *dtc5etaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoXiTrackCut           *tXiCut[numOfMultBins*numOfChTypes];
        AliFemtoXiTrackCut           *tAXiCut[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutPass1YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutFail1YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutPass1PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutFail1PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutPass2YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutFail2YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutPass2PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutFail2PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutPass3YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticleYPt   *cutFail3YPtetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutPass3PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorParticlePID   *cutFail3PIDetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorV0            *cutPass1V0[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorV0            *cutFail1V0[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorV0            *cutPass2V0[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorV0            *cutFail2V0[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorXi             *cutPass1Xi[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorXi             *cutFail1Xi[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorXi             *cutPass2Xi[numOfMultBins*numOfChTypes];
        AliFemtoCutMonitorXi             *cutFail2Xi[numOfMultBins*numOfChTypes];
        //       AliFemtoShareQualityTPCEntranceSepPairCut                      *sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
        // AliFemtoPairCutAntiGamma     *sqpcetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoPairCutRadialDistance                   *sqpcetaphitpc[numOfMultBins*numOfChTypes];
        //AliFemtoShareQualityPairCut                   *sqpcetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0PairCut               *sqp1cetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0TrackPairCut          *sqp2cetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0TrackPairCut          *sqp3cetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoV0TrackPairCut          *sqp4cetaphitpc[numOfMultBins*numOfChTypes];

        AliFemtoXiTrackPairCut           *tXiTrackPairCut[numOfMultBins*numOfChTypes];



        //      AliFemtoChi2CorrFctn                                    *cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
        AliFemtoPairCutPt               *ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
        AliFemtoQinvCorrFctn            *cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
        AliFemtoQinvCorrFctn            *cqinvtpc[numOfMultBins*numOfChTypes];
        AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections *cdedpetaphi[numOfMultBins*numOfChTypes];
        AliFemtoCorrFctnDEtaDPhiSimple  *cdedpetaphinocorr[numOfMultBins*numOfChTypes];
        //      AliFemtoCorrFctnDYDPhiSimpleWithCorrections     *cdydpyphinocorr[numOfMultBins*numOfChTypes];
        AliFemtoCorrFctnDEtaDPhiCorrections *cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];

        //tak bylo za ROOT5
        // AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[numOfMultBins**numOfkTbins];
        // AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[numOfMultBins**numOfkTbins];
        // AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[numOfMultBins**numOfkTbins];
        // AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[numOfMultBins**numOfkTbins];
        AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[500];
        AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[500];
        AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[500];
        AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[500];

        AliFemtoCorrFctnNonIdDR         *cnonidtpc[numOfMultBins*numOfChTypes];


        //Invarnaint MASS

        AliFemtoCorrFctnMezonPhi *IMassKK[numOfMultBins*numOfChTypes];

        // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
        // *** Begin pion-pion (positive) analysis ***
        int aniter = 0;
        for (int imult = 0; imult < numOfMultBins; imult++)
        {
                if (runmults[imult])
                {
                        for (int ichg = 0; ichg < numOfChTypes; ichg++)
                        {
                                if (runch[ichg])
                                {
                                        //std::cout<<"Krzyk 1"<<std::endl;
                                        aniter = ichg * numOfMultBins + imult;
                                        anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(zvertbino, -10.0, 10.0, multbino, multbins[imult], multbins[imult+1]);
                                        anetaphitpc[aniter]->SetNumEventsToMix(8);
                                        anetaphitpc[aniter]->SetMinSizePartCollection(0);
                                        anetaphitpc[aniter]->SetVerboseMode(kFALSE);//~~~~~~~~~~~~~~~~

                                        //*** Event cut ***
                                        mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
                                        mecetaphitpc[aniter]->SetEventMult(0.0,100000);
                                        mecetaphitpc[aniter]->SetVertZPos(-10,10);//cm

                                        //std::cout<<"Krzyk 2"<<std::endl;

                                        //****** event monitors **********
                                        cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                                        cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                                        mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);


                                        //std::cout<<"Krzyk 3"<<std::endl;
/*
                                        //Study the collection multiplicity distribution
                                        cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                                        cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                                        mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
*/
                                        // ***** single particle track cuts *********

                                        dtc1etaphitpc[aniter] = new AliFemtoKKTrackCutFull();
                                        dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(3.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTOF450_500(2.0);
                                        dtc1etaphitpc[aniter]->UseNsigmaTOF450_500(true);
                                        dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
                                        dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
                                        dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
                                        dtc1etaphitpc[aniter]->SetCharge(1.0);
                                        dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);

                                        cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                                        cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                                        mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

                                        /*
                                        dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
                                        dtc1etaphitpc[aniter]->SetCharge(1.0);
                                        dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
                                        dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal);
                                        dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
                                        dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
                                        */

                                        //dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal);
                                        //dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
                                        //dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
                                        //dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);

                                        //std::cout<<"Krzyk 4"<<std::endl;

                                        if (ichg == 0 ||ichg == 1 ||ichg == 2)//kaons 3-5
                                          {
                                                  //std::cout<<"Krzyk 5"<<std::endl;
                                            dtc1etaphitpc[aniter]->SetPt(0.2,5);
                                            dtc1etaphitpc[aniter]->SetMass(KaonMass);
                                            dtc1etaphitpc[aniter]->SetMostProbableKaon();//cut on Nsigma in pT not p
                                          }

                                        dtc2etaphitpc[aniter] = new AliFemtoKKTrackCutFull();
                                        dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(3.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTOF450_500(2.0);
                                        dtc2etaphitpc[aniter]->UseNsigmaTOF450_500(true);
                                        dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
                                        dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
                                        dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
                                        dtc2etaphitpc[aniter]->SetCharge(-1.0);
                                        dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);


                                        /*
                                        dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
                                        dtc2etaphitpc[aniter]->SetCharge(-1.0);
                                        dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
                                        dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal);
                                        dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
                                        dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
                                        */
                                        if (ichg == 0 ||ichg == 1 ||ichg == 2)//kaons 0-2
                                          {
                                                  //std::cout<<"Krzyk 5"<<std::endl;
                                            dtc2etaphitpc[aniter]->SetPt(0,5);
                                            dtc2etaphitpc[aniter]->SetMass(KaonMass);
                                            dtc2etaphitpc[aniter]->SetMostProbableKaon();//cut on Nsigma in pT not p
                                          }
/*
                                        //--------------V0 cuts-------------------
                                        //V0 first particle cut
                                        dtc4etaphitpc[aniter] = new AliFemtoV0TrackCut();
                                        dtc4etaphitpc[aniter]->SetMass(LambdaMass);
                                        dtc4etaphitpc[aniter]->SetEta(0.8); //0.8
                                        if(ichg>=20 && ichg<=26)
                                          dtc4etaphitpc[aniter]->SetPt(1.4,maxPt); //0.5,5.0
                                        if(ichg>=13 && ichg<=19)
                                          dtc4etaphitpc[aniter]->SetPt(0.6,1.4);
                                        if(ichg>=27 && ichg<=33)
                                          dtc4etaphitpc[aniter]->SetPt(0.6,maxPt);

                                        std::cout<<"Krzyk 6"<<std::endl;
                                        dtc4etaphitpc[aniter]->SetEtaDaughters(0.8);
                                        dtc4etaphitpc[aniter]->SetPtPosDaughter(0.3,4.0);
                                        dtc4etaphitpc[aniter]->SetPtNegDaughter(0.16,4.0);
                                        dtc4etaphitpc[aniter]->SetTPCnclsDaughters(70);
                                        dtc4etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
                                        dtc4etaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        dtc4etaphitpc[aniter]->SetParticleType(0);
                                        dtc4etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
                                        dtc4etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
                                        dtc4etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
                                        dtc4etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
                                        dtc4etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
                                        dtc4etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
                                        dtc4etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
                                        dtc4etaphitpc[aniter]->SetRadiusV0Min(0.5);
                                        dtc4etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
                                        dtc4etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
                                        dtc4etaphitpc[aniter]->SetRequireTOFPion(false);
                                        dtc4etaphitpc[aniter]->SetRequireTOFProton(false);

                                        //V0 second particle cut
                                        dtc5etaphitpc[aniter] = new AliFemtoV0TrackCut();
                                        dtc5etaphitpc[aniter]->SetMass(LambdaMass);
                                        dtc5etaphitpc[aniter]->SetEta(0.8);
                                        if(ichg>=20 && ichg<=26)
                                          dtc5etaphitpc[aniter]->SetPt(1.4,maxPt);
                                        if(ichg>=13 && ichg<=19)
                                          dtc5etaphitpc[aniter]->SetPt(0.6,1.4);
                                        if(ichg>=27 && ichg<=33)
                                          dtc5etaphitpc[aniter]->SetPt(0.6,maxPt);
                                        std::cout<<"Krzyk 7"<<std::endl;
                                        dtc5etaphitpc[aniter]->SetEtaDaughters(0.8);
                                        dtc5etaphitpc[aniter]->SetPtPosDaughter(0.16,4.0);
                                        dtc5etaphitpc[aniter]->SetPtNegDaughter(0.3,4.0);
                                        dtc5etaphitpc[aniter]->SetTPCnclsDaughters(70);
                                        dtc5etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
                                        dtc5etaphitpc[aniter]->SetParticleType(1);
          dtc5etaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        dtc5etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
                                        dtc5etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
                                        dtc5etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
                                        dtc5etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
                                        dtc5etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
                                        dtc5etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
                                        dtc5etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
                                        dtc5etaphitpc[aniter]->SetRadiusV0Min(0.5);
                                        dtc5etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
                                        dtc5etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
                                        dtc5etaphitpc[aniter]->SetRequireTOFPion(false);
                                        dtc5etaphitpc[aniter]->SetRequireTOFProton(false);
                                        std::cout<<"Krzyk 8"<<std::endl;
                                        //Cascade cuts
                                        //from J. Buxton
                                        //xi cut
                                        //NOTE: the SetMass call actually is important
                                        //      This should be set to the mass of the particle of interest, here the Xi
                                        //      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)

                                        //Xi -> Lam Pi-
                                        tXiCut[aniter] = new AliFemtoXiTrackCut();
                                        // %%%%%%%%%%%%%%%%%%%%%%%% Version 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        tXiCut[aniter]->SetChargeXi(-1);
                                        tXiCut[aniter]->SetParticleTypeXi(0);  //kXiMinus = 0
                                        tXiCut[aniter]->SetPtXi(0.8,100);
                                        tXiCut[aniter]->SetEtaXi(0.8);
                                        tXiCut[aniter]->SetMass(XiMass);
                                        tXiCut[aniter]->SetInvariantMassXi(XiMass-0.003,XiMass+0.003);
                                        tXiCut[aniter]->SetMaxDecayLengthXi(100.);
                                        tXiCut[aniter]->SetMinCosPointingAngleXi(0.9992);
                                        tXiCut[aniter]->SetMaxDcaXi(100);
                                        //XiDaughters
                                        tXiCut[aniter]->SetMaxDcaXiDaughters(0.3);


                                        //Bachelor cuts (here = PiM)
                                        tXiCut[aniter]->SetMinDcaXiBac(0.03);
                                        tXiCut[aniter]->SetEtaBac(0.8);
                                        tXiCut[aniter]->SetTPCnclsBac(70);
                                        tXiCut[aniter]->SetPtBac(0.,100.);
                                        tXiCut[aniter]->SetStatusBac(AliESDtrack::kTPCrefit);  //yes or no?


                                        //Lambda cuts (regular V0)
                                        tXiCut[aniter]->SetParticleType(0); //0=lambda
                                        tXiCut[aniter]->SetMinDcaV0(0.1);
                                        tXiCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
                                        tXiCut[aniter]->SetMinCosPointingAngle(0.998);
                                        tXiCut[aniter]->SetEta(0.8);
                                        tXiCut[aniter]->SetPt(0.0,100);
                                        tXiCut[aniter]->SetOnFlyStatus(kFALSE);
                                        tXiCut[aniter]->SetMaxV0DecayLength(100.);
                                        //Lambda daughter cuts
                                        tXiCut[aniter]->SetMinDaughtersToPrimVertex(0.1,0.1);
                                        tXiCut[aniter]->SetMaxDcaV0Daughters(0.8);
                                        tXiCut[aniter]->SetEtaDaughters(0.8);
                                        tXiCut[aniter]->SetPtPosDaughter(0.,99); //0.5 for protons
                                        tXiCut[aniter]->SetPtNegDaughter(0.,99); //0.16 for pions
                                        tXiCut[aniter]->SetTPCnclsDaughters(70);Vertex
                                        tXiCut[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?
                                        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                        tXiCut[aniter]->SetMinvPurityAidHistoXi("XiPurityAid","XiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
                                        tXiCut[aniter]->SetMinvPurityAidHistoV0("LambdaPurityAid","LambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);

                                        //antiXi cut
                                        //NOTE: the SetMass call actually is important
                                        //      This should be set to the mass of the particle of interest, here the Xi
                                        //      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)

                                        //AXi -> ALam Pi+

                                        tAXiCut[aniter] = new AliFemtoXiTrackCut();

                                        // %%%%%%%%%%%%%%%%%%%%%%%% Version 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                        tAXiCut[aniter]->SetChargeXi(1);
                                        tAXiCut[aniter]->SetParticleTypeXi(1); //kXiPlus = 1
                                        tAXiCut[aniter]->SetPtXi(0.8,100);
                                        tAXiCut[aniter]->SetEtaXi(0.8);
                                        tAXiCut[aniter]->SetMass(XiMass);
                                        tAXiCut[aniter]->SetInvariantMassXi(XiMass-0.003,XiMass+0.003);
                                        tAXiCut[aniter]->SetMaxDecayLengthXi(100.0);
                                        tAXiCut[aniter]->SetMinCosPointingAngleXi(0.9992);
                                        tAXiCut[aniter]->SetMaxDcaXi(100);
                                        //XiDaughters
                                        tAXiCut[aniter]->SetMaxDcaXiDaughters(0.3);


                                        //Bachelor cuts (here = PiP)
                                        tAXiCut[aniter]->SetMinDcaXiBac(0.03);
                                        tAXiCut[aniter]->SetEtaBac(0.8);
                                        tAXiCut[aniter]->SetTPCnclsBac(70);
                                        tAXiCut[aniter]->SetPtBac(0.,100);
                                        tAXiCut[aniter]->SetStatusBac(AliESDtrack::kTPCrefit);  //yes or no?


                                        //AntiLambda cuts (regular V0)
                                        tAXiCut[aniter]->SetParticleType(1); //1=anti-lambda
                                        tAXiCut[aniter]->SetMinDcaV0(0.1);
                                        tAXiCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
                                        tAXiCut[aniter]->SetMinCosPointingAngle(0.998);
                                        tAXiCut[aniter]->SetEta(0.8);
                                        tAXiCut[aniter]->SetPt(0.,100);
                                        tAXiCut[aniter]->SetOnFlyStatus(kTRUE);  //!!!!!!!!!!!!!!!!!!!!!!! Set to kTRUE!                                        tAXiCut[aniter]->SetMaxV0DecayLength(100.);
                                        //Lambda daughter cuts
                                        tAXiCut[aniter]->SetMinDaughtersToPrimVertex(0.1,0.1);
                                        tAXiCut[aniter]->SetMaxDcaV0Daughters(0.8);
                                        tAXiCut[aniter]->SetEtaDaughters(0.8);
                                        tAXiCut[aniter]->SetPtPosDaughter(0.,99); //0.16 for pions
                                        tAXiCut[aniter]->SetPtNegDaughter(0.,99); //0.5 for anti-protons
                                        tAXiCut[aniter]->SetTPCnclsDaughters(70);
                                        tAXiCut[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?

                                        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                        tAXiCut[aniter]->SetMinvPurityAidHistoXi("AXiPurityAid","AXiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
                                        tAXiCut[aniter]->SetMinvPurityAidHistoV0("AntiLambdaPurityAid","AntiLambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);

                                        */

                                        //****** DCA ******
                                        if(owndca){
                                                //std::cout<<"Krzyk 9"<<std::endl;
                                          dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01);      //     DCA xy
                                          //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
                                          dtc1etaphitpc[aniter]->SetMaxImpactZ(2);      //DCA Z
                                          dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01);      //     DCA xy
                                          dtc2etaphitpc[aniter]->SetMaxImpactZ(2);      //DCA Z
                                          if (ichg == 9){dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01);       //      DCA xy
                                          dtc3etaphitpc[aniter]->SetMaxImpactZ(2);}     //DCA Z
                                        }
                                        //****** Track quality cuts ******

                                        if(owncuts){
                                                //std::cout<<"Krzyk 9"<<std::endl;
                                          dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                          dtc1etaphitpc[aniter]->SetminTPCncls(70);
                                          dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                          dtc1etaphitpc[aniter]->SetLabel(kFALSE);
                                          //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                          dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0); // pisac
                                          //dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

                                          dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                          dtc2etaphitpc[aniter]->SetminTPCncls(70);
                                          dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                          dtc2etaphitpc[aniter]->SetLabel(kFALSE);
                                          //    dtc2etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                          dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                          //    dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                          if (ichg == 9){
                                                  //std::cout<<"Krzyk 10"<<std::endl;
                                            dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                            dtc3etaphitpc[aniter]->SetminTPCncls(70);
                                            dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                            dtc3etaphitpc[aniter]->SetLabel(kFALSE);
                                            //  dtc3etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                            dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                            //  dtc3etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                          }

                                        }
                                        //**************** track Monitors ***************

                                        if(ifMonitors)//ichg>8)
                                          {
                                                //std::cout<<"Krzyk 11"<<std::endl;
/*
                                            if(ichg>=17 && ichg<=33){
                                                        std::cout<<"Krzyk 12"<<std::endl;
                                              // //V0 monitors (memory leak problems?)
                                              cutPass1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass1%stpcM%i", chrgs[ichg], imult));
                                              cutFail1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail1%stpcM%i", chrgs[ichg], imult));
                                              dtc4etaphitpc[aniter]->AddCutMonitor(cutPass1V0[aniter], cutFail1V0[aniter]);

                                              cutPass2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass2%stpcM%i", chrgs[ichg], imult));
                                              cutFail2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail2%stpcM%i", chrgs[ichg], imult));
                                              dtc5etaphitpc[aniter]->AddCutMonitor(cutPass2V0[aniter], cutFail2V0[aniter]);

                                              anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                              anetaphitpc[aniter]->SetEnablePairMonitors(enablePairMonitors);

                                            }
*/
                                                        //Przemyśleć jak zrobić lepiej

                                                        //numOfMultBins
                                            //FULL
                                            //if(ichg<2 || ichg==0||ichg==1 || ichg == 2){
                                              cutPass3YPtetaphitpc[0] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", "Kp", imult),KaonMass);
                                              cutPass3YPtetaphitpc[1] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", "Km", imult),KaonMass);
                                              cutFail3YPtetaphitpc[0] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", "Kp", imult),KaonMass);
                                              cutFail3YPtetaphitpc[1] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", "Km", imult),KaonMass);
                                            //}
                                            //if(ichg==9) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
                                            //if(ichg==0)
                                            dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[0], cutFail3YPtetaphitpc[0]);
                                            //if(ichg==1)
                                            dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[1], cutFail3YPtetaphitpc[1]);




                                                        std::cout<<endl;
                                                        //std::cout<<"CHUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU"<<endl;
                                              cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
                                              cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),1);
                                                        dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);




                                          }

                                        //******** Two - track cuts ************
                                        sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
                                        sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);
                                        sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);
                 sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                                                                 sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);
                 sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
                 sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);
                                                                 sqpcetaphitpc[aniter]->SetPhiStarMin(kFALSE);
          sqpcetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        // sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
                                        //sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
                                        // sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);     // two track cuts on splitting and merging  //1- wylaczany 0 -wlaczany
                                        // sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);   //  ile moga miec wspolnych klastrow //1 - wylaczany, 0.05 - wlaczany
                                        // sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                                        // sqpcetaphitpc[aniter]->SetMaximumRadius(0.82);
                                        // sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);
                                        // sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
                                        // sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

                                        if (gammacut == 0)
                                          {
                                            sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
                                            sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
                                          }
                                        else if (gammacut == 1)
                                          {
                                            sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
                                            sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
                                          }
                                        // sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
                                        // sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
                                        // sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

                                        //V0 two-track cuts

/*
                                        sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
          sqp1cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                                        sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                                        sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                                        //sqp1cetaphitpc[aniter]->SetMinAvgSeparation(0,3);
                                        //sqp1cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
                                        //sqp1cetaphitpc[aniter]->SetMinAvgSeparation(2,0);
                                        //sqp1cetaphitpc[aniter]->SetMinAvgSeparation(3,3);

                                        sqp2cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-proton
          sqp2cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        sqp2cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track                                        sqp2cetaphitpc[aniter]->SetShareFractionMax(0.05);
                                        //sqp2cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                                        sqp2cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                                        sqp2cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                                        sqp2cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                                        sqp2cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
                                        //sqp2cetaphitpc[aniter]->SetMinAvgSeparation(0,11); //0 - track-pos, 1 - track-neg
                                        //sqp2cetaphitpc[aniter]->SetMinAvgSeparation(1,0);

                                        sqp3cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
          sqp3cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        sqp3cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track                                        sqp3cetaphitpc[aniter]->SetShareFractionMax(0.05);
                                        //sqp3cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                                        sqp3cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                                        sqp3cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                                        sqp3cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                                        sqp3cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
                                        //sqp3cetaphitpc[aniter]->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
                                        //sqp3cetaphitpc[aniter]->SetMinAvgSeparation(1,11);

                                        sqp4cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
          sqp4cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
                                        sqp4cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track                                        sqp4cetaphitpc[aniter]->SetShareFractionMax(0.05);
                                        //sqp4cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                                        sqp4cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                                        sqp4cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                                        sqp4cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                                        //SetMinAvgSeparation w if'ach ponizej


                                        tXiTrackPairCut[aniter] = new AliFemtoXiTrackPairCut(); //xi-proton, all combinations

*/

                                        //***** Setting cuts ***********
                                        // setting event cut
                                        anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                        //setting single track cuts
                                        if(ichg==0) //positive like-sign
                                        {
                                          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
                                        }
                                        if(ichg==1)//negative like-sign
                                        {
                                          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                                        }
                                        if(ichg==2)//unlike-sign
                                        {
                                          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                          anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                                        }

                                        IMassKK[aniter] = new AliFemtoCorrFctnMezonPhi(Form("cutPass%stpcM%i", chrgs[ichg], imult), 3000, 0.5, 2, KaonMass, KaonMass, 10, 0, 10);
                                        anetaphitpc[aniter]->AddCorrFctn(IMassKK[aniter]);

                                        Manager->AddAnalysis(anetaphitpc[aniter]);
                                }
                        }
                }
        }

        return Manager;
}
