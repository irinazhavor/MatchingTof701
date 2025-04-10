# MatchingTof701
1. Matching_Tofs.C - macro for carrying out the procedure of tracks from the inner tracking system to the TOF701 and matching them with TOF701 hits.
   
   input - dst from EOS.

   Info about the reconstucted vertices (PrimaryVertex.),  gem tracks (StsVector) and hits (BmnTof701Hit) is used.

   output - .root file with a tree, containing info about Tracks (extrapolated and matched in TOF701) and Tof701 matched hits. 
  
2. convertBmn_run8_Tofs.C

   RDataFrame interface is used to procces data.

   input - 1) .root file from the previous macro with independent matching;
           2) dst from EOS;
           3) digi from EOS (optional, not used for the matching examination, but one would need to remove .Define lines, utilising branches from digi .root file, i.e. 1408-1448);

   output - .root file, containing branches with info about both Global Tracks info and independent tracks, as well as TOF701 hits.
   
3. qa_converter_tofs.cc - plots final distributions, i.e. beta vs pq, m2 vs pq for both version of matching and saves them in the .root file
   
   Output histograms in .pdf are uploaded as well.
