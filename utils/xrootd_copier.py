import os 

din = "root://xrootd-cms.infn.it///eos/cms/store/group/phys_b2g/CASE/h5_files/UL/merged/"
odir  = "root://cmseosmgm01.fnal.gov:1094//store/user/oamram/case/sig_files/"


flist = [
        "QstarToQW_M_2000_mW_170_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX2000_MY400_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_2000_mW_25_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX2000_MY400_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_2000_mW_400_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX2000_MY80_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_2000_mW_80_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX2000_MY80_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_3000_mW_170_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX2000_MY80_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_3000_mW_25_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX2000_MY80_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_3000_mW_400_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX3000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_3000_mW_80_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX3000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_5000_mW_170_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX3000_MY170_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_5000_mW_25_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX3000_MY170_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_5000_mW_400_TuneCP2_13TeV-pythia8_TIMBER.h5",                                 "XToYYprimeTo4Q_MX3000_MY25_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "QstarToQW_M_5000_mW_80_TuneCP2_13TeV-pythia8_TIMBER.h5",                                  "XToYYprimeTo4Q_MX3000_MY25_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "RSGravitonToGluonGluon_kMpl01_M_1000_TuneCP5_13TeV_pythia8_TIMBER.h5",                    "XToYYprimeTo4Q_MX3000_MY25_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "RSGravitonToGluonGluon_kMpl01_M_2000_TuneCP5_13TeV_pythia8_TIMBER.h5",                    "XToYYprimeTo4Q_MX3000_MY25_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "RSGravitonToGluonGluon_kMpl01_M_3000_TuneCP5_13TeV_pythia8_TIMBER.h5",                    "XToYYprimeTo4Q_MX3000_MY400_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "RSGravitonToGluonGluon_kMpl01_M_5000_TuneCP5_13TeV_pythia8_TIMBER.h5",                    "XToYYprimeTo4Q_MX3000_MY400_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M2000_Mr170_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY400_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M2000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY400_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M3000_Mr170_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY80_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M3000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY80_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M5000_Mr170_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY80_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WkkToWRadionToWWW_M5000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",                  "XToYYprimeTo4Q_MX3000_MY80_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp2000_Bp170_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp2000_Bp25_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp2000_Bp400_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY170_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp2000_Bp80_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY170_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp3000_Bp170_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY25_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp3000_Bp25_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY25_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp3000_Bp400_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY25_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp3000_Bp80_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY25_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp5000_Bp170_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY400_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp5000_Bp25_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY400_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp5000_Bp400_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",             "XToYYprimeTo4Q_MX5000_MY400_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "WpToBpT_Wp5000_Bp80_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER.h5",              "XToYYprimeTo4Q_MX5000_MY400_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",  "XToYYprimeTo4Q_MX5000_MY80_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",   "XToYYprimeTo4Q_MX5000_MY80_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY170_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",  "XToYYprimeTo4Q_MX5000_MY80_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY170_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",   "XToYYprimeTo4Q_MX5000_MY80_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY25_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",   "YtoHH_Htott_Y2000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY25_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",    "YtoHH_Htott_Y3000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY25_MYprime400_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",   "YtoHH_Htott_Y5000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY25_MYprime80_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",    "ZpToTpTp_Zp2000_Tp400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY400_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",  "ZpToTpTp_Zp3000_Tp400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        "XToYYprimeTo4Q_MX2000_MY400_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",   "ZpToTpTp_Zp5000_Tp400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5",
        ]

for f in flist:
    print(f)

    os.system("xrdcp %s/%s ."  % (din, f))
    os.system("xrdcp -f %s %s" % (f, odir))
    os.system("rm %s" % f)
