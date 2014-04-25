
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "handler.h"
    
namespace gmx {

bool EnergyHandler::readFiles()
{
    gmx_bool     bCont;
    ener_file_t  fp = open_enx(fileName_.begin()->c_str(), "r");
    t_enxframe   frame;
    int          timecheck = 0;
    bool         bOK;

    printf("There are %d analysis tools registered in the energy handler.\n",
           (int) eap_.size());
    {
        int          nre;
        gmx_enxnm_t *enm = NULL;
        do_enxnms(fp, &nre, &enm);
        bOK = true;
        for(EnergyAnalysisPtrIterator eapi = eap_.begin();
            bOK && (eapi < eap_.end()); ++eapi)
        {
            EnergyAnalysisPtr eap = *eapi;
            bOK = bOK && eap->initAnalysis(nre, enm);
        }
        free_enxnms(nre, enm);
    }
    if (!bOK)
    {
        return false;
    }

    init_enxframe(&frame);
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &frame);
            if (bCont)
            {
                timecheck = check_times(frame.t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            for(EnergyAnalysisPtrIterator eapi = eap_.begin();
                bCont && (eapi < eap_.end()); ++eapi)
            {
                EnergyAnalysisPtr eap = *eapi;
                bCont = bCont && eap->addAnalysisFrame(&frame);
            }
        }
    }
    while (bCont);
    close_enx(fp);
    /* Printing a new line, just because the gromacs library prints step info
       while reading. */
    fprintf(stderr, "\n");

    for(EnergyAnalysisPtrIterator eapi = eap_.begin();
        bOK && (eapi < eap_.end()); ++eapi)
    {
        EnergyAnalysisPtr eap = *eapi;
        bOK = bOK && eap->finalizeAnalysis();
    }
    return true;
}

}
