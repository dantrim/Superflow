#ifndef SUPERFLOW_SUPERWEIGHT_H
#define SUPERFLOW_SUPERWEIGHT_H

//////////////////////////////////////////////////////////////////////////////
//
// Superweight
//
// structure for holding the components of the event weigth
//
// daniel.joseph.antrim@cern.ch
// July 2017
//
//////////////////////////////////////////////////////////////////////////////


namespace sflow {

    class Superweight
    {

        public :
            Superweight();

            double product();
            double product_multi();
            void reset();

            double susynt; // weight returned from SusyNtuple/MCWeighter::getMCWeight(...)
            double susynt_multi; // multi period susynt weight
            double pileup;
            double lepSf;
            double trigSf;
            double btagSf;
            double jvtSf;

    }; // class

} // namespace


#endif
