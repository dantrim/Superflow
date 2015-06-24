// Superflow.cxx
//

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>
#include "TGraphAsymmErrors.h"


#include "Superflow/Superflow.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"

using namespace std;

namespace sflow {

    ////////////////////////////////////////////
    // Constructor
    ////////////////////////////////////////////
    Superflow::Superflow()
    {
        m_dbg = true;
        setSelectTaus(true);

        m_runMode = SuperflowRunMode::null;
        
        /////////////////////////
        // Matrix Method
        /////////////////////////
        m_matrixFilename = "";
        m_use2dparametrization = false;
        m_allconfigured = false;
        m_matrix = nullptr;
        m_fake_region = "";
        selectBaseLineLeptons = false;

        /////////////////////////
        // Charge-flip
        /////////////////////////
        m_do_qflip = false;
        

        m_CutStore_Name_Exists = false;
        m_CutStoreUntitled     = 1;

        m_countWeights = true;

        m_outputFileName      = "";
        m_entry_list_FileName = "entrylist_";
        m_tree_name_auto      = "";
        m_outputFile          = nullptr;
        m_entryListFile       = nullptr;
        m_HFT                 = nullptr;

        m_output_array = nullptr;

        m_HFT_array = nullptr;

        m_entry_list_total       = nullptr;
        m_entry_list_single_tree = nullptr;

        m_period = ATLAS_period::null;
        m_stream = ATLAS_stream::null;

        m_varState = SupervarState::closed;
        m_superVar_hasFunction = false;
        m_superVar_hasNiceName = false;
        m_superVar_hasHFTName  = false;
        m_superVar_Untitled    = 1;

        m_mcWeighter = nullptr;

        m_nullExprFloat  = [](Superlink* sl, var_float*) -> double { return 0.0; };
        m_nullExprDouble = [](Superlink* sl, var_double*) -> double { return 0.0; };
        m_nullExprInt    = [](Superlink* sl, var_int*) -> int { return 0; };
        m_nullExprBool   = [](Superlink* sl, var_bool*) -> bool { return false; };
        m_nullExprVoid   = [](Superlink* sl, var_void*) {};

        m_varFloat  = nullptr;
        m_varDouble = nullptr;
        m_varInt    = nullptr;
        m_varBool   = nullptr;

        m_varFloat_array  = nullptr;
        m_varDouble_array = nullptr;
        m_varInt_array    = nullptr;
        m_varBool_array   = nullptr;

        m_sysState = SupersysState::closed;


        m_sysTemplate.reset();
        m_sys_hasNiceName   = false;
        m_sys_hasTreeName   = false;
        m_sys_hasType       = false;
        m_sys_hasSystematic = false;

        m_singleEventSyst = NtSys::NOM;
        m_RunSyst         = new Supersys(SupersysType::central);

        m_tree_leafs_size    = 0;
        m_weight_leaf_offset = 0;

        m_trigObj = nullptr;

        m_input_chain = nullptr;

        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::A] = "periodA.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::B] = "periodB.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::C] = "periodC.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::D] = "periodD.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::E] = "periodE.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::G] = "periodG.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::H] = "periodH.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::I] = "periodI.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::J] = "periodJ.physics_Egamma.PhysCont";
        m_data_periods[ATLAS_stream::Egamma][ATLAS_period::L] = "periodL.physics_Egamma.PhysCont";

        m_data_periods[ATLAS_stream::Muons][ATLAS_period::A] = "periodA.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::B] = "periodB.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::C] = "periodC.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::D] = "periodD.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::E] = "periodE.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::G] = "periodG.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::H] = "periodH.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::I] = "periodI.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::J] = "periodJ.physics_Muons.PhysCont";
        m_data_periods[ATLAS_stream::Muons][ATLAS_period::L] = "periodL.physics_Muons.PhysCont";

        m_NtSys_to_string[Susy::NtSys::EES_Z_UP]     = "EESZUP";
        m_NtSys_to_string[Susy::NtSys::EES_Z_DN]     = "EESZDOWN";
        m_NtSys_to_string[Susy::NtSys::EES_MAT_UP]   = "EESMATUP";
        m_NtSys_to_string[Susy::NtSys::EES_MAT_DN]   = "EESMATDOWN";
        m_NtSys_to_string[Susy::NtSys::EES_PS_UP]    = "EESPSUP";
        m_NtSys_to_string[Susy::NtSys::EES_PS_DN]    = "EESPSDOWN";
        m_NtSys_to_string[Susy::NtSys::EES_LOW_UP]   = "EESLOWUP";
        m_NtSys_to_string[Susy::NtSys::EES_LOW_DN]   = "EESLOWDOWN";
        m_NtSys_to_string[Susy::NtSys::EER_UP]       = "EERUP";
        m_NtSys_to_string[Susy::NtSys::EER_DN]       = "EERDOWN";
        m_NtSys_to_string[Susy::NtSys::MS_UP]        = "MSUP";
        m_NtSys_to_string[Susy::NtSys::MS_DN]        = "MSDOWN";
        m_NtSys_to_string[Susy::NtSys::ID_UP]        = "IDUP";
        m_NtSys_to_string[Susy::NtSys::ID_DN]        = "IDDOWN";
        m_NtSys_to_string[Susy::NtSys::JES_UP]       = "JESUP";
        m_NtSys_to_string[Susy::NtSys::JES_DN]       = "JESDOWN";
        m_NtSys_to_string[Susy::NtSys::JER]          = "JER";
        m_NtSys_to_string[Susy::NtSys::SCALEST_UP]   = "SCALESTUP";
        m_NtSys_to_string[Susy::NtSys::SCALEST_DN]   = "SCALESTDOWN";
        m_NtSys_to_string[Susy::NtSys::RESOST]       = "RESOST";
        m_NtSys_to_string[Susy::NtSys::TRIGSF_EL_UP] = "TRIGSFELUP";
        m_NtSys_to_string[Susy::NtSys::TRIGSF_EL_DN] = "TRIGSFELDN";
        m_NtSys_to_string[Susy::NtSys::TRIGSF_MU_UP] = "TRIGSFMUUP";
        m_NtSys_to_string[Susy::NtSys::TRIGSF_MU_DN] = "TRIGSFMUDN";
        m_NtSys_to_string[Susy::NtSys::TES_UP]       = "TESUP";
        m_NtSys_to_string[Susy::NtSys::TES_DN]       = "TESDOWN";
        m_NtSys_to_string[Susy::NtSys::JVF_UP]       = "JVFUP";
        m_NtSys_to_string[Susy::NtSys::JVF_DN]       = "JVFDOWN";
    }

    ////////////////////////////////////////////
    // Destructor
    ////////////////////////////////////////////
    Superflow::~Superflow()
    {}

    // --------------------------------------------------------------- //
    //  Superflow operator implementation
    //  Superflow operator implementation
    //  Superflow operator implementation
    // --------------------------------------------------------------- //

    ////////////////////////////////////////////
    // Cut operators
    ////////////////////////////////////////////

    /* ---------------------------
        CutName operator

        Feeds in a human-readable/descriptive name for a cut to be placed.
        --> Checks that the SupersysState and SupervarState are closed, meaning
            that any previous definition of a systematic or variable is finished
            and succesfully saved/stored.

        usage:
            <SuperflowObject> << CutName("lead lepton pT > 20 GeV");
       ---------------------------
    */
    Superflow& Superflow::operator<<(CutName cut_)
    {
        if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
            m_CutStore_Name_Exists = true;
            m_CutStoreNames.push_back(cut_.name);

            cout << app_name << "New cut: " << cut_.name << endl;
        }
        else {
            cout << app_name << "ERROR (Fatal): Cutflow operations are incorrectly ordered.";
            exit(1);
        }
        return *this;
    }
    
    /* ----------------------------
        Cut function operator

        Feeds in the (lambda) function expression to be evaluated for the cut.
        --> Checks that the SupersysState and SupervarState are closed, meaning
            that any previous definition of a systematic or variable is finished
            and succesfully saved/stored.
        --> Checks m_CutStore_Name_Exists variable to be sure that we have processed
            the CutName operator prior this call (i.e. so that this cut will have 
            a descriptive name). If this is not the case, a default name is given.
        --> Adds a place in the raw and weighted counters (m_RawCounter and m_WeightCounter).

         usage:
            <SuperflowObject> << CutName("at least 2 base leptons") << [](Superlink* sl) -> bool {
                return (sl->baseLeptons.size() >= 2);
            };
        --------------------------
    */
    Superflow& Superflow::operator<<(std::function<bool(Superlink*)> cut_)
    {
        if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
            m_CutStore.push_back(cut_);

            if (m_CutStore_Name_Exists) {
                m_CutStore_Name_Exists = false;
            }
            else {
                m_CutStoreNames.push_back("Untitled-" + to_string(m_CutStoreUntitled));
                cout << app_name << "New cut: " << "Untitled-" << m_CutStoreUntitled << endl;
                m_CutStoreUntitled++;
            }

            m_RawCounter.push_back(0.0);
            m_WeightCounter.push_back(0.0);
        }
        else {
            cout << app_name << "ERROR (Fatal): Cutflow operations are incorrectly ordered.";
            exit(1);
        }
        return *this;
    }

    ////////////////////////////////////////////
    // HFT operators
    //   These set what will be stored in the
    //   output ntuples
    ////////////////////////////////////////////

    /* ----------------------------
        NewVar operator
            
        Feeds in a human-readable/descriptive name for a variable to store in the output trees.
        --> This is NOT the name that will show up as the leaf name (i.e. when you do TBrowser/TTree::Draw())
        --> Checks that the SupersysState and SupervarState are closed, meaning
            that any previous definition of a systematic or variable is finished
            and succesfully saved/stored.
        --> Name is stored in vector m_varNiceName
        --> m_superVar_hasNiceName is set to true to ensure that the ordering of HFT var creation is correct.
            This is set to false once the full variable is saved/stored.
        --> Variable definition is to take place after the NewVar is fed in, so we set the 
            m_varState to SupervarState::open. This will be set to SupervarState::closed once
            the variable is saved/stored.

        usage:
            <SuperflowObject> << NewVar("stransverse mass");
 
       ----------------------------
    */
    Superflow& Superflow::operator<<(NewVar new_var_name) // this is the NiceName
    {
        if (m_varState == SupervarState::closed && m_sysState == SupersysState::closed) {
            m_varNiceName.push_back(new_var_name.name);
            m_superVar_hasNiceName = true;
            m_varState = SupervarState::open;
        }
        else {
            cout << app_name << "ERROR (Fatal): Close the Var using SaveVar().";
            exit(1);
        }
        return *this;
    }

    /* ----------------------------
        HFTname operator
        
        Feed in the name of the variable that will appear in the output ntuples (i.e. when you 
        do TTree::Draw() or look in the TBrowser).
        --> We must have already provided a "nice name" using the NewVar operator which means that
            m_varState == SupervarState::open.
        --> m_superVar_hasHFTName must be false, this operator provides the HFTname and sets this to true.
        --> The leaf name is stored in the vector m_varHFTName 

        usage:
            <SuperflowObject> << NewVar("stransverse mass") << HFTname("mT2");

       ----------------------------
    */
    Superflow& Superflow::operator<<(HFTname hft_name)
    {
        if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasHFTName) {
            m_varHFTName.push_back(hft_name.name);
            m_superVar_hasHFTName = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): Open a new Var using NewVar().";
            exit(1);
        }
        return *this;
    }

    /* ----------------------------
        Variable function operators
        
        Provide a (lambda) function expression using the Superlink object to 
        calculate/provide the values for the variables whose names we have
        already provided using the NewVar and HFTname operators.
        --> There are function expressions for each of the conceivable variable types:
                - float
                - double
                - int
                - bool
                - void
            and the expression is fed into the appropriate vector of functions. 
            E.G. m_varExprFloat is type vector<std::function<double(Superlink*, var_float*)>> 
        --> var_float, var_double, var_int, var_void, and var_bool are dummy classes that are
            used only for type identification in the signatures of these different operators
            (see Superflow/Supervar.h).
        --> A null expression is provided to all other function vectors that are not of
            the type of the current function expression being fed in. In this way all of the
            function vectors will be of the same length and the indices of the functions will
            unambiguously line up with the vectors for the NewVar and HFTname vectors.
        --> m_varState must == SupervarState::open, indicating that we are currently in the process of
            storing a variable.
        --> m_superVar_hasFunction must be false since we have, at this point, only just opened the
            the m_varState and have only provided names for this variable.
        --> m_varType holds values of the enum SupervarType -- one for each function expression 
            fed in and is essentially used as a way to determine how many variables we have
            to store. 

        --> NOTE: for var_void type function expressions the NewVar and HFTname operators do not need
                to be called as they are set automatically within void function. The checks for the m_varState, etc...
                are also not checked. The main purpose of the var_void type function expression operator
                is to, for example, clear a global object that is defined in the executable outside of 
                the function expressions but used across many of them.
                E.G.
                    JetVector central_light_jets;
                    <VAR OPERATOR> { central_light_jets.push_back(foo) }
                        ...
                    <VAR OPERATOR> { central_light_jets.at(0)->Pt() }
                        ...
                    <SuperflowObject> << [&](Superlink* sl, var_void*) { central_light_jets.clear(); };

        usage:
            <SuperflowObject> << NewVar("leading lepton transverse momenta");
            <SuperflowObject> << HFTname("lept1Pt");
            <SuperflowObject> << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Pt(); }

       ----------------------------
    */

    // float function
    Superflow& Superflow::operator<<(std::function<double(Superlink*, var_float*)> var_)
    {
        if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
            m_varExprFloat.push_back(var_); // fill
            m_varExprDouble.push_back(m_nullExprDouble);
            m_varExprInt.push_back(m_nullExprInt);
            m_varExprBool.push_back(m_nullExprBool);
            m_varExprVoid.push_back(m_nullExprVoid);

            m_varType.push_back(SupervarType::sv_float);
            m_superVar_hasFunction = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
            exit(1);
        }
        return *this;
    }
    // double function
    Superflow& Superflow::operator<<(std::function<double(Superlink*, var_double*)> var_)
    {
        if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
            m_varExprFloat.push_back(m_nullExprFloat);
            m_varExprDouble.push_back(var_); // fill
            m_varExprInt.push_back(m_nullExprInt);
            m_varExprBool.push_back(m_nullExprBool);
            m_varExprVoid.push_back(m_nullExprVoid);

            m_varType.push_back(SupervarType::sv_double);
            m_superVar_hasFunction = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
            exit(1);
        }
        return *this;
    }
    // int function
    Superflow& Superflow::operator<<(std::function<int(Superlink*, var_int*)> var_)
    {
        if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
            m_varExprFloat.push_back(m_nullExprFloat);
            m_varExprDouble.push_back(m_nullExprDouble);
            m_varExprInt.push_back(var_); // fill
            m_varExprBool.push_back(m_nullExprBool);
            m_varExprVoid.push_back(m_nullExprVoid);

            m_varType.push_back(SupervarType::sv_int);
            m_superVar_hasFunction = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
            exit(1);
        }
        return *this;
    }
    // bool function
    Superflow& Superflow::operator<<(std::function<bool(Superlink*, var_bool*)> var_)
    {
        if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
            m_varExprFloat.push_back(m_nullExprFloat);
            m_varExprDouble.push_back(m_nullExprDouble);
            m_varExprInt.push_back(m_nullExprInt);
            m_varExprBool.push_back(var_); // fill
            m_varExprVoid.push_back(m_nullExprVoid);

            m_varType.push_back(SupervarType::sv_bool);
            m_superVar_hasFunction = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
            exit(1);
        }
        return *this;
    }
    // void function
    Superflow& Superflow::operator<<(std::function<void(Superlink*, var_void*)> var_)
    {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(var_);// fill

        m_varNiceName.push_back("void");
        m_varHFTName.push_back("void");
        m_varType.push_back(SupervarType::sv_void);
        return *this;
    }

    /* ----------------------------
        SaveVar operator
        
        Check that we have all of the required pieces for storing a variable in the 
        output ntuples.
        
        --> If the var state is still open and we have defined a function (with a name) then we
            have all of the pieces. Set these flags back to false so that we can be ready to
            start defining the next variable to be stored.
        --> If there is no name, provide the default one.

        usage:
            <SuperflowObject> << NewVar("sub-leading jet eta");
            <SuperflowObject> << HFTname("jet2Eta");
            <SuperflowObject> << [](Superlink* sl, var_float*) -> double {
                return sl->jets->size() >= 2 ? sl->jets->at(1)->Eta() : 0.0;
            };
            <SuperflowObject> << SaveVar();
            
       ----------------------------
    */
    Superflow& Superflow::operator<<(SaveVar save_var)
    {
        if (m_varState == SupervarState::open && m_superVar_hasFunction && m_superVar_hasNiceName) {
            if (!m_superVar_hasHFTName) {
                m_varHFTName.push_back("Untitled_" + to_string(m_superVar_Untitled));
                m_superVar_Untitled++;
                m_superVar_hasHFTName = true;
            }

            m_varState = SupervarState::closed;
            m_superVar_hasFunction = false;
            m_superVar_hasNiceName = false;
            m_superVar_hasHFTName = false;

            cout << app_name << "New var: " << m_varNiceName.back() << endl;
            cout << app_name << "    HFT: " << m_varHFTName.back() << endl;
        }
        else {
            if (m_varState != SupervarState::open) {
                cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
            }
            else {
                cout << app_name << "ERROR (Fatal): A lambda-expression is required.";
            }
            exit(1);
        }
        return *this;
    }

    ////////////////////////////////////////////
    // Systematic operators
    ////////////////////////////////////////////

    Superflow& Superflow::operator<<(NewSystematic new_sys) // this is the NiceName
    {
        if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
            m_sysTemplate.name = new_sys.name;
            m_sys_hasNiceName = true;
            m_sysState = SupersysState::open;
        }
        else {
            cout << app_name << "ERROR (Fatal): Close using SaveSystematic().";
            exit(1);
        }
        return *this;
    }

    Superflow& Superflow::operator<<(TreeName tree_name)
    {
        if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasTreeName) {
            m_sysTemplate.tree_name = tree_name.name;
            m_sys_hasTreeName = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
            exit(1);
        }
        return *this;
    }

    Superflow& Superflow::operator<<(EventSystematic obj_)
    {
        if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasSystematic) {
            m_sysTemplate.event_syst = obj_.event_syst_;
            m_sysTemplate.weight_syst = SupersysWeight::null;
            m_sys_hasSystematic = true;

            m_sysTemplate.type = SupersysType::event;
            m_sys_hasType = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
            exit(1);
        }
        return *this;
    }

    Superflow& Superflow::operator<<(WeightSystematic obj_)
    {
        if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasSystematic) {
            m_sysTemplate.event_syst = Susy::NtSys::NOM;
            m_sysTemplate.weight_syst_up = obj_.weight_syst_up;
            m_sysTemplate.weight_syst_down = obj_.weight_syst_down;
            m_sys_hasSystematic = true;

            m_sysTemplate.type = SupersysType::weight;
            m_sys_hasType = true;
        }
        else {
            cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
            exit(1);
        }
        return *this;
    }

    Superflow& Superflow::operator<<(SaveSystematic save_var)
    {
        if (m_sysState == SupersysState::open
            && m_varState == SupervarState::closed
            && m_sys_hasNiceName
            && m_sys_hasTreeName
            && m_sys_hasType
            && m_sys_hasSystematic
            ) {
            if (m_sysTemplate.name == "") m_sysTemplate.name = m_sysTemplate.tree_name;

            m_sysState = SupersysState::closed;
            m_sys_hasNiceName = false;
            m_sys_hasTreeName = false;
            m_sys_hasSystematic = false;
            m_sys_hasType = false;

            m_sysStore.push_back(m_sysTemplate);
            m_sysTemplate.reset();

            cout << app_name << "New systematic: " << m_sysStore.back().name << endl;
            if (m_sysStore.back().type == SupersysType::event) {
                cout << app_name << "    event systematic: " << m_sysStore.back().tree_name << endl;
            }
            else if (m_sysStore.back().type == SupersysType::weight) {
                cout << app_name << "    weight systematics: " << m_sysStore.back().tree_name << weight_suffix_up
                    << "/" << m_sysStore.back().tree_name << weight_suffix_down << endl;
            }
            else {
                cout << app_name << "ERROR (Fatal): Impossible SupersysType.";
                exit(1);
            }

        }
        else {
            if (m_sysState != SupersysState::open) {
                cout << app_name << "ERROR (Fatal): First open a NewSystematic().";
            }
            else {
                cout << app_name << "ERROR (Fatal): Can't save incomplete systematic object.";
            }
            exit(1);
        }
        return *this;
    }

    // Superlink
    // Superlink
    // Superlink

    void Superflow::attach_superlink(Superlink* sl_)
    {
        //sl_->tools = this;
        sl_->tools = &m_nttools; // from SusyNtAna

        sl_->anaType = m_nttools.getAnaType();

        sl_->nt = &nt; // SusyNt
        sl_->weights = m_weights;
        sl_->nt_sys = m_RunSyst->event_syst;

        sl_->preElectrons = &m_preElectrons;
        sl_->preMuons = &m_preMuons;
        sl_->preJets = &m_preJets;

        sl_->baseLeptons = &m_baseLeptons;
        sl_->baseElectrons = &m_baseElectrons;
        sl_->baseMuons = &m_baseMuons;
        sl_->baseTaus = &m_baseTaus;
        sl_->baseJets = &m_baseJets;

        if (m_runMode == SuperflowRunMode::fakes) {
            sl_->leptons = &m_baseLeptons;
            sl_->electrons = &m_baseElectrons;
            sl_->muons = &m_baseMuons;
        }
        else{
            sl_->leptons = &m_signalLeptons;
            sl_->electrons = &m_signalElectrons;
            sl_->muons = &m_signalMuons;
        }
     //   sl_->leptons = selectBaseLineLeptons ? &m_baseLeptons : &m_signalLeptons;
     //   sl_->electrons = selectBaseLineLeptons ? &m_baseElectrons : &m_signalElectrons;
     //   sl_->muons = selectBaseLineLeptons ? &m_baseMuons : &m_signalMuons;

        sl_->taus = &m_signalTaus;
        sl_->jets = &m_signalJets;

        sl_->met = m_met;

        sl_->dileptonTrigger = m_trigObj;
    # warning not setting jvftool
      //  sl_->jvfTool = m_jvfTool;
    
        sl_->fakeMatrix = m_matrix;

        if (nt.evt()->isMC) {
            sl_->isMC = true;
        }
        else {
            sl_->isData = true;
        }
    }

    // TSelector States
    // TSelector States
    // TSelector States

    void Superflow::Begin(TTree* /*tree*/)
    {
        cout << app_name << m_sample << endl;
        cout << app_name << m_sample << endl;
        cout << app_name << m_sample << endl;

        if (m_runMode == SuperflowRunMode::null) {
            cout << app_name << "ERROR (Fatal): Missing call to Superflow::setRunMode()." << endl;
            exit(1);
        }

        if (m_varState != SupervarState::closed || m_sysState != SupersysState::closed) {
            cout << app_name << "ERROR (Fatal): Close the Var using SaveVar()." << endl;;
            exit(1);
        }

        if (m_runMode == SuperflowRunMode::single_event_syst && m_singleEventSyst == Susy::NtSys::NOM) {
            cout << app_name << "ERROR (Fatal): SuperflowRunMode::single_event_syst: Call setSingleEventSyst(SusyNtSys nt_syst_)." << endl;
            exit(1);
        }

        // determine number of weight systematics
        if (m_runMode == SuperflowRunMode::nominal_and_weight_syst || m_runMode == SuperflowRunMode::all_syst) {
            for (int i = 0; i < m_sysStore.size(); i++) {
                if (m_sysStore[i].type == SupersysType::weight) {
                    index_weight_sys.push_back(i);
                }
            }
        }

        // determine number of fake-weight systematics
        if (m_runMode == SuperflowRunMode::fakes) {
            for (int i = 0; i < m_sysStore.size(); i++) {
                if (m_sysStore[i].type == SupersysType::weight) {
                    index_weight_sys.push_back(i);
                }
            }
        }

        // determine number of event systematics
        if (m_runMode == SuperflowRunMode::all_syst) {
            for (int i = 0; i < m_sysStore.size(); i++) {
                if (m_sysStore[i].type == SupersysType::event) {
                    index_event_sys.push_back(i);
                }
            }
        }

        // Set single event systematic
        if (m_runMode == SuperflowRunMode::single_event_syst) {
            for (int i = 0; i < m_sysStore.size(); i++) {
                if (m_sysStore[i].type == SupersysType::event && m_sysStore[i].event_syst == m_singleEventSyst) {
                    delete m_RunSyst;

                    m_RunSyst = &m_sysStore[i]; // don't delete!!
                }
            }
        }

        cout << app_name << "in Begin()" << endl;
        SusyNtAna::Begin(0);

        string period = "Moriond";
        bool useReweightUtils = false;
        m_trigObj = new DilTrigLogic(period, useReweightUtils);

        // -------------- Configure ChargeFlip Tool [BEGIN] ---------------------//
        // Currently set-up to run ChargeFlip-00-00-11 which has the 
        // charge flip map chargeFlip.root and the functionality as used
        // here
        // TODO: re-configure for Emma's  new map (--Dantrim 11/18/2014)
        if(m_do_qflip) {  // global switch set in SuperflowQFlip executable
            //string chargeFlipInput = "../../ChargeFlip/data/chargeFlip.root";
        //    string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
        //    string chargeFlipInput = gSystem->ExpandPathName("$ROOTCOREBIN/data/ChargeFlip/chargeflip_map_12nov2014.root");
            string chargeFlipInput = gSystem->ExpandPathName("$ROOTCOREBIN/data/ChargeFlip/chargeflip_map_12nov2014_scale_with_mc_last_ptbin.root");
        //    string chargeFlipInput = "../../ChargeFlip/data/chargeflip_map_12nov2014.root";
            ifstream qflipFile(chargeFlipInput.data());
            if(qflipFile) {
                cout << "\n >>> Using chargeflip map: " << chargeFlipInput << ". " << endl;
            }
            else {
                cout << "ERROR (Fatal): chargeflip map file (" << chargeFlipInput << ") does not exist or cannot be opened. Exitting. " << endl;
                exit(1);
            }
            m_chargeFlip = new chargeFlip(chargeFlipInput);
        }
        // ---------------- Configure ChargeFlip Tool [END] ---------------------//
        // ----------------- Configure Matrix Tool [BEGIN] ------------------------//
        // cf https://github.com/gerbaudo/DileptonMatrixMethod/
        if (m_runMode == SuperflowRunMode::fakes) {   
           // m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/FakeMatrix_Oct_20.root");
           // m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/FakeMatrix_Nov_26.root");
           // m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/pass3_Summer2013.root");
           // m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/FakeMatrix_Dec_03.root");
           // m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/FakeMatrix_Dec_03_with_syst.root");
            m_matrixFilename = gSystem->ExpandPathName("$ROOTCOREBIN/data/DileptonMatrixMethod/FakeMatrix_Nov_26_with_syst.root");

       //     m_matrixFilename = "../../DileptonMatrixMethod/data/FakeMatrix_Oct_20.root";
            ifstream the_matrix(m_matrixFilename.data());
            if(!the_matrix){
                cout << "ERROR (Fatal): matrix method matrix at " << m_matrixFilename << " does not exist or cannot be opened. Exitting. " << endl;
                exit(1);
            }
            m_allconfigured = initMatrixTool();
        }
        // ----------------- Configure Matrix Tool [END] ------------------------//
    }

    void Superflow::Init(TTree* tree)
    {
        cout << app_name << "in Init()" << endl;
        SusyNtAna::Init(tree);

        if (!nt.evt()->isMC && !(m_runMode == SuperflowRunMode::fakes)) {
            m_runMode = SuperflowRunMode::data; // Changing the run mode!!
        }
        else if (m_runMode != SuperflowRunMode::fakes) {
            initMcWeighter(tree);
        }

        // output file name
        stringstream sfile_name_;

        // determine output file name
        if (m_runMode == SuperflowRunMode::data) {
            m_countWeights = false;

            // output file name
            stringstream sfile_name_;
            sfile_name_ << "CENTRAL_";

    # warning need to set up running over data15!!
            size_t find_period = m_sample.find("period");
            if (find_period != string::npos) {
                char data_per_ = m_sample.substr(find_period + 6, 1)[0];

                if (data_per_ >= 'A') {
                    ATLAS_period this_period = static_cast<ATLAS_period>(data_per_ - 'A');
                    cout << app_name << "Determined data period: " << data_per_ << endl;

                    m_period = this_period; // set the period

                    size_t find_Egamma = m_sample.find("Egamma");
                    size_t find_Muons = m_sample.find("Muons");

                    if (find_Egamma != string::npos) {
                        cout << app_name << "Determined stream: Egamma" << endl;
                        m_stream = ATLAS_stream::Egamma;
                    }
                    else if (find_Muons != string::npos) {
                        cout << app_name << "Determined stream: Muons" << endl;
                        m_stream = ATLAS_stream::Muons;
                    }
                    else {
                        cout << app_name << "ERROR (Fatal): Unknown data stream! Supported streams: Egamma, Muons." << endl;
                        exit(1);
                    }
                }

                sfile_name_ << m_data_periods[m_stream][m_period] << ".root";
                cout << app_name << "Setting output file name to: " << sfile_name_.str() << endl;

                m_outputFileName = sfile_name_.str();
                m_entry_list_FileName += m_data_periods[m_stream][m_period] + ".root";
            }
            else {
                cout << app_name << "ERROR (fatal): Failed to determine data period. Try to specify the full path." << endl;
                exit(1);
            }
        }
        else if (m_runMode == SuperflowRunMode::fakes) {

            if (m_fake_region.compare("") == 0) {
                cout << app_name << "ERROR (Fatal): No fake region! Example: /f \'razor\' " << endl;
                exit(1);
            }

            m_countWeights = false;

            // output file name
            stringstream sfile_name_;
            sfile_name_ << "CENTRAL_";

            size_t find_period = m_sample.find("period");
            if (find_period != string::npos) {
                char data_per_ = m_sample.substr(find_period + 6, 1)[0];

                if (data_per_ >= 'A') {
                    ATLAS_period this_period = static_cast<ATLAS_period>(data_per_ - 'A');
                    cout << app_name << "Determined data period: " << data_per_ << endl;

                    m_period = this_period; // set the period

                    size_t find_Egamma = m_sample.find("Egamma");
                    size_t find_Muons = m_sample.find("Muons");

                    if (find_Egamma != string::npos) {
                        cout << app_name << "Determined stream: Egamma" << endl;
                        m_stream = ATLAS_stream::Egamma;
                    }
                    else if (find_Muons != string::npos) {
                        cout << app_name << "Determined stream: Muons" << endl;
                        m_stream = ATLAS_stream::Muons;
                    }
                    else {
                        cout << app_name << "ERROR (Fatal): Unknown data stream! Supported streams: Egamma, Muons." << endl;
                        exit(1);
                    }
                }


                sfile_name_ << "fakes." << m_fake_region << "." << m_data_periods[m_stream][m_period] << ".root";
                cout << app_name << "Setting output file name to: " << sfile_name_.str() << endl;

                m_outputFileName = sfile_name_.str();
                m_entry_list_FileName += string("fakes.") + m_fake_region + "." + m_data_periods[m_stream][m_period] + ".root";
            }
            else {
                cout << app_name << "ERROR (fatal): Failed to determine data period. Try to specify the full path." << endl;
                exit(1);
            }
        }
        else if (m_runMode == SuperflowRunMode::single_event_syst) {
            if (m_NtSys_to_string.count(m_singleEventSyst) != 0) {
                stringstream sfile_name_; // output file name
                sfile_name_ << m_NtSys_to_string[m_singleEventSyst] << "_" << nt.evt()->mcChannel << ".root";
                cout << app_name << "Run mode: SuperflowRunMode::single_event_syst" << endl;
                cout << app_name << "Setting output file name to: " << sfile_name_.str() << endl;
                m_outputFileName = sfile_name_.str();
                m_entry_list_FileName += sfile_name_.str() + ".root";
            }
            else {
                cout << app_name << "ERROR (Fatal): Unknown event systematic! Code: " << static_cast<int>(m_singleEventSyst) << endl;
                exit(1);
            }
        }
        else if (nt.evt()->isMC) {
            stringstream sfile_name_; // output file name
            sfile_name_ << "CENTRAL_" << nt.evt()->mcChannel << ".root";
            if (m_runMode == SuperflowRunMode::nominal_and_weight_syst) {
                cout << app_name << "Run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
            }
            else {
                cout << app_name << "Run mode: SuperflowRunMode::nominal (weighted)" << endl;
            }
            cout << app_name << "Setting output file name to: " << sfile_name_.str() << endl;
            m_outputFileName = sfile_name_.str();
            m_entry_list_FileName += to_string(nt.evt()->mcChannel) + ".root";
        }
        else {
            cout << app_name << "ERROR (Fatal): Inconsistent setup." << endl;
            exit(1);
        }

        m_outputFile = new TFile(m_outputFileName.data(), "RECREATE");

        // output tree name
        stringstream tree_name;
        if (m_runMode == SuperflowRunMode::data) {
            tree_name << "id_" << m_data_periods[m_stream][m_period];
        }
        else if (m_runMode == SuperflowRunMode::fakes) {
            tree_name << "id_fakes." << m_data_periods[m_stream][m_period];
        }
        else {
            tree_name << "id_" << nt.evt()->mcChannel;
        }

        cout << app_name << "Tree name: " << tree_name.str() << endl;
        m_tree_name_auto = tree_name.str();

        // initialize total entry list (also see Notify();)
        m_entryListFile = new TFile(m_entry_list_FileName.data(), "RECREATE");
        m_entry_list_total = new TEntryList(m_tree_name_auto.data(), m_tree_name_auto.data());
        m_entry_list_total->SetDirectory(m_entryListFile);

        // initialize output tree
        m_HFT = new TTree(tree_name.str().data(), tree_name.str().data());
        m_HFT->SetDirectory(m_outputFile);
        m_HFT->SetAutoFlush(-16777216L);

        // define number of trees
        m_tree_leafs_size = m_varType.size() + 2 * index_weight_sys.size();// 2nd term may be zero
        m_weight_leaf_offset = m_varType.size();

        m_varFloat = new Float_t[m_tree_leafs_size]; // this one is larger to hold the syst_WEIGHT
        m_varDouble = new Double_t[m_varType.size()];
        m_varInt = new Int_t[m_varType.size()];
        m_varBool = new Bool_t[m_varType.size()];

        // initialize HFT
        for (int i = 0; i < m_tree_leafs_size; i++) m_varFloat[i] = 1.0;
        for (int i = 0; i < m_varType.size(); i++) m_varDouble[i] = 1.0;
        for (int i = 0; i < m_varType.size(); i++) m_varInt[i] = 0;
        for (int i = 0; i < m_varType.size(); i++) m_varBool[i] = false;

        for (int i = 0; i < m_varType.size(); i++) {
            switch (m_varType[i]) {
                case SupervarType::sv_void: break;
                case SupervarType::sv_float: {
                    string leaflist_ = m_varHFTName[i] + "/F";
                    m_HFT->Branch(m_varHFTName[i].data(), m_varFloat + i, leaflist_.data(), 65536);
                    break;
                }
                case SupervarType::sv_double: {
                    string leaflist_ = m_varHFTName[i] + "/D";
                    m_HFT->Branch(m_varHFTName[i].data(), m_varDouble + i, leaflist_.data(), 65536);
                    break;
                }
                case SupervarType::sv_int: {
                    string leaflist_ = m_varHFTName[i] + "/I";
                    m_HFT->Branch(m_varHFTName[i].data(), m_varInt + i, leaflist_.data(), 65536);
                    break;
                }
                case SupervarType::sv_bool: {
                    string leaflist_ = m_varHFTName[i] + "/O";
                    m_HFT->Branch(m_varHFTName[i].data(), m_varBool + i, leaflist_.data(), 65536);
                    break;
                }
            }
        }

        for (int i = 0; i < index_weight_sys.size(); i++) {
            string syst_var_name_up = weight_prefix + m_sysStore[index_weight_sys[i]].tree_name + weight_suffix_up;
            string syst_var_name_down = weight_prefix + m_sysStore[index_weight_sys[i]].tree_name + weight_suffix_down;

            string leaflist_up = syst_var_name_up + "/F";
            string leaflist_down = syst_var_name_down + "/F";

            cout << app_name << "Weight var trees: " << syst_var_name_up << ", " << syst_var_name_down << endl;

            m_HFT->Branch(syst_var_name_up.data(), m_varFloat + m_weight_leaf_offset + 2 * i, leaflist_up.data(), 65536);
            m_HFT->Branch(syst_var_name_down.data(), m_varFloat + m_weight_leaf_offset + 2 * i + 1, leaflist_down.data(), 65536);
        }

        // Make an output file for each event systematic.
        if (m_runMode == SuperflowRunMode::all_syst) {

            m_output_array = new TFile*[index_event_sys.size()];

            for (int i = 0; i < index_event_sys.size(); i++) {
                stringstream sfile_name_; // output file name
                sfile_name_ << m_NtSys_to_string[m_sysStore[index_event_sys[i]].event_syst] << "_" << nt.evt()->mcChannel << ".root";
                m_output_array[i] = new TFile(sfile_name_.str().data(), "RECREATE");
            }

            m_HFT_array = new TTree*[index_event_sys.size()];
            for (int i = 0; i < index_event_sys.size(); i++) {
                m_HFT_array[i] = new TTree(tree_name.str().data(), tree_name.str().data());
                m_HFT_array[i]->SetDirectory(m_output_array[i]);
                m_HFT_array[i]->SetAutoFlush(-16777216L);
            }

            m_varFloat_array = new Float_t*[index_event_sys.size()];
            m_varDouble_array = new Double_t*[index_event_sys.size()];
            m_varInt_array = new Int_t*[index_event_sys.size()];
            m_varBool_array = new Bool_t*[index_event_sys.size()];

            for (int i = 0; i < index_event_sys.size(); i++) {
                m_varFloat_array[i] = new Float_t[m_varType.size()];
                m_varDouble_array[i] = new Double_t[m_varType.size()];
                m_varInt_array[i] = new Int_t[m_varType.size()];
                m_varBool_array[i] = new Bool_t[m_varType.size()];

                for (int j = 0; j < m_varType.size(); j++) m_varFloat_array[i][j] = 1.0;
                for (int j = 0; j < m_varType.size(); j++) m_varDouble_array[i][j] = 1.0;
                for (int j = 0; j < m_varType.size(); j++) m_varInt_array[i][j] = 0;
                for (int j = 0; j < m_varType.size(); j++) m_varBool_array[i][j] = false;
            }

            for (int i = 0; i < index_event_sys.size(); i++) {
                for (int j = 0; j < m_varType.size(); j++) {
                    switch (m_varType[j]) {
                        case SupervarType::sv_void: break;
                        case SupervarType::sv_float: {
                            string leaflist_ = m_varHFTName[j] + "/F";
                            m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varFloat_array[i] + j, leaflist_.data(), 65536);
                            break;
                        }
                        case SupervarType::sv_double: {
                            string leaflist_ = m_varHFTName[j] + "/D";
                            m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varDouble_array[i] + j, leaflist_.data(), 65536);
                            break;
                        }
                        case SupervarType::sv_int: {
                            string leaflist_ = m_varHFTName[j] + "/I";
                            m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varInt_array[i] + j, leaflist_.data(), 65536);
                            break;
                        }
                        case SupervarType::sv_bool: {
                            string leaflist_ = m_varHFTName[j] + "/O";
                            m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varBool_array[i] + j, leaflist_.data(), 65536);
                            break;
                        }
                    }
                }
            }
        }
    }

    Bool_t Superflow::Notify()
    {
        static int tree_counter;

        if (m_entry_list_single_tree != nullptr) m_entry_list_total->Add(m_entry_list_single_tree);
        delete m_entry_list_single_tree;

        string new_list_name = m_tree_name_auto + "_" + to_string(tree_counter);
        tree_counter++;

        m_entry_list_single_tree = new TEntryList();
        m_entry_list_single_tree->SetTree(m_input_chain->GetTree());

        return kTRUE;
    }

    Bool_t Superflow::Process(Long64_t entry)
    {
        GetEntry(entry);
    #warning what is "save_entry_to_list"???
        bool save_entry_to_list = true;

        m_chainEntry++; // SusyNtAna counter

        if (m_chainEntry % 500 == 0) {
            cout << app_name << "**** Processing entry " << setw(6) << m_chainEntry
                << " run " << setw(6) << nt.evt()->run
                << " event " << setw(7) << nt.evt()->event << " ****" << endl;
        }

        // select baseline and signal objects
        bool removeLepsFromIso = false;

        // these are flags
        var_float* vf_ = nullptr;
        var_double* vd_ = nullptr;
        var_int* vi_ = nullptr;
        var_bool* vb_ = nullptr;
        var_void* vv_ = nullptr;

        switch (m_runMode) {
            case SuperflowRunMode::data: {
                clearObjects();
                selectObjects(m_RunSyst->event_syst, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)

                m_weights = new Superweight();
                Superlink* sl_ = new Superlink;
                attach_superlink(sl_);

                bool pass_cuts = true; // loop over and appply the cuts in m_CutStore.
                if (m_CutStore.size() > 0) {
                    for (int i = 0; i < m_CutStore.size(); i++) {
                        pass_cuts = m_CutStore[i](sl_); // run the cut function

                        if (pass_cuts) {
                            m_RawCounter[i]++;
                        }
                        else {
                            break;
                        }
                    }
                }

                if (pass_cuts) { // data passed cuts, so fill HFTs.
                    if (save_entry_to_list) {
                        m_entry_list_single_tree->Enter(entry);
                        save_entry_to_list = false;
                    }
                    for (int v_ = 0; v_ < m_varType.size(); v_++) {
                        switch (m_varType[v_]) {
                            case SupervarType::sv_float: {
                                m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_); break;
                            }
                            case SupervarType::sv_double: {
                                m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_); break;
                            }
                            case SupervarType::sv_int: {
                                m_varInt[v_] = m_varExprInt[v_](sl_, vi_); break;
                            }
                            case SupervarType::sv_bool: {
                                m_varBool[v_] = m_varExprBool[v_](sl_, vb_); break;
                            }
                            case SupervarType::sv_void: {
                                m_varExprVoid[v_](sl_, vv_); break;
                            }
                        }
                    }
                    m_HFT->Fill();
                }
                delete sl_;
                delete m_weights;

            } break;
            case SuperflowRunMode::nominal:
            case SuperflowRunMode::single_event_syst: {

                clearObjects();
                selectObjects(m_RunSyst->event_syst, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)

                m_weights = new Superweight();
                Superlink* sl_ = new Superlink;
                attach_superlink(sl_);

                bool pass_cuts = true;

                if (m_countWeights) {
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
                }

                if (m_CutStore.size() > 0) {
                    for (int i = 0; i < m_CutStore.size(); i++) {
                        pass_cuts = m_CutStore[i](sl_); // run the cut function
                        if (pass_cuts) {
                            m_RawCounter[i]++;
                            if (m_countWeights) m_WeightCounter[i] += m_weights->product();
                        }
                        else {
                            break;
                        }
                    }
                }

                if (pass_cuts) {
                    if (save_entry_to_list) {
                        m_entry_list_single_tree->Enter(entry);
                        save_entry_to_list = false;
                    }
                    if (!m_countWeights) {
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
                        m_WeightCounter[m_CutStore.size() - 1] += m_weights->product();
                    }

                    // FILL_HFTs
                    for (int v_ = 0; v_ < m_varType.size(); v_++) {
                        switch (m_varType[v_]) {
                            case SupervarType::sv_float: {
                                m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_); break;
                            }
                            case SupervarType::sv_double: {
                                m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_); break;
                            }
                            case SupervarType::sv_int: {
                                m_varInt[v_] = m_varExprInt[v_](sl_, vi_); break;
                            }
                            case SupervarType::sv_bool: {
                                m_varBool[v_] = m_varExprBool[v_](sl_, vb_); break;
                            }
                            case SupervarType::sv_void: {
                                m_varExprVoid[v_](sl_, vv_); break;
                            }
                        }
                    }
                    m_HFT->Fill();
                }
                delete sl_;
                delete m_weights;

            } break;
            case SuperflowRunMode::all_syst: {
                for (int i = 0; i < index_event_sys.size(); i++) { // loop over event systematics
                    clearObjects();
                    delete m_RunSyst;

                    m_RunSyst = &m_sysStore[index_event_sys[i]]; // don't delete!!
                    selectObjects(m_RunSyst->event_syst, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)

                    m_weights = new Superweight();
                    Superlink* sl_ = new Superlink;
                    attach_superlink(sl_);

                    bool pass_cuts = true;

                    if (m_CutStore.size() > 0) {
                        for (int i = 0; i < m_CutStore.size(); i++) {
                            pass_cuts = m_CutStore[i](sl_); // run the cut function
                            if (!pass_cuts) break;
                        }
                    }

                    if (pass_cuts) {
                        if (save_entry_to_list) {
                            m_entry_list_single_tree->Enter(entry);
                            save_entry_to_list = false;
                        }
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
                        // FILL_HFTs
                        for (int v_ = 0; v_ < m_varType.size(); v_++) {
                            switch (m_varType[v_]) {
                                case SupervarType::sv_float: {
                                    m_varFloat_array[i][v_] = m_varExprFloat[v_](sl_, vf_); break;
                                }
                                case SupervarType::sv_double: {
                                    m_varDouble_array[i][v_] = m_varExprDouble[v_](sl_, vd_); break;
                                }
                                case SupervarType::sv_int: {
                                    m_varInt_array[i][v_] = m_varExprInt[v_](sl_, vi_); break;
                                }
                                case SupervarType::sv_bool: {
                                    m_varBool_array[i][v_] = m_varExprBool[v_](sl_, vb_); break;
                                }
                                case SupervarType::sv_void: {
                                    m_varExprVoid[v_](sl_, vv_); break;
                                }
                            }
                        }
                        m_HFT_array[i]->Fill();
                    }
                    delete sl_;
                    delete m_weights;

                    m_RunSyst = nullptr;
                }
            } // we didn't break!!
            case SuperflowRunMode::nominal_and_weight_syst: {
                delete m_RunSyst;
                m_RunSyst = new Supersys(SupersysType::central);

                clearObjects();
                selectObjects(m_RunSyst->event_syst, removeLepsFromIso, TauID_medium);

                m_weights = new Superweight();
                Superlink* sl_ = new Superlink;
                attach_superlink(sl_);

                bool pass_cuts = true;

                if (m_CutStore.size() > 0) {
                    for (int i = 0; i < m_CutStore.size(); i++) {
                        pass_cuts = m_CutStore[i](sl_); // run the cut function

                        if (pass_cuts) {
                            m_RawCounter[i]++;
                        }
                        else {
                            break;
                        }
                    }
                }

                if (pass_cuts) {
                    if (save_entry_to_list) {
                        m_entry_list_single_tree->Enter(entry);
                        save_entry_to_list = false;
                    }
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);

                    double nom_eventweight = m_weights->product();
                    m_WeightCounter[m_CutStore.size() - 1] += m_weights->product();

                    // FILL HFTs
                    for (int v_ = 0; v_ < m_varType.size(); v_++) {
                        switch (m_varType[v_]) {
                            case SupervarType::sv_float: {
                                m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_); break;
                            }
                            case SupervarType::sv_double: {
                                m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_); break;
                            }
                            case SupervarType::sv_int: {
                                m_varInt[v_] = m_varExprInt[v_](sl_, vi_); break;
                            }
                            case SupervarType::sv_bool: {
                                m_varBool[v_] = m_varExprBool[v_](sl_, vb_); break;
                            }
                            case SupervarType::sv_void: {
                                m_varExprVoid[v_](sl_, vv_); break;
                            }
                        }
                    }

                    // FILL more HFTs
                    for (int w_ = 0; w_ < index_weight_sys.size(); w_++) {
                        Superweight* weightComponents_copy = new Superweight(*m_weights);

                        // Up variation
                        m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_up; // do up variation
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                        double up_weight = weightComponents_copy->product();
                        delete weightComponents_copy;

                        weightComponents_copy = new Superweight(*m_weights);
                        // Down variation
                        m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_down; // do down variation
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                        double down_weight = weightComponents_copy->product();

                        delete weightComponents_copy;
                        weightComponents_copy = nullptr;

                        if (up_weight < nom_eventweight && down_weight > nom_eventweight) {
                            double temp = up_weight; // reason to swap these
                            up_weight = down_weight;
                            down_weight = temp;
                        }
                        if (nom_eventweight > m_epsilon) {
                            up_weight = up_weight / nom_eventweight;
                            down_weight = down_weight / nom_eventweight;
                        }
                        else {
                            up_weight = 1.0;
                            down_weight = 1.0;
                        }

                        *(m_varFloat + m_weight_leaf_offset + 2 * w_) = up_weight; // put in TBranch
                        *(m_varFloat + m_weight_leaf_offset + 2 * w_ + 1) = down_weight;
                    }
                    m_RunSyst->weight_syst = SupersysWeight::null; // must reset this value!
                    m_HFT->Fill();
                }
                delete sl_;
                delete m_weights;
            } break;
            case SuperflowRunMode::fakes: {
                delete m_RunSyst;
                m_RunSyst = new Supersys(SupersysType::central);

                clearObjects();
                selectObjects(m_RunSyst->event_syst, removeLepsFromIso, TauID_medium);

                m_weights = new Superweight();
                Superlink* sl_ = new Superlink;
                attach_superlink(sl_);

                bool pass_cuts = true;

                if (m_CutStore.size() > 0) {
                    for (int i = 0; i < m_CutStore.size(); i++) {
                        pass_cuts = m_CutStore[i](sl_); // run the cut function

                        if (pass_cuts) {
                            m_RawCounter[i]++;
                        }
                        else {
                            break;
                        }
                    }
                }

                if (pass_cuts) {
                    if (save_entry_to_list) {
                        m_entry_list_single_tree->Enter(entry);
                        save_entry_to_list = false;
                    }
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);

                    double nom_eventweight = m_weights->product();
                    m_WeightCounter[m_CutStore.size() - 1] += m_weights->product();

                    // FILL HFTs
                    for (int v_ = 0; v_ < m_varType.size(); v_++) {
                        switch (m_varType[v_]) {
                            case SupervarType::sv_float: {
                                m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_); break;
                            }
                            case SupervarType::sv_double: {
                                m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_); break;
                            }
                            case SupervarType::sv_int: {
                                m_varInt[v_] = m_varExprInt[v_](sl_, vi_); break;
                            }
                            case SupervarType::sv_bool: {
                                m_varBool[v_] = m_varExprBool[v_](sl_, vb_); break;
                            }
                            case SupervarType::sv_void: {
                                m_varExprVoid[v_](sl_, vv_); break;
                            }
                        }
                    }

                    // FILL more HFTs
                    for (int w_ = 0; w_ < index_weight_sys.size(); w_++) {
                        Superweight* weightComponents_copy = new Superweight(*m_weights);

                        // Up variation
                        m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_up; // do up variation
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                        double up_weight = weightComponents_copy->product();
                        delete weightComponents_copy;

                        weightComponents_copy = new Superweight(*m_weights);
                        // Down variation
                        m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_down; // do down variation
                        computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                        double down_weight = weightComponents_copy->product();

                        delete weightComponents_copy;
                        weightComponents_copy = nullptr;

                        if (up_weight < nom_eventweight && down_weight > nom_eventweight) {
                            double temp = up_weight; // reason to swap these
                            up_weight = down_weight;
                            down_weight = temp;
                        }
                        if (nom_eventweight > m_epsilon) {
                            up_weight = up_weight / nom_eventweight;
                            down_weight = down_weight / nom_eventweight;
                        }
                        else {
                            up_weight = 1.0;
                            down_weight = 1.0;
                        }

                        *(m_varFloat + m_weight_leaf_offset + 2 * w_) = up_weight; // put in TBranch
                        *(m_varFloat + m_weight_leaf_offset + 2 * w_ + 1) = down_weight;
                    }
                    m_RunSyst->weight_syst = SupersysWeight::null; // must reset this value!
                    m_HFT->Fill();
                }
                delete sl_;
                delete m_weights;
            } break;
            default: break;
        }
        return kTRUE;
    }

    void Superflow::Terminate()
    {
        cout << app_name << "in Terminate()" << endl;

        cout << app_name << "Raw" << endl;
        cout << app_name << "Raw" << endl;
        cout << app_name << "Raw" << endl;
        cout << std::fixed;
        cout << std::setprecision(0);
        for (int i = 0; i < m_CutStore.size(); i++) {
            cout << app_name << "Cut " << pad_width(to_string(i), 2) << ": " << pad_width(m_CutStoreNames[i], 32) << ": " << m_RawCounter[i] << endl;
        }
        cout << app_name << endl << app_name << endl;

        cout << app_name << "Weighted" << endl;
        cout << app_name << "Weighted" << endl;
        cout << app_name << "Weighted" << endl;
        cout << std::resetiosflags(std::ios::floatfield);
        cout << std::resetiosflags(std::ios::adjustfield);
        cout << std::setprecision(6);
        for (int i = 0; i < m_CutStore.size(); i++) {
            cout << app_name << "Cut " << pad_width(to_string(i), 2) << ": " << pad_width(m_CutStoreNames[i], 32) << ": " << m_WeightCounter[i] << endl;
        }
        cout << app_name << endl << app_name << endl;

        m_outputFile->Write();
        for (int i = 0; i < index_event_sys.size(); i++) m_output_array[i]->Write();

        m_outputFile->Close();
        for (int i = 0; i < index_event_sys.size(); i++) m_output_array[i]->Close();

        cout << app_name << "Files OK." << endl;

        if (m_entry_list_single_tree != nullptr) {
            m_entry_list_total->Add(m_entry_list_single_tree); // last tree
        }

        m_entryListFile->Write();
        m_entryListFile->Close();

        delete m_outputFile;
        delete m_entryListFile;
        for (int i = 0; i < index_event_sys.size(); i++) delete m_output_array[i];
        delete[] m_output_array;

        if (m_entry_list_single_tree != nullptr) delete m_entry_list_single_tree;

        // Trees and entry-lists in the files are deleted with the file.

        SusyNtAna::Terminate();

        if (m_mcWeighter) delete m_mcWeighter;
        if (m_trigObj) delete m_trigObj;

        delete[] m_varFloat;
        delete[] m_varDouble;
        delete[] m_varInt;
        delete[] m_varBool;

        if (m_runMode == SuperflowRunMode::all_syst) {
            for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varFloat_array[i];
            delete m_varFloat_array;
            for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varDouble_array[i];
            delete m_varDouble_array;
            for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varInt_array[i];
            delete m_varInt_array;
            for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varBool_array[i];
            delete m_varBool_array;
        }

        if (m_runMode != SuperflowRunMode::single_event_syst) delete m_RunSyst;

        cout << app_name << "Done." << endl;
    }



    bool Superflow::initMcWeighter(TTree *tree)
    {
        bool success = false;
        if (tree) {
            string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
            m_mcWeighter = new MCWeighter(tree, xsecDir);

            bool isPmssmSample = false;
            if (sampleName().find("Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10") != string::npos) {
                isPmssmSample = true;
            }

         //   m_mcWeighter->parseAdditionalXsecFile("${ROOTCOREBIN}/data/Superflow/LFV.txt", /*m_dbg*/ false);

            if (isPmssmSample) {
                m_mcWeighter->setLabelBinCounter("Initial").clearAndRebuildSumwMap(m_tree);
            }
            if (m_dbg) {
                cout << app_name << "MCWeighter has been initialized." << endl;
            }
        }
        else {
            cout << app_name << "ERROR: Invalid input tree, cannot initialize MCWeighter." << endl;
        }
        return success;
    }

    // ----------------------- Initialize Matrix Tool [BEGIN] ------------------------ //
    bool Superflow::initMatrixTool()           
    {
        std::cout << "\n    ----- Initializating DileptonMatrixMethod    -----" << std::endl;
        m_matrix = new susy::fake::DileptonMatrixMethod();
        if(m_use2dparametrization){
            std::cout << " ----- Using 2D-Parametrized Fake Weight Calculation ----- \n" << std::endl;
        }
        susy::fake::Parametrization::Value p = (m_use2dparametrization ? susy::fake::Parametrization::PT_ETA : susy::fake::Parametrization::PT);
        std::vector<std::string> regions;
        regions.push_back(m_fake_region);            
        for(int i=0; i<regions.size(); ++i){
            std::cout<< ">> Available fake-region:  " << regions[i] << std::endl;
        }
        return m_matrix->configure(m_matrixFilename, regions, p, p, p, p);
    }
    // ----------------------- Initialize Matrix Tool [END] --------------------------- //

    bool Superflow::computeWeights(
        Susy::SusyNtObject &ntobj,
        MCWeighter &weighter,
        const LeptonVector& leptons,
        const JetVector& jets,
        Supersys* super_sys,
        Superweight* weights_
        )
    {

        if (m_runMode == SuperflowRunMode::fakes) {
            susy::fake::Systematic::Value fake_sys = susy::fake::Systematic::Value::SYS_NOM;
            bool do_fake_ = false;         // FAKE
            bool do_syst_ = false;
        
           switch (super_sys->weight_syst) {
               case SupersysWeight::ELFRUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_FR_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::ELFRDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_FR_DOWN;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::ELREUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_RE_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::ELREDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_RE_DOWN;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUFRUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_FR_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUFRDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_FR_DOWN;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUREUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_RE_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUREDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_RE_DOWN;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::ELFRACUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_FRAC_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::ELFRACDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_EL_FRAC_DO;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUFRACUP: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_FRAC_UP;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::MUFRACDOWN: {
                   fake_sys = susy::fake::Systematic::Value::SYS_MU_FRAC_DO;
                   do_fake_ = true;
                   do_syst_ = true;
               }   break;
               case SupersysWeight::null: {
                   do_fake_ = true;
                   do_syst_ = false;
               } break;
               default: break;
           }
                    if (do_fake_) {

                        Superlink* sl_ = new Superlink();
                        attach_superlink(sl_);

                        weights_->fake = PhysicsTools::getFakeWeight(sl_, m_fake_region, fake_sys, do_syst_);
                    }
         //   } // end if runMode fakes
        
        }

        // MCWeighter's susynt-weight calculation
        if (ntobj.evt()->isMC) {
            MCWeighter::WeightSys wSys = MCWeighter::Sys_NOM;


            bool do_susynt_w = false;

            switch (super_sys->weight_syst) {
                case SupersysWeight::XSUP: {
                    wSys = MCWeighter::Sys_XSEC_UP;
                    do_susynt_w = true;
                    break;
                }
                case SupersysWeight::XSDOWN: {
                    wSys = MCWeighter::Sys_XSEC_DN;
                    do_susynt_w = true;
                    break;
                }
                case SupersysWeight::PILEUPUP: {
                    wSys = MCWeighter::Sys_PILEUP_UP;
                    do_susynt_w = true;
                    break;
                }
                case SupersysWeight::PILEUPDOWN: {
                    wSys = MCWeighter::Sys_PILEUP_DN;
                    do_susynt_w = true;
                    break;
                }
                case SupersysWeight::null: {
                    do_susynt_w = true;
                }
                default: break;
            }
            if (do_susynt_w) {
                # warning setting susynt weight to 1 !!!!
                # warning setting susynt weight to 1 !!!!
                weights_->susynt = 1.0;
                //weights_->susynt = weighter.getMCWeight(ntobj.evt(), 10000, wSys);
            }

            // Other weight systematic variations
            if (leptons.size() > 1) {
                // vars.hasFiredTrig = m_trigObj->passDilEvtTrig(leptons, m_met->Et, nt.evt());
                // vars.hasTrigMatch = m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());
                const Lepton &l0 = *(leptons[0]);
                const Lepton &l1 = *(leptons[1]);

                bool do_lepSf_ = false; // do_lepSf_

                switch (super_sys->weight_syst) {
                    case SupersysWeight::ESFUP:
                    case SupersysWeight::ESFDOWN:
                    case SupersysWeight::MEFFUP:
                    case SupersysWeight::MEFFDOWN:
                    case SupersysWeight::null: {
                        do_lepSf_ = true;
                    } break;
                    default: break;
                }
                if (do_lepSf_) {
                    weights_->lepSf = (computeLeptonEfficiencySf(l0, super_sys->weight_syst) * computeLeptonEfficiencySf(l1, super_sys->weight_syst));
                }

                bool do_lep_triggers_ = false; // do_lep_triggers_
                SusyNtSys trig_sys = Susy::NtSys::NOM;

                switch (super_sys->weight_syst) {
                    case SupersysWeight::ETRIGREWUP: {
                        trig_sys = Susy::NtSys::TRIGSF_EL_UP;
                        do_lep_triggers_ = true;
                    } break;
                    case SupersysWeight::ETRIGREWDOWN: {
                        trig_sys = Susy::NtSys::TRIGSF_EL_DN;
                        do_lep_triggers_ = true;
                    } break;
                    case SupersysWeight::MTRIGREWUP: {
                        trig_sys = Susy::NtSys::TRIGSF_MU_UP;
                        do_lep_triggers_ = true;
                    } break;
                    case SupersysWeight::MTRIGREWDOWN: {
                        trig_sys = Susy::NtSys::TRIGSF_MU_DN;
                        do_lep_triggers_ = true;
                    } break;
                    case SupersysWeight::null: {
                        do_lep_triggers_ = true;
                    } break;
                    default: break;
                }
                if (do_lep_triggers_) {
                    weights_->trigger = computeDileptonTriggerWeight(leptons, trig_sys);
                }

                bool do_btag_ = false; // do_btag_

                switch (super_sys->weight_syst) {
                    case SupersysWeight::BJETUP:
                    case SupersysWeight::BJETDOWN:
                    case SupersysWeight::CJETUP:
                    case SupersysWeight::CJETDOWN:
                    case SupersysWeight::BMISTAGUP:
                    case SupersysWeight::BMISTAGDOWN:
                    case SupersysWeight::null: {
                        do_btag_ = true;
                    } break;
                    default: break;
                }
                if (do_btag_) {
                    weights_->btag = 1.0;
                    #warning fixing btag weight to 1.0!
                    //weights_->btag = computeBtagWeight(jets, nt.evt(), super_sys->weight_syst);
                }

                if(m_do_qflip){ // global switch set in SuperflowQFlip executable
                    bool do_qflip_ = false;
            
                    switch (super_sys->weight_syst) {
                        case SupersysWeight::BKGMETHODUP:
                        case SupersysWeight::BKGMETHODDOWN:
                        case SupersysWeight::null: {
                            do_qflip_ = true;
                        } break;
                        default: break;
                    } // switch
                    bool  isMuMu(isMM(leptons));
                    bool  isGenSS(isGenuineSS(leptons));
                    if(do_qflip_ && !isMuMu && !isGenSS){
                        weights_->qflip = computeChargeFlipWeight(leptons, super_sys->weight_syst);
                    }
                } // end if(m_do_qflip)

                // ISR uncertainty
                bool do_isr_ = false;
                switch (super_sys->weight_syst) {
                    case SupersysWeight :: ISRUP : 
                    case SupersysWeight :: ISRDOWN :
                    case SupersysWeight :: null : {
                        do_isr_ = true;
                    } break;
                    default: break;
                } // switch
                if(do_isr_){
                    Superlink* sl_ = new Superlink;
                    attach_superlink(sl_);
                    weights_->isr = get_isr_weights(sl_, super_sys->weight_syst);
                    delete sl_;
                }
            }
        }

        return true;
    }


    double Superflow::computeDileptonTriggerWeight(const LeptonVector &leptons, const SusyNtSys sys)
    {
        double trigW = 1.0;
        if (leptons.size() == 2) {

            trigW = m_trigObj->getTriggerWeight(leptons, nt.evt()->isMC, m_met->Et, m_signalJets.size(), nt.evt()->nVtx, sys);
            bool twIsInvalid = !(trigW >= 0) || trigW < 0.0;
            trigW = twIsInvalid ? 0.0 : trigW;
        }
        return trigW;
    }

    double Superflow::computeBtagWeight(const JetVector& jets, const Susy::Event* evt, SupersysWeight sys)
    {
        // cout << app_name << "in computeBtagWeight" << endl;
        //JetVector taggableJets = SusyNtTools::getBTagSFJets2Lep(jets);
        //return SusyNtTools::bTagSF(evt, taggableJets, evt->mcChannel, supersys_to_btagsys(sys));
        Superlink* sl_ = new Superlink();
        attach_superlink(sl_);
        JetVector taggableJets = sl_->tools->getBTagSFJets2Lep(jets);
        double btagsf = sl_->tools->bTagSF(evt, taggableJets, evt->mcChannel, supersys_to_btagsys(sys));
        delete sl_;
        return btagsf;
    }

    double Superflow::computeLeptonEfficiencySf(const Susy::Lepton &lep, const SupersysWeight sys)
    {
        double effFactor = 1.0;
        double sf = lep.effSF;
        double delta = 0.0;

        if (lep.isEle()) {
            if (sys == SupersysWeight::ESFUP) {
                delta = (+lep.errEffSF);
            }
            else if (sys == SupersysWeight::ESFDOWN) {
                delta = (-lep.errEffSF);
            }
        }
        else if (lep.isMu()) {
            if (sys == SupersysWeight::MEFFUP) {
                delta = (+lep.errEffSF);
            }
            else if (sys == SupersysWeight::MEFFDOWN) {
                delta = (-lep.errEffSF);
            }
        }

        effFactor = (sf + delta); // ?? Seems odd.
        return effFactor;
    }

    float Superflow::computeChargeFlipWeight(const LeptonVector &leptons, const SupersysWeight sys)
    {
        float qflipProb = 1.0;
    
     //   if(!isGenSS && isMuMu ) return 0.0;
     //   if(!isGenSS && !isMuMu){
            const LeptonVector &ls = leptons;
            Susy::Lepton* l0 = ls[0];
            Susy::Lepton* l1 = ls[1];
            int pdg0 = l0->isEle() ? 11 : 13;
            int pdg1 = l1->isEle() ? 11 : 13;
            TVector2 met(m_met->lv().Px(), m_met->lv().Py());
            int _sys = (sys==SupersysWeight::BKGMETHODUP ? +1 :
                        (sys==SupersysWeight::BKGMETHODDOWN ? -1 : 0)); // convert to convention used in ChargeFlip (sys==0 --> nominal)
            bool isData=false; // this is used only in new verisions of ChargeFlip package
            m_chargeFlip->setSeed(nt.evt()->event);
            float flipProb(m_chargeFlip->OS2SS(pdg0, l0, pdg1, l1, _sys, isData, chargeFlip::dataonly)); // newer version
            //float flipProb(m_chargeFlip->OS2SS(pdg0, l0, pdg1, l1, &met, syst));
            float overlapFrac(m_chargeFlip->overlapFrac().first);
            
            qflipProb = flipProb*overlapFrac;
//            if(sys==SupersysWeight::BKGMETHODUP) qflipProb *= 1.5;
//            if(sys==SupersysWeight::BKGMETHODDOWN) qflipProb *= 0.5;

//            return qflipProb;
            // check if final qflipProb > 0.0 && < 1.0 before returning?
//        }

        return qflipProb;
//        else{
//            return 1.0;
//        }
    }
    bool Superflow::isGenuineSS(const LeptonVector& leps)
    {
        bool lessThanTwo        = (leps.size()<2) ? true : false;
        bool qFlipPresent       = hasQFlip(leps)  ? true : false;
        float qq                = leps.at(0)->q * leps.at(1)->q;
        bool ntOS               = (qq<0) ? true : false;
        if(!lessThanTwo && !qFlipPresent && !ntOS ) return true;
        else{
            return false;
        }
    }
    bool Superflow::hasQFlip(const LeptonVector& leps)
    {
        bool lessThanTwo        = (leps.size()<2) ? true : false;
        const Susy::Lepton* l0  = leps.at(0);
        const Susy::Lepton* l1  = leps.at(1);
        bool l0_isQFlip         = l0->isEle() ? (static_cast<const Susy::Electron*>(l0))->isChargeFlip : false;
        bool l1_isQFlip         = l1->isEle() ? (static_cast<const Susy::Electron*>(l1))->isChargeFlip : false;
        if(!lessThanTwo && (l0_isQFlip || l1_isQFlip) ) {
            return true;
        }
        else{
            return false;
        }
    }
    bool Superflow::isMM(const LeptonVector& leps)
    {
        return (leps.at(0)->isMu() && leps.at(1)->isMu());
    }

    double Superflow::get_isr_weights(Superlink* sl, const SupersysWeight sys)
    {
        double isr_weight = 1.;
        // grab the weights file 
        string weights_filename = gSystem->ExpandPathName("$ROOTCOREBIN/data/Superflow/FinalPtC1C1Weights.root");
        TFile* weights_file = new TFile(weights_filename.c_str(), "READ");
        TGraphAsymmErrors* weights = (TGraphAsymmErrors*) weights_file->Get("Weight");
        double* EXhigh = weights->GetEXhigh();
        double* EXlow  = weights->GetEXlow();
        double* EYhigh = weights->GetEYhigh();
        double* EYlow  = weights->GetEYlow();
        double c1c1Pt = 0.0;
        std::vector<Susy::TruthParticle> c1c1Parts;
        for(int itp=0; itp < sl->nt->tpr()->size(); itp++){
            if(fabs(sl->nt->tpr()->at(itp).pdgId)==1000024){
                c1c1Parts.push_back(sl->nt->tpr()->at(itp));
            } // if itp
        } // tpr loop
        if(c1c1Parts.size()==2){
            c1c1Pt = (c1c1Parts[0] + c1c1Parts[1]).Pt();
            for(int i=0; i < weights->GetN(); ++i){
                double x=0.;
                double y=0.;
                weights->GetPoint(i, x, y);
                if( (c1c1Pt > (x-EXlow[i])) && 
                    ( c1c1Pt < (x+EXhigh[i])) ) {
                    if     (sys==SupersysWeight::null)    { isr_weight = y; }
                    else if(sys==SupersysWeight::ISRUP)   { isr_weight = y + EYhigh[i]; }
                    else if(sys==SupersysWeight::ISRDOWN) { isr_weight = y - EYlow[i]; }
                 //   if     (sys==SupersysWeight::null)    { isr_weight = y; cout << "ptC1C1: " << c1c1Pt << "  nom: " << isr_weight << endl; }
                 //   else if(sys==SupersysWeight::ISRUP)   { isr_weight = y + EYhigh[i]; cout << "    up: " << isr_weight << endl; }
                 //   else if(sys==SupersysWeight::ISRDOWN) { isr_weight = y - EYlow[i]; cout << "     down: " << isr_weight << endl; }
                } // in pt range
            } // i
        } // if c1c1 found
        weights_file->Close();
        delete weights_file;
        delete weights;
        return isr_weight;
    }
            

    void Superflow::setCountWeights(bool value) ///> public function, if set true it prints the weighted cutflow
    {
        m_countWeights = value;
    }

    void Superflow::setRunMode(SuperflowRunMode run_mode_) ///> public function
    {
        m_runMode = run_mode_;
    }

    void Superflow::setSingleEventSyst(SusyNtSys nt_syst_)
    {
        m_singleEventSyst = nt_syst_;
    }

    void Superflow::setChain(TChain* input_chain_)
    {
        m_input_chain = input_chain_;
    }

    void Superflow::setFakeRegion(string fk_reg)
    {
        m_fake_region = fk_reg;
        selectBaseLineLeptons = true;
    }
    void Superflow::setFake2dParam(bool use2d)
    {
        m_use2dparametrization = use2d;
    }
    void Superflow::setQFlip()
    {
        m_do_qflip = true;
    }
}

/* snippets

//debug// cout << app_name << endl;
//debug// cout << app_name << m_NtSys_to_string[m_RunSyst->event_syst] << endl;
//debug// cout << app_name << m_NtSys_to_string[m_RunSyst->event_syst] << endl;
//debug// cout << app_name << m_NtSys_to_string[m_RunSyst->event_syst] << endl;
//debug// cout << app_name << endl;

// cout << app_name << "Weight variation: " << m_sysStore[index_weight_sys[w_]].tree_name << endl;
// cout << "    nom: " << nom_eventweight << endl;
// cout << "    up : " << up_weight << endl;
// cout << "    dwn: " << down_weight << endl;
// cout << app_name << endl << app_name << endl;


*/
