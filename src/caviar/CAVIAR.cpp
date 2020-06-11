
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// This file is part of the CAVIAR package.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#include "caviar/CAVIAR.h"

#include <string>

#include "caviar/interpreter/all.h"

namespace caviar {

#if defined(CAVIAR_WITH_MPI)
CAVIAR::CAVIAR (int argc, char **argv, MPI_Comm mpi_comm) :
  mpi_comm {mpi_comm},
#else
CAVIAR::CAVIAR (int argc, char **argv) :
#endif

  comm {new interpreter::Communicator {this}},
  error {new interpreter::Error {this}},
  output {new interpreter::Output {this}},
  input {new interpreter::Input {this}},
  object_handler {new interpreter::Object_handler {this}},
  object_container {new interpreter::Object_container {this}},
  object_creator {new interpreter::Object_creator {this}},  
  in {std::cin.rdbuf()},
  out {std::cout.rdbuf()},
  err {std::cerr.rdbuf()},
  log_flag {true},
  out_flag {true},
  err_flag {true},
  argc{argc},
  argv{argv} {
    if (comm->me == 0) log.open ("log");

    // this value is set to '1' because one input class is already constructed
    interpreter_num_Input_class = 1;

    interpreter_break_called = false;
    interpreter_continue_called = false;
}

#if defined(CAVIAR_WITH_MPI)

#else
CAVIAR::CAVIAR (std::string str) :

  comm {new interpreter::Communicator {this}},
  error {new interpreter::Error {this}},
  output {new interpreter::Output {this}},
  input {new interpreter::Input {this}},
  object_handler {new interpreter::Object_handler {this}},
  object_container {new interpreter::Object_container {this}},
  object_creator {new interpreter::Object_creator {this}},  
  in {std::cin.rdbuf()},
  out {std::cout.rdbuf()},
  err {std::cerr.rdbuf()},
  log_flag {true},
  out_flag {true},
  err_flag {true},
  construct_str{str}
  {
    if (comm->me == 0) log.open ("log");

    // this value is set to '1' because one input class is already constructed
    interpreter_num_Input_class = 1;

    interpreter_break_called = false;
    interpreter_continue_called = false;
}
#endif

CAVIAR::~CAVIAR () {
  delete input;
  delete output;
  delete error;
  delete comm;
  delete object_handler;  
  delete object_container;    
  delete object_creator;    
}

void CAVIAR::execute () {
  std::string greeting = "CAVIAR-";
  greeting += std::to_string (CAVIAR_MAJOR_VERSION);
  greeting += ".";
  greeting += std::to_string (CAVIAR_MINOR_VERSION);
  greeting += ".";
  greeting += std::to_string (CAVIAR_PATCH_VERSION);
    
  output->info(greeting);
    
  input->read ();
}


} // namespace caviar

#define BOOST_BIND_GLOBAL_PLACEHOLDERS 1 // To ignore the warning: note: #pragma message: The practice of declaring the Bind placeholders (_1, _2, ...) in the global namespace is deprecated. Please use <boost/bind/bind.hpp> + using namespace boost::placeholders, or define BOOST_BIND_GLOBAL_PLACEHOLDERS to retain the current behavior.


#include <boost/python.hpp>
//#include <boost/python/list.hpp>
//#include <boost/python/extract.hpp>

using namespace boost::python;
#include "caviar/objects/unique.h"

#include "caviar/objects/unique/all.h"
#include "caviar/objects/atom_data/all.h"
#include "caviar/objects/domain/all.h"
#include "caviar/objects/neighborlist/all.h"
#include "caviar/objects/neighborlist/all.h"
#include "caviar/objects/integrator/all.h"
#include "caviar/objects/writer/all.h"
#include "caviar/objects/md_simulator/all.h"
#include "caviar/objects/force_field/all.h"

void  export_unique();
void  export_atom_data();
void  export_force_field();
void  export_domain();
void  export_writer();
void  export_md_simulator();
void  export_neighborlist();

BOOST_PYTHON_MODULE(caviarmd)
{
    class_<caviar::CAVIAR,boost::noncopyable>("caviar", init<std::string>())
        .def("execute", &caviar::CAVIAR::execute)

    ;


    //class_<caviar::objects::Unique>("unique",init<caviar::CAVIAR*>())
    //    .def("verify", &caviar::objects::Unique::verify_settings)
    //;

    export_unique();
    export_atom_data();
    export_force_field();
    export_writer();
    export_neighborlist();
    export_md_simulator();
    export_domain();

};


void export_unique()
{
    namespace bp = boost::python;
    // map the unique namespace to a sub-module
    // make "from mypackage.unique import <whatever>" work
    bp::object uniqueModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.unique"))));
    // make "from mypackage import unique" work
    bp::scope().attr("unique") = uniqueModule;
    // set the current scope to the new sub-module
    bp::scope unique_scope = uniqueModule;
    // export stuff in the unique namespace

    class_<caviar::objects::unique::Atom>("Atom",init<caviar::CAVIAR*>())
     
    ;
    // etc.
}


void export_atom_data()
{

    namespace bp = boost::python;
    // map the atom_data namespace to a sub-module
    // make "from mypackage.atom_data import <whatever>" work
    bp::object atom_dataModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.atom_data"))));
    // make "from mypackage import atom_data" work
    bp::scope().attr("atom_data") = atom_dataModule;
    // set the current scope to the new sub-module
    bp::scope atom_data_scope = atom_dataModule;
    // export stuff in the atom_data namespace

    class_<caviar::objects::atom_data::Basic>("Basic",init<caviar::CAVIAR*>())
     
    ;
    // etc.

}


void export_force_field()
{
    namespace bp = boost::python;
    // map the force_field namespace to a sub-module
    // make "from mypackage.force_field import <whatever>" work
    bp::object force_fieldModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.force_field"))));
    // make "from mypackage import force_field" work
    bp::scope().attr("force_field") = force_fieldModule;
    // set the current scope to the new sub-module
    bp::scope force_field_scope = force_fieldModule;
    // export stuff in the force_field namespace

    class_<caviar::objects::force_field::Lj>("Lj",init<caviar::CAVIAR*>())
     
    ;
    // etc.
}


void export_neighborlist()
{
    namespace bp = boost::python;
    // map the neighborlist namespace to a sub-module
    // make "from mypackage.neighborlist import <whatever>" work
    bp::object neighborlistModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.neighborlist"))));
    // make "from mypackage import neighborlist" work
    bp::scope().attr("neighborlist") = neighborlistModule;
    // set the current scope to the new sub-module
    bp::scope neighborlist_scope = neighborlistModule;
    // export stuff in the neighborlist namespace

    class_<caviar::objects::neighborlist::Verlet_list>("Verlet_list",init<caviar::CAVIAR*>())
     
    ;

    class_<caviar::objects::neighborlist::Cell_list>("Cell_list",init<caviar::CAVIAR*>())
    ;
    // etc.
}


void export_md_simulator()
{
    namespace bp = boost::python;
    // map the md_simulator namespace to a sub-module
    // make "from mypackage.md_simulator import <whatever>" work
    bp::object md_simulatorModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.md_simulator"))));
    // make "from mypackage import md_simulator" work
    bp::scope().attr("md_simulator") = md_simulatorModule;
    // set the current scope to the new sub-module
    bp::scope md_simulator_scope = md_simulatorModule;
    // export stuff in the md_simulator namespace

    class_<caviar::objects::md_simulator::Basic>("Basic",init<caviar::CAVIAR*>())
     
    ;
    // etc.
}


void export_domain()
{
    namespace bp = boost::python;
    // map the domain namespace to a sub-module
    // make "from mypackage.domain import <whatever>" work
    bp::object domainModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.domain"))));
    // make "from mypackage import domain" work
    bp::scope().attr("domain") = domainModule;
    // set the current scope to the new sub-module
    bp::scope domain_scope = domainModule;
    // export stuff in the domain namespace

    class_<caviar::objects::domain::Box>("Box",init<caviar::CAVIAR*>())
     
    ;
    // etc.
}


void export_writer()
{

    namespace bp = boost::python;
    // map the writer namespace to a sub-module
    // make "from mypackage.writer import <whatever>" work
    bp::object writerModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.writer"))));
    // make "from mypackage import writer" work
    bp::scope().attr("writer") = writerModule;
    // set the current scope to the new sub-module
    bp::scope writer_scope = writerModule;
    // export stuff in the writer namespace

    //class_<caviar::objects::writer::Atom_dataW>("Atom_data",init<caviar::CAVIAR*>())     
    //;

    class_<caviar::objects::writer::Force_field>("Force_field",init<caviar::CAVIAR*>())     
    ;
    // etc.

}
