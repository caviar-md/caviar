
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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#include <boost/python/list.hpp>
//#include <boost/python/extract.hpp>

/*
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter
{
  /// @note Registers converter from a python interable type to the
  ///       provided type.
  template <typename Container>
  iterable_converter&
  from_python()
  {
    boost::python::converter::registry::push_back(
      &iterable_converter::convertible,
      &iterable_converter::construct<Container>,
      boost::python::type_id<Container>());

    // Support chaining.
    return *this;
  }

  /// @brief Check if PyObject is iterable.
  static void* convertible(PyObject* object)
  {
    return PyObject_GetIter(object) ? object : NULL;
  }

  /// @brief Convert iterable PyObject to C++ container type.
  ///
  /// Container Concept requirements:
  ///
  ///   * Container::value_type is CopyConstructable.
  ///   * Container can be constructed and populated with two iterators.
  ///     I.e. Container(begin, end)
  template <typename Container>
  static void construct(
    PyObject* object,
    boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    namespace python = boost::python;
    // Object is a borrowed reference, so create a handle indicting it is
    // borrowed for proper reference counting.
    python::handle<> handle(python::borrowed(object));

    // Obtain a handle to the memory block that the converter has allocated
    // for the C++ type.
    typedef python::converter::rvalue_from_python_storage<Container>
                                                                storage_type;
    void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

    typedef python::stl_input_iterator<typename Container::value_type>
                                                                    iterator;

    // Allocate the C++ type into the converter's memory block, and assign
    // its handle to the converter's convertible variable.  The C++
    // container is populated by passing the begin and end iterators of
    // the python object to the container's constructor.
    new (storage) Container(
      iterator(python::object(handle)), // begin
      iterator());                      // end
    data->convertible = storage;
  }
};


*/

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
void  export_domain();
void  export_atom_data();
/*

void  export_force_field();

void  export_writer();
void  export_md_simulator();
void  export_neighborlist();
*/

BOOST_PYTHON_MODULE(caviarmd)
{

    class_<caviar::CAVIAR,boost::noncopyable>("caviar", init<std::string>())
        .def("execute", &caviar::CAVIAR::execute)

    ;
    /*
    // used for std::vector get() ?? 
    class_<std::vector<double> >("double_vector")
        .def(vector_indexing_suite<std::vector<double> >())
    ;
    class_<std::vector<int> >("int_vector")
        .def(vector_indexing_suite<std::vector<int> >())
    ;
    class_<std::vector<unsigned int> >("uint_vector")
        .def(vector_indexing_suite<std::vector<unsigned int> >())
    ;


 // Register interable conversions. Used in setters
  iterable_converter()
    // Build-in type.
    .from_python<std::vector<double> >()
    .from_python<std::vector<int> >()
    .from_python<std::vector<unsigned int> >()
    // Each dimension needs to be convertable.
    .from_python<std::vector<std::string> >()
    .from_python<std::vector<std::vector<std::string> > >()
    // User type.
    //.from_python<std::list<foo> >()
    ;


     */
    export_unique();
    export_domain();
    export_atom_data();
    //export_force_field();
    //export_writer();
    //export_neighborlist();
    //export_md_simulator();


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

    caviar::objects::unique::export_py_Atom ();

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

    caviar::objects::domain::export_py_Box ();

}




void export_atom_data()
{

    namespace bp = boost::python;


 class_<caviar::objects::Atom_data,std::shared_ptr<caviar::objects::Atom_data>, boost::noncopyable>("Atom_data",no_init); // necessary line

    //class_<Derived,std::shared_ptr<Derived>,bases<AbstractBase>,boost::noncopyable>("Derived");


    // map the atom_data namespace to a sub-module
    // make "from mypackage.atom_data import <whatever>" work
    bp::object atom_dataModule(bp::handle<>(bp::borrowed(PyImport_AddModule("caviarmd.atom_data"))));
    // make "from mypackage import atom_data" work
    bp::scope().attr("atom_data") = atom_dataModule;
    // set the current scope to the new sub-module
    bp::scope atom_data_scope = atom_dataModule;
    // export stuff in the atom_data namespace

    caviar::objects::atom_data::export_py_Basic ();

}

/*

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

    caviar::objects::force_field::export_py_Lj ();


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

    caviar::objects::neighborlist::export_py_Verlet_list ();
    caviar::objects::neighborlist::export_py_Cell_list ();

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

    caviar::objects::md_simulator::export_py_Basic ();
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



}
*/
