
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

#include "caviar/interpreter/object_creator.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
bool Object_creator::string_variable (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);

  std::string NAME = "";
  bool in_file = true;
  std::string r = "";


// there are four different cases of creation of a real varable:
// 'string '
// 'string NAME'
// 'string NAME SOMETHING_INT '
// 'string NAME=SOMETHING_INT '

  auto t = parser->get_val_token();
  // 'string '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'string NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // string NAME = SOMETHING_STRING
    r = parser -> get_string();
  } else {
    // string NAME SOMETHING_STRING
    parser -> keep_current_token();
    r = parser -> get_string();
  }

  int index = object_container->string_variable.size ();
  object_container->string_variable.emplace_back (r);
  object_handler::Dictionary dict (object_handler::gdst("string_variable"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  return in_file;

}

bool Object_creator::boolean_variable (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);

  std::string NAME = "";
  bool in_file = true;
  bool r = 0;

// there are four different cases of creation of a boolean varable:
// 'bool '
// 'bool NAME'
// 'bool NAME SOMETHING_bool '
// 'bool NAME=SOMETHING_bool '

  auto t = parser->get_val_token();
  // 'bool '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'bool NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // bool NAME = SOMETHING_bool
    r = parser -> get_bool();
  } else {
    // bool NAME SOMETHING_bool
    parser -> keep_current_token();
    r = parser -> get_bool();
  }

  int index = object_container->boolean_variable.size ();
  object_container->boolean_variable.emplace_back (r);
  object_handler::Dictionary dict (object_handler::gdst("boolean_variable"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  return in_file;

}

bool Object_creator::int_variable (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);

  std::string NAME = "";
  bool in_file = true;
  int r = 0;

// there are four different cases of creation of a real varable:
// 'int '
// 'int NAME'
// 'int NAME SOMETHING_INT '
// 'int NAME=SOMETHING_INT '

  auto t = parser->get_val_token();
  // 'int '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'int NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // int NAME = SOMETHING_INT
    r = parser -> get_int();
  } else {
    // int NAME SOMETHING_INT
    parser -> keep_current_token();
    r = parser -> get_int();
  }

  int index = object_container->int_variable.size ();
  object_container->int_variable.emplace_back (r);
  object_handler::Dictionary dict (object_handler::gdst("int_variable"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  return in_file;

}

bool Object_creator::real_variable (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);
  std::string NAME = "";
  bool in_file = true;
  double r = 0;
// there are four different cases of creation of a real varable:
// 'real '
// 'real NAME'
// 'real NAME SOMETHING_REAL '
// 'real NAME=SOMETHING_REAL '

  auto t = parser->get_val_token();
  // 'real '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'real NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // real NAME = SOMETHING_REAL
    r = parser -> get_real();
  } else {
    // real NAME SOMETHING_REAL
    parser -> keep_current_token();
    r = parser -> get_real();
  }

  int index = object_container->real_variable.size ();
  object_container->real_variable.emplace_back (r);

  object_handler::Dictionary dict (object_handler::gdst("real_variable"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));
  
  return in_file; 
}


bool Object_creator::int_2d_vector (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);
  std::string NAME = "";

  bool in_file = true;
  Vector2D<int> r{0,0};

// there are four different cases of creation of a int_3d_vector:
// 'int3d '
// 'int3d NAME'
// 'int3d NAME SOMETHING.x SOMETHING.y SOMETHING.z '
// 'int3d NAME=SOMETHING.x SOMETHING.y SOMETHING.z  '

  auto t = parser->get_val_token();
  // 'int3d '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'int3d NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // int3d NAME = SOMETHING.x SOMETHING.y SOMETHING.z 
    r.x = parser -> get_int();
    r.y = parser -> get_int();
  } else {
    // int3d NAME SOMETHING.x SOMETHING.y SOMETHING.z 
    parser -> keep_current_token();
    r.x = parser -> get_int();
    r.y = parser -> get_int();
  }

  int index = object_container->int_2d_vector.size ();
  object_container->int_2d_vector.emplace_back (r);

  object_handler::Dictionary dict (object_handler::gdst("int_2d_vector"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  //std::cout <<"i3d: " << v << std::endl;
  return in_file; //WARNING
}


bool Object_creator::real_2d_vector (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);
  std::string NAME = "";

  bool in_file = true;
  Vector2D<double> r{0,0};

// there are four different cases of creation of a real_3d_vector:
// 'real3d '
// 'real3d NAME'
// 'real3d NAME SOMETHING.x SOMETHING.y SOMETHING.z '
// 'real3d NAME=SOMETHING.x SOMETHING.y SOMETHING.z  '

  auto t = parser->get_val_token();
  // 'real3d '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'real3d NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // real3d NAME = SOMETHING.x SOMETHING.y SOMETHING.z 
    r.x = parser -> get_real();
    r.y = parser -> get_real();
  } else {
    // real3d NAME SOMETHING.x SOMETHING.y SOMETHING.z 
    parser -> keep_current_token();
    r.x = parser -> get_real();
    r.y = parser -> get_real();
  }


  int index = object_container->real_2d_vector.size ();
  object_container->real_2d_vector.emplace_back (r);

  object_handler::Dictionary dict (object_handler::gdst("real_2d_vector"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  //std::cout <<"r3d: " << v << std::endl;
  return in_file; //WARNING
}


bool Object_creator::int_3d_vector (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);
  std::string NAME = "";

  bool in_file = true;
  Vector<int> r{0,0,0};

// there are four different cases of creation of a int_3d_vector:
// 'int3d '
// 'int3d NAME'
// 'int3d NAME SOMETHING.x SOMETHING.y SOMETHING.z '
// 'int3d NAME=SOMETHING.x SOMETHING.y SOMETHING.z  '

  auto t = parser->get_val_token();
  // 'int3d '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'int3d NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // int3d NAME = SOMETHING.x SOMETHING.y SOMETHING.z 
    r.x = parser -> get_int();
    r.y = parser -> get_int();
    r.z = parser -> get_int();
  } else {
    // int3d NAME SOMETHING.x SOMETHING.y SOMETHING.z 
    parser -> keep_current_token();
    r.x = parser -> get_int();
    r.y = parser -> get_int();
    r.z = parser -> get_int();
  }

  int index = object_container->int_3d_vector.size ();
  object_container->int_3d_vector.emplace_back (r);

  object_handler::Dictionary dict (object_handler::gdst("int_3d_vector"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  //std::cout <<"i3d: " << v << std::endl;
  return in_file; //WARNING
}


bool Object_creator::real_3d_vector (Parser * parser) {
  std::string this_object_full_type = __func__;
  output->info_create (this_object_full_type);
  std::string NAME = "";

  bool in_file = true;
  Vector<double> r{0,0,0};

// there are four different cases of creation of a real_3d_vector:
// 'real3d '
// 'real3d NAME'
// 'real3d NAME SOMETHING.x SOMETHING.y SOMETHING.z '
// 'real3d NAME=SOMETHING.x SOMETHING.y SOMETHING.z  '

  auto t = parser->get_val_token();
  // 'real3d '
  if (t.kind == caviar::interpreter::Kind::eol) {return true;}
  if (t.kind == caviar::interpreter::Kind::eof) {return false;}

  // 'real3d NAME'
  if (t.kind == caviar::interpreter::Kind::identifier) {
    NAME_ASSIGN_CHECK(t)
    NAME = t.string_value;
    object_container->all_names.insert(NAME);
  }

  t = parser->get_raw_token();
  if (t.kind == caviar::interpreter::Kind::eol) {

  } else if (t.kind == caviar::interpreter::Kind::eof) {
    in_file = false;
  } else if (t.kind == caviar::interpreter::Kind::assign) {
    // real3d NAME = SOMETHING.x SOMETHING.y SOMETHING.z 
    r.x = parser -> get_real();
    r.y = parser -> get_real();
    r.z = parser -> get_real();
  } else {
    // real3d NAME SOMETHING.x SOMETHING.y SOMETHING.z 
    parser -> keep_current_token();
    r.x = parser -> get_real();
    r.y = parser -> get_real();
    r.z = parser -> get_real();
  }


  int index = object_container->real_3d_vector.size ();
  object_container->real_3d_vector.emplace_back (r);

  object_handler::Dictionary dict (object_handler::gdst("real_3d_vector"), index);  
  object_container->dictionary.insert (std::make_pair(NAME,dict));

  //std::cout <<"r3d: " << v << std::endl;
  return in_file; //WARNING
}
} //interpreter
CAVIAR_NAMESPACE_CLOSE


