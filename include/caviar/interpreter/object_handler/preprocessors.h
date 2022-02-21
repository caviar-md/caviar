
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_PREPROCESSORS_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_PREPROCESSORS_H

// ==================== 
// ==================== MAKE_ERROR_MASSAGE
// ====================

#define MAKE_MASSAGE_CUSTOM(ERR,MASSAGE_STRING_START,OBJECT_TYPE,MASSAGE_STRING_END) \
  std::string ERR = MASSAGE_STRING_START;\
  ERR.append(#OBJECT_TYPE);  \
  ERR.append(MASSAGE_STRING_END);

#define MAKE_ERROR_MASSAGE_EXPECTED(ERR,OBJECT_TYPE,MASSAGE_STRING_END) \
  std::string ERR = "expected a (an)";\
  ERR.append(#OBJECT_TYPE);  \
  ERR.append(MASSAGE_STRING_END);
// ==================== 
// ==================== GET_A_TOKEN
// ====================

#define GET_A_TOKEN_FOR_CREATION \
  auto token = parser->get_val_token(); \
  if (token.kind == caviar::interpreter::Kind::eol) break; \
  if (token.kind == caviar::interpreter::Kind::eof) {in_file = false; break;} 

#define GET_A_TOKEN_FOR_CREATION_NAME(VARIABLE) \
  auto VARIABLE = parser->get_val_token(); \
  if (VARIABLE.kind == caviar::interpreter::Kind::eol) break; \
  if (VARIABLE.kind == caviar::interpreter::Kind::eof) {in_file = false; break;} 
  
#define GET_A_TOKEN_FOR_CREATION_NAME_RETURN(VARIABLE) \
  auto VARIABLE = parser->get_val_token(); \
  if (VARIABLE.kind == caviar::interpreter::Kind::eol) return true; \
  if (VARIABLE.kind == caviar::interpreter::Kind::eof) return true; 
// ==================== 
// ==================== NAME ASSIGN
// ====================
//  error->all (FC_FILE_LINE_FUNC_PARSE, "This NAME is reserved as object_creator command");
//  if (fptr->input->commands_map.count (name) == 1) 

#define NAME_ASSIGN_CHECK(VARIABLE) \
{\
  if (VARIABLE.kind != caviar::interpreter::Kind::identifier)  \
    error->all (FC_FILE_LINE_FUNC_PARSE, "expected string for NAME");\
  auto name = VARIABLE.string_value;\
  if (commands_map.count (name) == 1) \
      error->all (FC_FILE_LINE_FUNC_PARSE, "This NAME is reserved as input command");\
  if (object_container->all_names.count (name) == 1) \
    error->all (FC_FILE_LINE_FUNC_PARSE, "This NAME is used before and cannot be re-assigned to other objects");\
}


#define ASSIGN_NAME \
  if (!name_called) {\
    name_called = true;\
    NAME_ASSIGN_CHECK(token)\
    NAME = token.string_value;\
    object_container->all_names.insert(NAME);\
  }


#define ASSIGN_NAME_WITHOUT_KEYWORD \
  if (name_called==false) { \
    NAME_ASSIGN_CHECK(token)  \
    NAME = token.string_value; \
    name_called = true; \
  } else { \
    error->all (FC_FILE_LINE_FUNC_PARSE, " Unknown variable or command"); \
  } 

// ===========================
// ===========================
// ===========================

#define GET_A_BOOL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::string st = token.string_value; \
    if (string_cmp(st,"off") || string_cmp(st,"false")) \
      VARIABLE = false; \
    if (string_cmp(st,"on") || string_cmp(st,"true")) \
      VARIABLE = true; \
  } else if (token.kind == caviar::interpreter::Kind::int_number) {\
    if (token.int_value == 0) \
      VARIABLE = false; \
    if (token.int_value == 1) \
      VARIABLE = true; \
  } else { \
    std::string error_massage = " ";\
    error_massage.append(FUNCTION_NAME);  \
    error_massage.append(ERROR_MASSAGE); \
    error_massage.append(#VARIABLE);  \
    error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
  }

#define GET_A_BOOL(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
{\
  auto token = parser->get_val_token(); \
  GET_A_BOOL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
}

// ===========================
// ===========================
// ===========================


#define GET_A_STRING_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.kind != caviar::interpreter::Kind::identifier) {\
    std::string error_massage = " ";\
    error_massage.append(FUNCTION_NAME);  \
    error_massage.append(ERROR_MASSAGE); \
    error_massage.append(#VARIABLE);  \
    error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
  }\
  VARIABLE = token.string_value; \

#define GET_A_STRING(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
{\
  auto token = parser->get_val_token(); \
  GET_A_STRING_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
}


#define CHECK_NAME_EXISTANCE(VARIABLE,ITERATOR,FUNCTION_NAME,ERROR_MASSAGE) \
{ \
  ITERATOR = object_container->dictionary.find(VARIABLE); \
  if (ITERATOR == object_container->dictionary.end())   {\
    std::string error_massage = " ";\
    error_massage.append(FUNCTION_NAME);  \
    error_massage.append(ERROR_MASSAGE); \
    error_massage.append("Undefined object NAME: ");  \
    error_massage.append(VARIABLE);  \
    error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
  }\
}
// NNT : NO_NEW_TOKEN

#define GET_A_REAL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
if (token.kind == caviar::interpreter::Kind::real_number)\
  VARIABLE = token.real_value;\
else if (token.kind == caviar::interpreter::Kind::int_number)\
  VARIABLE = token.int_value;\
else {\
    std::string error_massage = " ";\
    error_massage.append(FUNCTION_NAME);  \
    error_massage.append(ERROR_MASSAGE); \
    error_massage.append("Expected a real number as ");  \
    error_massage.append(#VARIABLE);  \
    error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
  }\



#define GET_A_REAL(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
{\
  auto token = parser->get_val_token(); \
  GET_A_REAL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}

#define GET_OR_CHOOSE_A_REAL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable")) { \
        VARIABLE = object_container->int_variable[it->second.index];\
      } else if (it->second.type == caviar::interpreter::object_handler::gdst("real_variable")) { \
        VARIABLE = object_container->real_variable[it->second.index];\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(token.string_value);\
        error_massage.append(" is not referred to a real or int number ");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
  GET_A_REAL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  }

#define GET_OR_CHOOSE_A_REAL(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_REAL_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}


#define GET_A_INT_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind != caviar::interpreter::Kind::int_number) {\
    std::string error_massage = " ";\
    error_massage.append(FUNCTION_NAME);  \
    error_massage.append(ERROR_MASSAGE); \
    error_massage.append("Expected an int number as ");  \
    error_massage.append(#VARIABLE);  \
    error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
  }\
  VARIABLE = token.int_value; 


#define GET_A_INT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_A_INT_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}



#define GET_OR_CHOOSE_A_INT_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable")) { \
        VARIABLE = object_container->int_variable[it->second.index];\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(token.string_value);\
        error_massage.append(" is not referred to a int number ");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
    GET_A_INT_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  }

#define GET_OR_CHOOSE_A_INT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_INT_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}


// ======================================
// ======================================  2D VECTORS
// ======================================

#define GET_OR_CHOOSE_A_REAL_2D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("real_2d_vector")) { \
        VARIABLE = object_container->real_2d_vector[it->second.index];\
      } else if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable") || it->second.type == caviar::interpreter::object_handler::gdst("real_variable")){ \
        auto &vx = VARIABLE.x;\
        auto &vy = VARIABLE.y;\
        GET_OR_CHOOSE_A_REAL_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_REAL(vy,FUNCTION_NAME,ERROR_MASSAGE)\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(" expected real number or a 2D vector NAME");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
    auto &vx = VARIABLE.x;\
    auto &vy = VARIABLE.y;\
    GET_OR_CHOOSE_A_REAL_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_REAL(vy,FUNCTION_NAME,ERROR_MASSAGE)\
  }



#define GET_OR_CHOOSE_A_REAL_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_REAL_2D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}


#define GET_OR_CHOOSE_A_INT_2D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("int_2d_vector")) { \
        VARIABLE = object_container->int_2d_vector[it->second.index];\
      } else if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable") ){ \
        auto &vx = VARIABLE.x;\
        auto &vy = VARIABLE.y;\
        GET_OR_CHOOSE_A_INT_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_INT(vy,FUNCTION_NAME,ERROR_MASSAGE)\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(" expected int number or a 2D vector NAME");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
    auto &vx = VARIABLE.x;\
    auto &vy = VARIABLE.y;\
    GET_OR_CHOOSE_A_INT_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_INT(vy,FUNCTION_NAME,ERROR_MASSAGE)\
  }


#define GET_OR_CHOOSE_A_INT_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_INT_2D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}


// ======================================
// ======================================  3D VECTORS
// ======================================

#define GET_OR_CHOOSE_A_REAL_3D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("real_3d_vector")) { \
        VARIABLE = object_container->real_3d_vector[it->second.index];\
      } else if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable") || it->second.type == caviar::interpreter::object_handler::gdst("real_variable")){ \
        auto &vx = VARIABLE.x;\
        auto &vy = VARIABLE.y;\
        auto &vz = VARIABLE.z;\
        GET_OR_CHOOSE_A_REAL_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_REAL(vy,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_REAL(vz,FUNCTION_NAME,ERROR_MASSAGE)\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(" expected real number or a 3D vector NAME");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
    auto &vx = VARIABLE.x;\
    auto &vy = VARIABLE.y;\
    auto &vz = VARIABLE.z;\
    GET_OR_CHOOSE_A_REAL_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_REAL(vy,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_REAL(vz,FUNCTION_NAME,ERROR_MASSAGE)\
  }



#define GET_OR_CHOOSE_A_REAL_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_REAL_3D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}


#define GET_OR_CHOOSE_A_INT_3D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
  if (token.kind == caviar::interpreter::Kind::identifier) {\
    std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;\
    it = object_container->dictionary.find(token.string_value);\
    if (it == object_container->dictionary.end()) { \
      std::string error_massage = " ";\
      error_massage.append(FUNCTION_NAME);  \
      error_massage.append(ERROR_MASSAGE); \
      error_massage.append("Invalid object NAME ");  \
      error_massage.append(#VARIABLE);  \
      error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
    } else {\
      if (it->second.type == caviar::interpreter::object_handler::gdst("int_3d_vector")) { \
        VARIABLE = object_container->int_3d_vector[it->second.index];\
      } else if (it->second.type == caviar::interpreter::object_handler::gdst("int_variable") ){ \
        auto &vx = VARIABLE.x;\
        auto &vy = VARIABLE.y;\
        auto &vz = VARIABLE.z;\
        GET_OR_CHOOSE_A_INT_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_INT(vy,FUNCTION_NAME,ERROR_MASSAGE)\
        GET_OR_CHOOSE_A_INT(vz,FUNCTION_NAME,ERROR_MASSAGE)\
      } else {\
        std::string error_massage = " ";\
        error_massage.append(FUNCTION_NAME);  \
        error_massage.append(ERROR_MASSAGE); \
        error_massage.append(" expected int number or a 3D vector NAME");  \
        error->all (FC_FILE_LINE_FUNC_PARSE, error_massage );\
      }\
    }\
  } else {\
    auto &vx = VARIABLE.x;\
    auto &vy = VARIABLE.y;\
    auto &vz = VARIABLE.z;\
    GET_OR_CHOOSE_A_INT_NNT(vx,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_INT(vy,FUNCTION_NAME,ERROR_MASSAGE)\
    GET_OR_CHOOSE_A_INT(vz,FUNCTION_NAME,ERROR_MASSAGE)\
  }


#define GET_OR_CHOOSE_A_INT_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
{\
  auto token = parser->get_val_token(); \
  GET_OR_CHOOSE_A_INT_3D_VECTOR_NNT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE)\
}

//===========
//=========== VARIABLE ASSIGN
//===========
#define ASSIGN_STRING(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_A_STRING(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 
  
#define ASSIGN_REAL(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_REAL(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 

#define ASSIGN_INT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_INT(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 

#define ASSIGN_REAL_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_REAL_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 

#define ASSIGN_INT_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_INT_2D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 



#define ASSIGN_REAL_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_REAL_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 

#define ASSIGN_INT_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  if (token.string_value == #VARIABLE) { \
    GET_OR_CHOOSE_A_INT_3D_VECTOR(VARIABLE,FUNCTION_NAME,ERROR_MASSAGE) \
  } 
  
//===========
//=========== MATRIX ASSIGN
//===========

#define GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(vector_name)\
  unsigned int i = 0 , j = 0;\
  double vector_value = 0;\
  GET_OR_CHOOSE_A_INT(i,"","")\
  GET_OR_CHOOSE_A_INT(j,"","")\
  GET_OR_CHOOSE_A_REAL(vector_value,"","")\
  if (vector_name.size() < i+1) vector_name.resize(i+1);\
  if (vector_name[i].size() < j+1) vector_name[i].resize(j+1);\
  vector_name[i][j] = vector_value;
  
#define GET_A_STDVECTOR_REAL_ELEMENT(vector_name)\
  unsigned int i = 0;\
  double vector_value = 0;\
  GET_OR_CHOOSE_A_INT(i,"","")\
  GET_OR_CHOOSE_A_REAL(vector_value,"","")\
  if (vector_name.size() < i+1) vector_name.resize(i+1);\
  vector_name[i] = vector_value;

  
#define GET_A_STDVECTOR_INT_ELEMENT(vector_name)\
  unsigned int i = 0;\
  int vector_value = 0;\
  GET_OR_CHOOSE_A_INT(i,"","")\
  GET_OR_CHOOSE_A_INT(vector_value,"","")\
  if (vector_name.size() < i+1) vector_name.resize(i+1);\
  vector_name[i] = vector_value;

//===========
//=========== FIND OBJECT POINTER
//===========
// NO ITERATOR WILL BE CREATED
#define FIND_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME) \
  std::string name_to_find_;\
  MAKE_ERROR_MASSAGE_EXPECTED(err,OBJECT_TYPE,"name.")\
  GET_A_STRING(name_to_find_,"", err)\
  CHECK_NAME_EXISTANCE(name_to_find_, ITERATOR_NAME, "","")\
  if (ITERATOR_NAME->second.type != caviar::interpreter::object_handler::gdst( #OBJECT_TYPE ))\
    error->all(FC_FILE_LINE_FUNC_PARSE,": undefined object. ");

#define FIND_OBJECT_BY_NAME(OBJECT_TYPE,ITERATOR_NAME) \
  std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator ITERATOR_NAME;\
  FIND_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME)
/*
  std::string name_to_find_;
  MAKE_ERROR_MASSAGE_EXPECTED(err,OBJECT_TYPE,"name.")
  GET_A_STRING(name_to_find_,"", err)
  CHECK_NAME_EXISTANCE(name_to_find_, ITERATOR_NAME, "","")
  if (ITERATOR_NAME->second.type != caviar::interpreter::object_handler::gdst( #OBJECT_TYPE ))
    error->all(FC_FILE_LINE_FUNC_PARSE,": undefined object. ");
*/

        
#endif
