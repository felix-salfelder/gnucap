/*                               -*- C++ -*-
 * Copyright (C) 2007 Albert Davis
 * Author: Albert Davis <aldavis@gnu.org>
 *
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 */
#include "e_cardlist.h"
#include "c_comand.h"
#include "constant.h"
#include "globals.h"
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <sys/stat.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/wait.h>

// #include "boost/filesystem/operations.hpp"
// #include "boost/filesystem/fstream.hpp"
// namespace fs = boost::filesystem;
/*--------------------------------------------------------------------------*/
namespace {

#if 0
  bool search_file( std::string &name ){
    const char* h ="HOME";
    const char* home= getenv(h);

    if (name[0] == '/') return true; // fixme

    std::string pathlist[4] = { OPT::libpath , "./",  LIBDIR,
      std::string(home) + std::string("/.gnucap/lib/") };

    // FIXME. use libpath 

    for(int i=1; i<4 ; i++) {
      if ( FILE* tmp = fopen( (pathlist[i] + "/" + name).c_str(), "r" ) ) {
        fclose(tmp);
        name = pathlist[i]+"/"+name;
        return true;
      }
      trace0( (" not found " + pathlist[i] + "/" + name).c_str());
    }
    return false;

  }
#endif
/*--------------------------------------------------------------------------*/
using std::string;
/*--------------------------------------------------------------------------*/
std::map<const std::string, void*> attach_list;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class CMD_ATTACH : public CMD {
  static void compile(string& filename, string source, string make);
  static void* do_attach(string filename, int flags, bool force=false);
public:
  void do_it(CS& cmd, CARD_LIST*)
  {
    unsigned here = cmd.cursor();
    int dl_scope = RTLD_LOCAL;
    int check = RTLD_NOW;
    string make = OS::getenv("GNUCAP_MAKE", GNUCAP_MAKE);
    // RTLD_NOW means to resolve symbols on loading
    // RTLD_LOCAL means symbols defined in a plugin are local
    do {
      if (cmd.umatch("public ")) {
	dl_scope = RTLD_GLOBAL;
	// RTLD_GLOBAL means symbols defined in a plugin are global
	// Use this when a plugin depends on another.
      }else if (cmd.umatch("lazy|force")) {
	check = RTLD_LAZY;
	// RTLD_LAZY means to defer resolving symbols until needed
	// Use when a plugin will not load because of unresolved symbols,
	// but it may work without it.
      }else{
	Get(cmd,"make{file}", &make);
      }
    } while (cmd.more() && !cmd.stuck(&here));
    trace1("attach::do_it", make);

    string file_name;
    cmd >> file_name;

    OMSTREAM _out = IO::mstdout;
    _out.outset(cmd);

    if(file_name==""){
      string comma;
      for (std::map<std::string, void*>::iterator
	  ii = attach_list.begin(); ii != attach_list.end(); ++ii) {
	if (ii->second) {
	  _out << comma << ii->first;
	  comma = ",\n";
	}else{
	}
      }
      _out << "\n";
      return;
    }

    void*& handle = attach_list[file_name];
    trace2("...", file_name, handle);
    if (handle) {
      if (CARD_LIST::card_list.is_empty()) {
	cmd.warn(bDANGER, here, "\"" + file_name + "\": already loaded, replacing");
	dlclose(handle);
	handle = NULL;
      }else{untested();
	cmd.reset(here);
	throw Exception_CS("already loaded, cannot replace when there is a circuit", cmd);
      }
    }else{
    }

    string source_filename(file_name);
    // FIXME: incomplete... some more control...
    // global list of supported suffixes?
    if (file_name.size()>3 && !strcmp(file_name.c_str()+file_name.size()-3,".so")) {
      source_filename = "";
    }else if (file_name.size()>3 && file_name.c_str()[file_name.size()-3] == '.') {
      file_name[file_name.size()-2]='s';
      file_name[file_name.size()-1]='o';

      if(file_name[0]=='/') { itested();
      } else {
	char* cwd = get_current_dir_name(); // POSIX, no C++ implementation available
	source_filename = string(cwd) + "/" + source_filename;
	free(cwd);
      }
    } else {
      source_filename = "";
    }

    if (source_filename!="") {
      trace1("attach", source_filename);
      assert(source_filename[0]=='/');
      try {
	compile(file_name, source_filename, make);
      }catch(Exception& e){
	cmd.reset(here);
	throw Exception_CS(e.message(), cmd);
      }
    }else{
    }

    handle = dlopen(file_name.c_str(), check | dl_scope);
    const char* e = dlerror();
    if (check == RTLD_LAZY) {
    }else if (handle) {
      const char* (*name)() = (const char*(*)()) dlsym(handle, "interface_name");
      if (name){
      }else{
	dlclose(handle);
	handle = NULL;
	throw Exception_CS("missing interface", cmd);
      }
    }
    if (e){
      cmd.reset(here);
      throw Exception_CS(e, cmd);
    }
#if 0
    try {
      assert(!handle);
      handle = do_attach(file_name, check | dl_scope, force);
    } catch (Exception& e) {
      trace0("do_attach threw");
      cmd.reset(here);
      throw Exception_CS(e.message(), cmd);
    }
#endif
    trace0("done attach");
  }
} p1;
DISPATCHER<CMD>::INSTALL d1(&command_dispatcher, "attach|load", &p1);
/*--------------------------------------------------------------------------*/
#if 0
// overengineered gncap-uf approach
void* CMD_ATTACH::do_attach(string file_name, int flags, bool force)
{ untested();
  void* handle = dlopen(file_name.c_str(), flags);
  char* e = dlerror();
  if(e || !handle ) { untested();
    throw Exception("cannot attach (" + string(e) + ")");
  }else if (handle) { untested();
    const char* (*name)() = (const char*(*)()) dlsym(handle, "interface_name");
    e = dlerror();
    if (force) { untested();
      trace1("forced load...", file_name);
    } else if (e && !force) { untested();
      dlclose(handle);
      throw Exception(file_name + " lacks interface information");
    } else if ((e || !name) && !force) { untested();
      dlclose(handle); untested();
      throw Exception("lacks interface name");
    } else { untested();
    }

    unsigned (*version)() = (unsigned(*)()) dlsym(handle, "interface_version");
    e = dlerror();
    if (force) {
    } else if ((e || !version) && !force) { untested();
      dlclose(handle);
      throw Exception("lacks interface version");
    } else if (strcmp(name(), interface_name())) {
      string n(name());
      dlclose(handle);
      throw Exception(file_name + ": wrong interface ("+ n +
				", not " + string(interface_name()) + ")");
    } else if (interface_version() == version()) {
    } else if (HAVE_GIT_REPO) { untested();
      throw Exception("loading " + file_name + ": plugin (" + to_string(version()) +
           ") doesnt match git revision (" + to_string(interface_version()) + ")");
    } else { untested();
      throw Exception(string("plugin too ") + ((interface_version() < version())?"new":"old"));
    }
    assert(handle);
    return handle;
  }else{itested();
    throw Exception(dlerror());
  }
}
#endif
/*--------------------------------------------------------------------------*/
void CMD_ATTACH::compile(string &filename, string source_filename, string make)
{
  struct stat ccattrib;
  int ccstat = stat(source_filename.c_str(), &ccattrib);
  if (ccstat) {
    throw Exception("cannot compile: " + source_filename +
		    " does not exist (" + to_string(ccstat) + ")\n");
  } else {
  }

  struct stat soattrib;
  int sostat = stat(filename.c_str(), &soattrib);
  string gnucap_cppflags = OS::getenv("GNUCAP_CPPFLAGS", GNUCAP_CPPFLAGS);

  FILE* f = NULL;
  char* sf2 = strdup(filename.c_str());
  string f1(basename(sf2)); free(sf2);
  char* tmplt = NULL;
  char* t = NULL;

  if (filename[0]!='/') {
    filename = "./" + filename;
  }

  if (!sostat && ccattrib.st_mtime < soattrib.st_mtime) {
    trace0("so exists and is newer");
    return;
  } else if(!sostat) {
    trace0("so needs to be refreshed");
    f = fopen(filename.c_str(),"a");
    if (f) {
      trace0("so is writable");
      fclose(f);
    } else { untested();
      tmplt = strdup("/tmp/soXXXXXX");
      t = mkdtemp(tmplt);
      if (!t) throw Exception("cannot create temporary file");
      filename = string(t) + "/" + f1;
    }
  } else {
  }
  free(tmplt);
  t = NULL;

  char* sf1 = strdup(source_filename.c_str());
  char* fn1 = strdup(filename.c_str());
  string d1(dirname(fn1));
  string d2(dirname(sf1));

  int childret;
  pid_t p = vfork();
  if (p) {
    waitpid(p, &childret, 0);
  } else {
    error(bDEBUG, "calling " + make + " -C" + d1 + " VPATH=" + d2 + " " + f1 + "\n");
    int ret = execlp( make.c_str(), make.c_str(),
		      "-C", d1.c_str(),
		      ("GNUCAP_CPPFLAGS="+gnucap_cppflags+"").c_str(),
		      ("VPATH=" + d2).c_str(),
		      f1.c_str(), (char*) NULL);
    _exit(ret);
  }
  free(sf1);
  free(fn1);
//  if (t) { incomplete();
//    // rm -r t;
//  }

  if(childret){
    throw Exception("cannot make " + filename + "(" + to_string(childret) + ")");
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class CMD_DETACH : public CMD {
public:
  void do_it(CS& cmd, CARD_LIST*)
  {
    if (CARD_LIST::card_list.is_empty()) {
      unsigned here = cmd.cursor();
      std::string file_name;
      cmd >> file_name;
      
      void* handle = attach_list[file_name];
      if (handle) {
	dlclose(handle);
	attach_list[file_name] = NULL;
      }else{
	cmd.reset(here);
	throw Exception_CS("plugin not attached", cmd);
      }
    }else{itested();
      throw Exception_CS("detach prohibited when there is a circuit", cmd);
    }
  }
} p2;
DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "detach|unload", &p2);
/*--------------------------------------------------------------------------*/
class CMD_DETACH_ALL : public CMD {
public:
  void do_it(CS& cmd, CARD_LIST*)
  {
    if (CARD_LIST::card_list.is_empty()) {
      for (std::map<std::string, void*>::iterator
	     ii = attach_list.begin(); ii != attach_list.end(); ++ii) {
	void* handle = ii->second;
	if (handle) {
	  dlclose(handle);
	  ii->second = NULL;
	}else{
	}
      }
    } else {
      throw Exception_CS("detach prohibited when there is a circuit", cmd);
    }
  }
} p3;
DISPATCHER<CMD>::INSTALL d3(&command_dispatcher, "detach_all", &p3);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
