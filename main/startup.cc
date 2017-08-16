
// read system startup file
// read cwd startup
//   else read user startup
//
//   BUG: relative paths in include files do not work.

#include "io_trace.h"
#include "c_comand.h"
#include "io_.h"
#include "e_cardlist.h"

// hack...
#define CWDSTARTFILE   "gnucap.rc"


using std::string;

bool startup_recursive()
{
  char buf[PATH_MAX];
  char buf2[PATH_MAX];

  char* cwd = getcwd(buf, PATH_MAX);
  char* start = getcwd(buf2, PATH_MAX);
  size_t cur = strlen(cwd);

  for (;;) {
    trace2("..", cur, cwd);
    std::string f = findfile(CWDSTARTFILE, cwd, R_OK);

    if ("" != f) { itested();
      trace2("get", cwd, f);
      try{
        CMD::command("get " + f, &CARD_LIST::card_list);
      } catch(Exception e){
        error(bDANGER, "%s\n",e.message().c_str());
      }
      chdir(start);
      return true;
    }

    while (--cur > 0 && cwd[cur] != *ENDDIR);
    if (cur <= 0) {
      chdir(start);
      return false;
    }
    if (chdir("..")) { untested();
      chdir(start);
      throw Exception("Cannot change to '%s/..'", cwd);
    }
    cwd[cur] = '\0';

  }
}

/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
