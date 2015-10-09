#define ADD_VERSION

#include "io_trace.h"
// fake build against newer kernel revision
#include "patchlev_other.h"
#include <stdio.h>


static class versioncheck{
	public:
	versioncheck(){ untested();
		if(version_revision() < VERSION_REVISION){ untested();
			printf("bogus4: built against %d revision %d\n",
					VERSION_CURRENT, VERSION_REVISION);
	   }else{
		}
	}
} a;
