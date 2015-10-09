#define ADD_VERSION

#include <stdio.h>
#include <patchlev_fake_0.h>

static class versioncheck{
	public:
	versioncheck(){ untested();
		if(version_revision() > VERSION_REVISION){
			printf("bogus3: loading, although kernel revision higher\n");
	   }	
	}
} B;
