// do nothing but define bogus version
// need to define a function, otherwise dlopen fails, whether or not RTLD_LAZY

extern "C" {
	extern const char* VERY_DIFFERENT_INTERFACE();
	void verschk(){ VERY_DIFFERENT_INTERFACE(); }
};

#include "stdio.h"
static class versioncheck{
	public:
	versioncheck(){
		printf("trying to call undefined symbol\n");
		printf("bogus1: kernel %s undefined\n", VERY_DIFFERENT_INTERFACE());
		printf("unreachable\n");
	}
} a;

#include <stdio.h>
